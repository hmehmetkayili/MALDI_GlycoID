"""Spectrum peak preprocessing, filtering, and profile signal processing.

Strictly enforces the scientific boundary between centroid peak lists and
profile mass spectra:
- Profile-specific operations (baseline subtraction, local-maxima peak picking,
  numerical trapezoidal peak area integration, FWHM, and local noise S/N) are
  performed ONLY when the input spectrum is classified as `SpectrumDataType.PROFILE`.
- If a centroid or unclassified peak list is supplied, profile-only operations
  are skipped and an explicit warning is recorded.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np

from app.parsers.base_models import (
    SpectrumData,
    SpectrumDataType,
    SpectrumPeak,
)


# NumPy >= 2.0 renamed `trapz` to `trapezoid` (and later removed `trapz`)
_trapezoid = getattr(np, "trapezoid", None) or getattr(np, "trapz")


@dataclass(frozen=True)
class PeakPreprocessingConfig:
    """Configuration for optional peak filtering and profile peak detection."""

    min_intensity: float = 0.0
    min_relative_intensity_pct: float = 0.0
    min_mz: Optional[float] = None
    max_mz: Optional[float] = None
    min_snr: Optional[float] = None
    apply_profile_baseline_subtraction: bool = False
    apply_profile_peak_detection: bool = False
    integrate_profile_peak_area: bool = True
    baseline_window_points: int = 51
    baseline_percentile: float = 10.0


@dataclass(frozen=True)
class ProcessedSpectrumResult:
    """Output of `preprocess_spectrum` containing eligible peaks, relative intensities, and warnings."""

    spectrum: SpectrumData
    eligible_peaks: Tuple[SpectrumPeak, ...]
    relative_intensities_pct: Tuple[float, ...]  # 100 * peak.intensity / max_eligible_intensity
    max_eligible_intensity: float
    was_profile_processed: bool
    warnings: Tuple[str, ...] = field(default_factory=tuple)
    isotope_pool_peaks: Tuple[SpectrumPeak, ...] = field(default_factory=tuple)


def estimate_profile_baseline_and_noise(
    intensities: np.ndarray,
    window_points: int = 51,
    percentile: float = 10.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Estimate local baseline and noise floor (MAD) on a 1D profile intensity array."""
    n = len(intensities)
    if n < 5:
        return np.zeros_like(intensities), np.ones_like(intensities)

    half_w = max(2, window_points // 2)
    baseline = np.zeros(n, dtype=np.float64)
    noise = np.zeros(n, dtype=np.float64)

    # Compute local percentile baseline and MAD noise on strided control points, then interpolate
    step = max(1, half_w // 2)
    ctrl_indices = list(range(0, n, step))
    if ctrl_indices[-1] != n - 1:
        ctrl_indices.append(n - 1)

    ctrl_base: List[float] = []
    ctrl_noise: List[float] = []
    for idx in ctrl_indices:
        lo = max(0, idx - half_w)
        hi = min(n, idx + half_w + 1)
        window = intensities[lo:hi]
        b_val = float(np.percentile(window, percentile))
        mad = float(np.median(np.abs(window - np.median(window)))) * 1.4826
        ctrl_base.append(b_val)
        ctrl_noise.append(max(mad, 1e-6))

    x_all = np.arange(n, dtype=np.float64)
    x_ctrl = np.array(ctrl_indices, dtype=np.float64)
    baseline = np.interp(x_all, x_ctrl, np.array(ctrl_base, dtype=np.float64))
    noise = np.interp(x_all, x_ctrl, np.array(ctrl_noise, dtype=np.float64))
    return baseline, noise


def detect_and_integrate_profile_peaks(
    mzs: np.ndarray,
    intensities: np.ndarray,
    subtract_baseline: bool = True,
    integrate_area: bool = True,
    min_snr: float = 3.0,
    window_points: int = 51,
    percentile: float = 10.0,
    noise_override: Optional[np.ndarray] = None,
    min_fwhm_points: int = 0,
) -> List[SpectrumPeak]:
    """Detect local maxima on a profile spectrum and compute baseline-corrected intensity, S/N, FWHM, and area.

    `noise_override` supplies a per-point noise estimate (e.g. measured before baseline clipping);
    `min_fwhm_points` rejects single-point spikes by requiring at least that many points at or above
    half-maximum (0 disables the check).
    """
    n = len(mzs)
    if n < 5:
        return []

    baseline, noise = estimate_profile_baseline_and_noise(
        intensities,
        window_points=window_points,
        percentile=percentile,
    )
    if noise_override is not None and len(noise_override) == n:
        noise = np.maximum(np.asarray(noise_override, dtype=np.float64), 1e-6)
    corr_int = np.maximum(0.0, intensities - baseline) if subtract_baseline else np.maximum(0.0, intensities)

    detected: List[SpectrumPeak] = []
    # Vectorised candidate search: strict local maxima (or left edge of flat top) above the S/N floor.
    # Only these candidates enter the (Python) refinement loop, which makes profile peak picking fast.
    inner = corr_int[1:-1]
    cand_mask = (inner > 0.0) & (inner > corr_int[:-2]) & (inner >= corr_int[2:])
    with np.errstate(divide="ignore", invalid="ignore"):
        snr_all = np.where(noise[1:-1] > 0, inner / noise[1:-1], np.inf)
    cand_mask &= snr_all >= min_snr
    for i in (np.flatnonzero(cand_mask) + 1).tolist():
        y_curr = corr_int[i]
        if y_curr > 0.0:
            snr_val = float(y_curr / noise[i]) if noise[i] > 0 else None
            if snr_val is not None and snr_val < min_snr:
                continue

            # Walk left and right to find local valleys / half-maximum bounds
            left = i
            while left > 0 and corr_int[left - 1] <= corr_int[left] and corr_int[left - 1] > 0.05 * y_curr:
                left -= 1
            right = i
            while right < n - 1 and corr_int[right + 1] <= corr_int[right] and corr_int[right + 1] > 0.05 * y_curr:
                right += 1

            # Parabolic 3-point centroid refinement around index i
            y_l, y_c, y_r = corr_int[i - 1], corr_int[i], corr_int[i + 1]
            denom = y_l - 2.0 * y_c + y_r
            if abs(denom) > 1e-12:
                delta_idx = 0.5 * (y_l - y_r) / denom
                delta_idx = max(-0.5, min(0.5, float(delta_idx)))
                step_mz = (mzs[i + 1] - mzs[i - 1]) * 0.5
                centroid_mz = float(mzs[i] + delta_idx * step_mz)
            else:
                centroid_mz = float(mzs[i])

            area_val: Optional[float] = None
            if integrate_area and right > left:
                area_val = float(_trapezoid(corr_int[left : right + 1], mzs[left : right + 1]))

            # Estimate FWHM at 50% height
            half_h = 0.5 * y_curr
            l_hm = left
            while l_hm < i and corr_int[l_hm] < half_h:
                l_hm += 1
            r_hm = right
            while r_hm > i and corr_int[r_hm] < half_h:
                r_hm -= 1
            fwhm_val = float(mzs[r_hm] - mzs[l_hm]) if r_hm > l_hm else None
            if min_fwhm_points > 0 and (r_hm - l_hm + 1) < min_fwhm_points:
                continue

            detected.append(
                SpectrumPeak(
                    mz=centroid_mz,
                    intensity=float(y_curr),
                    area=area_val,
                    snr=round(snr_val, 3) if snr_val is not None else None,
                    fwhm=round(fwhm_val, 5) if fwhm_val is not None else None,
                    row_index=i + 1,
                )
            )

    return detected


def preprocess_spectrum(
    spectrum: SpectrumData,
    config: Optional[PeakPreprocessingConfig] = None,
) -> ProcessedSpectrumResult:
    """Filter and optionally peak-pick a `SpectrumData` object according to scientific rules."""
    if config is None:
        config = PeakPreprocessingConfig()

    warnings: List[str] = list(spectrum.warnings)
    was_profile_processed = False
    working_peaks: List[SpectrumPeak] = list(spectrum.peaks)

    if config.apply_profile_baseline_subtraction or config.apply_profile_peak_detection:
        if spectrum.data_type != SpectrumDataType.PROFILE:
            warnings.append(
                "Profile-specific signal processing (baseline subtraction / peak detection / numerical area integration) "
                f"was requested, but input spectrum is classified as '{spectrum.data_type.value}'. "
                "Insufficient profile data for peak-area calculation; proceeding with existing centroid peaks."
            )
        else:
            mzs = np.array([p.mz for p in working_peaks], dtype=np.float64)
            ints = np.array([p.intensity for p in working_peaks], dtype=np.float64)
            working_peaks = detect_and_integrate_profile_peaks(
                mzs=mzs,
                intensities=ints,
                subtract_baseline=config.apply_profile_baseline_subtraction,
                integrate_area=config.integrate_profile_peak_area,
                min_snr=config.min_snr if config.min_snr is not None else 1.5,
                window_points=config.baseline_window_points,
                percentile=config.baseline_percentile,
            )
            was_profile_processed = True
    elif spectrum.data_type == SpectrumDataType.PROFILE and len(working_peaks) >= 100:
        # When a raw profile spectrum is passed directly to MS1 search without explicit peak picking,
        # extract local maxima so baseline points are not searched as thousands of centroids.
        mzs = np.array([p.mz for p in working_peaks], dtype=np.float64)
        ints = np.array([p.intensity for p in working_peaks], dtype=np.float64)
        working_peaks = detect_and_integrate_profile_peaks(
            mzs=mzs,
            intensities=ints,
            subtract_baseline=True,
            integrate_area=config.integrate_profile_peak_area,
            min_snr=config.min_snr if config.min_snr is not None else 3.0,
            window_points=config.baseline_window_points,
            percentile=config.baseline_percentile,
        )
        was_profile_processed = True
        warnings.append(
            f"Automatically detected and integrated {len(working_peaks)} centroid peaks from profile spectrum."
        )

    # When an S/N threshold is requested but peaks carry no S/N (typical for centroid lists),
    # estimate a relative local S/N so the filter is still meaningful (and say so explicitly).
    if config.min_snr is not None and config.min_snr > 0 and any(p.snr is None for p in working_peaks):
        from app.analysis.spectrum_corrections import fill_missing_snr

        working_peaks, n_est = fill_missing_snr(working_peaks)
        warnings.append(
            f"S/N filter requested but {n_est} peaks had no instrument S/N; "
            "a relative local-median S/N estimate was used for them."
        )

    # Apply m/z range, absolute intensity, and S/N filters
    filtered: List[SpectrumPeak] = []
    for p in working_peaks:
        if p.intensity <= 0.0:
            continue
        if config.min_mz is not None and p.mz < config.min_mz:
            continue
        if config.max_mz is not None and p.mz > config.max_mz:
            continue
        if p.intensity < config.min_intensity:
            continue
        if config.min_snr is not None and p.snr is not None and p.snr < config.min_snr:
            continue
        filtered.append(p)

    if not filtered:
        return ProcessedSpectrumResult(
            spectrum=spectrum,
            eligible_peaks=(),
            relative_intensities_pct=(),
            max_eligible_intensity=0.0,
            was_profile_processed=was_profile_processed,
            warnings=tuple(warnings + ["No peaks remained after applying preprocessing filters."]),
        )

    max_int = max(p.intensity for p in filtered)
    # Peaks passing the m/z / S/N filters but BEFORE the relative-intensity cut: kept for isotope
    # envelope scoring so that weak M+1/M+2 satellites are not discarded by the intensity threshold.
    isotope_pool = tuple(filtered)

    # Apply relative intensity % threshold if specified
    if config.min_relative_intensity_pct > 0.0 and max_int > 0.0:
        filtered = [
            p for p in filtered if (p.intensity / max_int) * 100.0 >= config.min_relative_intensity_pct
        ]
        max_int = max((p.intensity for p in filtered), default=0.0)

    rel_pcts = tuple((p.intensity / max_int) * 100.0 if max_int > 0 else 0.0 for p in filtered)

    return ProcessedSpectrumResult(
        spectrum=spectrum,
        eligible_peaks=tuple(filtered),
        relative_intensities_pct=rel_pcts,
        max_eligible_intensity=max_int,
        was_profile_processed=was_profile_processed,
        warnings=tuple(warnings),
        isotope_pool_peaks=isotope_pool,
    )


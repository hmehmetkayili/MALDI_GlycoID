"""Optional MALDI-TOF spectrum corrections, S/N estimation, and import-time filtering.

Every correction in this module is OPT-IN (nothing is applied silently) and every
applied step is logged in `CorrectionResult.steps` so that the processing history
is fully transparent and reproducible.

Available corrections (see `list_available_corrections()` for the UI catalogue):

Profile spectra only (they require continuous signal and are skipped with a warning
on centroid peak lists):
    - Savitzky-Golay smoothing
    - SNIP baseline subtraction (chemical-noise / matrix-background removal)
    - Local-maxima peak picking with MAD-based S/N, FWHM and area

Centroid or profile:
    - m/z range trimming
    - S/N estimation for peaks without an instrument S/N value (flagged as *estimated*)
    - S/N threshold filtering
    - Relative-intensity threshold filtering
    - MALDI matrix / contaminant ion exclusion (DHB, CHCA, SA; positive & negative mode)
    - Merge of near-duplicate centroids (centroiding artefacts)
    - Intensity normalisation (base peak / TIC / median)
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from app.parsers.base_models import SpectrumData, SpectrumDataType, SpectrumPeak

NEUTRON_SPACING_DA = 1.0033548378  # mean 13C-12C / generic isotope spacing

# Singly-charged ion masses (Da). Used only for contaminant annotation.
_PROTON = 1.00727646688
_NA_ION = 22.98921677
_K_ION = 38.96315790

# Neutral monoisotopic masses of common MALDI matrices.
_MATRICES: Dict[str, float] = {
    "DHB": 154.026609,
    "CHCA": 189.042593,
    "SA": 224.068473,
    "DHAP": 152.047344,
}


def matrix_ion_list(polarity: str = "positive") -> List[Tuple[str, float]]:
    """Return (label, m/z) of common singly-charged matrix ions for the given polarity."""
    ions: List[Tuple[str, float]] = []
    water = 18.010565
    for name, m in _MATRICES.items():
        if polarity.lower().startswith("neg"):
            ions.append((f"{name} [M-H]-", m - _PROTON))
            ions.append((f"{name} [2M-H]-", 2 * m - _PROTON))
        else:
            ions.append((f"{name} [M+H]+", m + _PROTON))
            ions.append((f"{name} [M+Na]+", m + _NA_ION))
            ions.append((f"{name} [M+K]+", m + _K_ION))
            ions.append((f"{name} [M-H2O+H]+", m - water + _PROTON))
            ions.append((f"{name} [2M+H]+", 2 * m + _PROTON))
            ions.append((f"{name} [2M-H2O+H]+", 2 * m - water + _PROTON))
            ions.append((f"{name} [2M+Na]+", 2 * m + _NA_ION))
    return ions


@dataclass(frozen=True)
class CorrectionConfig:
    """User-selected optional corrections. All defaults are 'off' (no-op)."""

    min_mz: Optional[float] = None
    max_mz: Optional[float] = None
    # profile-only
    smoothing: bool = False
    smoothing_window_points: int = 9
    smoothing_polyorder: int = 3
    baseline_correction: bool = False
    baseline_iterations: int = 60
    pick_peaks: bool = False
    pick_min_snr: float = 3.0
    # centroid or profile
    estimate_missing_snr: bool = True
    min_snr: Optional[float] = None
    min_relative_intensity_pct: float = 0.0
    remove_matrix_peaks: bool = False
    matrix_polarity: str = "positive"
    matrix_tolerance_ppm: float = 100.0
    merge_close_peaks: bool = False
    merge_tolerance_ppm: float = 20.0
    normalize: Optional[str] = None  # None | "base_peak" | "tic" | "median"

    @staticmethod
    def from_dict(d: Dict[str, Any]) -> "CorrectionConfig":
        def _f(key: str) -> Optional[float]:
            v = d.get(key)
            if v is None or v == "":
                return None
            return float(v)

        norm = d.get("normalize")
        norm = str(norm).strip().lower() if norm else None
        if norm in ("", "none"):
            norm = None
        if norm not in (None, "base_peak", "tic", "median"):
            raise ValueError(f"Unsupported normalization '{norm}'. Use base_peak, tic, or median.")
        min_snr = _f("min_snr")
        return CorrectionConfig(
            min_mz=_f("min_mz"),
            max_mz=_f("max_mz"),
            smoothing=bool(d.get("smoothing", False)),
            smoothing_window_points=int(d.get("smoothing_window_points", 9)),
            smoothing_polyorder=int(d.get("smoothing_polyorder", 3)),
            baseline_correction=bool(d.get("baseline_correction", False)),
            baseline_iterations=int(d.get("baseline_iterations", 60)),
            pick_peaks=bool(d.get("pick_peaks", False)),
            pick_min_snr=float(d.get("pick_min_snr", 3.0)),
            estimate_missing_snr=bool(d.get("estimate_missing_snr", True)),
            min_snr=min_snr if (min_snr is not None and min_snr > 0) else None,
            min_relative_intensity_pct=float(d.get("min_relative_intensity_pct", 0.0) or 0.0),
            remove_matrix_peaks=bool(d.get("remove_matrix_peaks", False)),
            matrix_polarity=str(d.get("matrix_polarity", "positive")),
            matrix_tolerance_ppm=float(d.get("matrix_tolerance_ppm", 100.0)),
            merge_close_peaks=bool(d.get("merge_close_peaks", False)),
            merge_tolerance_ppm=float(d.get("merge_tolerance_ppm", 20.0)),
            normalize=norm,
        )

    def is_noop(self) -> bool:
        return self == CorrectionConfig()


@dataclass(frozen=True)
class CorrectionResult:
    """Corrected spectrum with a transparent step-by-step processing log."""

    spectrum: SpectrumData
    steps: Tuple[Dict[str, Any], ...]
    warnings: Tuple[str, ...]
    peaks_before: int
    peaks_after: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            "peaks_before": self.peaks_before,
            "peaks_after": self.peaks_after,
            "steps": [dict(s) for s in self.steps],
            "warnings": list(self.warnings),
        }


def list_available_corrections() -> List[Dict[str, Any]]:
    """Catalogue of optional corrections for UI presentation."""
    return [
        {"key": "smoothing", "label": "Savitzky-Golay smoothing", "applies_to": "profile",
         "purpose": "Reduces high-frequency detector noise before peak picking; preserves peak width."},
        {"key": "baseline_correction", "label": "SNIP baseline subtraction", "applies_to": "profile",
         "purpose": "Removes chemical-noise / matrix background that rises toward low m/z."},
        {"key": "pick_peaks", "label": "Peak picking (centroiding)", "applies_to": "profile",
         "purpose": "Converts profile data to centroids with S/N, FWHM and area; makes search much faster."},
        {"key": "min_snr", "label": "S/N threshold", "applies_to": "both",
         "purpose": "Drops noise peaks. S/N is estimated locally when the file has no S/N column."},
        {"key": "min_relative_intensity_pct", "label": "Relative-intensity threshold", "applies_to": "both",
         "purpose": "Drops peaks below a % of the base peak."},
        {"key": "min_mz/max_mz", "label": "m/z range trimming", "applies_to": "both",
         "purpose": "Removes matrix/peptide region below the glycan range."},
        {"key": "remove_matrix_peaks", "label": "Matrix / contaminant exclusion", "applies_to": "both",
         "purpose": "Removes DHB, CHCA, SA, DHAP matrix ions (and clusters) from the peak list."},
        {"key": "merge_close_peaks", "label": "Merge near-duplicate centroids", "applies_to": "both",
         "purpose": "Merges split centroids of one peak (intensity-weighted m/z)."},
        {"key": "normalize", "label": "Intensity normalisation", "applies_to": "both",
         "purpose": "Scales intensities to base peak=100, TIC=100 or median=1 for sample comparison."},
        {"key": "calibration", "label": "Mass calibration (before search)", "applies_to": "both",
         "purpose": "Corrects systematic m/z drift (linear/quadratic) using internal glycan or peptide calibrants."},
        {"key": "exclude_isotope_satellites", "label": "Exclude isotope satellites from annotation", "applies_to": "search",
         "purpose": "Peaks that are M+1/M+2 of a lower peak are not annotated separately (still used for idotp)."},
    ]


# --------------------------------------------------------------------------- #
# Signal-processing primitives
# --------------------------------------------------------------------------- #

def savitzky_golay(y: np.ndarray, window: int = 9, polyorder: int = 3) -> np.ndarray:
    """Savitzky-Golay smoothing (scipy if available; clipped to non-negative)."""
    n = len(y)
    window = int(window)
    if window % 2 == 0:
        window += 1
    window = min(window, n if n % 2 == 1 else n - 1)
    polyorder = min(int(polyorder), window - 1)
    if window < 3 or polyorder < 1:
        return y.copy()
    from scipy.signal import savgol_filter

    return np.maximum(0.0, savgol_filter(y, window_length=window, polyorder=polyorder, mode="interp"))


def snip_baseline(y: np.ndarray, iterations: int = 60) -> np.ndarray:
    """SNIP (Statistics-sensitive Non-linear Iterative Peak-clipping) baseline estimate."""
    n = len(y)
    if n < 5:
        return np.zeros_like(y)
    iterations = max(1, min(int(iterations), (n - 1) // 2))
    # LLS transform compresses dynamic range
    v = np.log(np.log(np.sqrt(np.maximum(y, 0.0) + 1.0) + 1.0) + 1.0)
    for p in range(iterations, 0, -1):
        avg = 0.5 * (v[: n - 2 * p] + v[2 * p:])
        v[p: n - p] = np.minimum(v[p: n - p], avg)
    base = (np.exp(np.exp(v) - 1.0) - 1.0) ** 2 - 1.0
    return np.clip(base, 0.0, y)


def estimate_centroid_snr(
    mzs: np.ndarray,
    intensities: np.ndarray,
    window_da: float = 300.0,
) -> np.ndarray:
    """Estimate a relative S/N for centroid lists lacking an instrument-reported S/N.

    The noise level is the local median peak intensity (window +-`window_da`), which for a
    glycan peak list is dominated by low-level chemical-noise centroids. When fewer than
    5 neighbours are available the global median is used. The value is a *relative
    estimate* and is flagged as such in the processing log.
    """
    n = len(mzs)
    if n == 0:
        return np.array([], dtype=float)
    order = np.argsort(mzs)
    m = mzs[order]
    y = intensities[order]
    global_med = float(np.median(y)) if n else 1.0
    out = np.zeros(n, dtype=float)
    lo = np.searchsorted(m, m - window_da, side="left")
    hi = np.searchsorted(m, m + window_da, side="right")
    for i in range(n):
        seg = y[lo[i]: hi[i]]
        noise = float(np.median(seg)) if len(seg) >= 5 else global_med
        out[i] = y[i] / noise if noise > 0 else 0.0
    res = np.zeros(n, dtype=float)
    res[order] = out
    return res


def fill_missing_snr(peaks: Sequence[SpectrumPeak]) -> Tuple[List[SpectrumPeak], int]:
    """Return peaks with estimated S/N where missing, plus the number of estimated values."""
    missing = [i for i, p in enumerate(peaks) if p.snr is None]
    if not missing:
        return list(peaks), 0
    mzs = np.array([p.mz for p in peaks], dtype=float)
    ys = np.array([p.intensity for p in peaks], dtype=float)
    est = estimate_centroid_snr(mzs, ys)
    out = list(peaks)
    for i in missing:
        out[i] = replace(peaks[i], snr=round(float(est[i]), 3))
    return out, len(missing)


def find_isotope_satellite_indices(
    mzs: np.ndarray,
    intensities: np.ndarray,
    tolerance_ppm: float = 25.0,
    charges: Sequence[int] = (1,),
    max_offset: int = 3,
    max_ratio: float = 1.6,
) -> set:
    """Indices of peaks that look like M+k isotopologues of a *lower* m/z peak.

    A peak j is a satellite of i when mz_j ~ mz_i + k*1.00335/z (k=1..max_offset) and
    I_j <= max_ratio * I_i (a satellite far more intense than its parent is not a satellite).
    Peaks are never removed - the caller decides what to do with the returned indices.
    """
    n = len(mzs)
    if n < 2:
        return set()
    order = np.argsort(mzs)
    m = mzs[order]
    y = intensities[order]
    sat: set = set()
    for i in range(n):
        for z in charges:
            for k in range(1, max_offset + 1):
                target = m[i] + k * NEUTRON_SPACING_DA / max(1, abs(int(z)))
                tol = target * tolerance_ppm * 1e-6
                lo = int(np.searchsorted(m, target - tol, side="left"))
                hi = int(np.searchsorted(m, target + tol, side="right"))
                for j in range(lo, hi):
                    if j != i and y[j] <= max_ratio * y[i]:
                        sat.add(int(order[j]))
    return sat


# --------------------------------------------------------------------------- #
# Spectrum rebuild helper
# --------------------------------------------------------------------------- #

def _rebuild(spectrum: SpectrumData, peaks: Sequence[SpectrumPeak], **overrides: Any) -> SpectrumData:
    peaks_t = tuple(sorted(peaks, key=lambda p: p.mz))
    return replace(
        spectrum,
        peaks=peaks_t,
        valid_points_count=len(peaks_t),
        min_mz=min((p.mz for p in peaks_t), default=0.0),
        max_mz=max((p.mz for p in peaks_t), default=0.0),
        max_intensity=max((p.intensity for p in peaks_t), default=0.0),
        **overrides,
    )


# --------------------------------------------------------------------------- #
# Main entry point
# --------------------------------------------------------------------------- #

def apply_corrections(spectrum: SpectrumData, config: Optional[CorrectionConfig] = None) -> CorrectionResult:
    """Apply the user-selected optional corrections and return the log."""
    cfg = config or CorrectionConfig()
    warnings: List[str] = []
    steps: List[Dict[str, Any]] = []
    spec = spectrum
    n_before = len(spectrum.peaks)
    is_profile = spec.data_type == SpectrumDataType.PROFILE

    def log(name: str, **info: Any) -> None:
        steps.append({"step": name, "peaks_after": len(spec.peaks), **info})

    # 1. m/z range trimming
    if cfg.min_mz is not None or cfg.max_mz is not None:
        keep = [
            p for p in spec.peaks
            if (cfg.min_mz is None or p.mz >= cfg.min_mz) and (cfg.max_mz is None or p.mz <= cfg.max_mz)
        ]
        removed = len(spec.peaks) - len(keep)
        spec = _rebuild(spec, keep)
        log("mz_range", min_mz=cfg.min_mz, max_mz=cfg.max_mz, removed=removed)

    # 2-4. Profile-only steps
    profile_steps = cfg.smoothing or cfg.baseline_correction or cfg.pick_peaks
    if profile_steps and not is_profile:
        warnings.append(
            "Smoothing / baseline correction / peak picking require a profile spectrum; "
            f"input is '{spec.data_type.value}', these steps were skipped."
        )
    elif profile_steps and len(spec.peaks) >= 5:
        from app.analysis.peak_processing import (
            detect_and_integrate_profile_peaks,
            estimate_profile_baseline_and_noise,
        )

        # Work on arrays and rebuild the (large) peak tuple only once at the end of the profile block.
        mzs = np.array([p.mz for p in spec.peaks], dtype=float)
        ys = np.array([p.intensity for p in spec.peaks], dtype=float)
        noise_est: Optional[np.ndarray] = None
        profile_logs: List[Tuple[str, Dict[str, Any]]] = []
        if cfg.smoothing:
            ys = savitzky_golay(ys, cfg.smoothing_window_points, cfg.smoothing_polyorder)
            profile_logs.append(("smoothing", {"method": "savitzky_golay",
                                               "window": cfg.smoothing_window_points,
                                               "polyorder": cfg.smoothing_polyorder}))
        if cfg.pick_peaks:
            # Noise is measured on the smoothed signal BEFORE baseline clipping; after subtracting a baseline
            # the clipped residual no longer has a meaningful MAD (it collapses and inflates every S/N).
            _, noise_est = estimate_profile_baseline_and_noise(ys)
        if cfg.baseline_correction:
            base = snip_baseline(ys, cfg.baseline_iterations)
            ys = np.maximum(0.0, ys - base)
            profile_logs.append(("baseline_correction", {"method": "snip", "iterations": cfg.baseline_iterations,
                                                         "median_baseline": round(float(np.median(base)), 4)}))
        if cfg.smoothing or cfg.baseline_correction:
            spec = _rebuild(spec, [
                SpectrumPeak(mz=p.mz, intensity=float(v), area=p.area, snr=p.snr, fwhm=p.fwhm,
                             charge=p.charge, row_index=p.row_index)
                for p, v in zip(spec.peaks, ys)
            ])
            for name, info in profile_logs:
                log(name, **info)
        if cfg.pick_peaks:
            picked = detect_and_integrate_profile_peaks(
                mzs=mzs,
                intensities=ys,
                # baseline already removed when requested; otherwise let the detector estimate it
                subtract_baseline=not cfg.baseline_correction,
                integrate_area=True,
                min_snr=cfg.pick_min_snr,
                noise_override=noise_est,
                min_fwhm_points=3,
            )
            spec = _rebuild(spec, picked, data_type=SpectrumDataType.CENTROID)
            is_profile = False
            log("peak_picking", min_snr=cfg.pick_min_snr, picked=len(picked))
            if not picked:
                warnings.append("Peak picking returned no peaks; relax the S/N threshold.")

    # 5. S/N estimation + filtering
    if (cfg.estimate_missing_snr and cfg.min_snr is not None) and spec.peaks:
        filled, n_est = fill_missing_snr(spec.peaks)
        if n_est:
            spec = _rebuild(spec, filled)
            warnings.append(
                f"S/N was not present for {n_est} peaks; a relative local-median S/N estimate was used."
            )
            log("snr_estimation", estimated=n_est, method="local_median_noise")
    if cfg.min_snr is not None:
        keep = [p for p in spec.peaks if p.snr is None or p.snr >= cfg.min_snr]
        removed = len(spec.peaks) - len(keep)
        spec = _rebuild(spec, keep)
        log("snr_filter", min_snr=cfg.min_snr, removed=removed)

    # 6. Relative intensity filter
    if cfg.min_relative_intensity_pct > 0 and spec.peaks:
        mx = max(p.intensity for p in spec.peaks)
        keep = [p for p in spec.peaks if p.intensity >= mx * cfg.min_relative_intensity_pct / 100.0]
        removed = len(spec.peaks) - len(keep)
        spec = _rebuild(spec, keep)
        log("relative_intensity_filter", threshold_pct=cfg.min_relative_intensity_pct, removed=removed)

    # 7. Matrix / contaminant exclusion
    if cfg.remove_matrix_peaks and spec.peaks:
        ions = matrix_ion_list(cfg.matrix_polarity)
        mz_arr = np.array([m for _, m in ions])
        removed_labels: List[str] = []
        keep = []
        for p in spec.peaks:
            hit = np.where(np.abs(mz_arr - p.mz) / mz_arr * 1e6 <= cfg.matrix_tolerance_ppm)[0]
            if len(hit):
                removed_labels.append(f"{ions[int(hit[0])][0]} @ {p.mz:.4f}")
            else:
                keep.append(p)
        spec = _rebuild(spec, keep)
        log("matrix_exclusion", polarity=cfg.matrix_polarity, removed=len(removed_labels),
            removed_ions=removed_labels[:20])

    # 8. Merge near-duplicate centroids (only meaningful for centroid lists)
    if cfg.merge_close_peaks and spec.peaks:
        if is_profile:
            warnings.append("Merging close centroids was skipped: spectrum is profile data.")
        else:
            merged: List[SpectrumPeak] = []
            group: List[SpectrumPeak] = []

            def flush() -> None:
                if not group:
                    return
                tot = sum(g.intensity for g in group)
                mz = sum(g.mz * g.intensity for g in group) / tot if tot > 0 else group[0].mz
                top = max(group, key=lambda g: g.intensity)
                area = sum(g.area for g in group if g.area is not None) if any(
                    g.area is not None for g in group) else None
                merged.append(replace(top, mz=float(mz), intensity=float(tot), area=area))

            for p in spec.peaks:
                if group and (p.mz - group[-1].mz) / group[-1].mz * 1e6 > cfg.merge_tolerance_ppm:
                    flush()
                    group = []
                group.append(p)
            flush()
            removed = len(spec.peaks) - len(merged)
            spec = _rebuild(spec, merged)
            log("merge_close_peaks", tolerance_ppm=cfg.merge_tolerance_ppm, merged=removed)

    # 9. Normalisation
    if cfg.normalize and spec.peaks:
        ints = np.array([p.intensity for p in spec.peaks], dtype=float)
        if cfg.normalize == "base_peak":
            denom, scale = float(ints.max()), 100.0
        elif cfg.normalize == "tic":
            denom, scale = float(ints.sum()), 100.0
        else:
            denom, scale = float(np.median(ints)), 1.0
        if denom > 0:
            f = scale / denom
            spec = _rebuild(
                spec,
                [replace(p, intensity=p.intensity * f, area=p.area * f if p.area is not None else None)
                 for p in spec.peaks],
            )
            log("normalization", mode=cfg.normalize, scale_factor=f)
            warnings.append(
                f"Intensities normalised ({cfg.normalize}); absolute intensity thresholds now refer to normalised units."
            )

    if steps:
        meta = dict(spec.metadata)
        meta["corrections_applied"] = [s["step"] for s in steps]
        spec = replace(spec, metadata=meta, warnings=spec.warnings + tuple(warnings))
    return CorrectionResult(
        spectrum=spec,
        steps=tuple(steps),
        warnings=tuple(warnings),
        peaks_before=n_before,
        peaks_after=len(spec.peaks),
    )

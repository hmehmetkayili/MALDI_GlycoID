"""Internal and external mass calibration engine for MALDI-GlycoID (Stage 8).

Implements transparent linear and quadratic least-squares mass calibration using
known reference calibrant ions (internal N-glycan calibrants, peptide calibration
standards, or user-defined reference lists).

Scientific Safeguards:
1. Never silently calibrates a spectrum: calibration requires explicit execution,
   logs all matched/unmatched calibrants, and attaches provenance metadata to
   the calibrated `SpectrumData`.
2. Minimum calibrant count safeguards:
   - Linear model (`mz_cal = a * mz_obs + b`) requires >= 2 matched calibrants
     (warns if < 3 since 2 points yield zero degrees of freedom).
   - Quadratic model (`mz_cal = a * mz_obs^2 + b * mz_obs + c`) requires >= 3
     matched calibrants (warns if < 4 since 3 points yield zero degrees of freedom).
3. Extrapolation & span warnings: warns when matched calibrants cover less than
   30% of the observed spectrum m/z range or when post-calibration residuals
   exceed `max_residual_ppm_threshold`.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from enum import Enum
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from app.core.ion_calculator import calculate_mz_from_neutral_mass
from app.core.mass_calculator import calculate_peptide_neutral_mass, calculate_released_glycan_mass
from app.core.mass_error import (
    MassToleranceUnit,
    calculate_abs_error_ppm,
    calculate_signed_error_da,
    calculate_signed_error_ppm,
    is_within_mass_tolerance,
)
from app.parsers.base_models import SpectrumData, SpectrumDataType, SpectrumPeak


class CalibrationModelType(str, Enum):
    """Supported polynomial mass calibration models."""

    LINEAR = "linear"
    QUADRATIC = "quadratic"


@dataclass(frozen=True)
class CalibrantReference:
    """Reference ion used for internal or external mass calibration."""

    name: str
    theoretical_mz: float
    composition_or_formula: Optional[str] = None
    adduct: Optional[str] = None


@dataclass(frozen=True)
class MatchedCalibrantPoint:
    """Diagnostic record for a single matched reference calibrant peak."""

    name: str
    theoretical_mz: float
    observed_mz: float
    calibrated_mz: float
    intensity: float
    error_da_before: float
    error_ppm_before: float
    residual_da_after: float
    residual_ppm_after: float


@dataclass(frozen=True)
class CalibrationSettings:
    """Configuration for reference peak matching and polynomial calibration."""

    model_type: CalibrationModelType = CalibrationModelType.LINEAR
    match_tolerance: float = 150.0
    match_tolerance_unit: MassToleranceUnit = MassToleranceUnit.PPM
    min_rel_intensity_pct: float = 0.5
    max_residual_ppm_threshold: float = 25.0
    require_calibration_success: bool = False
    reject_outliers: bool = True


@dataclass(frozen=True)
class MassCalibrationResult:
    """Complete result and residual diagnostics of a mass calibration run."""

    sample_id: str
    applied: bool
    model_type: CalibrationModelType
    coefficients: Tuple[float, ...]
    matched_calibrants: Tuple[MatchedCalibrantPoint, ...]
    unmatched_calibrants: Tuple[CalibrantReference, ...]
    rms_ppm_before: Optional[float]
    rms_ppm_after: Optional[float]
    median_abs_ppm_before: Optional[float]
    median_abs_ppm_after: Optional[float]
    mean_signed_ppm_before: Optional[float]
    mean_signed_ppm_after: Optional[float]
    max_abs_ppm_before: Optional[float]
    max_abs_ppm_after: Optional[float]
    calibrant_mz_span_da: float
    calibrated_spectrum: SpectrumData
    warnings: Tuple[str, ...] = field(default_factory=tuple)


def get_default_nglycan_calibrants(
    adduct: str = "[M+Na]+",
    labeling: str = "NONE",
) -> List[CalibrantReference]:
    """Return exact theoretical reference calibrants for canonical human N-glycans.

    All theoretical m/z values are computed directly from the Stage 1 mass engine.
    """
    reference_panel = [
        ("Man5 (Hex5HexNAc2)", "Hex5HexNAc2"),
        ("FA2 / G0F (Hex3HexNAc4dHex1)", "Hex3HexNAc4dHex1"),
        ("FA2G1 / G1F (Hex4HexNAc4dHex1)", "Hex4HexNAc4dHex1"),
        ("A2G2 / G2 (Hex5HexNAc4)", "Hex5HexNAc4"),
        ("FA2G2 / G2F (Hex5HexNAc4dHex1)", "Hex5HexNAc4dHex1"),
        ("A2G2S1 (Hex5HexNAc4Neu5Ac1)", "Hex5HexNAc4Neu5Ac1"),
        ("FA2G2S1 (Hex5HexNAc4dHex1Neu5Ac1)", "Hex5HexNAc4dHex1Neu5Ac1"),
        ("A2G2S2 (Hex5HexNAc4Neu5Ac2)", "Hex5HexNAc4Neu5Ac2"),
    ]
    calibrants: List[CalibrantReference] = []
    for label, comp_str in reference_panel:
        neutral_res = calculate_released_glycan_mass(comp_str, labeling=labeling)
        ion_res = calculate_mz_from_neutral_mass(neutral_res.neutral_monoisotopic_mass, adduct_code=adduct)
        calibrants.append(
            CalibrantReference(
                name=f"{label} {adduct}",
                theoretical_mz=ion_res.theoretical_mz,
                composition_or_formula=comp_str,
                adduct=adduct,
            )
        )
    return calibrants


def get_default_peptide_calibrants() -> List[CalibrantReference]:
    """Return standard MALDI-TOF peptide calibration reference ions ([M+H]+)."""
    peptides = [
        ("Bradykinin (1-7)", "RPPGFSP"),
        ("Angiotensin II", "DRVYIHPF"),
        ("Angiotensin I", "DRVYIHPFHL"),
        ("Substance P (free acid)", "RPKPQQFFGLM"),
        ("ACTH clip (1-17)", "SYSMEHFRWGKPVGKKR"),
        ("ACTH clip (18-39)", "RPVKVYPNGAEDESAEAFPLEF"),
    ]
    calibrants: List[CalibrantReference] = []
    for name, seq in peptides:
        pep_res = calculate_peptide_neutral_mass(seq)
        ion_res = calculate_mz_from_neutral_mass(pep_res.neutral_monoisotopic_mass, adduct_code="[M+H]+")
        calibrants.append(
            CalibrantReference(
                name=f"{name} [M+H]+",
                theoretical_mz=ion_res.theoretical_mz,
                composition_or_formula=seq,
                adduct="[M+H]+",
            )
        )
    return calibrants


CALIBRANT_PANELS: Dict[str, str] = {
    "nglycan_na": "Canonical human N-glycans [M+Na]+ (positive)",
    "nglycan_h": "Canonical human N-glycans [M+H]+ (positive)",
    "nglycan_neg": "Canonical human N-glycans [M-H]- (negative)",
    "peptide_h": "MALDI peptide calibration standard [M+H]+ (positive)",
}


def get_calibrants_for_panel(panel: str, labeling: str = "NONE") -> List[CalibrantReference]:
    """Return calibrants for a named reference panel (see `CALIBRANT_PANELS`)."""
    key = panel.strip().lower()
    if key == "peptide_h":
        return get_default_peptide_calibrants()
    if key == "nglycan_h":
        return get_default_nglycan_calibrants(adduct="[M+H]+", labeling=labeling)
    if key == "nglycan_neg":
        return get_default_nglycan_calibrants(adduct="[M-H]-", labeling=labeling)
    if key == "nglycan_na":
        return get_default_nglycan_calibrants(adduct="[M+Na]+", labeling=labeling)
    raise ValueError(f"Unknown calibrant panel '{panel}'. Choose from {sorted(CALIBRANT_PANELS)} or 'auto'.")


def select_best_calibrant_panel(
    spectrum: SpectrumData,
    settings: Optional["CalibrationSettings"] = None,
    labeling: str = "NONE",
    candidate_panels: Optional[Sequence[str]] = None,
) -> Tuple[str, Dict[str, int]]:
    """Pick the panel with the most matched calibrants (ties: panel order). Returns (panel, counts)."""
    panels = list(candidate_panels or ("nglycan_na", "nglycan_neg", "nglycan_h", "peptide_h"))
    counts: Dict[str, int] = {}
    for name in panels:
        matched, _ = match_calibrants_to_spectrum(spectrum, get_calibrants_for_panel(name, labeling), settings)
        counts[name] = len(matched)
    best = max(panels, key=lambda n: (counts[n], -panels.index(n)))
    return best, counts


def _refine_profile_apex(peaks: Sequence[SpectrumPeak], idx: int) -> SpectrumPeak:
    """Refine the m/z of a profile maximum with an intensity-weighted centroid over its upper half.

    On profile data the highest raw point is only accurate to half the point spacing; the weighted
    centroid over the points above 50 % of the apex height is far more accurate (sub-ppm to a few ppm).
    """
    apex = peaks[idx]
    thr = 0.5 * apex.intensity
    lo = idx
    while lo > 0 and peaks[lo - 1].intensity >= thr and peaks[lo - 1].intensity <= peaks[lo].intensity * 1.0000001:
        lo -= 1
    hi = idx
    while hi < len(peaks) - 1 and peaks[hi + 1].intensity >= thr and peaks[hi + 1].intensity <= peaks[hi].intensity * 1.0000001:
        hi += 1
    seg = peaks[lo: hi + 1]
    tot = sum(p.intensity for p in seg)
    if len(seg) < 3 or tot <= 0:
        return apex
    mz = sum(p.mz * p.intensity for p in seg) / tot
    return replace(apex, mz=float(mz))


def match_calibrants_to_spectrum(
    spectrum: SpectrumData,
    calibrants: Sequence[CalibrantReference],
    settings: Optional[CalibrationSettings] = None,
) -> Tuple[List[Tuple[CalibrantReference, SpectrumPeak]], List[CalibrantReference]]:
    """Match reference calibrants to the highest-intensity eligible peak within tolerance.

    Uses binary search on the m/z-sorted peak array (fast on large profile spectra) and refines the
    apex m/z of profile data with an intensity-weighted centroid.
    """
    if settings is None:
        settings = CalibrationSettings()

    if not spectrum.peaks or not calibrants:
        return [], list(calibrants)

    base_peak_intensity = max(p.intensity for p in spectrum.peaks)
    min_abs_intensity = (
        base_peak_intensity * (settings.min_rel_intensity_pct / 100.0)
        if base_peak_intensity > 0
        else 0.0
    )

    sorted_peaks = sorted(spectrum.peaks, key=lambda p: p.mz)
    is_profile = spectrum.data_type == SpectrumDataType.PROFILE
    mz_arr = np.array([p.mz for p in sorted_peaks], dtype=np.float64)

    matched_pairs: List[Tuple[CalibrantReference, SpectrumPeak]] = []
    unmatched: List[CalibrantReference] = []
    used_peak_indices: set[int] = set()

    for cal in calibrants:
        if settings.match_tolerance_unit == MassToleranceUnit.PPM:
            half_w = cal.theoretical_mz * settings.match_tolerance * 1e-6 * 1.5 + 1e-9
        else:
            half_w = settings.match_tolerance * 1.5 + 1e-9
        lo = int(np.searchsorted(mz_arr, cal.theoretical_mz - half_w, side="left"))
        hi = int(np.searchsorted(mz_arr, cal.theoretical_mz + half_w, side="right"))

        candidates: List[Tuple[int, SpectrumPeak, float]] = []
        for idx in range(lo, hi):
            pk = sorted_peaks[idx]
            if idx in used_peak_indices or pk.intensity < min_abs_intensity:
                continue
            if is_within_mass_tolerance(
                observed=pk.mz,
                theoretical=cal.theoretical_mz,
                tolerance=settings.match_tolerance,
                unit=settings.match_tolerance_unit,
            ):
                candidates.append((idx, pk, calculate_abs_error_ppm(pk.mz, cal.theoretical_mz)))

        if not candidates:
            unmatched.append(cal)
            continue

        # Select candidate with highest intensity (tie-break by smallest ppm error)
        candidates.sort(key=lambda item: (-item[1].intensity, item[2]))
        best_idx, best_peak, _ = candidates[0]
        used_peak_indices.add(best_idx)
        if is_profile:
            best_peak = _refine_profile_apex(sorted_peaks, best_idx)
        matched_pairs.append((cal, best_peak))

    return matched_pairs, unmatched


def apply_calibration_coefficients(
    mz_values: Sequence[float],
    model_type: CalibrationModelType,
    coefficients: Sequence[float],
) -> List[float]:
    """Evaluate polynomial calibration model on a sequence of observed m/z values."""
    arr = np.asarray(mz_values, dtype=np.float64)
    if model_type == CalibrationModelType.LINEAR:
        if len(coefficients) != 2:
            raise ValueError("Linear calibration requires exactly 2 coefficients: (slope_a, intercept_b).")
        a, b = coefficients
        return (a * arr + b).tolist()
    elif model_type == CalibrationModelType.QUADRATIC:
        if len(coefficients) != 3:
            raise ValueError(
                "Quadratic calibration requires exactly 3 coefficients: (quad_a, slope_b, intercept_c)."
            )
        a, b, c = coefficients
        return (a * (arr ** 2) + b * arr + c).tolist()
    else:
        raise ValueError(f"Unsupported calibration model type: {model_type}")


def calibrate_spectrum(
    spectrum: SpectrumData,
    calibrants: Sequence[CalibrantReference],
    settings: Optional[CalibrationSettings] = None,
) -> MassCalibrationResult:
    """Perform internal or external mass calibration on a `SpectrumData` object.

    Never silently modifies a spectrum: if insufficient calibrants match, returns
    `applied=False` with the original uncalibrated `SpectrumData` and explicit warnings
    (or raises `ValueError` when `settings.require_calibration_success=True`).
    """
    if settings is None:
        settings = CalibrationSettings()

    if settings.match_tolerance <= 0:
        raise ValueError(f"Calibration match tolerance must be positive (got {settings.match_tolerance}).")

    warnings: List[str] = []
    matched_pairs, unmatched = match_calibrants_to_spectrum(spectrum, calibrants, settings)

    min_required = 2 if settings.model_type == CalibrationModelType.LINEAR else 3
    recommended_min = 3 if settings.model_type == CalibrationModelType.LINEAR else 4

    # Robust fit: iteratively drop the calibrant with the largest post-fit residual when it exceeds
    # the residual threshold (typically a mis-assigned peak), while keeping the minimum count.
    if settings.reject_outliers and len(matched_pairs) > min_required:
        degree_ = 1 if settings.model_type == CalibrationModelType.LINEAR else 2
        while len(matched_pairs) > min_required:
            o = np.array([pk.mz for _, pk in matched_pairs], dtype=np.float64)
            t = np.array([c.theoretical_mz for c, _ in matched_pairs], dtype=np.float64)
            fit = np.polyfit(o, t, deg=degree_)
            res_ppm = (np.polyval(fit, o) - t) / t * 1e6
            worst = int(np.argmax(np.abs(res_ppm)))
            if abs(res_ppm[worst]) <= settings.max_residual_ppm_threshold:
                break
            cal_rej, _ = matched_pairs.pop(worst)
            unmatched.append(cal_rej)
            warnings.append(
                f"Calibrant '{cal_rej.name}' rejected as outlier (residual {res_ppm[worst]:+.1f} ppm "
                f"> {settings.max_residual_ppm_threshold:.1f} ppm)."
            )

    if len(matched_pairs) < min_required:
        msg = (
            f"Insufficient matched calibrants for {settings.model_type.value} calibration: "
            f"matched {len(matched_pairs)} of {len(calibrants)} reference ions "
            f"(minimum {min_required} required); spectrum remains uncalibrated."
        )
        if settings.require_calibration_success:
            raise ValueError(msg)
        warnings.append(msg)

        # Build diagnostic points before calibration (calibrated_mz = observed_mz)
        diag_points = [
            MatchedCalibrantPoint(
                name=cal.name,
                theoretical_mz=cal.theoretical_mz,
                observed_mz=pk.mz,
                calibrated_mz=pk.mz,
                intensity=pk.intensity,
                error_da_before=calculate_signed_error_da(pk.mz, cal.theoretical_mz),
                error_ppm_before=calculate_signed_error_ppm(pk.mz, cal.theoretical_mz),
                residual_da_after=calculate_signed_error_da(pk.mz, cal.theoretical_mz),
                residual_ppm_after=calculate_signed_error_ppm(pk.mz, cal.theoretical_mz),
            )
            for cal, pk in matched_pairs
        ]

        ppm_before = [abs(p.error_ppm_before) for p in diag_points]
        signed_before = [p.error_ppm_before for p in diag_points]
        rms_before = float(np.sqrt(np.mean(np.square(signed_before)))) if signed_before else None
        med_before = float(np.median(ppm_before)) if ppm_before else None
        mean_before = float(np.mean(signed_before)) if signed_before else None
        max_before = float(np.max(ppm_before)) if ppm_before else None

        identity_coeffs = (1.0, 0.0) if settings.model_type == CalibrationModelType.LINEAR else (0.0, 1.0, 0.0)
        return MassCalibrationResult(
            sample_id=spectrum.sample_id,
            applied=False,
            model_type=settings.model_type,
            coefficients=identity_coeffs,
            matched_calibrants=tuple(diag_points),
            unmatched_calibrants=tuple(unmatched),
            rms_ppm_before=rms_before,
            rms_ppm_after=rms_before,
            median_abs_ppm_before=med_before,
            median_abs_ppm_after=med_before,
            mean_signed_ppm_before=mean_before,
            mean_signed_ppm_after=mean_before,
            max_abs_ppm_before=max_before,
            max_abs_ppm_after=max_before,
            calibrant_mz_span_da=0.0,
            calibrated_spectrum=spectrum,
            warnings=tuple(warnings),
        )

    if len(matched_pairs) < recommended_min:
        warnings.append(
            f"Low degrees of freedom: {len(matched_pairs)} calibrants used for {settings.model_type.value} "
            f"fit (>= {recommended_min} recommended to evaluate non-zero residuals)."
        )

    obs_mz = np.array([pk.mz for _, pk in matched_pairs], dtype=np.float64)
    theo_mz = np.array([cal.theoretical_mz for cal, _ in matched_pairs], dtype=np.float64)

    degree = 1 if settings.model_type == CalibrationModelType.LINEAR else 2
    coeffs_arr = np.polyfit(obs_mz, theo_mz, deg=degree)
    coeffs_tuple = tuple(float(c) for c in coeffs_arr)

    cal_obs_mz = apply_calibration_coefficients(obs_mz.tolist(), settings.model_type, coeffs_tuple)

    diag_points: List[MatchedCalibrantPoint] = []
    for (cal, pk), cal_mz_val in zip(matched_pairs, cal_obs_mz):
        err_da_pre = calculate_signed_error_da(pk.mz, cal.theoretical_mz)
        err_ppm_pre = calculate_signed_error_ppm(pk.mz, cal.theoretical_mz)
        res_da_post = calculate_signed_error_da(cal_mz_val, cal.theoretical_mz)
        res_ppm_post = calculate_signed_error_ppm(cal_mz_val, cal.theoretical_mz)
        diag_points.append(
            MatchedCalibrantPoint(
                name=cal.name,
                theoretical_mz=cal.theoretical_mz,
                observed_mz=pk.mz,
                calibrated_mz=cal_mz_val,
                intensity=pk.intensity,
                error_da_before=err_da_pre,
                error_ppm_before=err_ppm_pre,
                residual_da_after=res_da_post,
                residual_ppm_after=res_ppm_post,
            )
        )

    abs_ppm_before = [abs(p.error_ppm_before) for p in diag_points]
    signed_ppm_before = [p.error_ppm_before for p in diag_points]
    abs_ppm_after = [abs(p.residual_ppm_after) for p in diag_points]
    signed_ppm_after = [p.residual_ppm_after for p in diag_points]

    rms_before = float(np.sqrt(np.mean(np.square(signed_ppm_before))))
    rms_after = float(np.sqrt(np.mean(np.square(signed_ppm_after))))
    med_before = float(np.median(abs_ppm_before))
    med_after = float(np.median(abs_ppm_after))
    mean_before = float(np.mean(signed_ppm_before))
    mean_after = float(np.mean(signed_ppm_after))
    max_before = float(np.max(abs_ppm_before))
    max_after = float(np.max(abs_ppm_after))

    calibrant_span = float(np.max(obs_mz) - np.min(obs_mz))
    spec_span = (
        (spectrum.max_mz - spectrum.min_mz)
        if (spectrum.max_mz is not None and spectrum.min_mz is not None)
        else 0.0
    )
    if spec_span > 100.0 and calibrant_span < 0.30 * spec_span:
        warnings.append(
            f"Narrow calibrant m/z span ({calibrant_span:.1f} Da vs spectrum span {spec_span:.1f} Da); "
            f"extrapolating outside [{ float(np.min(obs_mz)):.1f}, {float(np.max(obs_mz)):.1f}] may increase mass error."
        )

    if max_after > settings.max_residual_ppm_threshold:
        warnings.append(
            f"Post-calibration maximum residual ({max_after:.2f} ppm) exceeds threshold "
            f"({settings.max_residual_ppm_threshold:.2f} ppm); inspect calibrant assignments for misidentification."
        )

    # Apply calibration to all peaks in the spectrum
    all_obs_mz = [p.mz for p in spectrum.peaks]
    all_cal_mz = apply_calibration_coefficients(all_obs_mz, settings.model_type, coeffs_tuple)

    calibrated_peaks = tuple(
        SpectrumPeak(
            mz=float(new_mz),
            intensity=p.intensity,
            area=p.area,
            snr=p.snr,
            fwhm=p.fwhm,
            charge=p.charge,
            row_index=p.row_index,
        )
        for p, new_mz in zip(spectrum.peaks, all_cal_mz)
    )

    new_metadata: Dict[str, object] = dict(spectrum.metadata)
    new_metadata.update(
        {
            "is_calibrated": True,
            "calibration_model": settings.model_type.value,
            "calibration_coefficients": coeffs_tuple,
            "calibrants_matched_count": len(diag_points),
            "rms_ppm_before_calibration": round(rms_before, 4),
            "rms_ppm_after_calibration": round(rms_after, 4),
        }
    )

    calibrated_spec = replace(
        spectrum,
        peaks=calibrated_peaks,
        min_mz=min((p.mz for p in calibrated_peaks), default=0.0),
        max_mz=max((p.mz for p in calibrated_peaks), default=0.0),
        warnings=spectrum.warnings + tuple(warnings),
        metadata=new_metadata,
    )

    return MassCalibrationResult(
        sample_id=spectrum.sample_id,
        applied=True,
        model_type=settings.model_type,
        coefficients=coeffs_tuple,
        matched_calibrants=tuple(diag_points),
        unmatched_calibrants=tuple(unmatched),
        rms_ppm_before=rms_before,
        rms_ppm_after=rms_after,
        median_abs_ppm_before=med_before,
        median_abs_ppm_after=med_after,
        mean_signed_ppm_before=mean_before,
        mean_signed_ppm_after=mean_after,
        max_abs_ppm_before=max_before,
        max_abs_ppm_after=max_after,
        calibrant_mz_span_da=calibrant_span,
        calibrated_spectrum=calibrated_spec,
        warnings=tuple(warnings),
    )

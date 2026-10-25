"""Quality control (QC) evaluation engine for MALDI-GlycoID (Stage 8).

Implements MassyTools-inspired single-spectrum and multi-sample batch quality
control metrics with explicit threshold provenance classification:
- `USER_DEFINED`: Explicitly configured by the analyst for the study.
- `INSTRUMENT_SPECIFIC`: Derived from MALDI-TOF instrument mode (reflectron vs. linear).
- `LITERATURE_DERIVED`: Grounded in published glycomics QC conventions (e.g., MassyTools S/N >= 9,
  annotated signal fraction, and mass accuracy stability).
- `PROVISIONAL_DEFAULT`: Conservative software default where universal consensus does not exist.

Scientific Safeguards:
1. Profile-dependent metrics (median S/N, fraction S/N >= 3 or >= 10, FWHM, resolving power)
   are computed ONLY when profile spectra or S/N-annotated peaks are provided; they are
   strictly `None` (`profile_metrics_available = False`) on centroid-only peak lists.
2. Every `QCCriterionCheck` records its numerical observation, pass/warn thresholds,
   resulting `QCStatus` (`PASS`, `WARN`, `FAIL`), and `QCThresholdSource`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from app.analysis.batch_processing import BatchAnalysisResult, BatchItemStatus
from app.analysis.glycan_search import MS1SearchResult
from app.analysis.quantification import MultiSampleQuantificationMatrix, SampleQuantificationResult
from app.core.calibration import MassCalibrationResult
from app.parsers.base_models import SpectrumData, SpectrumDataType


class QCStatus(str, Enum):
    """Quality control evaluation outcome."""

    PASS = "PASS"
    WARN = "WARN"
    FAIL = "FAIL"


class QCThresholdSource(str, Enum):
    """Provenance category of a quality control threshold."""

    USER_DEFINED = "user_defined"
    INSTRUMENT_SPECIFIC = "instrument_specific"
    LITERATURE_DERIVED = "literature_derived"
    PROVISIONAL_DEFAULT = "provisional_default"


@dataclass(frozen=True)
class QCCriterionCheck:
    """Single evaluated QC rule with explicit threshold provenance."""

    metric_name: str
    observed_value: Optional[float]
    pass_threshold: str
    warn_threshold: str
    status: QCStatus
    threshold_source: QCThresholdSource
    rationale: str


@dataclass(frozen=True)
class QCToleranceThresholds:
    """Configurable QC thresholds and their provenance sources."""

    min_matched_peaks_pass: int = 5
    min_matched_peaks_warn: int = 2
    min_matched_peaks_source: QCThresholdSource = QCThresholdSource.PROVISIONAL_DEFAULT

    min_match_rate_pct_pass: float = 20.0
    min_match_rate_pct_warn: float = 10.0
    min_match_rate_source: QCThresholdSource = QCThresholdSource.PROVISIONAL_DEFAULT

    max_median_abs_ppm_pass: float = 25.0
    max_median_abs_ppm_warn: float = 50.0
    max_median_abs_ppm_source: QCThresholdSource = QCThresholdSource.INSTRUMENT_SPECIFIC

    max_abs_mean_signed_ppm_pass: float = 20.0
    max_abs_mean_signed_ppm_warn: float = 40.0
    max_abs_mean_signed_ppm_source: QCThresholdSource = QCThresholdSource.INSTRUMENT_SPECIFIC

    min_annotated_intensity_fraction_pct_pass: float = 30.0
    min_annotated_intensity_fraction_pct_warn: float = 15.0
    min_annotated_intensity_source: QCThresholdSource = QCThresholdSource.LITERATURE_DERIVED

    min_median_sn_pass: float = 9.0
    min_median_sn_warn: float = 3.0
    min_median_sn_source: QCThresholdSource = QCThresholdSource.LITERATURE_DERIVED

    max_calibration_rms_ppm_pass: float = 15.0
    max_calibration_rms_ppm_warn: float = 30.0
    max_calibration_rms_source: QCThresholdSource = QCThresholdSource.INSTRUMENT_SPECIFIC

    max_median_replicate_cv_pct_pass: float = 20.0
    max_median_replicate_cv_pct_warn: float = 35.0
    max_replicate_cv_source: QCThresholdSource = QCThresholdSource.LITERATURE_DERIVED

    max_batch_missingness_pct_pass: float = 30.0
    max_batch_missingness_pct_warn: float = 50.0
    max_batch_missingness_source: QCThresholdSource = QCThresholdSource.PROVISIONAL_DEFAULT


@dataclass(frozen=True)
class SpectrumQCReport:
    """Comprehensive quality control report for a single spectrum."""

    sample_id: str
    overall_status: QCStatus
    data_type: str
    imported_peaks_count: int
    eligible_peaks_count: int
    matched_peaks_count: int
    match_rate_pct: float
    quantified_glycans_count: int
    median_abs_ppm_error: Optional[float]
    mean_signed_ppm_error: Optional[float]
    iqr_ppm_error: Optional[float]
    rms_ppm_error: Optional[float]
    p95_abs_ppm_error: Optional[float]
    total_spectrum_intensity: float
    annotated_glycan_intensity: float
    annotated_intensity_fraction_pct: float
    base_peak_intensity: float
    dynamic_range_log10: Optional[float]
    profile_metrics_available: bool
    median_matched_sn: Optional[float]
    fraction_matched_sn_ge_3: Optional[float]
    fraction_matched_sn_ge_10: Optional[float]
    median_fwhm_da: Optional[float]
    mean_resolving_power: Optional[float]
    is_calibrated: bool
    calibration_rms_ppm_before: Optional[float]
    calibration_rms_ppm_after: Optional[float]
    criteria_checks: Tuple[QCCriterionCheck, ...]
    warnings: Tuple[str, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class BatchQCReport:
    """Multi-sample batch quality control summary and outlier diagnostics."""

    batch_id: str
    overall_status: QCStatus
    total_samples: int
    passed_samples_count: int
    warned_samples_count: int
    failed_samples_count: int
    per_sample_reports: Tuple[SpectrumQCReport, ...]
    batch_missingness_pct: Optional[float]
    median_replicate_cv_pct: Optional[float]
    max_replicate_cv_pct: Optional[float]
    flagged_outlier_sample_ids: Tuple[str, ...]
    batch_criteria_checks: Tuple[QCCriterionCheck, ...]
    warnings: Tuple[str, ...] = field(default_factory=tuple)


def _combine_statuses(statuses: Sequence[QCStatus]) -> QCStatus:
    if any(s == QCStatus.FAIL for s in statuses):
        return QCStatus.FAIL
    if any(s == QCStatus.WARN for s in statuses):
        return QCStatus.WARN
    return QCStatus.PASS


def evaluate_spectrum_qc(
    spectrum: SpectrumData,
    search_result: MS1SearchResult,
    quantification_result: Optional[SampleQuantificationResult] = None,
    calibration_result: Optional[MassCalibrationResult] = None,
    thresholds: Optional[QCToleranceThresholds] = None,
) -> SpectrumQCReport:
    """Evaluate MassyTools-inspired quality control metrics for a single spectrum."""
    if thresholds is None:
        thresholds = QCToleranceThresholds()

    warnings: List[str] = []
    checks: List[QCCriterionCheck] = []

    imported_cnt = spectrum.valid_points_count
    eligible_cnt = search_result.eligible_peaks_searched
    matched_cnt = search_result.matched_peaks_count
    match_rate = (matched_cnt / eligible_cnt * 100.0) if eligible_cnt > 0 else 0.0
    quant_cnt = (
        len(quantification_result.glycan_abundances)
        if quantification_result is not None
        else len({m.glycan_composition for m in search_result.top_matches_per_peak})
    )

    # 1. Matched peaks count check
    if matched_cnt >= thresholds.min_matched_peaks_pass:
        mp_status = QCStatus.PASS
    elif matched_cnt >= thresholds.min_matched_peaks_warn:
        mp_status = QCStatus.WARN
    else:
        mp_status = QCStatus.FAIL
    checks.append(
        QCCriterionCheck(
            metric_name="matched_peaks_count",
            observed_value=float(matched_cnt),
            pass_threshold=f">= {thresholds.min_matched_peaks_pass}",
            warn_threshold=f">= {thresholds.min_matched_peaks_warn}",
            status=mp_status,
            threshold_source=thresholds.min_matched_peaks_source,
            rationale="Minimum number of annotated N-glycan peaks required for reliable compositional profiling.",
        )
    )

    # 2. Match rate % check
    if match_rate >= thresholds.min_match_rate_pct_pass:
        mr_status = QCStatus.PASS
    elif match_rate >= thresholds.min_match_rate_pct_warn:
        mr_status = QCStatus.WARN
    else:
        mr_status = QCStatus.FAIL
    checks.append(
        QCCriterionCheck(
            metric_name="match_rate_pct",
            observed_value=round(match_rate, 2),
            pass_threshold=f">= {thresholds.min_match_rate_pct_pass:.1f}%",
            warn_threshold=f">= {thresholds.min_match_rate_pct_warn:.1f}%",
            status=mr_status,
            threshold_source=thresholds.min_match_rate_source,
            rationale="Percentage of eligible MS1 peaks assigned to candidate glycan compositions.",
        )
    )

    # 3. Mass error metrics from top-ranked match per peak
    top_matches = search_result.top_matches_per_peak
    if top_matches:
        abs_ppms = np.array([m.abs_error_ppm for m in top_matches], dtype=np.float64)
        signed_ppms = np.array([m.signed_error_ppm for m in top_matches], dtype=np.float64)
        med_abs_ppm = float(np.median(abs_ppms))
        mean_signed_ppm = float(np.mean(signed_ppms))
        q75, q25 = np.percentile(abs_ppms, [75, 25])
        iqr_ppm = float(q75 - q25)
        rms_ppm = float(np.sqrt(np.mean(np.square(signed_ppms))))
        p95_abs_ppm = float(np.percentile(abs_ppms, 95))

        if med_abs_ppm <= thresholds.max_median_abs_ppm_pass:
            ppm_status = QCStatus.PASS
        elif med_abs_ppm <= thresholds.max_median_abs_ppm_warn:
            ppm_status = QCStatus.WARN
        else:
            ppm_status = QCStatus.FAIL
        checks.append(
            QCCriterionCheck(
                metric_name="median_abs_ppm_error",
                observed_value=round(med_abs_ppm, 4),
                pass_threshold=f"<= {thresholds.max_median_abs_ppm_pass:.1f} ppm",
                warn_threshold=f"<= {thresholds.max_median_abs_ppm_warn:.1f} ppm",
                status=ppm_status,
                threshold_source=thresholds.max_median_abs_ppm_source,
                rationale="Reflectron MALDI-TOF mass accuracy stability across matched glycan ions.",
            )
        )

        abs_bias = abs(mean_signed_ppm)
        if abs_bias <= thresholds.max_abs_mean_signed_ppm_pass:
            bias_status = QCStatus.PASS
        elif abs_bias <= thresholds.max_abs_mean_signed_ppm_warn:
            bias_status = QCStatus.WARN
        else:
            bias_status = QCStatus.FAIL
        checks.append(
            QCCriterionCheck(
                metric_name="mean_signed_ppm_error",
                observed_value=round(mean_signed_ppm, 4),
                pass_threshold=f"|mean| <= {thresholds.max_abs_mean_signed_ppm_pass:.1f} ppm",
                warn_threshold=f"|mean| <= {thresholds.max_abs_mean_signed_ppm_warn:.1f} ppm",
                status=bias_status,
                threshold_source=thresholds.max_abs_mean_signed_ppm_source,
                rationale="Detects systematic mass calibration drift (uncalibrated TOF shift).",
            )
        )
    else:
        med_abs_ppm = None
        mean_signed_ppm = None
        iqr_ppm = None
        rms_ppm = None
        p95_abs_ppm = None
        warnings.append("No matched glycan peaks available to compute mass error QC statistics.")

    # 4. Intensity distribution metrics
    total_spec_intensity = float(sum(p.intensity for p in spectrum.peaks))
    base_pk_intensity = float(max((p.intensity for p in spectrum.peaks), default=0.0))
    annotated_intensity = float(sum(m.observed_intensity for m in top_matches))
    annotated_fraction_pct = (
        (annotated_intensity / total_spec_intensity * 100.0) if total_spec_intensity > 0 else 0.0
    )

    matched_intensities = [m.observed_intensity for m in top_matches if m.observed_intensity > 0]
    if len(matched_intensities) >= 2:
        dyn_range_log10 = float(np.log10(max(matched_intensities) / min(matched_intensities)))
    else:
        dyn_range_log10 = None

    if annotated_fraction_pct >= thresholds.min_annotated_intensity_fraction_pct_pass:
        int_status = QCStatus.PASS
    elif annotated_fraction_pct >= thresholds.min_annotated_intensity_fraction_pct_warn:
        int_status = QCStatus.WARN
    else:
        int_status = QCStatus.FAIL
    checks.append(
        QCCriterionCheck(
            metric_name="annotated_intensity_fraction_pct",
            observed_value=round(annotated_fraction_pct, 2),
            pass_threshold=f">= {thresholds.min_annotated_intensity_fraction_pct_pass:.1f}%",
            warn_threshold=f">= {thresholds.min_annotated_intensity_fraction_pct_warn:.1f}%",
            status=int_status,
            threshold_source=thresholds.min_annotated_intensity_source,
            rationale="MassyTools-inspired total analyte signal fraction explained by annotated N-glycans.",
        )
    )

    # 5. Profile / S/N / FWHM metrics (ONLY if profile data or S/N attributes are present)
    # Match observed m/z back to SpectrumPeak objects for S/N, FWHM, and resolving power
    matched_mz_set = {round(m.observed_mz, 4) for m in top_matches}
    matched_spec_peaks = [p for p in spectrum.peaks if round(p.mz, 4) in matched_mz_set]
    sn_values = [p.snr for p in matched_spec_peaks if p.snr is not None]
    fwhm_values = [p.fwhm for p in matched_spec_peaks if p.fwhm is not None and p.fwhm > 0]
    rp_values = [
        p.mz / p.fwhm
        for p in matched_spec_peaks
        if p.fwhm is not None and p.fwhm > 0
    ]

    if sn_values:
        profile_metrics_available = True
        med_sn = float(np.median(sn_values))
        frac_sn_3 = float(sum(1 for v in sn_values if v >= 3.0) / len(sn_values))
        frac_sn_10 = float(sum(1 for v in sn_values if v >= 10.0) / len(sn_values))
        med_fwhm = float(np.median(fwhm_values)) if fwhm_values else None
        mean_rp = float(np.mean(rp_values)) if rp_values else None

        if med_sn >= thresholds.min_median_sn_pass:
            sn_status = QCStatus.PASS
        elif med_sn >= thresholds.min_median_sn_warn:
            sn_status = QCStatus.WARN
        else:
            sn_status = QCStatus.FAIL
        checks.append(
            QCCriterionCheck(
                metric_name="median_matched_sn",
                observed_value=round(med_sn, 2),
                pass_threshold=f">= {thresholds.min_median_sn_pass:.1f}",
                warn_threshold=f">= {thresholds.min_median_sn_warn:.1f}",
                status=sn_status,
                threshold_source=thresholds.min_median_sn_source,
                rationale="MassyTools-inspired signal-to-noise criterion (S/N >= 9 for quantitation, >= 3 for LOD).",
            )
        )
    else:
        profile_metrics_available = False
        med_sn = None
        frac_sn_3 = None
        frac_sn_10 = None
        med_fwhm = None
        mean_rp = None
        if spectrum.data_type == SpectrumDataType.CENTROID:
            warnings.append(
                "Profile S/N and FWHM QC metrics are unavailable on centroid peak lists without S/N columns."
            )

    # 6. Calibration residual check (if calibration result provided or tagged in metadata)
    is_calibrated = bool(
        (calibration_result is not None and calibration_result.applied)
        or spectrum.metadata.get("is_calibrated", False)
    )
    cal_rms_before: Optional[float] = None
    cal_rms_after: Optional[float] = None

    if calibration_result is not None:
        cal_rms_before = calibration_result.rms_ppm_before
        cal_rms_after = calibration_result.rms_ppm_after
        warnings.extend(calibration_result.warnings)
        if calibration_result.applied and cal_rms_after is not None:
            if cal_rms_after <= thresholds.max_calibration_rms_ppm_pass:
                cal_status = QCStatus.PASS
            elif cal_rms_after <= thresholds.max_calibration_rms_ppm_warn:
                cal_status = QCStatus.WARN
            else:
                cal_status = QCStatus.FAIL
            checks.append(
                QCCriterionCheck(
                    metric_name="calibration_rms_ppm_after",
                    observed_value=round(cal_rms_after, 4),
                    pass_threshold=f"<= {thresholds.max_calibration_rms_ppm_pass:.1f} ppm",
                    warn_threshold=f"<= {thresholds.max_calibration_rms_ppm_warn:.1f} ppm",
                    status=cal_status,
                    threshold_source=thresholds.max_calibration_rms_source,
                    rationale="Root-mean-square mass residual of reference calibrants after polynomial fit.",
                )
            )

    overall_status = _combine_statuses([c.status for c in checks])

    return SpectrumQCReport(
        sample_id=spectrum.sample_id,
        overall_status=overall_status,
        data_type=spectrum.data_type.value,
        imported_peaks_count=imported_cnt,
        eligible_peaks_count=eligible_cnt,
        matched_peaks_count=matched_cnt,
        match_rate_pct=round(match_rate, 2),
        quantified_glycans_count=quant_cnt,
        median_abs_ppm_error=round(med_abs_ppm, 4) if med_abs_ppm is not None else None,
        mean_signed_ppm_error=round(mean_signed_ppm, 4) if mean_signed_ppm is not None else None,
        iqr_ppm_error=round(iqr_ppm, 4) if iqr_ppm is not None else None,
        rms_ppm_error=round(rms_ppm, 4) if rms_ppm is not None else None,
        p95_abs_ppm_error=round(p95_abs_ppm, 4) if p95_abs_ppm is not None else None,
        total_spectrum_intensity=total_spec_intensity,
        annotated_glycan_intensity=annotated_intensity,
        annotated_intensity_fraction_pct=round(annotated_fraction_pct, 2),
        base_peak_intensity=base_pk_intensity,
        dynamic_range_log10=round(dyn_range_log10, 3) if dyn_range_log10 is not None else None,
        profile_metrics_available=profile_metrics_available,
        median_matched_sn=round(med_sn, 2) if med_sn is not None else None,
        fraction_matched_sn_ge_3=round(frac_sn_3, 4) if frac_sn_3 is not None else None,
        fraction_matched_sn_ge_10=round(frac_sn_10, 4) if frac_sn_10 is not None else None,
        median_fwhm_da=round(med_fwhm, 4) if med_fwhm is not None else None,
        mean_resolving_power=round(mean_rp, 1) if mean_rp is not None else None,
        is_calibrated=is_calibrated,
        calibration_rms_ppm_before=round(cal_rms_before, 4) if cal_rms_before is not None else None,
        calibration_rms_ppm_after=round(cal_rms_after, 4) if cal_rms_after is not None else None,
        criteria_checks=tuple(checks),
        warnings=tuple(warnings),
    )


def evaluate_batch_qc(
    batch_result: BatchAnalysisResult,
    thresholds: Optional[QCToleranceThresholds] = None,
    calibration_results_by_sample: Optional[Dict[str, MassCalibrationResult]] = None,
) -> BatchQCReport:
    """Evaluate cohort-wide batch QC metrics, missingness, replicate CVs, and outlier samples."""
    if thresholds is None:
        thresholds = QCToleranceThresholds()
    if calibration_results_by_sample is None:
        calibration_results_by_sample = {}

    sample_reports: List[SpectrumQCReport] = []
    outliers: List[str] = []
    warnings: List[str] = list(batch_result.warnings)
    batch_checks: List[QCCriterionCheck] = []

    for outcome in batch_result.sample_outcomes:
        if outcome.status == BatchItemStatus.COMPLETED and outcome.spectrum and outcome.ms1_result:
            cal_res = calibration_results_by_sample.get(outcome.sample_id)
            rep = evaluate_spectrum_qc(
                spectrum=outcome.spectrum,
                search_result=outcome.ms1_result,
                quantification_result=outcome.quantification,
                calibration_result=cal_res,
                thresholds=thresholds,
            )
            sample_reports.append(rep)
            if rep.overall_status == QCStatus.FAIL:
                outliers.append(outcome.sample_id)
        elif outcome.status == BatchItemStatus.FAILED:
            outliers.append(outcome.sample_id)
            warnings.append(f"Sample '{outcome.sample_id}' failed processing and is flagged as a QC failure.")

    passed_cnt = sum(1 for r in sample_reports if r.overall_status == QCStatus.PASS)
    warned_cnt = sum(1 for r in sample_reports if r.overall_status == QCStatus.WARN)
    failed_cnt = sum(1 for r in sample_reports if r.overall_status == QCStatus.FAIL) + batch_result.failed_files_count

    # Evaluate batch missingness and replicate CV from the abundance matrix
    matrix: Optional[MultiSampleQuantificationMatrix] = batch_result.abundance_matrix
    batch_missingness_pct: Optional[float] = None
    med_cv_pct: Optional[float] = None
    max_cv_pct: Optional[float] = None

    if matrix is not None and matrix.replicate_summaries:
        missing_vals = [s.missingness_pct for s in matrix.replicate_summaries]
        cv_vals = [
            s.cv_pct
            for s in matrix.replicate_summaries
            if s.cv_pct is not None
        ]

        if missing_vals:
            batch_missingness_pct = float(np.mean(missing_vals))
            if batch_missingness_pct <= thresholds.max_batch_missingness_pct_pass:
                miss_status = QCStatus.PASS
            elif batch_missingness_pct <= thresholds.max_batch_missingness_pct_warn:
                miss_status = QCStatus.WARN
            else:
                miss_status = QCStatus.FAIL
            batch_checks.append(
                QCCriterionCheck(
                    metric_name="batch_missingness_pct",
                    observed_value=round(batch_missingness_pct, 2),
                    pass_threshold=f"<= {thresholds.max_batch_missingness_pct_pass:.1f}%",
                    warn_threshold=f"<= {thresholds.max_batch_missingness_pct_warn:.1f}%",
                    status=miss_status,
                    threshold_source=thresholds.max_batch_missingness_source,
                    rationale="Mean percentage of missing glycan observations across batch samples (unimputed).",
                )
            )

        if cv_vals:
            med_cv_pct = float(np.median(cv_vals))
            max_cv_pct = float(np.max(cv_vals))
            if med_cv_pct <= thresholds.max_median_replicate_cv_pct_pass:
                cv_status = QCStatus.PASS
            elif med_cv_pct <= thresholds.max_median_replicate_cv_pct_warn:
                cv_status = QCStatus.WARN
            else:
                cv_status = QCStatus.FAIL
            batch_checks.append(
                QCCriterionCheck(
                    metric_name="median_replicate_cv_pct",
                    observed_value=round(med_cv_pct, 2),
                    pass_threshold=f"<= {thresholds.max_median_replicate_cv_pct_pass:.1f}%",
                    warn_threshold=f"<= {thresholds.max_median_replicate_cv_pct_warn:.1f}%",
                    status=cv_status,
                    threshold_source=thresholds.max_replicate_cv_source,
                    rationale="Median coefficient of variation (%) of relative abundance across replicates.",
                )
            )

    all_statuses = [r.overall_status for r in sample_reports] + [c.status for c in batch_checks]
    if batch_result.failed_files_count > 0:
        all_statuses.append(QCStatus.FAIL)
    overall_batch_status = _combine_statuses(all_statuses) if all_statuses else QCStatus.FAIL

    return BatchQCReport(
        batch_id=batch_result.batch_id,
        overall_status=overall_batch_status,
        total_samples=batch_result.total_files,
        passed_samples_count=passed_cnt,
        warned_samples_count=warned_cnt,
        failed_samples_count=failed_cnt,
        per_sample_reports=tuple(sample_reports),
        batch_missingness_pct=round(batch_missingness_pct, 2) if batch_missingness_pct is not None else None,
        median_replicate_cv_pct=round(med_cv_pct, 2) if med_cv_pct is not None else None,
        max_replicate_cv_pct=round(max_cv_pct, 2) if max_cv_pct is not None else None,
        flagged_outlier_sample_ids=tuple(outliers),
        batch_criteria_checks=tuple(batch_checks),
        warnings=tuple(warnings),
    )

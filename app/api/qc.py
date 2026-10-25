"""Mass calibration, quality control evaluation, and statistical/ML API endpoints (Stage 12)."""

from __future__ import annotations

import base64
from typing import Any, Dict, List

from app.analysis.quality_control import evaluate_spectrum_qc
from app.analysis.quantification import (
    DenominatorSetPolicy,
    GlycanAbundanceRecord,
    GlycanAggregationPolicy,
    SampleQuantificationResult,
    summarize_replicates_and_build_matrix,
)
from app.analysis.statistics import (
    MissingValueImputationMethod,
    render_pca_plot_png,
    run_group_comparisons,
    run_leakage_free_classification,
    run_pca_analysis,
)
from app.api.compat import APIRouter, HTTPException, Request
from app.api.state import session_state
from app.core.calibration import (
    CalibrationModelType,
    CalibrationSettings,
    calibrate_spectrum,
    get_calibrants_for_panel,
    select_best_calibrant_panel,
)
from app.core.mass_error import MassToleranceUnit

router = APIRouter(prefix="/api/qc", tags=["qc"])


@router.post("/calibrate")
def calibrate_spectrum_endpoint(request: Request) -> Dict[str, Any]:
    """Run mass calibration on the imported spectrum (see `run_calibration`)."""
    return run_calibration(request.json_sync())


def run_calibration(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Run linear or quadratic internal/external mass calibration on the imported spectrum.

    Calibration is always computed from the uncalibrated (imported / corrected) spectrum so that
    repeated runs never compound. Run it BEFORE the MS1 search; the search then uses the calibrated
    m/z axis. Any previous search/quantification results are invalidated.

    `reference_set` may be `nglycan_na`, `nglycan_h`, `nglycan_neg`, `peptide_h` or `auto`
    (the panel matching the most peaks is chosen).
    """
    spectrum = session_state.get_base_spectrum()

    try:
        model_str = str(payload.get("model_type", "linear")).strip().lower()
        model_type = CalibrationModelType.QUADRATIC if model_str == "quadratic" else CalibrationModelType.LINEAR
        ref_set = str(payload.get("reference_set", "nglycan_na")).strip().lower()
        tol_ppm = float(payload.get("match_tolerance_ppm", 150.0))
        labeling = str(payload.get("labeling", "NONE")).strip().upper() or "NONE"

        cal_settings = CalibrationSettings(
            model_type=model_type,
            match_tolerance=tol_ppm,
            match_tolerance_unit=MassToleranceUnit.PPM,
            require_calibration_success=False,
        )
        panel_counts: Dict[str, int] = {}
        if ref_set == "auto":
            ref_set, panel_counts = select_best_calibrant_panel(spectrum, cal_settings, labeling=labeling)
        calibrants = get_calibrants_for_panel(ref_set, labeling=labeling)
        cal_res = calibrate_spectrum(
            spectrum=spectrum,
            calibrants=calibrants,
            settings=cal_settings,
        )
        session_state.set_calibration(cal_res)

        return {
            "sample_id": cal_res.sample_id,
            "applied": cal_res.applied,
            "reference_set": ref_set,
            "panel_match_counts": panel_counts,
            "model_type": cal_res.model_type.value,
            "coefficients": list(cal_res.coefficients),
            "matched_calibrants_count": len(cal_res.matched_calibrants),
            "unmatched_calibrants_count": len(cal_res.unmatched_calibrants),
            "rms_ppm_before": cal_res.rms_ppm_before,
            "rms_ppm_after": cal_res.rms_ppm_after,
            "median_abs_ppm_before": cal_res.median_abs_ppm_before,
            "median_abs_ppm_after": cal_res.median_abs_ppm_after,
            "calibrant_mz_span_da": cal_res.calibrant_mz_span_da,
            "matched_calibrants": [
                {
                    "name": pt.name,
                    "theoretical_mz": round(pt.theoretical_mz, 6),
                    "observed_mz": round(pt.observed_mz, 6),
                    "calibrated_mz": round(pt.calibrated_mz, 6),
                    "error_ppm_before": round(pt.error_ppm_before, 4),
                    "residual_ppm_after": round(pt.residual_ppm_after, 4),
                }
                for pt in cal_res.matched_calibrants
            ],
            "warnings": list(cal_res.warnings),
        }
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Mass calibration error: {exc}") from exc


@router.post("/calibrate/reset")
def reset_calibration_endpoint() -> Dict[str, Any]:
    """Discard the calibration and return to the uncalibrated (imported / corrected) spectrum."""
    session_state.get_or_load_spectrum()
    session_state.reset_calibration()
    return {"is_calibrated": False, "message": "Calibration removed; spectrum is uncalibrated."}


@router.post("/evaluate")
def evaluate_qc_endpoint() -> Dict[str, Any]:
    """Evaluate MassyTools-inspired quality control criteria on the active spectrum."""
    session_state.ensure_baseline_analysis()
    spectrum = session_state.get_or_load_spectrum()
    if session_state.ms1_search_result is None:
        raise HTTPException(status_code=400, detail="MS1 search result unavailable for QC evaluation.")

    qc_rep = evaluate_spectrum_qc(
        spectrum=spectrum,
        search_result=session_state.ms1_search_result,
        quantification_result=session_state.quantification_result,
        calibration_result=session_state.calibration_result,
    )
    session_state.spectrum_qc_report = qc_rep
    return {
        "sample_id": qc_rep.sample_id,
        "overall_status": qc_rep.overall_status.value,
        "data_type": qc_rep.data_type,
        "imported_peaks_count": qc_rep.imported_peaks_count,
        "eligible_peaks_count": qc_rep.eligible_peaks_count,
        "matched_peaks_count": qc_rep.matched_peaks_count,
        "match_rate_pct": round(qc_rep.match_rate_pct, 2),
        "quantified_glycans_count": qc_rep.quantified_glycans_count,
        "median_abs_ppm_error": (
            round(qc_rep.median_abs_ppm_error, 4) if qc_rep.median_abs_ppm_error is not None else None
        ),
        "mean_signed_ppm_error": (
            round(qc_rep.mean_signed_ppm_error, 4) if qc_rep.mean_signed_ppm_error is not None else None
        ),
        "annotated_intensity_fraction_pct": round(qc_rep.annotated_intensity_fraction_pct, 2),
        "profile_metrics_available": qc_rep.profile_metrics_available,
        "median_matched_sn": qc_rep.median_matched_sn,
        "is_calibrated": qc_rep.is_calibrated,
        "criteria_checks": [
            {
                "metric_name": c.metric_name,
                "observed_value": c.observed_value,
                "pass_threshold": c.pass_threshold,
                "warn_threshold": c.warn_threshold,
                "status": c.status.value,
                "threshold_source": c.threshold_source.value,
                "rationale": c.rationale,
            }
            for c in qc_rep.criteria_checks
        ],
        "warnings": list(qc_rep.warnings),
    }


def _build_synthetic_cohort_quantifications() -> tuple[List[SampleQuantificationResult], Dict[str, str]]:
    """Build a 6-sample (3 Control, 3 Disease) synthetic N-glycan cohort for statistical validation."""
    cohort_profiles = [
        (
            "CTRL_1",
            "Control",
            [
                ("Hex5HexNAc2", "high_mannose", 14.0),
                ("Hex3HexNAc4dHex1", "complex", 32.0),
                ("Hex4HexNAc4dHex1", "complex", 28.0),
                ("Hex5HexNAc4dHex1", "complex", 18.0),
                ("Hex5HexNAc4Neu5Ac2", "complex", 8.0),
            ],
        ),
        (
            "CTRL_2",
            "Control",
            [
                ("Hex5HexNAc2", "high_mannose", 15.0),
                ("Hex3HexNAc4dHex1", "complex", 31.0),
                ("Hex4HexNAc4dHex1", "complex", 27.5),
                ("Hex5HexNAc4dHex1", "complex", 18.5),
                ("Hex5HexNAc4Neu5Ac2", "complex", 8.0),
            ],
        ),
        (
            "CTRL_3",
            "Control",
            [
                ("Hex5HexNAc2", "high_mannose", 13.5),
                ("Hex3HexNAc4dHex1", "complex", 33.0),
                ("Hex4HexNAc4dHex1", "complex", 27.0),
                ("Hex5HexNAc4dHex1", "complex", 17.5),
                ("Hex5HexNAc4Neu5Ac2", "complex", 9.0),
            ],
        ),
        (
            "DIS_1",
            "Disease",
            [
                ("Hex5HexNAc2", "high_mannose", 8.0),
                ("Hex3HexNAc4dHex1", "complex", 16.0),
                ("Hex4HexNAc4dHex1", "complex", 22.0),
                ("Hex5HexNAc4dHex1", "complex", 26.0),
                ("Hex5HexNAc4Neu5Ac2", "complex", 28.0),
            ],
        ),
        (
            "DIS_2",
            "Disease",
            [
                ("Hex5HexNAc2", "high_mannose", 7.5),
                ("Hex3HexNAc4dHex1", "complex", 15.5),
                ("Hex4HexNAc4dHex1", "complex", 23.0),
                ("Hex5HexNAc4dHex1", "complex", 25.0),
                ("Hex5HexNAc4Neu5Ac2", "complex", 29.0),
            ],
        ),
        (
            "DIS_3",
            "Disease",
            [
                ("Hex5HexNAc2", "high_mannose", 8.5),
                ("Hex3HexNAc4dHex1", "complex", 17.0),
                ("Hex4HexNAc4dHex1", "complex", 21.5),
                ("Hex5HexNAc4dHex1", "complex", 25.5),
                ("Hex5HexNAc4Neu5Ac2", "complex", 27.5),
            ],
        ),
    ]

    sample_results: List[SampleQuantificationResult] = []
    sample_groups: Dict[str, str] = {}
    for sid, grp, comp_list in cohort_profiles:
        sample_groups[sid] = grp
        abunds = tuple(
            GlycanAbundanceRecord(
                sample_id=sid,
                glycan_composition=comp,
                candidate_class=g_cls,
                aggregation_policy_used=GlycanAggregationPolicy.PRIMARY_ADDUCT_ONLY,
                primary_adduct_filter="[M+Na]+",
                aggregated_intensity=rel_pct * 1000.0,
                aggregated_area=None,
                base_peak_relative_intensity_pct=rel_pct * 3.0,
                relative_abundance_intensity_pct=rel_pct,
                relative_abundance_area_pct=None,
                contributing_peaks_count=1,
                contributing_adducts=("[M+Na]+",),
                contributing_mzs=(1500.0,),
                peak_provenance_history=(f"synthetic_{comp}",),
                contains_ambiguous_peak=False,
                warnings=(),
            )
            for comp, g_cls, rel_pct in comp_list
        )
        sample_results.append(
            SampleQuantificationResult(
                sample_id=sid,
                denominator_policy=DenominatorSetPolicy.ALL_ELIGIBLE_SPECTRUM_PEAKS,
                aggregation_policy=GlycanAggregationPolicy.PRIMARY_ADDUCT_ONLY,
                primary_adduct="[M+Na]+",
                exclude_ambiguous_from_aggregation=False,
                denominator_max_intensity=35000.0,
                total_annotated_intensity=100000.0,
                total_annotated_area=None,
                has_peak_area_data=False,
                peak_records=(),
                glycan_abundances=abunds,
                warnings=(),
            )
        )
    return sample_results, sample_groups


@router.post("/statistics-demo")
def run_statistics_demo_endpoint(request: Request) -> Dict[str, Any]:
    """Run PCA, BH-FDR group comparisons, and leakage-free cross-validation on a multi-sample cohort."""
    payload = request.json_sync()
    imp_str = str(payload.get("imputation_method", "half_min_observed")).strip().lower()
    clf_str = str(payload.get("classifier_type", "logistic_regression")).strip().lower()

    try:
        imp_method = MissingValueImputationMethod(imp_str)
    except ValueError:
        imp_method = MissingValueImputationMethod.HALF_MIN_OBSERVED

    try:
        sample_results, sample_groups = _build_synthetic_cohort_quantifications()
        matrix = summarize_replicates_and_build_matrix(sample_results)
        pca_res = run_pca_analysis(
            matrix=matrix,
            sample_groups=sample_groups,
            imputation_method=imp_method,
        )
        group_comps = run_group_comparisons(
            matrix=matrix,
            sample_groups=sample_groups,
            imputation_method=imp_method,
        )
        cv_res = run_leakage_free_classification(
            matrix=matrix,
            sample_groups=sample_groups,
            classifier_type=clf_str,
            n_splits=3,
        )
        pca_png = render_pca_plot_png(pca_res)
        pca_b64 = base64.b64encode(pca_png).decode("ascii")

        session_state.abundance_matrix = matrix
        session_state.pca_result = pca_res
        session_state.group_comparisons = group_comps
        session_state.cv_classification_result = cv_res

        return {
            "pca": {
                "n_components": pca_res.n_components,
                "explained_variance_ratio": list(pca_res.explained_variance_ratio),
                "cumulative_variance_ratio": list(pca_res.cumulative_variance_ratio),
                "imputation_method_used": pca_res.imputation_method_used,
                "sample_scores": {k: list(v) for k, v in pca_res.sample_scores.items()},
                "warnings": list(pca_res.warnings),
            },
            "classification": {
                "classifier_name": cv_res.classifier_name,
                "classes": list(cv_res.classes),
                "n_samples": cv_res.n_samples,
                "n_splits": cv_res.n_splits,
                "oof_accuracy": cv_res.oof_accuracy,
                "oof_balanced_accuracy": cv_res.oof_balanced_accuracy,
                "oof_roc_auc": cv_res.oof_roc_auc,
                "top_features_mean_importance": [
                    {"glycan_composition": name, "importance": val}
                    for name, val in cv_res.top_features_mean_importance
                ],
                "leakage_prevention_statement": cv_res.leakage_prevention_statement,
                "warnings": list(cv_res.warnings),
            },
            "group_comparisons": [
                {
                    "glycan_composition": g.glycan_composition,
                    "comparison_type": g.comparison_type,
                    "groups_compared": list(g.groups_compared),
                    "group_means": g.group_means,
                    "log2_fold_change": g.log2_fold_change,
                    "parametric_test_name": g.parametric_test_name,
                    "parametric_p_value": g.parametric_p_value,
                    "parametric_fdr_q_value": g.parametric_fdr_q_value,
                    "nonparametric_test_name": g.nonparametric_test_name,
                    "nonparametric_p_value": g.nonparametric_p_value,
                    "nonparametric_fdr_q_value": g.nonparametric_fdr_q_value,
                    "significant_fdr_0_05": g.significant_fdr_0_05,
                }
                for g in group_comps
            ],
            "pca_plot_base64": pca_b64,
        }
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Statistical suite error: {exc}") from exc

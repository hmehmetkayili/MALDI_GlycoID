"""Unified export and reporting orchestrator for MALDI-GlycoID (Stage 10).

Builds all 12 standardized export categories and packages them into:
- Individual CSV and TSV files
- Multi-sheet `.xlsx` Excel workbooks
- Structured `.json` analytical manifests
- Multi-page `.pdf` analytical summary reports
- Static `.png` spectrum and PCA plots
- Complete downloadable `.zip` export bundles
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import io
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple
import zipfile

from app.analysis.batch_processing import BatchAnalysisResult
from app.analysis.glycan_search import MS1SearchResult
from app.analysis.glycopeptide_search import GlycopeptideScreeningResult
from app.analysis.quality_control import BatchQCReport, SpectrumQCReport
from app.analysis.quantification import MultiSampleQuantificationMatrix, SampleQuantificationResult
from app.analysis.statistics import (
    CorrelationPairRecord,
    CrossValidatedClassificationResult,
    DescriptiveGlycanStat,
    GroupComparisonRecord,
    PCAResult,
    SampleDerivedTraitsReport,
    UnivariateROCResult,
    compute_derived_glycan_traits,
    render_annotated_spectrum_png,
    render_pca_plot_png,
)
from app import __version__ as SOFTWARE_VERSION
from app.core.calibration import MassCalibrationResult
from app.exports.csv_export import export_rows_to_csv_bytes, export_rows_to_tsv_bytes
from app.exports.excel_export import export_tables_to_xlsx_bytes
from app.exports.json_export import export_analysis_to_json_bytes
from app.exports.pdf_export import export_summary_report_to_pdf_bytes
from app.libraries.composition_generator import GlycanLibrary
from app.parsers.base_models import SpectrumData

REQUIRED_12_EXPORT_CATEGORIES: Tuple[str, ...] = (
    "annotated_peaks",
    "all_glycan_candidates",
    "quantified_glycans",
    "replicate_abundance_matrix",
    "derived_glycan_traits",
    "glycopeptide_candidates",
    "calibration_diagnostics",
    "quality_control_summary",
    "statistical_analysis",
    "glycan_reference_library",
    "analysis_parameters_and_reproducibility",
    "processing_and_error_log",
)


@dataclass(frozen=True)
class AnalysisExportBundleInput:
    """Container of single-sample or batch outputs for unified export generation."""

    sample_id: str = "Sample_1"
    spectrum: Optional[SpectrumData] = None
    ms1_search_result: Optional[MS1SearchResult] = None
    quantification_result: Optional[SampleQuantificationResult] = None
    abundance_matrix: Optional[MultiSampleQuantificationMatrix] = None
    derived_traits_reports: Sequence[SampleDerivedTraitsReport] = field(default_factory=tuple)
    glycopeptide_result: Optional[GlycopeptideScreeningResult] = None
    calibration_result: Optional[MassCalibrationResult] = None
    spectrum_qc_report: Optional[SpectrumQCReport] = None
    batch_result: Optional[BatchAnalysisResult] = None
    batch_qc_report: Optional[BatchQCReport] = None
    descriptive_stats: Sequence[DescriptiveGlycanStat] = field(default_factory=tuple)
    group_comparisons: Sequence[GroupComparisonRecord] = field(default_factory=tuple)
    correlations: Sequence[CorrelationPairRecord] = field(default_factory=tuple)
    roc_results: Sequence[UnivariateROCResult] = field(default_factory=tuple)
    pca_result: Optional[PCAResult] = None
    cv_classification_result: Optional[CrossValidatedClassificationResult] = None
    reference_library: Optional[GlycanLibrary] = None


def build_standardized_export_tables(
    bundle: AnalysisExportBundleInput,
) -> Dict[str, List[Dict[str, Any]]]:
    """Build all 12 standardized export categories as lists of row dictionaries."""
    tables: Dict[str, List[Dict[str, Any]]] = {cat: [] for cat in REQUIRED_12_EXPORT_CATEGORIES}

    # 1. annotated_peaks (top match per peak)
    if bundle.ms1_search_result is not None:
        tables["annotated_peaks"] = [m.to_dict() for m in bundle.ms1_search_result.top_matches_per_peak]
    elif bundle.batch_result is not None:
        tables["annotated_peaks"] = [m.to_dict() for m in bundle.batch_result.combined_annotated_peaks]

    # 2. all_glycan_candidates (one-to-many candidate table)
    if bundle.ms1_search_result is not None:
        tables["all_glycan_candidates"] = [m.to_dict() for m in bundle.ms1_search_result.matches]
    elif bundle.batch_result is not None:
        tables["all_glycan_candidates"] = [m.to_dict() for m in bundle.batch_result.combined_glycan_candidates]

    # 3. quantified_glycans
    if bundle.quantification_result is not None:
        tables["quantified_glycans"] = [g.to_dict() for g in bundle.quantification_result.glycan_abundances]
    elif bundle.batch_result is not None:
        rows_q: List[Dict[str, Any]] = []
        for sq in bundle.batch_result.per_sample_quantifications:
            rows_q.extend(g.to_dict() for g in sq.glycan_abundances)
        tables["quantified_glycans"] = rows_q

    # 4. replicate_abundance_matrix
    eff_matrix = bundle.abundance_matrix or (
        bundle.batch_result.abundance_matrix if bundle.batch_result is not None else None
    )
    if eff_matrix is not None:
        matrix_rows: List[Dict[str, Any]] = []
        for s in eff_matrix.replicate_summaries:
            row = {
                "glycan_composition": s.glycan_composition,
                "candidate_class": s.candidate_class,
                "total_replicates": s.total_replicates,
                "observed_replicates": s.observed_replicates,
                "missing_replicates": s.missing_replicates,
                "missingness_pct": round(s.missingness_pct, 2),
                "mean_relative_abundance_pct": s.mean_relative_abundance_pct,
                "sd_relative_abundance_pct": s.sd_relative_abundance_pct,
                "cv_pct": s.cv_pct,
                "imputation_applied": s.imputation_applied,
            }
            for sid, val in s.sample_values_relative_abundance_pct.items():
                row[f"sample_{sid}_rel_abundance_pct"] = round(val, 4) if val is not None else None
            matrix_rows.append(row)
        tables["replicate_abundance_matrix"] = matrix_rows

    # 5. derived_glycan_traits
    trait_reports = list(bundle.derived_traits_reports)
    if not trait_reports:
        if bundle.quantification_result is not None:
            trait_reports.append(compute_derived_glycan_traits(bundle.quantification_result))
        elif bundle.batch_result is not None:
            for sq in bundle.batch_result.per_sample_quantifications:
                trait_reports.append(compute_derived_glycan_traits(sq))

    trait_rows: List[Dict[str, Any]] = []
    for tr in trait_reports:
        for t in tr.traits:
            trait_rows.append(
                {
                    "sample_id": tr.sample_id,
                    "trait_code": t.trait_code,
                    "trait_name": t.trait_name,
                    "value": t.value,
                    "unit": t.unit,
                    "is_composition_inferred_hypothesis": t.is_composition_inferred_hypothesis,
                    "structural_disclaimer": t.structural_disclaimer,
                }
            )
    tables["derived_glycan_traits"] = trait_rows

    # 6. glycopeptide_candidates
    if bundle.glycopeptide_result is not None:
        tables["glycopeptide_candidates"] = [m.to_dict() for m in bundle.glycopeptide_result.matches]
    elif bundle.batch_result is not None and bundle.batch_result.combined_glycopeptide_candidates:
        tables["glycopeptide_candidates"] = [
            m.to_dict() for m in bundle.batch_result.combined_glycopeptide_candidates
        ]

    # 7. calibration_diagnostics
    if bundle.calibration_result is not None:
        cal = bundle.calibration_result
        cal_rows: List[Dict[str, Any]] = []
        for pt in cal.matched_calibrants:
            cal_rows.append(
                {
                    "sample_id": cal.sample_id,
                    "calibration_applied": cal.applied,
                    "model_type": cal.model_type.value,
                    "coefficients": list(cal.coefficients),
                    "calibrant_name": pt.name,
                    "theoretical_mz": round(pt.theoretical_mz, 6),
                    "observed_mz": round(pt.observed_mz, 6),
                    "calibrated_mz": round(pt.calibrated_mz, 6),
                    "intensity": round(pt.intensity, 4),
                    "error_da_before": round(pt.error_da_before, 6),
                    "error_ppm_before": round(pt.error_ppm_before, 4),
                    "residual_da_after": round(pt.residual_da_after, 6),
                    "residual_ppm_after": round(pt.residual_ppm_after, 4),
                    "rms_ppm_before": cal.rms_ppm_before,
                    "rms_ppm_after": cal.rms_ppm_after,
                }
            )
        if not cal_rows:
            cal_rows.append(
                {
                    "sample_id": cal.sample_id,
                    "calibration_applied": cal.applied,
                    "model_type": cal.model_type.value,
                    "warnings": "; ".join(cal.warnings),
                }
            )
        tables["calibration_diagnostics"] = cal_rows

    # 8. quality_control_summary
    qc_rows: List[Dict[str, Any]] = []
    if bundle.spectrum_qc_report is not None:
        r = bundle.spectrum_qc_report
        qc_rows.append(
            {
                "sample_id": r.sample_id,
                "overall_status": r.overall_status.value,
                "data_type": r.data_type,
                "imported_peaks_count": r.imported_peaks_count,
                "eligible_peaks_count": r.eligible_peaks_count,
                "matched_peaks_count": r.matched_peaks_count,
                "match_rate_pct": r.match_rate_pct,
                "quantified_glycans_count": r.quantified_glycans_count,
                "median_abs_ppm_error": r.median_abs_ppm_error,
                "mean_signed_ppm_error": r.mean_signed_ppm_error,
                "annotated_intensity_fraction_pct": r.annotated_intensity_fraction_pct,
                "profile_metrics_available": r.profile_metrics_available,
                "median_matched_sn": r.median_matched_sn,
                "is_calibrated": r.is_calibrated,
                "threshold_sources": "; ".join(
                    f"{c.metric_name}:{c.threshold_source.value}({c.status.value})"
                    for c in r.criteria_checks
                ),
            }
        )
    elif bundle.batch_qc_report is not None:
        for r in bundle.batch_qc_report.per_sample_reports:
            qc_rows.append(
                {
                    "sample_id": r.sample_id,
                    "overall_status": r.overall_status.value,
                    "data_type": r.data_type,
                    "imported_peaks_count": r.imported_peaks_count,
                    "eligible_peaks_count": r.eligible_peaks_count,
                    "matched_peaks_count": r.matched_peaks_count,
                    "match_rate_pct": r.match_rate_pct,
                    "quantified_glycans_count": r.quantified_glycans_count,
                    "median_abs_ppm_error": r.median_abs_ppm_error,
                    "mean_signed_ppm_error": r.mean_signed_ppm_error,
                    "annotated_intensity_fraction_pct": r.annotated_intensity_fraction_pct,
                    "profile_metrics_available": r.profile_metrics_available,
                    "median_matched_sn": r.median_matched_sn,
                    "is_calibrated": r.is_calibrated,
                }
            )
    elif bundle.batch_result is not None:
        qc_rows = [dict(row) for row in bundle.batch_result.qc_summary_table]
    tables["quality_control_summary"] = qc_rows

    # 9. statistical_analysis
    stat_rows: List[Dict[str, Any]] = []
    for gc in bundle.group_comparisons:
        stat_rows.append({"record_type": "group_comparison", **asdict(gc)})
    for ds in bundle.descriptive_stats:
        stat_rows.append({"record_type": "descriptive_stat", **asdict(ds)})
    for cp in bundle.correlations:
        stat_rows.append({"record_type": "correlation", **asdict(cp)})
    for rc in bundle.roc_results:
        d = asdict(rc)
        d.pop("fpr_points", None)
        d.pop("tpr_points", None)
        stat_rows.append({"record_type": "univariate_roc", **d})
    if bundle.cv_classification_result is not None:
        stat_rows.append(
            {
                "record_type": "cross_validated_classification",
                **asdict(bundle.cv_classification_result),
            }
        )
    tables["statistical_analysis"] = stat_rows

    # 10. glycan_reference_library
    if bundle.reference_library is not None:
        tables["glycan_reference_library"] = [e.to_dict() for e in bundle.reference_library.entries]

    # 11. analysis_parameters_and_reproducibility
    now_iso = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    if bundle.batch_result is not None:
        tables["analysis_parameters_and_reproducibility"] = [asdict(bundle.batch_result.reproducibility)]
    elif bundle.ms1_search_result is not None:
        s = bundle.ms1_search_result.settings
        tables["analysis_parameters_and_reproducibility"] = [
            {
                "software_version": SOFTWARE_VERSION,
                "sample_id": bundle.ms1_search_result.sample_id,
                "search_mode": s.search_mode.value,
                "mass_tolerance": s.tolerance,
                "mass_tolerance_unit": s.tolerance_unit.value,
                "ionization_polarity": s.polarity.value,
                "adduct_codes": list(s.adduct_codes),
                "charge_states": list(s.charge_states),
                "labeling_method": s.labeling,
                "enable_isotope_scoring": s.enable_isotope_scoring,
                "library_name": bundle.reference_library.library_name if bundle.reference_library else "Custom/Targeted",
                "library_version": bundle.reference_library.library_version if bundle.reference_library else "1.0.0",
                "exported_at_utc": now_iso,
            }
        ]
    else:
        tables["analysis_parameters_and_reproducibility"] = [
            {
                "software_version": SOFTWARE_VERSION,
                "sample_id": bundle.sample_id,
                "exported_at_utc": now_iso,
            }
        ]

    # 12. processing_and_error_log
    log_rows: List[Dict[str, Any]] = []
    if bundle.batch_result is not None:
        for l_entry in bundle.batch_result.processing_log:
            log_rows.append(asdict(l_entry))
        for f_row in bundle.batch_result.failed_files_report:
            log_rows.append(
                {
                    "timestamp_utc": f_row.get("timestamp_utc", now_iso),
                    "level": "FAILED_FILE",
                    "sample_id": f_row.get("sample_id"),
                    "filename": f_row.get("filename"),
                    "message": f_row.get("error_message"),
                }
            )
    else:
        warnings_list: List[str] = []
        if bundle.spectrum is not None:
            warnings_list.extend(bundle.spectrum.warnings)
        if bundle.ms1_search_result is not None:
            warnings_list.extend(bundle.ms1_search_result.warnings)
        if bundle.quantification_result is not None:
            warnings_list.extend(bundle.quantification_result.warnings)
        if bundle.calibration_result is not None:
            warnings_list.extend(bundle.calibration_result.warnings)
        if bundle.spectrum_qc_report is not None:
            warnings_list.extend(bundle.spectrum_qc_report.warnings)

        log_rows.append(
            {
                "timestamp_utc": now_iso,
                "level": "INFO",
                "sample_id": bundle.sample_id,
                "message": "Export bundle generated successfully.",
            }
        )
        for w in warnings_list:
            log_rows.append(
                {
                    "timestamp_utc": now_iso,
                    "level": "WARNING",
                    "sample_id": bundle.sample_id,
                    "message": w,
                }
            )
    tables["processing_and_error_log"] = log_rows

    return tables


def generate_complete_export_zip_bundle(
    bundle: AnalysisExportBundleInput,
    include_library_csv: bool = True,
) -> bytes:
    """Package all 12 export categories across CSV, TSV, multi-sheet XLSX, JSON, PDF, and PNG plots into a ZIP."""
    tables = build_standardized_export_tables(bundle)
    repro_rows = tables.get("analysis_parameters_and_reproducibility", [{}])
    repro_meta = repro_rows[0] if repro_rows else {}

    zip_buf = io.BytesIO()
    with zipfile.ZipFile(zip_buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        # 1. Individual CSV and TSV files for each of the 12 categories
        for idx, cat in enumerate(REQUIRED_12_EXPORT_CATEGORIES, start=1):
            if cat == "glycan_reference_library" and not include_library_csv:
                continue
            rows = tables.get(cat, [])
            zf.writestr(f"csv/{idx:02d}_{cat}.csv", export_rows_to_csv_bytes(rows))
            zf.writestr(f"tsv/{idx:02d}_{cat}.tsv", export_rows_to_tsv_bytes(rows))

        # 2. Multi-sheet Excel (.xlsx) workbook containing all 12 categories
        xlsx_bytes = export_tables_to_xlsx_bytes(tables)
        zf.writestr(f"MALDI_GlycoID_{bundle.sample_id}_Workbook.xlsx", xlsx_bytes)

        # 3. Complete JSON manifest
        json_bytes = export_analysis_to_json_bytes(tables, reproducibility_metadata=repro_meta)
        zf.writestr(f"MALDI_GlycoID_{bundle.sample_id}_Results.json", json_bytes)

        # 4. Multi-page PDF Summary Report
        pdf_bytes = export_summary_report_to_pdf_bytes(
            standardized_tables=tables,
            spectrum=bundle.spectrum,
            search_result=bundle.ms1_search_result,
            quantification_result=bundle.quantification_result,
            qc_report=bundle.spectrum_qc_report,
        )
        zf.writestr(f"MALDI_GlycoID_{bundle.sample_id}_Summary_Report.pdf", pdf_bytes)

        # 5. Static analytical plots (PNG)
        if bundle.spectrum is not None:
            spec_png = render_annotated_spectrum_png(bundle.spectrum, bundle.ms1_search_result)
            zf.writestr("plots/annotated_ms1_spectrum.png", spec_png)
        if bundle.pca_result is not None:
            pca_png = render_pca_plot_png(bundle.pca_result)
            zf.writestr("plots/pca_scores_plot.png", pca_png)

    return zip_buf.getvalue()

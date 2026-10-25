"""Export and reporting download API endpoints for MALDI-GlycoID (Stage 12)."""

from __future__ import annotations

from app.api.compat import APIRouter, HTTPException, Request, Response
from app.api.state import session_state
from app.exports import (
    REQUIRED_12_EXPORT_CATEGORIES,
    AnalysisExportBundleInput,
    build_standardized_export_tables,
    export_analysis_to_json_bytes,
    export_rows_to_csv_bytes,
    export_rows_to_tsv_bytes,
    export_summary_report_to_pdf_bytes,
    export_tables_to_xlsx_bytes,
    generate_complete_export_zip_bundle,
)

router = APIRouter(prefix="/api/exports", tags=["exports"])


def _build_current_export_bundle() -> AnalysisExportBundleInput:
    """Assemble an `AnalysisExportBundleInput` from the active session state."""
    session_state.ensure_baseline_analysis()
    spectrum = session_state.get_or_load_spectrum()
    library = session_state.get_or_load_library()
    traits_list = (
        [session_state.derived_traits_report]
        if session_state.derived_traits_report is not None
        else []
    )
    return AnalysisExportBundleInput(
        sample_id=spectrum.sample_id,
        spectrum=spectrum,
        ms1_search_result=session_state.ms1_search_result,
        quantification_result=session_state.quantification_result,
        abundance_matrix=session_state.abundance_matrix,
        derived_traits_reports=traits_list,
        glycopeptide_result=session_state.glycopeptide_result,
        calibration_result=session_state.calibration_result,
        spectrum_qc_report=session_state.spectrum_qc_report,
        batch_result=session_state.batch_result,
        batch_qc_report=session_state.batch_qc_report,
        group_comparisons=session_state.group_comparisons,
        pca_result=session_state.pca_result,
        cv_classification_result=session_state.cv_classification_result,
        reference_library=library,
    )


@router.get("/download")
def download_export_endpoint(request: Request) -> Response:
    """Download complete ZIP bundle, XLSX workbook, PDF report, JSON manifest, or individual CSV/TSV category."""
    fmt = request.query_params.get("format", "zip").strip().lower()
    category = request.query_params.get("category", "").strip().lower()

    bundle = _build_current_export_bundle()
    tables = build_standardized_export_tables(bundle)

    if fmt == "zip":
        zip_bytes = generate_complete_export_zip_bundle(bundle)
        return Response(
            content=zip_bytes,
            media_type="application/zip",
            headers={
                "Content-Disposition": f'attachment; filename="MALDI_GlycoID_{bundle.sample_id}_Export_Bundle.zip"'
            },
        )

    if fmt == "xlsx":
        xlsx_bytes = export_tables_to_xlsx_bytes(tables)
        return Response(
            content=xlsx_bytes,
            media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            headers={
                "Content-Disposition": f'attachment; filename="MALDI_GlycoID_{bundle.sample_id}_Results.xlsx"'
            },
        )

    if fmt == "pdf":
        pdf_bytes = export_summary_report_to_pdf_bytes(
            standardized_tables=tables,
            spectrum=bundle.spectrum,
            search_result=bundle.ms1_search_result,
            quantification_result=bundle.quantification_result,
            qc_report=bundle.spectrum_qc_report,
        )
        return Response(
            content=pdf_bytes,
            media_type="application/pdf",
            headers={
                "Content-Disposition": f'attachment; filename="MALDI_GlycoID_{bundle.sample_id}_Report.pdf"'
            },
        )

    if fmt == "json":
        repro_rows = tables.get("analysis_parameters_and_reproducibility", [{}])
        repro_meta = repro_rows[0] if repro_rows else {"sample_id": bundle.sample_id}
        json_bytes = export_analysis_to_json_bytes(
            tables,
            reproducibility_metadata=repro_meta,
        )
        return Response(
            content=json_bytes,
            media_type="application/json",
            headers={
                "Content-Disposition": f'attachment; filename="MALDI_GlycoID_{bundle.sample_id}_Manifest.json"'
            },
        )

    if fmt in ("csv", "tsv"):
        target_cat = category or "annotated_peaks"
        if target_cat not in REQUIRED_12_EXPORT_CATEGORIES:
            raise HTTPException(
                status_code=400,
                detail=f"Unknown export category '{target_cat}'. Valid categories: {REQUIRED_12_EXPORT_CATEGORIES}",
            )
        rows = tables[target_cat]
        if fmt == "csv":
            data_bytes = export_rows_to_csv_bytes(rows)
            mime = "text/csv; charset=utf-8"
        else:
            data_bytes = export_rows_to_tsv_bytes(rows)
            mime = "text/tab-separated-values; charset=utf-8"
        return Response(
            content=data_bytes,
            media_type=mime,
            headers={
                "Content-Disposition": f'attachment; filename="{bundle.sample_id}_{target_cat}.{fmt}"'
            },
        )

    raise HTTPException(
        status_code=400,
        detail=f"Unsupported export format '{fmt}'. Supported formats: zip, xlsx, pdf, json, csv, tsv.",
    )

"""Multi-sample batch processing and safe ZIP archive upload API endpoints (Stage 12)."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

from app.analysis.batch_processing import (
    BatchSampleInput,
    execute_batch_analysis,
    extract_batch_zip_safely,
    retry_failed_batch_items,
    sanitize_filename,
)
from app.analysis.glycan_search import MS1SearchSettings, SearchMode
from app.analysis.quality_control import evaluate_batch_qc
from app.api.compat import APIRouter, HTTPException, Request
from app.api.state import session_state
from app.config import EXAMPLES_DIR, settings
from app.core.mass_error import MassToleranceUnit

router = APIRouter(prefix="/api/batch", tags=["batch"])


@router.post("/demo")
def run_demo_batch_endpoint() -> Dict[str, Any]:
    """Run a 3-sample synthetic batch workflow across CSV, XY, and TXT example spectra."""
    library = session_state.get_or_load_library()
    samples = [
        BatchSampleInput(
            sample_id="SYNTH_REP_1_CSV",
            filename="synthetic_native_nglycans.csv",
            content=EXAMPLES_DIR / "synthetic_native_nglycans.csv",
            group="Control",
            replicate_id="Rep1",
        ),
        BatchSampleInput(
            sample_id="SYNTH_REP_2_XY",
            filename="synthetic_native_nglycans.xy",
            content=EXAMPLES_DIR / "synthetic_native_nglycans.xy",
            group="Control",
            replicate_id="Rep2",
        ),
        BatchSampleInput(
            sample_id="SYNTH_REP_3_TXT",
            filename="synthetic_native_nglycans_centroid.txt",
            content=EXAMPLES_DIR / "synthetic_native_nglycans_centroid.txt",
            group="Treated",
            replicate_id="Rep1",
        ),
    ]
    search_settings = MS1SearchSettings(
        search_mode=SearchMode.LIBRARY,
        tolerance=20.0,
        tolerance_unit=MassToleranceUnit.PPM,
        adduct_codes=("[M+Na]+", "[M+H]+", "[M+K]+"),
        charge_states=(1,),
    )
    batch_res = execute_batch_analysis(
        samples=samples,
        library=library,
        search_settings=search_settings,
        batch_id="BATCH_SYNTH_DEMO",
    )
    batch_qc = evaluate_batch_qc(batch_res)
    session_state.batch_result = batch_res
    session_state.abundance_matrix = batch_res.abundance_matrix
    session_state.batch_qc_report = batch_qc
    return batch_res.to_summary_dict()


@router.post("/run")
def run_batch_upload_endpoint(request: Request) -> Dict[str, Any]:
    """Execute fault-isolated batch analysis on uploaded spectrum files or a `.zip` archive."""
    form = request.form_sync()
    uploads = form.getlist("files")
    if not uploads:
        single = form.get("file")
        if single is not None:
            uploads = [single]

    if not uploads:
        raise HTTPException(
            status_code=400,
            detail="No files provided for batch processing under field 'files'.",
        )

    default_group = str(form.get("default_group", "Cohort_A") or "Cohort_A").strip()
    tol_ppm = float(form.get("tolerance_ppm", 20.0) or 20.0)

    batch_inputs: List[BatchSampleInput] = []
    try:
        for idx, up in enumerate(uploads, start=1):
            if not hasattr(up, "filename"):
                continue
            raw_fname = str(up.filename or "")
            if ".." in raw_fname or raw_fname.startswith("/") or raw_fname.startswith("\\"):
                raise HTTPException(
                    status_code=400,
                    detail=f"Unsafe path traversal filename rejected in batch upload: '{raw_fname}'",
                )
            safe_name = sanitize_filename(raw_fname)
            ext = Path(safe_name).suffix.lower()
            if ext == ".zip":
                extracted = extract_batch_zip_safely(up.content)
                for ex in extracted:
                    batch_inputs.append(
                        BatchSampleInput(
                            sample_id=ex.sample_id,
                            filename=ex.filename,
                            content=ex.content,
                            group=default_group,
                            replicate_id=ex.replicate_id,
                        )
                    )
            else:
                batch_inputs.append(
                    BatchSampleInput(
                        sample_id=Path(safe_name).stem or f"Sample_{idx}",
                        filename=safe_name,
                        content=up.content,
                        group=default_group,
                        replicate_id=f"Rep_{idx}",
                    )
                )

        if len(batch_inputs) > settings.max_batch_files:
            raise HTTPException(
                status_code=400,
                detail=f"Batch upload exceeds maximum allowed file count ({settings.max_batch_files}).",
            )

        library = session_state.get_or_load_library()
        search_settings = MS1SearchSettings(
            search_mode=SearchMode.LIBRARY,
            tolerance=tol_ppm,
            tolerance_unit=MassToleranceUnit.PPM,
            adduct_codes=("[M+Na]+", "[M+H]+", "[M+K]+"),
            charge_states=(1,),
        )
        batch_res = execute_batch_analysis(
            samples=batch_inputs,
            library=library,
            search_settings=search_settings,
        )
        session_state.batch_result = batch_res
        session_state.abundance_matrix = batch_res.abundance_matrix
        session_state.batch_qc_report = evaluate_batch_qc(batch_res)
        return batch_res.to_summary_dict()
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Batch processing error: {exc}") from exc


@router.post("/retry")
def retry_batch_endpoint(request: Request) -> Dict[str, Any]:
    """Retry failed batch samples using updated replacement inputs."""
    if session_state.batch_result is None:
        raise HTTPException(status_code=400, detail="No previous batch run available to retry.")

    library = session_state.get_or_load_library()
    updated = retry_failed_batch_items(
        previous_result=session_state.batch_result,
        replacement_samples={},
        library=library,
    )
    session_state.batch_result = updated
    return updated.to_summary_dict()

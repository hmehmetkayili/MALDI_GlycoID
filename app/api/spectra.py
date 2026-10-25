"""Spectrum file upload, validation, and synthetic example loading API endpoints (Stage 12)."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from app.analysis.batch_processing import sanitize_filename
from app.api.compat import APIRouter, HTTPException, Request
from app.api.state import session_state
from app.config import EXAMPLES_DIR, settings
from app.parsers import (
    AmbiguousColumnError,
    BrukerFormatConversionRequiredError,
    DuplicateMzPolicy,
    SpectrumDataType,
    SpectrumParseError,
    parse_spectrum_file,
)

router = APIRouter(prefix="/api/spectra", tags=["spectra"])


@router.get("/example")
def load_example_spectrum(request: Request) -> Dict[str, Any]:
    """Load one of the built-in synthetic validation spectra (`csv`, `xy`, or `txt`)."""
    fmt = request.query_params.get("format", "csv").strip().lower()
    if fmt == "csv":
        example_path = EXAMPLES_DIR / "synthetic_native_nglycans.csv"
        sample_id = "SYNTH_NATIVE_CSV"
    elif fmt == "xy":
        example_path = EXAMPLES_DIR / "synthetic_native_nglycans.xy"
        sample_id = "SYNTH_NATIVE_XY"
    elif fmt == "txt":
        example_path = EXAMPLES_DIR / "synthetic_native_nglycans_centroid.txt"
        sample_id = "SYNTH_NATIVE_TXT"
    else:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown synthetic example format '{fmt}'. Choose from 'csv', 'xy', or 'txt'.",
        )

    if not example_path.is_file():
        raise HTTPException(status_code=404, detail=f"Example file not found: {example_path.name}")

    spectrum = parse_spectrum_file(example_path, sample_id=sample_id)
    session_state.active_spectrum = spectrum
    session_state.ms1_search_result = None
    session_state.quantification_result = None
    session_state.derived_traits_report = None
    session_state.calibration_result = None
    session_state.spectrum_qc_report = None
    summary = spectrum.to_summary_dict(preview_limit=100)
    summary["raw_text_lines"] = example_path.read_text(encoding="utf-8", errors="replace").splitlines()[:25]
    return summary


@router.post("/upload")
def upload_spectrum_endpoint(request: Request) -> Dict[str, Any]:
    """Validate and parse an uploaded mass spectrum or peak list file."""
    form = request.form_sync()
    upload = form.get("file")
    if upload is None or not hasattr(upload, "filename"):
        raise HTTPException(status_code=400, detail="No spectrum file provided under form field 'file'.")

    raw_filename = str(upload.filename or "").strip()
    if not raw_filename:
        raise HTTPException(status_code=400, detail="Uploaded spectrum file must have a non-empty filename.")

    # Reject explicit path traversal attempts
    if ".." in raw_filename or raw_filename.startswith("/") or raw_filename.startswith("\\"):
        raise HTTPException(
            status_code=400,
            detail="Invalid filename: path traversal sequences ('..', leading slashes) are rejected.",
        )

    safe_name = sanitize_filename(raw_filename)
    ext = Path(safe_name).suffix.lower()

    if len(upload.content) == 0:
        raise HTTPException(status_code=400, detail="Uploaded spectrum file is empty (0 bytes).")

    if len(upload.content) > settings.max_upload_size_bytes:
        raise HTTPException(
            status_code=413,
            detail=f"Uploaded file exceeds maximum allowed size ({settings.max_upload_size_bytes} bytes).",
        )

    sample_id_raw = str(form.get("sample_id", "") or "").strip()
    sample_id = sample_id_raw if sample_id_raw else Path(safe_name).stem

    dtype_str = str(form.get("user_specified_data_type", "") or "").strip().lower()
    user_dtype: Optional[SpectrumDataType] = None
    if dtype_str in ("centroid", "profile", "unknown"):
        user_dtype = SpectrumDataType(dtype_str)

    dup_str = str(form.get("duplicate_policy", "keep_max_intensity") or "keep_max_intensity").strip().lower()
    try:
        dup_policy = DuplicateMzPolicy(dup_str)
    except ValueError:
        dup_policy = DuplicateMzPolicy.KEEP_MAX_INTENSITY

    delimiter = str(form.get("delimiter", "") or "") or None
    if delimiter == "\\t":
        delimiter = "\t"
    decimal_sep = str(form.get("decimal_separator", "") or "") or None
    mz_col = str(form.get("mz_column", "") or "").strip() or None
    int_col = str(form.get("intensity_column", "") or "").strip() or None
    area_col = str(form.get("area_column", "") or "").strip() or None
    snr_col = str(form.get("snr_column", "") or "").strip() or None

    try:
        spectrum = parse_spectrum_file(
            content=upload.content,
            filename=safe_name,
            sample_id=sample_id,
            delimiter=delimiter,
            decimal_separator=decimal_sep,
            mz_column=mz_col,
            intensity_column=int_col,
            area_column=area_col,
            snr_column=snr_col,
            user_specified_data_type=user_dtype,
            duplicate_policy=dup_policy,
        )
    except BrukerFormatConversionRequiredError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except AmbiguousColumnError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except SpectrumParseError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Spectrum parsing failed: {exc}") from exc

    session_state.active_spectrum = spectrum
    session_state.ms1_search_result = None
    session_state.quantification_result = None
    session_state.derived_traits_report = None
    session_state.calibration_result = None
    session_state.spectrum_qc_report = None
    summary = spectrum.to_summary_dict(preview_limit=100)
    if ext in (".txt", ".xy", ".csv", ".tsv"):
        summary["raw_text_lines"] = upload.content.decode("utf-8-sig", errors="replace").splitlines()[:25]
    return summary

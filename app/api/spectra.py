"""Spectrum file upload, validation, and synthetic example loading API endpoints (Stage 12)."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from app.analysis.batch_processing import sanitize_filename
from app.analysis.spectrum_corrections import (
    CorrectionConfig,
    apply_corrections,
    list_available_corrections,
)
from app.api.compat import APIRouter, HTTPException, Request
from app.api.state import session_state
from app.config import EXAMPLES_DIR, settings
from app.core.calibration import CALIBRANT_PANELS
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
    """Load one of the built-in synthetic validation spectra (`csv`, `xy`, `txt`, `complex_pos`, or `complex_neg`)."""
    fmt = request.query_params.get("format", "csv").strip().lower()
    if fmt in ("complex_pos", "complex_neg"):
        from app.analysis.synthetic_spectra import generate_complex_nglycan_spectrum

        pol = "negative" if fmt == "complex_neg" else "positive"
        syn = generate_complex_nglycan_spectrum(polarity=pol)
        session_state.set_imported_spectrum(syn.spectrum)
        resp = _spectrum_response(syn.spectrum)
        resp["synthetic_ground_truth_count"] = len({t.composition for t in syn.truth})
        resp["synthetic_injected_ppm_offset"] = syn.ppm_offset
        return resp
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
            detail=(
                f"Unknown synthetic example format '{fmt}'. "
                "Choose from 'csv', 'xy', 'txt', 'complex_pos', or 'complex_neg'."
            ),
        )

    if not example_path.is_file():
        raise HTTPException(status_code=404, detail=f"Example file not found: {example_path.name}")

    spectrum = parse_spectrum_file(example_path, sample_id=sample_id)
    session_state.set_imported_spectrum(spectrum)
    return _spectrum_response(spectrum, raw_lines=example_path.read_text(encoding="utf-8", errors="replace").splitlines()[:25])


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

    # Optional import-time S/N filter (instrument S/N where present, otherwise local-median estimate)
    import_min_snr_raw = str(form.get("min_snr", "") or "").strip()
    import_correction = None
    if import_min_snr_raw:
        try:
            import_min_snr = float(import_min_snr_raw.replace(",", "."))
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=f"Invalid min_snr value '{import_min_snr_raw}'.") from exc
        if import_min_snr < 0:
            raise HTTPException(status_code=400, detail="min_snr must be >= 0.")
        if import_min_snr > 0:
            import_correction = apply_corrections(spectrum, CorrectionConfig(min_snr=import_min_snr))
            spectrum = import_correction.spectrum
            if not spectrum.peaks:
                raise HTTPException(
                    status_code=400,
                    detail=f"No peaks remain after applying S/N >= {import_min_snr}; lower the threshold.",
                )

    session_state.set_imported_spectrum(spectrum)
    if import_correction is not None:
        session_state.correction_result = import_correction
    raw_lines = None
    if ext in (".txt", ".xy", ".csv", ".tsv"):
        raw_lines = upload.content.decode("utf-8-sig", errors="replace").splitlines()[:25]
    return _spectrum_response(spectrum, raw_lines=raw_lines)


def _spectrum_response(spectrum: Any, raw_lines: Optional[list] = None) -> Dict[str, Any]:
    """Standard spectrum summary payload incl. calibration / correction provenance."""
    summary = spectrum.to_summary_dict(preview_limit=100)
    if raw_lines is not None:
        summary["raw_text_lines"] = raw_lines
    summary["is_calibrated"] = bool(spectrum.metadata.get("is_calibrated", False))
    corr = session_state.correction_result
    summary["corrections"] = corr.to_dict() if corr is not None else None
    summary["imported_peak_count"] = (
        len(session_state.original_spectrum.peaks) if session_state.original_spectrum is not None else None
    )
    return summary


@router.get("/corrections")
def list_corrections_endpoint() -> Dict[str, Any]:
    """Catalogue of optional MALDI spectrum corrections and available calibrant panels."""
    return {
        "corrections": list_available_corrections(),
        "calibrant_panels": CALIBRANT_PANELS,
    }


@router.post("/process")
def process_spectrum_endpoint(request: Request) -> Dict[str, Any]:
    """Apply optional corrections / filters to the ORIGINAL imported spectrum (before calibration)."""
    return run_spectrum_processing(request.json_sync())


def run_spectrum_processing(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Apply optional corrections / filters to the ORIGINAL imported spectrum (before calibration).

    Re-running always starts from the imported data, so corrections never compound. Any previous
    calibration and search results are invalidated and must be re-run on the corrected spectrum.
    """
    session_state.get_or_load_spectrum()
    original = session_state.original_spectrum
    try:
        cfg = CorrectionConfig.from_dict(payload or {})
        result = apply_corrections(original, cfg)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Spectrum processing error: {exc}") from exc
    if not result.spectrum.peaks:
        raise HTTPException(
            status_code=400,
            detail="No peaks remain after the selected corrections/filters; relax the thresholds.",
        )
    session_state.set_processed_spectrum(result.spectrum, result)
    return _spectrum_response(result.spectrum)

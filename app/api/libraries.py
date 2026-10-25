"""N-Glycan reference and composition library API endpoints for MALDI-GlycoID (Stage 12)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from app.api.compat import APIRouter, HTTPException, Request, Response
from app.api.state import session_state
from app.config import settings
from app.core.mass_error import MassToleranceUnit
from app.libraries.composition_generator import (
    CompositionConstraints,
    GlycanClass,
    generate_human_n_glycan_library,
)
from app.libraries.library_io import (
    export_library_to_delimited,
    export_library_to_json,
    import_library_from_delimited,
    import_library_from_json,
)
from app.libraries.library_search import LibrarySearchFilter, get_entry_adduct_mz, query_glycan_library
from app.libraries.library_validation import validate_glycan_library

router = APIRouter(prefix="/api/libraries", tags=["libraries"])

# Negative-mode ions shown next to the positive-mode ones in the library table
NEGATIVE_DISPLAY_ADDUCTS = ("[M-H]-", "[M+Cl]-", "[M-2H]2-")


def _entry_row(entry: Any) -> Dict[str, Any]:
    """Serialize a library entry for the UI, always including negative-mode m/z values."""
    row = entry.to_dict()
    adducts = dict(row.get("adduct_mzs", {}))
    negative: Dict[str, float] = {}
    for code in NEGATIVE_DISPLAY_ADDUCTS:
        mz = adducts.get(code)
        if mz is None or mz <= 0:
            mz = round(get_entry_adduct_mz(entry, code), 8)
            adducts[code] = mz
        negative[code] = mz
    row["adduct_mzs"] = adducts
    row["negative_mode_mz"] = negative
    return row


@router.get("/default")
def get_default_library() -> Dict[str, Any]:
    """Reset and return the default 333-entry human N-glycan library."""
    library = session_state.reset_to_default_library()
    return {
        "library_name": library.library_name,
        "library_version": library.library_version,
        "entry_count": library.entry_count,
        "entries": [_entry_row(e) for e in library.entries],
    }


@router.post("/generate")
def generate_library_endpoint(request: Request) -> Dict[str, Any]:
    """Generate and validate a custom human N-glycan composition library."""
    payload = request.json_sync()
    try:
        include_neu5gc = bool(payload.get("include_neu5gc", False))
        constraints = CompositionConstraints(
            min_hex=int(payload.get("min_hex", 3)),
            max_hex=int(payload.get("max_hex", 12)),
            min_hexnac=int(payload.get("min_hexnac", 2)),
            max_hexnac=int(payload.get("max_hexnac", 8)),
            min_dhex=int(payload.get("min_dhex", 0)),
            max_dhex=int(payload.get("max_dhex", 4)),
            min_neu5ac=int(payload.get("min_neu5ac", 0)),
            max_neu5ac=int(payload.get("max_neu5ac", 4)),
            include_neu5gc=include_neu5gc,
            max_neu5gc=int(payload.get("max_neu5gc", 2 if include_neu5gc else 0)),
            labeling=str(payload.get("labeling", "NONE")).strip().upper(),
            sialic_linkage_expansion=bool(payload.get("expand_sialic_linkages", False)),
        )
        library = generate_human_n_glycan_library(
            constraints=constraints,
            library_name=str(payload.get("library_name", f"Custom_Human_NGlycans_{constraints.labeling}")),
            library_version=str(payload.get("library_version", "1.0.0")),
        )
        val_res = validate_glycan_library(library)
        if not val_res.is_valid:
            raise HTTPException(
                status_code=400,
                detail=f"Generated library failed validation: {'; '.join(val_res.errors)}",
            )
        session_state.active_library = library
        return {
            "library_name": library.library_name,
            "library_version": library.library_version,
            "entry_count": library.entry_count,
            "validation_valid": val_res.is_valid,
            "entries": [_entry_row(e) for e in library.entries],
        }
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Library generation error: {exc}") from exc


@router.post("/query")
def query_library_endpoint(request: Request) -> Dict[str, Any]:
    """Filter and search the active N-glycan library offline."""
    payload = request.json_sync()
    library = session_state.get_or_load_library()
    try:
        g_class_str = payload.get("glycan_class")
        classes: Optional[List[GlycanClass]] = None
        if g_class_str:
            classes = [GlycanClass(str(g_class_str).strip().lower())]

        tol_unit_str = str(payload.get("tolerance_unit", "ppm")).strip().lower()
        tol_unit = MassToleranceUnit.PPM if tol_unit_str == "ppm" else MassToleranceUnit.DA

        search_filter = LibrarySearchFilter(
            composition_query=payload.get("composition_query"),
            min_neutral_mass=float(payload["min_neutral_mass"]) if payload.get("min_neutral_mass") is not None else None,
            max_neutral_mass=float(payload["max_neutral_mass"]) if payload.get("max_neutral_mass") is not None else None,
            target_mass=float(payload["target_mass"]) if payload.get("target_mass") is not None else None,
            target_mass_is_mz=bool(payload.get("target_mass_is_mz", False)),
            adduct_code=str(payload.get("adduct_code", "[M+Na]+")),
            tolerance=float(payload.get("tolerance", 20.0)),
            tolerance_unit=tol_unit,
            glycan_classes=classes,
            labeling=payload.get("labeling"),
            fucosylated_only=payload.get("fucosylated_only"),
            sialylated_only=payload.get("sialylated_only"),
            exclude_neu5gc=bool(payload.get("exclude_neu5gc", False)),
        )
        matches = query_glycan_library(library=library, search_filter=search_filter)
        return {
            "library_name": library.library_name,
            "total_in_library": library.entry_count,
            "matched_count": len(matches),
            "entries": [_entry_row(m.entry) for m in matches],
        }
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Library query error: {exc}") from exc


@router.get("/export")
def export_library_endpoint(request: Request) -> Response:
    """Export the active N-glycan library in JSON, CSV, or TSV format."""
    fmt = request.query_params.get("format", "json").strip().lower()
    library = session_state.get_or_load_library()
    if fmt == "json":
        text_out = json.dumps(library.to_dict(), indent=2, ensure_ascii=False)
        return Response(
            content=text_out.encode("utf-8"),
            media_type="application/json",
            headers={"Content-Disposition": f'attachment; filename="{library.library_name}.json"'},
        )
    elif fmt in ("csv", "tsv"):
        delim = "," if fmt == "csv" else "\t"
        text_out = export_library_to_delimited(library, delimiter=delim)
        mime = "text/csv; charset=utf-8" if fmt == "csv" else "text/tab-separated-values; charset=utf-8"
        return Response(
            content=text_out.encode("utf-8"),
            media_type=mime,
            headers={"Content-Disposition": f'attachment; filename="{library.library_name}.{fmt}"'},
        )
    raise HTTPException(status_code=400, detail=f"Unsupported library export format: '{fmt}'")


@router.post("/import")
def import_library_endpoint(request: Request) -> Dict[str, Any]:
    """Import and validate a custom N-glycan library from JSON, CSV, or TSV."""
    form = request.form_sync()
    upload = form.get("file")
    if upload is None or not hasattr(upload, "filename"):
        raise HTTPException(status_code=400, detail="No library file provided under form field 'file'.")

    ext = Path(upload.filename).suffix.lower()
    if ext not in settings.allowed_library_extensions:
        raise HTTPException(
            status_code=400,
            detail=f"Unsupported library file extension '{ext}'. Allowed: {settings.allowed_library_extensions}",
        )

    raw_text = upload.content.decode("utf-8-sig", errors="replace")
    try:
        if ext == ".json":
            library = import_library_from_json(raw_text)
        else:
            delim = "\t" if ext == ".tsv" else ","
            library = import_library_from_delimited(
                raw_text,
                delimiter=delim,
                library_name=Path(upload.filename).stem,
            )
        val_res = validate_glycan_library(library)
        if not val_res.is_valid:
            raise HTTPException(
                status_code=400,
                detail=f"Imported library failed validation: {'; '.join(val_res.errors)}",
            )
        session_state.active_library = library
        return {
            "library_name": library.library_name,
            "library_version": library.library_version,
            "entry_count": library.entry_count,
            "validation_valid": val_res.is_valid,
            "entries": [_entry_row(e) for e in library.entries],
        }
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Failed to import library: {exc}") from exc

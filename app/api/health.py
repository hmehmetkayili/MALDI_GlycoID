"""Health and platform status API endpoint for MALDI-GlycoID (Stage 12)."""

from __future__ import annotations

from typing import Any, Dict

from app import __version__ as SOFTWARE_VERSION
from app.api.compat import APIRouter
from app.api.state import session_state
from app.config import settings

router = APIRouter(prefix="/api", tags=["health"])


@router.get("/health")
def health_check() -> Dict[str, Any]:
    """Return service health, version, and scientific configuration status."""
    library = session_state.get_or_load_library()
    return {
        "status": "ok",
        "app_name": settings.app_name,
        "version": SOFTWARE_VERSION,
        "active_library_name": library.library_name,
        "active_library_version": library.library_version,
        "active_library_entries_count": library.entry_count,
        "max_upload_size_bytes": settings.max_upload_size_bytes,
        "allowed_spectrum_extensions": list(settings.allowed_spectrum_extensions),
    }

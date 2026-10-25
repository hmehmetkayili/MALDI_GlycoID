"""Main web application entrypoint and router registration for MALDI-GlycoID (Stage 12)."""

from __future__ import annotations

from pathlib import Path

from app import __version__ as SOFTWARE_VERSION
from app.api import batch, exports, glycopeptides, health, libraries, qc, search, spectra
from app.api.compat import FastAPI, HTMLResponse, run_http_server
from app.config import BASE_DIR, settings

TEMPLATES_DIR = BASE_DIR / "app" / "templates"
STATIC_DIR = BASE_DIR / "app" / "static"


def create_app() -> FastAPI:
    """Create and configure the MALDI-GlycoID web application instance."""
    application = FastAPI(
        title=settings.app_name,
        version=SOFTWARE_VERSION,
        description=(
            "Scientifically transparent MALDI-TOF mass spectrometry platform for "
            "human N-glycan profiling, relative quantification, and glycopeptide screening."
        ),
    )

    application.mount_static("/static", STATIC_DIR)

    @application.get("/")
    def index_page() -> HTMLResponse:
        index_path = TEMPLATES_DIR / "index.html"
        html_content = index_path.read_text(encoding="utf-8")
        return HTMLResponse(content=html_content, status_code=200)

    application.include_router(health.router)
    application.include_router(libraries.router)
    application.include_router(spectra.router)
    application.include_router(search.router)
    application.include_router(glycopeptides.router)
    application.include_router(batch.router)
    application.include_router(qc.router)
    application.include_router(exports.router)

    return application


app = create_app()


if __name__ == "__main__":
    print(f"Starting {settings.app_name} v{SOFTWARE_VERSION} on http://{settings.host}:{settings.port}")
    run_http_server(app, host=settings.host, port=settings.port)

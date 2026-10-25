"""Global configuration settings and constants for MALDI-GlycoID."""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Tuple

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
LIBRARIES_DIR = DATA_DIR / "libraries"
EXAMPLES_DIR = DATA_DIR / "examples"
SCHEMAS_DIR = DATA_DIR / "schemas"
TEMP_UPLOAD_DIR = Path(os.getenv("MALDI_GLYCOID_TMP_DIR", str( Path(tempfile.gettempdir()) / "maldi_glycoid_uploads" )))

# Ensure matplotlib has a writable config cache directory in cloud/container environments (e.g. Render)
if "MPLCONFIGDIR" not in os.environ:
    _mpl_dir = Path(tempfile.gettempdir()) / "maldi_glycoid_mpl"
    try:
        _mpl_dir.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(_mpl_dir)
    except OSError:
        pass


def _default_host() -> str:
    if "MALDI_GLYCOID_HOST" in os.environ:
        return os.environ["MALDI_GLYCOID_HOST"]
    if "RENDER" in os.environ or "PORT" in os.environ:
        return "0.0.0.0"
    return "127.0.0.1"


def _default_port() -> int:
    raw = os.getenv("PORT", os.getenv("MALDI_GLYCOID_PORT", "8000"))
    try:
        return int(raw)
    except ValueError:
        return 8000


@dataclass(frozen=True)
class AppConfig:
    """Centralized configuration for MALDI-GlycoID."""

    app_name: str = "MALDI-GlycoID"
    version: str = "1.0.0"
    host: str = field(default_factory=_default_host)
    port: int = field(default_factory=_default_port)
    debug: bool = field(default_factory=lambda: os.getenv("MALDI_GLYCOID_DEBUG", "false").lower() == "true")

    # Security and upload limits
    max_upload_size_bytes: int = field(
        default_factory=lambda: int(os.getenv("MALDI_GLYCOID_MAX_UPLOAD_BYTES", str(50 * 1024 * 1024)))
    )
    max_batch_files: int = field(
        default_factory=lambda: int(os.getenv("MALDI_GLYCOID_MAX_BATCH_FILES", "100"))
    )
    allowed_spectrum_extensions: Tuple[str, ...] = (
        ".txt",
        ".xy",
        ".csv",
        ".tsv",
        ".mzml",
        ".mzxml",
    )
    allowed_library_extensions: Tuple[str, ...] = (
        ".csv",
        ".tsv",
        ".json",
    )
    allowed_fasta_extensions: Tuple[str, ...] = (
        ".fasta",
        ".fa",
        ".faa",
        ".txt",
    )

    # Default scientific search presets
    default_ppm_tolerance: float = 20.0
    preset_ppm_tolerances: Tuple[float, ...] = (5.0, 10.0, 20.0, 50.0, 100.0)
    default_da_tolerance: float = 0.05
    max_charge_state: int = 5
    min_charge_state: int = 1


settings = AppConfig()

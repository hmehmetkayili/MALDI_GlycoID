"""Package the complete MALDI-GlycoID v1.0.0 repository into a clean release ZIP archive (Stage 13)."""

from __future__ import annotations

from pathlib import Path
import zipfile

ROOT_DIR = Path(__file__).resolve().parent.parent
OUTPUT_ZIP = ROOT_DIR / "MALDI-GlycoID-v1.0.0.zip"

EXCLUDED_DIRS = {
    "__pycache__",
    ".git",
    ".venv",
    "venv",
    "env",
    "tmp_uploads",
    ".pytest_cache",
    ".mypy_cache",
    ".ruff_cache",
}
EXCLUDED_SUFFIXES = {".pyc", ".pyo", ".pyd"}


def build_release_zip(output_path: Path = OUTPUT_ZIP) -> Path:
    """Create `MALDI-GlycoID-v1.0.0.zip` containing all source, data, tests, docs, and manifests."""
    files_added = 0
    with zipfile.ZipFile(output_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for item in sorted(ROOT_DIR.rglob("*")):
            if item.is_dir():
                continue
            rel_parts = item.relative_to(ROOT_DIR).parts
            if any(part in EXCLUDED_DIRS for part in rel_parts):
                continue
            if item.suffix.lower() in EXCLUDED_SUFFIXES:
                continue
            if item.resolve() == output_path.resolve():
                continue
            arcname = f"MALDI-GlycoID-v1.0.0/{item.relative_to(ROOT_DIR).as_posix()}"
            zf.write(item, arcname=arcname)
            files_added += 1

    print(
        f"[OK] Packaged {files_added} files into {output_path.name} "
        f"({output_path.stat().st_size:,} bytes)."
    )
    return output_path


if __name__ == "__main__":
    build_release_zip()

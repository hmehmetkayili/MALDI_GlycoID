"""CLI utility to generate and validate MALDI-GlycoID N-glycan libraries offline."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from app.config import LIBRARIES_DIR
from app.libraries import (
    CompositionConstraints,
    export_library_to_delimited,
    export_library_to_json,
    generate_human_n_glycan_library,
    import_library_from_delimited,
    import_library_from_json,
    validate_glycan_library,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate or validate a MALDI-GlycoID N-glycan library.")
    parser.add_argument(
        "--file",
        type=str,
        default=None,
        help="Path to an existing JSON/CSV/TSV library to validate. If omitted, generates and saves the default human N-glycan library.",
    )
    args = parser.parse_args()

    if args.file:
        p = Path(args.file)
        if not p.is_file():
            print(f"[ERROR] File not found: {p}", file=sys.stderr)
            return 1
        if p.suffix.lower() == ".json":
            lib = import_library_from_json(p)
        else:
            lib = import_library_from_delimited(p)
    else:
        constraints = CompositionConstraints(
            min_hex=3,
            max_hex=10,
            min_hexnac=2,
            max_hexnac=7,
            min_dhex=0,
            max_dhex=3,
            min_neu5ac=0,
            max_neu5ac=4,
            include_neu5gc=False,
            labeling="NONE",
        )
        lib = generate_human_n_glycan_library(
            constraints=constraints,
            library_name="Default Human N-Glycan Biosynthetic Composition Library (Native)",
            library_version="1.0.0",
        )
        LIBRARIES_DIR.mkdir(parents=True, exist_ok=True)
        export_library_to_json(lib, output_path=LIBRARIES_DIR / "human_nglycan_native_v1.json")
        export_library_to_delimited(lib, delimiter=",", output_path=LIBRARIES_DIR / "human_nglycan_native_v1.csv")

    report = validate_glycan_library(lib)
    print(f"Library: {lib.library_name} (v{lib.library_version})")
    print(f"Entries: {report.entry_count}")
    print(f"Valid:   {report.is_valid}")
    if report.errors:
        print("Errors:")
        for err in report.errors:
            print(f"  - {err}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

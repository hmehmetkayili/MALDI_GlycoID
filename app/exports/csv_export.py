"""CSV and TSV table builders and exporters for MALDI-GlycoID (Stage 10).

Provides standardized English column schemas and UTF-8 / UTF-8-BOM delimited
exporters for all 12 analytical result categories:
1. `annotated_peaks`
2. `all_glycan_candidates`
3. `quantified_glycans`
4. `replicate_abundance_matrix`
5. `derived_glycan_traits`
6. `glycopeptide_candidates`
7. `calibration_diagnostics`
8. `quality_control_summary`
9. `statistical_analysis`
10. `glycan_reference_library`
11. `analysis_parameters_and_reproducibility`
12. `processing_and_error_log`
"""

from __future__ import annotations

import csv
import io
from typing import Any, Dict, List, Mapping, Optional, Sequence


def export_rows_to_delimited_bytes(
    rows: Sequence[Mapping[str, Any]],
    delimiter: str = ",",
    fieldnames: Optional[Sequence[str]] = None,
    include_utf8_bom: bool = False,
) -> bytes:
    """Serialize a list of dictionaries to CSV (`,`) or TSV (`\\t`) UTF-8 bytes."""
    if fieldnames is None:
        seen_keys: List[str] = []
        seen_set: set[str] = set()
        for row in rows:
            for k in row.keys():
                if k not in seen_set:
                    seen_set.add(k)
                    seen_keys.append(k)
        headers = seen_keys
    else:
        headers = list(fieldnames)

    buf = io.StringIO()
    if headers:
        writer = csv.DictWriter(
            buf,
            fieldnames=headers,
            delimiter=delimiter,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            clean_row: Dict[str, Any] = {}
            for k in headers:
                val = row.get(k)
                if val is None:
                    clean_row[k] = ""
                elif isinstance(val, (list, tuple)):
                    clean_row[k] = "; ".join(str(x) for x in val)
                elif isinstance(val, dict):
                    clean_row[k] = "; ".join(f"{dk}={dv}" for dk, dv in val.items())
                else:
                    clean_row[k] = val
            writer.writerow(clean_row)

    encoding = "utf-8-sig" if include_utf8_bom else "utf-8"
    return buf.getvalue().encode(encoding)


def export_rows_to_csv_bytes(
    rows: Sequence[Mapping[str, Any]],
    fieldnames: Optional[Sequence[str]] = None,
    include_utf8_bom: bool = False,
) -> bytes:
    """Serialize table rows to standard comma-separated values (CSV) bytes."""
    return export_rows_to_delimited_bytes(
        rows=rows,
        delimiter=",",
        fieldnames=fieldnames,
        include_utf8_bom=include_utf8_bom,
    )


def export_rows_to_tsv_bytes(
    rows: Sequence[Mapping[str, Any]],
    fieldnames: Optional[Sequence[str]] = None,
) -> bytes:
    """Serialize table rows to standard tab-separated values (TSV) bytes."""
    return export_rows_to_delimited_bytes(
        rows=rows,
        delimiter="\t",
        fieldnames=fieldnames,
        include_utf8_bom=False,
    )

"""Structured JSON exporter for MALDI-GlycoID (Stage 10).

Serializes all 12 standardized export categories, reproducibility parameters,
and scientific interpretation disclaimers into a deterministic UTF-8 JSON payload.
"""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
from enum import Enum
import json
from typing import Any, Dict, Mapping


def _json_default_serializer(obj: Any) -> Any:
    """Convert dataclasses, Enums, tuples, and NumPy scalars to JSON-serializable primitives."""
    if isinstance(obj, Enum):
        return obj.value
    if is_dataclass(obj) and not isinstance(obj, type):
        return asdict(obj)
    if hasattr(obj, "tolist"):
        return obj.tolist()
    if hasattr(obj, "item"):
        return obj.item()
    return str(obj)


def export_analysis_to_json_bytes(
    standardized_tables: Mapping[str, Any],
    reproducibility_metadata: Optional[Mapping[str, Any]] = None,
    pretty: bool = True,
) -> bytes:
    """Serialize the complete set of standardized analytical tables and metadata to UTF-8 JSON bytes."""
    payload: Dict[str, Any] = {
        "schema_name": "MALDI-GlycoID Analytical Export Bundle",
        "schema_version": "1.0.0",
        "exported_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "scientific_disclaimers": [
            "MS1 mass matches identify monosaccharide compositions and adduct candidates; structural isomers (branching, core vs. antennary fucosylation, and glycosidic linkages without linkage-specific derivatization) cannot be distinguished from MS1 mass alone.",
            "Glycan relative intensities and relative abundances are compositional semi-quantitative metrics and do not represent absolute molar concentrations without spiked internal standards and validated response factors.",
            "Heuristic candidate ranking scores combine mass accuracy, relative intensity, and optional isotope similarity (idotp) and are not calibrated false-discovery-rate probabilities.",
        ],
        "reproducibility_metadata": dict(reproducibility_metadata) if reproducibility_metadata else {},
        "tables": dict(standardized_tables),
    }

    indent = 2 if pretty else None
    text = json.dumps(payload, default=_json_default_serializer, indent=indent, ensure_ascii=False)
    return text.encode("utf-8")

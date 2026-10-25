"""Import and export of N-glycan libraries in JSON, CSV, and TSV formats.

Preserves library versioning, provenance, license metadata, and structural annotation
status across round-trip imports and exports. Never implies that generated theoretical
compositions are experimentally verified reference standards.
"""

from __future__ import annotations

import csv
import io
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

from app.core.ion_calculator import calculate_mz_from_neutral_mass
from app.core.isotope_scoring import generate_theoretical_isotope_envelope
from app.core.mass_calculator import (
    calculate_released_glycan_mass,
    format_glycan_composition,
    parse_glycan_composition,
)
from app.core.scientific_constants import ValidationStatus
from app.libraries.composition_generator import (
    CompositionConstraints,
    GlycanClass,
    GlycanLibrary,
    GlycanLibraryEntry,
    StructuralAnnotationStatus,
    classify_n_glycan_composition,
)
from app.parsers.base_models import decode_text_bytes

DEFAULT_EXPORT_ADDUCTS: Tuple[str, ...] = (
    "[M+H]+",
    "[M+Na]+",
    "[M+K]+",
    "[M+2H]2+",
    "[M+2Na]2+",
    "[M-H]-",
)


def _neg_mz(entry: GlycanLibraryEntry, adduct_code: str) -> float:
    """Negative-mode m/z of an entry: precomputed value if present, else computed from the neutral mass."""
    mz = entry.adduct_mzs.get(adduct_code)
    if mz is not None and mz > 0:
        return float(mz)
    return float(calculate_mz_from_neutral_mass(entry.theoretical_neutral_mass, adduct_code).theoretical_mz)


def export_library_to_json(
    library: GlycanLibrary,
    output_path: Optional[Path] = None,
    indent: int = 2,
) -> str:
    """Serialize a `GlycanLibrary` to a JSON string and optionally write to `output_path`."""
    payload = {
        "library_name": library.library_name,
        "library_version": library.library_version,
        "created_at_utc": library.created_at_utc,
        "labeling": library.labeling,
        "entry_count": library.entry_count,
        "provenance_summary": library.provenance_summary,
        "constraints_used": library.constraints_used.to_dict() if library.constraints_used else None,
        "warnings": list(library.warnings),
        "entries": [e.to_dict() for e in library.entries],
    }
    json_str = json.dumps(payload, indent=indent)
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json_str, encoding="utf-8")
    return json_str


def export_library_to_delimited(
    library: GlycanLibrary,
    delimiter: str = ",",
    output_path: Optional[Path] = None,
) -> str:
    """Serialize a `GlycanLibrary` to CSV (`delimiter=','`) or TSV (`delimiter='\\t'`)."""
    buf = io.StringIO()
    fieldnames = [
        "library_id",
        "composition",
        "Hex",
        "HexNAc",
        "dHex",
        "Neu5Ac",
        "Neu5Gc",
        "glycan_class",
        "antennary_hypothesis",
        "is_fucosylated",
        "is_sialylated",
        "contains_neu5gc",
        "theoretical_neutral_mass",
        "hill_formula",
        "mz_M_plus_H_1pos",
        "mz_M_plus_Na_1pos",
        "mz_M_plus_K_1pos",
        "mz_M_plus_2H_2pos",
        "mz_M_plus_2Na_2pos",
        "mz_M_minus_H_1neg",
        "mz_M_plus_Cl_1neg",
        "mz_M_minus_2H_2neg",
        "labeling",
        "default_charge",
        "structural_annotation_status",
        "isomer_group_note",
        "source_provenance",
        "reference_id",
        "license_info",
        "validation_status",
        "library_version",
        "alpha23_sialic_count",
        "alpha26_sialic_count",
    ]
    writer = csv.DictWriter(buf, fieldnames=fieldnames, delimiter=delimiter, lineterminator="\n")
    writer.writeheader()

    for e in library.entries:
        writer.writerow(
            {
                "library_id": e.library_id,
                "composition": e.composition,
                "Hex": e.residue_counts.get("Hex", 0),
                "HexNAc": e.residue_counts.get("HexNAc", 0),
                "dHex": e.residue_counts.get("dHex", 0),
                "Neu5Ac": e.residue_counts.get("Neu5Ac", 0),
                "Neu5Gc": e.residue_counts.get("Neu5Gc", 0),
                "glycan_class": e.glycan_class.value,
                "antennary_hypothesis": e.antennary_hypothesis,
                "is_fucosylated": e.is_fucosylated,
                "is_sialylated": e.is_sialylated,
                "contains_neu5gc": e.contains_neu5gc,
                "theoretical_neutral_mass": f"{e.theoretical_neutral_mass:.8f}",
                "hill_formula": e.hill_formula,
                "mz_M_plus_H_1pos": f"{e.adduct_mzs.get('[M+H]+', 0.0):.8f}",
                "mz_M_plus_Na_1pos": f"{e.adduct_mzs.get('[M+Na]+', 0.0):.8f}",
                "mz_M_plus_K_1pos": f"{e.adduct_mzs.get('[M+K]+', 0.0):.8f}",
                "mz_M_plus_2H_2pos": f"{e.adduct_mzs.get('[M+2H]2+', 0.0):.8f}",
                "mz_M_plus_2Na_2pos": f"{e.adduct_mzs.get('[M+2Na]2+', 0.0):.8f}",
                "mz_M_minus_H_1neg": f"{_neg_mz(e, '[M-H]-'):.8f}",
                "mz_M_plus_Cl_1neg": f"{_neg_mz(e, '[M+Cl]-'):.8f}",
                "mz_M_minus_2H_2neg": f"{_neg_mz(e, '[M-2H]2-'):.8f}",
                "labeling": e.labeling,
                "default_charge": e.default_charge,
                "structural_annotation_status": e.structural_annotation_status.value,
                "isomer_group_note": e.isomer_group_note,
                "source_provenance": e.source_provenance,
                "reference_id": e.reference_id or "",
                "license_info": e.license_info,
                "validation_status": e.validation_status.value,
                "library_version": e.library_version,
                "alpha23_sialic_count": "" if e.alpha23_sialic_count is None else e.alpha23_sialic_count,
                "alpha26_sialic_count": "" if e.alpha26_sialic_count is None else e.alpha26_sialic_count,
            }
        )

    out_str = buf.getvalue()
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(out_str, encoding="utf-8")
    return out_str


def import_library_from_json(content: Union[str, bytes, Path]) -> GlycanLibrary:
    """Load a `GlycanLibrary` from JSON content or file path."""
    if isinstance(content, Path):
        text = content.read_text(encoding="utf-8")
    elif isinstance(content, bytes):
        text, _ = decode_text_bytes(content)
    else:
        text = content

    data = json.loads(text)
    if not isinstance(data, dict) or "entries" not in data:
        raise ValueError("Invalid library JSON format: root object must contain 'entries'")

    lib_version = str(data.get("library_version", "1.0.0"))
    lib_name = str(data.get("library_name", "Imported Glycan Library"))
    labeling_default = str(data.get("labeling", "NONE")).upper()

    constraints_dict = data.get("constraints_used")
    constraints_obj: Optional[CompositionConstraints] = None
    if isinstance(constraints_dict, dict):
        allowed_cls = constraints_dict.get("allowed_classes")
        constraints_obj = CompositionConstraints(
            min_hex=int(constraints_dict.get("min_hex", 3)),
            max_hex=int(constraints_dict.get("max_hex", 12)),
            min_hexnac=int(constraints_dict.get("min_hexnac", 2)),
            max_hexnac=int(constraints_dict.get("max_hexnac", 8)),
            min_dhex=int(constraints_dict.get("min_dhex", 0)),
            max_dhex=int(constraints_dict.get("max_dhex", 4)),
            min_neu5ac=int(constraints_dict.get("min_neu5ac", 0)),
            max_neu5ac=int(constraints_dict.get("max_neu5ac", 4)),
            include_neu5gc=bool(constraints_dict.get("include_neu5gc", False)),
            max_neu5gc=int(constraints_dict.get("max_neu5gc", 0)),
            include_paucimannose=bool(constraints_dict.get("include_paucimannose", False)),
            enforce_biosynthetic_rules=bool(constraints_dict.get("enforce_biosynthetic_rules", True)),
            labeling=str(constraints_dict.get("labeling", labeling_default)),
            sialic_linkage_expansion=bool(constraints_dict.get("sialic_linkage_expansion", False)),
            adduct_codes=tuple(constraints_dict.get("adduct_codes", DEFAULT_EXPORT_ADDUCTS)),
            allowed_classes=tuple(GlycanClass(c) for c in allowed_cls) if allowed_cls else None,
        )

    entries: List[GlycanLibraryEntry] = []
    for idx, row in enumerate(data["entries"], start=1):
        res_counts = row.get("residue_counts")
        if not res_counts:
            comp_raw = str(row.get("composition", "")).split("[")[0]
            res_counts = parse_glycan_composition(comp_raw)
        else:
            res_counts = {str(k): int(v) for k, v in res_counts.items() if int(v) > 0}

        label_code = str(row.get("labeling", labeling_default)).upper().strip()
        a23 = int(row["alpha23_sialic_count"]) if row.get("alpha23_sialic_count") not in (None, "") else None
        a26 = int(row["alpha26_sialic_count"]) if row.get("alpha26_sialic_count") not in (None, "") else None

        # Always compute authoritative Stage 1 mass & formula to verify or populate missing fields
        mass_res = calculate_released_glycan_mass(
            composition=res_counts,
            labeling=label_code,
            alpha23_count=a23,
            alpha26_count=a26,
        )
        neutral_mass = float(row.get("theoretical_neutral_mass", mass_res.neutral_monoisotopic_mass))

        stored_adducts = row.get("adduct_mzs")
        if isinstance(stored_adducts, dict) and stored_adducts:
            adduct_mzs = {str(k): float(v) for k, v in stored_adducts.items()}
        else:
            adduct_mzs = {}
            for code in DEFAULT_EXPORT_ADDUCTS:
                adduct_mzs[code] = calculate_mz_from_neutral_mass(
                    neutral_mass=neutral_mass,
                    adduct_code=code,
                ).theoretical_mz

        na_ion = calculate_mz_from_neutral_mass(
            neutral_mass=neutral_mass,
            adduct_code="[M+Na]+",
            neutral_formula=mass_res.elemental_formula,
        )
        assert na_ion.ion_formula is not None
        iso_env = generate_theoretical_isotope_envelope(
            ion_formula=na_ion.ion_formula,
            monoisotopic_mz=na_ion.theoretical_mz,
            charge=1,
            num_peaks=3,
        )

        _, auto_class, auto_hyp, auto_note, _ = classify_n_glycan_composition(res_counts)
        g_class_str = row.get("glycan_class", auto_class.value)
        g_class = GlycanClass(g_class_str) if g_class_str in GlycanClass._value2member_map_ else auto_class

        ann_str = row.get("structural_annotation_status", StructuralAnnotationStatus.COMPOSITION_CANDIDATE.value)
        ann_status = (
            StructuralAnnotationStatus(ann_str)
            if ann_str in StructuralAnnotationStatus._value2member_map_
            else StructuralAnnotationStatus.COMPOSITION_CANDIDATE
        )

        val_str = row.get("validation_status", ValidationStatus.VALIDATED_ENABLED.value)
        val_status = (
            ValidationStatus(val_str)
            if val_str in ValidationStatus._value2member_map_
            else ValidationStatus.VALIDATED_ENABLED
        )

        entries.append(
            GlycanLibraryEntry(
                library_id=str(row.get("library_id", f"MGID_IMP_{idx:04d}")),
                composition=str(row.get("composition", format_glycan_composition(res_counts))),
                residue_counts=res_counts,
                glycan_class=g_class,
                antennary_hypothesis=str(row.get("antennary_hypothesis", auto_hyp)),
                is_fucosylated=res_counts.get("dHex", 0) > 0,
                is_sialylated=(res_counts.get("Neu5Ac", 0) + res_counts.get("Neu5Gc", 0)) > 0,
                contains_neu5gc=res_counts.get("Neu5Gc", 0) > 0,
                theoretical_neutral_mass=neutral_mass,
                elemental_formula=mass_res.elemental_formula,
                hill_formula=str(row.get("hill_formula", mass_res.hill_formula)),
                adduct_mzs=adduct_mzs,
                labeling=label_code,
                default_charge=int(row.get("default_charge", 1)),
                expected_isotope_mzs_na=tuple(
                    row.get("expected_isotope_mzs_na", [p.mz for p in iso_env])
                ),
                expected_isotope_abundances_na=tuple(
                    row.get("expected_isotope_abundances_na", [p.relative_abundance for p in iso_env])
                ),
                structural_annotation_status=ann_status,
                isomer_group_note=str(row.get("isomer_group_note", auto_note)),
                source_provenance=str(row.get("source_provenance", "User-Imported Reference Library")),
                reference_id=str(row["reference_id"]) if row.get("reference_id") else None,
                license_info=str(row.get("license_info", "Unspecified User Import License")),
                validation_status=val_status,
                library_version=str(row.get("library_version", lib_version)),
                alpha23_sialic_count=a23,
                alpha26_sialic_count=a26,
                warnings=tuple(row.get("warnings", [])),
            )
        )

    return GlycanLibrary(
        library_name=lib_name,
        library_version=lib_version,
        created_at_utc=str(data.get("created_at_utc", datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"))),
        labeling=labeling_default,
        constraints_used=constraints_obj,
        entries=tuple(entries),
        provenance_summary=str(data.get("provenance_summary", "Imported from JSON library file.")),
        warnings=tuple(data.get("warnings", [])),
    )


def import_library_from_delimited(
    content: Union[str, bytes, Path],
    delimiter: Optional[str] = None,
    library_name: str = "Imported Delimited Glycan Library",
    library_version: str = "1.0.0",
) -> GlycanLibrary:
    """Import a `GlycanLibrary` from a CSV or TSV file."""
    if isinstance(content, Path):
        library_name = content.stem
        text, _ = decode_text_bytes(content.read_bytes())
        if delimiter is None:
            delimiter = "\t" if content.suffix.lower() == ".tsv" else ","
    elif isinstance(content, bytes):
        text, _ = decode_text_bytes(content)
    else:
        text = content

    if delimiter is None:
        first_line = text.splitlines()[0] if text.splitlines() else ""
        delimiter = "\t" if "\t" in first_line else ","

    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    entries: List[GlycanLibraryEntry] = []

    for idx, row in enumerate(reader, start=1):
        # Extract residue counts either from explicit columns (Hex, HexNAc, ...) or 'composition' column
        counts: Dict[str, int] = {}
        for col_sym in ("Hex", "HexNAc", "dHex", "Neu5Ac", "Neu5Gc", "Pent", "HexA"):
            if col_sym in row and str(row[col_sym]).strip().isdigit():
                val = int(str(row[col_sym]).strip())
                if val > 0:
                    counts[col_sym] = val

        comp_raw = str(row.get("composition", "")).strip()
        if not counts and comp_raw:
            base_comp_str = comp_raw.split("[")[0]
            counts = parse_glycan_composition(base_comp_str)

        if not counts:
            raise ValueError(f"Row {idx} in imported library lacks valid monosaccharide composition or counts.")

        label_code = str(row.get("labeling", "NONE") or "NONE").upper().strip()
        a23_raw = str(row.get("alpha23_sialic_count", "")).strip()
        a26_raw = str(row.get("alpha26_sialic_count", "")).strip()
        a23 = int(a23_raw) if a23_raw.isdigit() else None
        a26 = int(a26_raw) if a26_raw.isdigit() else None

        mass_res = calculate_released_glycan_mass(
            composition=counts,
            labeling=label_code,
            alpha23_count=a23,
            alpha26_count=a26,
        )

        if row.get("theoretical_neutral_mass") and str(row["theoretical_neutral_mass"]).strip():
            neutral_mass = float(str(row["theoretical_neutral_mass"]).strip())
        else:
            neutral_mass = mass_res.neutral_monoisotopic_mass

        adduct_mzs: Dict[str, float] = {}
        for code in DEFAULT_EXPORT_ADDUCTS:
            adduct_mzs[code] = calculate_mz_from_neutral_mass(
                neutral_mass=neutral_mass,
                adduct_code=code,
            ).theoretical_mz

        na_ion = calculate_mz_from_neutral_mass(
            neutral_mass=neutral_mass,
            adduct_code="[M+Na]+",
            neutral_formula=mass_res.elemental_formula,
        )
        assert na_ion.ion_formula is not None
        iso_env = generate_theoretical_isotope_envelope(
            ion_formula=na_ion.ion_formula,
            monoisotopic_mz=na_ion.theoretical_mz,
            charge=1,
            num_peaks=3,
        )

        _, auto_class, auto_hyp, auto_note, auto_warns = classify_n_glycan_composition(counts)
        g_class_str = str(row.get("glycan_class", "")).strip()
        g_class = GlycanClass(g_class_str) if g_class_str in GlycanClass._value2member_map_ else auto_class

        ann_str = str(row.get("structural_annotation_status", "")).strip()
        ann_status = (
            StructuralAnnotationStatus(ann_str)
            if ann_str in StructuralAnnotationStatus._value2member_map_
            else StructuralAnnotationStatus.COMPOSITION_CANDIDATE
        )

        ref_id = str(row.get("reference_id", "")).strip() or None
        source_prov = str(row.get("source_provenance", "")).strip() or "User-Imported Delimited Reference Library"
        lic_info = str(row.get("license_info", "")).strip() or "User-Supplied Reference Entry"

        entries.append(
            GlycanLibraryEntry(
                library_id=str(row.get("library_id", "")).strip() or f"MGID_CSV_{idx:04d}",
                composition=comp_raw or format_glycan_composition(counts),
                residue_counts=counts,
                glycan_class=g_class,
                antennary_hypothesis=str(row.get("antennary_hypothesis", "")).strip() or auto_hyp,
                is_fucosylated=counts.get("dHex", 0) > 0,
                is_sialylated=(counts.get("Neu5Ac", 0) + counts.get("Neu5Gc", 0)) > 0,
                contains_neu5gc=counts.get("Neu5Gc", 0) > 0,
                theoretical_neutral_mass=neutral_mass,
                elemental_formula=mass_res.elemental_formula,
                hill_formula=str(row.get("hill_formula", "")).strip() or mass_res.hill_formula,
                adduct_mzs=adduct_mzs,
                labeling=label_code,
                default_charge=int(str(row.get("default_charge", "1")).strip() or "1"),
                expected_isotope_mzs_na=tuple(p.mz for p in iso_env),
                expected_isotope_abundances_na=tuple(p.relative_abundance for p in iso_env),
                structural_annotation_status=ann_status,
                isomer_group_note=str(row.get("isomer_group_note", "")).strip() or auto_note,
                source_provenance=source_prov,
                reference_id=ref_id,
                license_info=lic_info,
                validation_status=ValidationStatus.VALIDATED_ENABLED,
                library_version=str(row.get("library_version", "")).strip() or library_version,
                alpha23_sialic_count=a23,
                alpha26_sialic_count=a26,
                warnings=tuple(auto_warns),
            )
        )

    if not entries:
        raise ValueError("Imported delimited library file contains no valid glycan entries.")

    return GlycanLibrary(
        library_name=library_name,
        library_version=library_version,
        created_at_utc=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        labeling=entries[0].labeling,
        constraints_used=None,
        entries=tuple(entries),
        provenance_summary=f"Imported {len(entries)} entries from delimited library source.",
        warnings=(
            "Imported entries retain user-supplied source and license provenance; composition-level entries do not imply unique isomer identification.",
        ),
    )

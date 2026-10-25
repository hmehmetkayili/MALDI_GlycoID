"""Scientific parameter and biosynthetic validation utilities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Mapping, Tuple

from app.core.mass_error import MassToleranceUnit
from app.core.scientific_constants import (
    ADDUCT_DEFINITIONS,
    IonPolarity,
    LABELING_REGISTRY,
    MONOSACCHARIDE_RESIDUES,
)


@dataclass(frozen=True)
class ValidationReport:
    """Structured validation outcome with errors and scientific warnings."""

    is_valid: bool
    errors: Tuple[str, ...]
    warnings: Tuple[str, ...]


def validate_search_parameters(
    tolerance: float,
    tolerance_unit: MassToleranceUnit,
    polarity: IonPolarity,
    adduct_codes: List[str],
    charge_states: List[int],
    labeling: str,
) -> ValidationReport:
    """Validate user-supplied MS1 search parameters for chemical and analytical consistency."""
    errors: List[str] = []
    warnings: List[str] = []

    if tolerance <= 0:
        errors.append(f"Mass tolerance must be strictly positive (got {tolerance}).")
    elif tolerance_unit == MassToleranceUnit.PPM and tolerance > 500:
        warnings.append(
            f"Large mass tolerance ({tolerance} ppm) increases false-positive and isobaric candidate overlap."
        )
    elif tolerance_unit == MassToleranceUnit.DA and tolerance > 1.0:
        warnings.append(
            f"Large mass tolerance ({tolerance} Da) may cause severe adduct and composition ambiguity."
        )

    if not adduct_codes:
        errors.append("At least one ion adduct must be selected.")
    else:
        for code in adduct_codes:
            if code not in ADDUCT_DEFINITIONS:
                errors.append(f"Unsupported adduct code: '{code}'.")
            else:
                adduct = ADDUCT_DEFINITIONS[code]
                if adduct.polarity != polarity:
                    errors.append(
                        f"Adduct '{code}' has polarity '{adduct.polarity.value}', which conflicts with selected polarity '{polarity.value}'."
                    )

    for z in charge_states:
        if z < 1 or z > 5:
            errors.append(f"Charge state magnitude must be between 1 and 5 (got {z}).")

    label_clean = labeling.upper().strip()
    if label_clean not in LABELING_REGISTRY:
        errors.append(f"Unsupported labeling/derivatization method: '{labeling}'.")

    return ValidationReport(
        is_valid=len(errors) == 0,
        errors=tuple(errors),
        warnings=tuple(warnings),
    )


def validate_human_n_glycan_composition(composition: Mapping[str, int]) -> ValidationReport:
    """Check whether a monosaccharide composition satisfies canonical N-glycan rules.

    Checks:
    1. All residue symbols are recognized and non-negative.
    2. Contains the chitobiose / trimannosyl N-glycan core (Hex >= 3, HexNAc >= 2) unless flagged as truncated/paucimannose.
    3. Sialic acids (Neu5Ac + Neu5Gc) do not exceed available non-reducing Hex/HexNAc antennae.
    4. Warns if Neu5Gc is present (non-human dietary/trace or non-human mammalian origin).
    """
    errors: List[str] = []
    warnings: List[str] = []

    if not composition:
        errors.append("Composition is empty.")
        return ValidationReport(False, tuple(errors), tuple(warnings))

    for sym, cnt in composition.items():
        if sym not in MONOSACCHARIDE_RESIDUES:
            errors.append(f"Unknown residue symbol: '{sym}'.")
        elif cnt < 0:
            errors.append(f"Negative count ({cnt}) for residue '{sym}'.")

    if errors:
        return ValidationReport(False, tuple(errors), tuple(warnings))

    hex_cnt = composition.get("Hex", 0)
    hexnac_cnt = composition.get("HexNAc", 0)
    fuc_cnt = composition.get("dHex", 0)
    neu5ac_cnt = composition.get("Neu5Ac", 0)
    neu5gc_cnt = composition.get("Neu5Gc", 0)
    total_sia = neu5ac_cnt + neu5gc_cnt

    if hexnac_cnt < 2 or hex_cnt < 3:
        warnings.append(
            f"Composition (Hex={hex_cnt}, HexNAc={hexnac_cnt}) lacks the full canonical N-glycan trimannosyl core (Hex3HexNAc2); may represent a paucimannose or truncated species."
        )

    if total_sia > 0 and hexnac_cnt <= 2:
        warnings.append(
            "Sialylated N-glycan with HexNAc <= 2 is biosynthetically unusual (lacks GlcNAc antennae)."
        )

    # In standard N-glycans, each antenna requires at least one GlcNAc beyond the 2 core GlcNAc residues
    max_antennae = max(0, hexnac_cnt - 2)
    if total_sia > max_antennae and max_antennae > 0:
        warnings.append(
            f"Total sialic acid count ({total_sia}) exceeds estimated GlcNAc antennae ({max_antennae}); requires polysialylation or di-sialylated motifs."
        )

    if fuc_cnt > (hexnac_cnt - 1) and hexnac_cnt >= 2:
        warnings.append(
            f"Deoxyhexose/Fucose count ({fuc_cnt}) is high relative to HexNAc count ({hexnac_cnt})."
        )

    if neu5gc_cnt > 0:
        warnings.append(
            "Contains Neu5Gc: humans lack functional CMAH enzyme; Neu5Gc in human samples typically reflects dietary incorporation, biotherapeutic origin, or isobaric oxidations."
        )

    return ValidationReport(
        is_valid=len(errors) == 0,
        errors=tuple(errors),
        warnings=tuple(warnings),
    )

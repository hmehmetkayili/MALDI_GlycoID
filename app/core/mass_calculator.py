"""Glycan, derivatization, and peptide neutral mass calculation engine.

Enforces strict separation between:
1. Released native reducing glycans (residue sum + terminal H2O).
2. Released glycans with reducing-end labels (residue sum + H2O + net label reaction delta).
3. Permethylated glycans (site-specific CH2 methylation across residues + reducing end).
4. Linkage-specific sialic-acid derivatized glycans (alpha2,3 lactonization vs alpha2,6 esterification/amidation).
5. Glycan moieties attached to peptides (dehydratively linked residue sum WITHOUT terminal H2O).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Tuple, Union

from app.core.formula_calculator import (
    add_formulas,
    calculate_formula_monoisotopic_mass,
    format_hill_formula,
)
from app.core.scientific_constants import (
    AMINO_ACID_RESIDUES,
    ELEMENT_MONOISOTOPIC_MASSES,
    GlycanForm,
    LABELING_REGISTRY,
    METHYLENE_FORMULA,
    METHYLENE_MASS_DA,
    MONOSACCHARIDE_RESIDUES,
    PEPTIDE_MODIFICATIONS,
    RESIDUE_ALIAS_MAP,
    SIALIC_DMA_ALPHA23_SHIFT_DA,
    SIALIC_DMA_ALPHA26_SHIFT_DA,
    SIALIC_EEC_ALPHA23_SHIFT_DA,
    SIALIC_EEC_ALPHA26_SHIFT_DA,
    ValidationStatus,
    WATER_FORMULA,
    WATER_MASS_DA,
)

# Canonical display order for monosaccharide compositions
CANONICAL_RESIDUE_ORDER: Tuple[str, ...] = (
    "Hex",
    "HexNAc",
    "dHex",
    "Neu5Ac",
    "Neu5Gc",
    "Pent",
    "HexA",
)

# Regex patterns for flexible composition parsing:
# Supports "Hex5HexNAc4dHex1Neu5Ac2", "Hex(5)HexNAc(4)Fuc(1)NeuAc(2)", "H5N4F1S2", etc.
_COMP_PAREN_RE = re.compile(r"([A-Za-z0-9_]+)\s*\(\s*(\d+)\s*\)")
_COMP_INLINE_RE = re.compile(
    r"(HexNAc|Neu5Ac|Neu5Gc|NeuAc|NeuGc|GlcNAc|GalNAc|dHex|HexA|Pent|Fuc|Hex|Man|Gal|Glc|Xyl|Sia|H|N|F|S|G|X)(\d+)"
)


def parse_glycan_composition(composition: Union[str, Mapping[str, int]]) -> Dict[str, int]:
    """Parse a glycan composition string or dict into canonical monosaccharide counts.

    Examples of supported string formats:
      - 'Hex5HexNAc4dHex1Neu5Ac2'
      - 'Hex(5)HexNAc(4)Fuc(1)Neu5Ac(2)'
      - 'H5N4F1S2'
      - 'Hex:5, HexNAc:4, Fuc:1'

    Raises ValueError if any residue is unrecognized or has a negative count.
    """
    if isinstance(composition, Mapping):
        normalized: Dict[str, int] = {}
        for raw_key, count in composition.items():
            key_clean = str(raw_key).strip()
            if key_clean not in RESIDUE_ALIAS_MAP:
                raise ValueError(f"Unknown monosaccharide symbol: '{raw_key}'")
            if not isinstance(count, int) or isinstance(count, bool) or count < 0:
                raise ValueError(f"Monosaccharide count for '{raw_key}' must be a non-negative integer, got {count!r}")
            if count > 0:
                canon = RESIDUE_ALIAS_MAP[key_clean]
                normalized[canon] = normalized.get(canon, 0) + count
        return normalized

    if not isinstance(composition, str):
        raise ValueError(f"Composition must be a string or mapping, got {type(composition).__name__}")

    text = composition.strip()
    if not text:
        raise ValueError("Glycan composition string cannot be empty")

    # Case 1: Key:Value separated by commas/spaces (e.g., "Hex:5, HexNAc:4")
    if ":" in text:
        counts: Dict[str, int] = {}
        for part in re.split(r"[,;]+", text):
            part = part.strip()
            if not part:
                continue
            if ":" not in part:
                raise ValueError(f"Invalid key:value pair '{part}' in composition '{composition}'")
            k, v = part.split(":", 1)
            k = k.strip()
            v_str = v.strip()
            if k not in RESIDUE_ALIAS_MAP:
                raise ValueError(f"Unknown monosaccharide symbol '{k}' in '{composition}'")
            if not v_str.isdigit():
                raise ValueError(f"Invalid residue count '{v_str}' for '{k}' in '{composition}'")
            cnt = int(v_str)
            if cnt > 0:
                canon = RESIDUE_ALIAS_MAP[k]
                counts[canon] = counts.get(canon, 0) + cnt
        if not counts:
            raise ValueError(f"No positive monosaccharide counts found in '{composition}'")
        return counts

    # Case 2: Parenthesized format e.g. Hex(5)HexNAc(4)
    if "(" in text:
        counts = {}
        cleaned = re.sub(r"[\s_\-]+", "", text)
        pos = 0
        for m in _COMP_PAREN_RE.finditer(cleaned):
            if m.start() != pos:
                raise ValueError(f"Unrecognized syntax in composition '{composition}' at '{cleaned[pos:m.start()]}'")
            pos = m.end()
            sym, cnt_str = m.group(1), m.group(2)
            if sym not in RESIDUE_ALIAS_MAP:
                raise ValueError(f"Unknown monosaccharide symbol '{sym}' in '{composition}'")
            cnt = int(cnt_str)
            if cnt > 0:
                canon = RESIDUE_ALIAS_MAP[sym]
                counts[canon] = counts.get(canon, 0) + cnt
        if pos != len(cleaned) or not counts:
            raise ValueError(f"Could not parse composition '{composition}'")
        return counts

    # Case 3: Inline format e.g. Hex5HexNAc4dHex1Neu5Ac2 or H5N4F1S2
    cleaned = re.sub(r"[\s_\-]+", "", text)
    pos = 0
    counts = {}
    for m in _COMP_INLINE_RE.finditer(cleaned):
        if m.start() != pos:
            raise ValueError(f"Unrecognized residue token in composition '{composition}' near '{cleaned[pos:]}'")
        pos = m.end()
        sym, cnt_str = m.group(1), m.group(2)
        canon = RESIDUE_ALIAS_MAP[sym]
        cnt = int(cnt_str)
        if cnt > 0:
            counts[canon] = counts.get(canon, 0) + cnt
    if pos != len(cleaned) or not counts:
        raise ValueError(f"Invalid or empty glycan composition string: '{composition}'")
    return counts


def format_glycan_composition(
    composition: Mapping[str, int],
    use_fuc_alias: bool = False,
) -> str:
    """Format canonical monosaccharide dictionary into a deterministic composition string."""
    parts: List[str] = []
    seen = set()
    for sym in CANONICAL_RESIDUE_ORDER:
        cnt = composition.get(sym, 0)
        if cnt > 0:
            out_sym = "Fuc" if (sym == "dHex" and use_fuc_alias) else sym
            parts.append(f"{out_sym}{cnt}")
            seen.add(sym)
    for sym in sorted(composition.keys()):
        if sym not in seen and composition[sym] > 0:
            parts.append(f"{sym}{composition[sym]}")
    return "".join(parts)


@dataclass(frozen=True)
class GlycanMassResult:
    """Complete, traceable result of a glycan neutral mass calculation."""

    canonical_composition: str
    residue_counts: Mapping[str, int]
    glycan_form: GlycanForm
    labeling_code: str
    residue_sum_mass: float
    terminal_mass_contribution: float
    derivatization_mass_shift: float
    neutral_monoisotopic_mass: float
    elemental_formula: Mapping[str, int]
    hill_formula: str
    permethylation_sites_total: int = 0
    alpha23_sialic_count: Optional[int] = None
    alpha26_sialic_count: Optional[int] = None
    notes: Tuple[str, ...] = field(default_factory=tuple)


def calculate_residue_sum_mass_and_formula(
    composition: Union[str, Mapping[str, int]]
) -> Tuple[float, Dict[str, int], Dict[str, int]]:
    """Calculate the dehydratively linked residue sum mass and elemental formula (NO terminal H2O).

    Returns:
      (residue_sum_mass_da, elemental_formula_dict, canonical_residue_counts)
    """
    counts = parse_glycan_composition(composition)
    if not counts:
        raise ValueError("Glycan composition must contain at least one monosaccharide residue")

    total_mass = 0.0
    formula: Dict[str, int] = {}
    for sym, count in counts.items():
        res = MONOSACCHARIDE_RESIDUES[sym]
        total_mass += res.residue_mass * count
        formula = add_formulas(formula, res.residue_formula, scale=count)

    return total_mass, formula, counts


def calculate_attached_glycan_moiety_mass(
    composition: Union[str, Mapping[str, int]],
    sialic_derivatization: str = "NONE",
    alpha23_count: Optional[int] = None,
    alpha26_count: Optional[int] = None,
) -> GlycanMassResult:
    """Calculate the mass of an N-glycan moiety covalently attached to a peptide side-chain.

    CRITICAL SCIENTIFIC CONVENTION:
    When an N-glycan is attached to an asparagine residue in a glycopeptide, the glycosidic
    bond between GlcNAc-1 and the Asn amide nitrogen is formed with loss of H2O relative to
    the free reducing glycan and free peptide. Because the peptide neutral mass already
    includes all atoms of the unmodified Asn residue (-NH2), the attached glycan moiety mass
    is strictly the sum of the anhydro-monosaccharide residue masses (NO terminal reducing-end H2O).
    """
    residue_sum_mass, formula, counts = calculate_residue_sum_mass_and_formula(composition)
    sialic_code = sialic_derivatization.upper().strip()

    deriv_shift = 0.0
    notes: List[str] = [
        "Peptide-attached glycosyl moiety: calculated as sum of anhydro-residue masses without reducing-end H2O."
    ]

    if sialic_code != "NONE":
        if sialic_code not in ("EEC_SIALIC", "DMA_SIALIC"):
            raise ValueError(
                f"Labeling method '{sialic_code}' reduces or modifies a free reducing end and is not applicable to peptide-attached glycans."
            )
        total_sialic = counts.get("Neu5Ac", 0) + counts.get("Neu5Gc", 0)
        a23 = 0 if alpha23_count is None else alpha23_count
        a26 = total_sialic if alpha26_count is None and alpha23_count is None else (0 if alpha26_count is None else alpha26_count)
        if a23 < 0 or a26 < 0 or (a23 + a26) != total_sialic:
            raise ValueError(
                f"Linkage-specific sialic acid counts (alpha2,3={a23}, alpha2,6={a26}) must sum to total sialic acids ({total_sialic})."
            )
        if sialic_code == "EEC_SIALIC":
            deriv_shift = a26 * SIALIC_EEC_ALPHA26_SHIFT_DA + a23 * SIALIC_EEC_ALPHA23_SHIFT_DA
            formula = add_formulas(formula, {"C": 2, "H": 4}, scale=a26)
            formula = add_formulas(formula, WATER_FORMULA, scale=-a23)
        elif sialic_code == "DMA_SIALIC":
            deriv_shift = a26 * SIALIC_DMA_ALPHA26_SHIFT_DA + a23 * SIALIC_DMA_ALPHA23_SHIFT_DA
            formula = add_formulas(formula, {"C": 2, "H": 5, "N": 1}, scale=a26)
            formula = add_formulas(formula, {"O": 1}, scale=-a26)
            formula = add_formulas(formula, WATER_FORMULA, scale=-a23)
        notes.append(f"Applied linkage-specific sialic acid derivatization ({sialic_code}): alpha2,3={a23}, alpha2,6={a26}.")
    else:
        a23 = None
        a26 = None

    neutral_mass = residue_sum_mass + deriv_shift
    return GlycanMassResult(
        canonical_composition=format_glycan_composition(counts),
        residue_counts=counts,
        glycan_form=GlycanForm.PEPTIDE_ATTACHED,
        labeling_code=sialic_code,
        residue_sum_mass=residue_sum_mass,
        terminal_mass_contribution=0.0,
        derivatization_mass_shift=deriv_shift,
        neutral_monoisotopic_mass=neutral_mass,
        elemental_formula=formula,
        hill_formula=format_hill_formula(formula),
        permethylation_sites_total=0,
        alpha23_sialic_count=a23,
        alpha26_sialic_count=a26,
        notes=tuple(notes),
    )


def calculate_released_glycan_mass(
    composition: Union[str, Mapping[str, int]],
    labeling: str = "NONE",
    alpha23_count: Optional[int] = None,
    alpha26_count: Optional[int] = None,
    reducing_end_tag_with_sialic: str = "NONE",
) -> GlycanMassResult:
    """Calculate the neutral monoisotopic mass and elemental formula of a released glycan.

    Supports:
    - Native free reducing glycans (`labeling='NONE'`):
        neutral_mass = sum(residue_masses) + m(H2O)
    - Reduced alditols (`labeling='REDUCED'`):
        neutral_mass = sum(residue_masses) + m(H2O) + m(H2)
    - Reducing-end fluorescent/ionizable labels (`2-AB`, `2-AA`, `PROCAINAMIDE`, `PA`):
        neutral_mass = sum(residue_masses) + m(H2O) + net_label_shift
    - Permethylated non-reduced glycans (`PERMETHYLATED`):
        neutral_mass = sum(residue_masses) + m(H2O) + (sum(sites_i * count_i) + 2) * m(CH2)
    - Permethylated reduced alditols (`PERMETHYLATED_REDUCED`):
        neutral_mass = sum(residue_masses) + m(H2O) + m(H2) + (sum(sites_i * count_i) + 3) * m(CH2)
    - Linkage-specific sialic acid derivatization (`EEC_SIALIC`, `DMA_SIALIC`), optionally
      combined with a reducing-end label (`reducing_end_tag_with_sialic`).
    """
    label_code = labeling.upper().strip()
    if label_code not in LABELING_REGISTRY:
        raise ValueError(
            f"Unsupported labeling/derivatization method: '{labeling}'. "
            f"Supported methods: {', '.join(LABELING_REGISTRY.keys())}"
        )

    label_def = LABELING_REGISTRY[label_code]
    if label_def.validation_status == ValidationStatus.NOT_IMPLEMENTED:
        raise NotImplementedError(f"Labeling method '{label_code}' is marked NOT_IMPLEMENTED")

    residue_sum_mass, formula, counts = calculate_residue_sum_mass_and_formula(composition)

    # Every released glycan starts from the free reducing glycan (residue sum + terminal H2O)
    terminal_mass = WATER_MASS_DA
    formula = add_formulas(formula, WATER_FORMULA, scale=1)

    notes: List[str] = ["Released glycan includes terminal reducing-end H2O (+18.0105646865 Da)."]
    permethylation_sites = 0
    deriv_shift = 0.0
    glycan_form = GlycanForm.RELEASED_FREE

    if label_code == "NONE":
        glycan_form = GlycanForm.RELEASED_FREE

    elif label_def.category == "reducing_end_label":
        glycan_form = GlycanForm.RELEASED_LABELED
        assert label_def.net_mass_shift_from_free_reducing is not None
        deriv_shift = label_def.net_mass_shift_from_free_reducing
        if label_def.net_formula_added:
            formula = add_formulas(formula, label_def.net_formula_added, scale=1)
        if label_def.net_formula_removed:
            formula = add_formulas(formula, label_def.net_formula_removed, scale=-1)
        notes.append(f"Applied reducing-end label {label_def.name} ({deriv_shift:+.6f} Da relative to free reducing glycan).")

    elif label_def.category == "permethylation":
        glycan_form = GlycanForm.PERNETHYLATED
        residue_sites = sum(
            MONOSACCHARIDE_RESIDUES[sym].permethylation_sites * cnt for sym, cnt in counts.items()
        )
        if label_code == "PERMETHYLATED":
            # Native reducing end H2O contributes 2 free hydroxyls (1 anomeric -OH + 1 C-terminal -OH)
            terminal_sites = 2
            permethylation_sites = residue_sites + terminal_sites
            deriv_shift = permethylation_sites * METHYLENE_MASS_DA
            formula = add_formulas(formula, METHYLENE_FORMULA, scale=permethylation_sites)
            notes.append(
                f"Permethylated (non-reduced): {residue_sites} residue sites + 2 terminal sites = {permethylation_sites} CH2 groups (+{deriv_shift:.6f} Da)."
            )
        elif label_code == "PERMETHYLATED_REDUCED":
            # Reduction adds H2 (+2.015650 Da) and opens the ring, creating 3 terminal -OH sites instead of 2
            terminal_sites = 3
            permethylation_sites = residue_sites + terminal_sites
            reduction_mass = 2 * ELEMENT_MONOISOTOPIC_MASSES["H"]
            deriv_shift = reduction_mass + permethylation_sites * METHYLENE_MASS_DA
            formula = add_formulas(formula, {"H": 2}, scale=1)
            formula = add_formulas(formula, METHYLENE_FORMULA, scale=permethylation_sites)
            notes.append(
                f"Permethylated (reduced alditol): +H2 reduction + ({residue_sites} residue sites + 3 alditol terminal sites = {permethylation_sites} CH2 groups, +{deriv_shift:.6f} Da)."
            )

    elif label_def.category == "linkage_specific_sialic":
        glycan_form = GlycanForm.RELEASED_LABELED
        total_sialic = counts.get("Neu5Ac", 0) + counts.get("Neu5Gc", 0)
        a23 = 0 if alpha23_count is None else alpha23_count
        a26 = (
            total_sialic
            if (alpha26_count is None and alpha23_count is None)
            else (0 if alpha26_count is None else alpha26_count)
        )
        if a23 < 0 or a26 < 0 or (a23 + a26) != total_sialic:
            raise ValueError(
                f"Linkage-specific sialic acid counts (alpha2,3={a23}, alpha2,6={a26}) must sum to total sialic acids ({total_sialic})."
            )
        if label_code == "EEC_SIALIC":
            sialic_shift = a26 * SIALIC_EEC_ALPHA26_SHIFT_DA + a23 * SIALIC_EEC_ALPHA23_SHIFT_DA
            formula = add_formulas(formula, {"C": 2, "H": 4}, scale=a26)
            formula = add_formulas(formula, WATER_FORMULA, scale=-a23)
        else:  # DMA_SIALIC
            sialic_shift = a26 * SIALIC_DMA_ALPHA26_SHIFT_DA + a23 * SIALIC_DMA_ALPHA23_SHIFT_DA
            formula = add_formulas(formula, {"C": 2, "H": 5, "N": 1}, scale=a26)
            formula = add_formulas(formula, {"O": 1}, scale=-a26)
            formula = add_formulas(formula, WATER_FORMULA, scale=-a23)

        deriv_shift += sialic_shift
        notes.append(f"Linkage-specific sialic acid derivatization ({label_code}): alpha2,3={a23}, alpha2,6={a26} ({sialic_shift:+.6f} Da).")

        tag_code = reducing_end_tag_with_sialic.upper().strip()
        if tag_code != "NONE":
            if tag_code not in LABELING_REGISTRY or LABELING_REGISTRY[tag_code].category != "reducing_end_label":
                raise ValueError(f"Invalid reducing-end tag '{reducing_end_tag_with_sialic}' combined with {label_code}")
            tag_def = LABELING_REGISTRY[tag_code]
            assert tag_def.net_mass_shift_from_free_reducing is not None
            deriv_shift += tag_def.net_mass_shift_from_free_reducing
            if tag_def.net_formula_added:
                formula = add_formulas(formula, tag_def.net_formula_added, scale=1)
            if tag_def.net_formula_removed:
                formula = add_formulas(formula, tag_def.net_formula_removed, scale=-1)
            notes.append(f"Combined with reducing-end tag {tag_def.name} ({tag_def.net_mass_shift_from_free_reducing:+.6f} Da).")

        alpha23_count = a23
        alpha26_count = a26

    neutral_mass = residue_sum_mass + terminal_mass + deriv_shift

    # Cross-check against elemental formula monoisotopic mass
    formula_mass = calculate_formula_monoisotopic_mass(formula)
    if abs(neutral_mass - formula_mass) > 1e-7:
        raise RuntimeError(
            f"Internal mass consistency check failed: residue/shift sum ({neutral_mass:.8f}) != elemental formula mass ({formula_mass:.8f})"
        )

    return GlycanMassResult(
        canonical_composition=format_glycan_composition(counts),
        residue_counts=counts,
        glycan_form=glycan_form,
        labeling_code=label_code,
        residue_sum_mass=residue_sum_mass,
        terminal_mass_contribution=terminal_mass,
        derivatization_mass_shift=deriv_shift,
        neutral_monoisotopic_mass=neutral_mass,
        elemental_formula=formula,
        hill_formula=format_hill_formula(formula),
        permethylation_sites_total=permethylation_sites,
        alpha23_sialic_count=alpha23_count,
        alpha26_sialic_count=alpha26_count,
        notes=tuple(notes),
    )


@dataclass(frozen=True)
class PeptideMassResult:
    """Traceable neutral mass result for an unmodified or modified peptide sequence."""

    sequence: str
    length: int
    unmodified_neutral_mass: float
    modifications_applied: Tuple[str, ...]
    modification_mass_shift: float
    neutral_monoisotopic_mass: float
    elemental_formula: Mapping[str, int]
    hill_formula: str


def calculate_peptide_neutral_mass(
    sequence: str,
    modifications: Optional[Mapping[str, int]] = None,
    fixed_carbamidomethyl_cys: bool = False,
) -> PeptideMassResult:
    """Calculate the neutral monoisotopic mass of a peptide sequence.

    Convention:
      unmodified_mass = sum(amino_acid_residue_masses) + m(H2O) (for intact H- and -OH termini).
    """
    seq = sequence.strip().upper()
    if not seq:
        raise ValueError("Peptide sequence cannot be empty")

    formula: Dict[str, int] = dict(WATER_FORMULA)
    unmod_mass = WATER_MASS_DA

    for aa in seq:
        if aa not in AMINO_ACID_RESIDUES:
            raise ValueError(f"Unsupported or non-canonical amino acid symbol '{aa}' in sequence '{sequence}'")
        aa_def = AMINO_ACID_RESIDUES[aa]
        unmod_mass += aa_def.monoisotopic_mass
        formula = add_formulas(formula, aa_def.formula, scale=1)

    mod_shift = 0.0
    mods_applied: List[str] = []

    if fixed_carbamidomethyl_cys:
        cys_count = seq.count("C")
        if cys_count > 0:
            c_mod = PEPTIDE_MODIFICATIONS["Carbamidomethyl_C"]
            mod_shift += c_mod.mass_delta * cys_count
            formula = add_formulas(formula, c_mod.formula_delta, scale=cys_count)
            mods_applied.append(f"Carbamidomethyl_C x{cys_count}")

    if modifications:
        for mod_code, mod_count in modifications.items():
            if mod_count == 0:
                continue
            if mod_code not in PEPTIDE_MODIFICATIONS:
                raise ValueError(f"Unknown peptide modification code: '{mod_code}'")
            if mod_count < 0:
                raise ValueError(f"Modification count for '{mod_code}' cannot be negative")
            if mod_code == "Carbamidomethyl_C" and fixed_carbamidomethyl_cys:
                continue  # Already applied
            mod_def = PEPTIDE_MODIFICATIONS[mod_code]
            mod_shift += mod_def.mass_delta * mod_count
            formula = add_formulas(formula, mod_def.formula_delta, scale=mod_count)
            mods_applied.append(f"{mod_code} x{mod_count}")

    total_mass = unmod_mass + mod_shift
    return PeptideMassResult(
        sequence=seq,
        length=len(seq),
        unmodified_neutral_mass=unmod_mass,
        modifications_applied=tuple(mods_applied),
        modification_mass_shift=mod_shift,
        neutral_monoisotopic_mass=total_mass,
        elemental_formula=formula,
        hill_formula=format_hill_formula(formula),
    )


def calculate_glycopeptide_neutral_mass(
    peptide_sequence: str,
    glycan_composition: Union[str, Mapping[str, int]],
    peptide_modifications: Optional[Mapping[str, int]] = None,
    fixed_carbamidomethyl_cys: bool = False,
    sialic_derivatization: str = "NONE",
    alpha23_count: Optional[int] = None,
    alpha26_count: Optional[int] = None,
) -> Tuple[float, PeptideMassResult, GlycanMassResult, Dict[str, int]]:
    """Calculate the neutral mass of a glycopeptide using the peptide-attached glycan convention.

    Neutral glycopeptide mass = Peptide neutral mass (with terminal H2O)
                              + Attached glycan moiety mass (WITHOUT terminal H2O).
    """
    pep_res = calculate_peptide_neutral_mass(
        sequence=peptide_sequence,
        modifications=peptide_modifications,
        fixed_carbamidomethyl_cys=fixed_carbamidomethyl_cys,
    )
    gly_res = calculate_attached_glycan_moiety_mass(
        composition=glycan_composition,
        sialic_derivatization=sialic_derivatization,
        alpha23_count=alpha23_count,
        alpha26_count=alpha26_count,
    )
    combined_formula = add_formulas(pep_res.elemental_formula, gly_res.elemental_formula, scale=1)
    total_neutral_mass = pep_res.neutral_monoisotopic_mass + gly_res.neutral_monoisotopic_mass
    return total_neutral_mass, pep_res, gly_res, combined_formula

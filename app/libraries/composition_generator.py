"""Automated human N-glycan composition library generator.

Implements biosynthetically constrained human N-glycan composition generation
using the Stage 1 authoritative mass and ion calculators. Strictly distinguishes
composition-level candidates from topology hypotheses and structurally validated
reference entries.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from app.core.ion_calculator import calculate_mz_from_neutral_mass
from app.core.isotope_scoring import generate_theoretical_isotope_envelope
from app.core.mass_calculator import (
    calculate_released_glycan_mass,
    format_glycan_composition,
)
from app.core.scientific_constants import ValidationStatus


class GlycanClass(str, Enum):
    """Biosynthetic class of an N-glycan composition."""

    PAUCIMANNOSE = "paucimannose"
    HIGH_MANNOSE = "high_mannose"
    HYBRID = "hybrid"
    COMPLEX = "complex"
    UNCLASSIFIED = "unclassified"


class StructuralAnnotationStatus(str, Enum):
    """Distinguishes MS1 composition candidates from topological and reference structures."""

    COMPOSITION_CANDIDATE = "COMPOSITION_CANDIDATE"
    TOPOLOGY_HYPOTHESIS = "TOPOLOGY_HYPOTHESIS"
    LINKAGE_SPECIFIC_DERIVATIZED = "LINKAGE_SPECIFIC_DERIVATIZED"
    STRUCTURALLY_VALIDATED_REFERENCE = "STRUCTURALLY_VALIDATED_REFERENCE"


@dataclass(frozen=True)
class CompositionConstraints:
    """Configurable residue-count limits and biosynthetic rules for N-glycan generation."""

    min_hex: int = 3
    max_hex: int = 12
    min_hexnac: int = 2
    max_hexnac: int = 8
    min_dhex: int = 0
    max_dhex: int = 4
    min_neu5ac: int = 0
    max_neu5ac: int = 4
    include_neu5gc: bool = False
    max_neu5gc: int = 0
    include_paucimannose: bool = False
    enforce_biosynthetic_rules: bool = True
    max_antennae_hypothesis: int = 5
    allow_bisecting_glcnac: bool = True
    allow_polylacnac: bool = True
    labeling: str = "NONE"
    sialic_linkage_expansion: bool = False  # Expand alpha2,3 / alpha2,6 isomers for EEC_SIALIC / DMA_SIALIC
    adduct_codes: Tuple[str, ...] = ("[M+H]+", "[M+Na]+", "[M+K]+", "[M+2H]2+", "[M+2Na]2+", "[M-H]-")
    allowed_classes: Optional[Tuple[GlycanClass, ...]] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "min_hex": self.min_hex,
            "max_hex": self.max_hex,
            "min_hexnac": self.min_hexnac,
            "max_hexnac": self.max_hexnac,
            "min_dhex": self.min_dhex,
            "max_dhex": self.max_dhex,
            "min_neu5ac": self.min_neu5ac,
            "max_neu5ac": self.max_neu5ac,
            "include_neu5gc": self.include_neu5gc,
            "max_neu5gc": self.max_neu5gc,
            "include_paucimannose": self.include_paucimannose,
            "enforce_biosynthetic_rules": self.enforce_biosynthetic_rules,
            "max_antennae_hypothesis": self.max_antennae_hypothesis,
            "allow_bisecting_glcnac": self.allow_bisecting_glcnac,
            "allow_polylacnac": self.allow_polylacnac,
            "labeling": self.labeling,
            "sialic_linkage_expansion": self.sialic_linkage_expansion,
            "adduct_codes": list(self.adduct_codes),
            "allowed_classes": [c.value for c in self.allowed_classes] if self.allowed_classes else None,
        }


@dataclass(frozen=True)
class GlycanLibraryEntry:
    """Single record in a MALDI-GlycoID reference or theoretical composition library."""

    library_id: str
    composition: str
    residue_counts: Mapping[str, int]
    glycan_class: GlycanClass
    antennary_hypothesis: str  # e.g., "Core/Oligomannose", "Biantennary hypothesis", "Triantennary / Bisecting isomeric group"
    is_fucosylated: bool
    is_sialylated: bool
    contains_neu5gc: bool
    theoretical_neutral_mass: float
    elemental_formula: Mapping[str, int]
    hill_formula: str
    adduct_mzs: Mapping[str, float]  # e.g. {"[M+H]+": ..., "[M+Na]+": ...}
    labeling: str
    default_charge: int
    expected_isotope_mzs_na: Tuple[float, ...]  # Top 3 theoretical isotope m/z for [M+Na]+
    expected_isotope_abundances_na: Tuple[float, ...]  # Top 3 relative abundances for [M+Na]+
    structural_annotation_status: StructuralAnnotationStatus
    isomer_group_note: str
    source_provenance: str
    reference_id: Optional[str] = None
    license_info: str = "CC0-1.0 / Theoretical Calculation"
    validation_status: ValidationStatus = ValidationStatus.VALIDATED_ENABLED
    library_version: str = "1.0.0"
    alpha23_sialic_count: Optional[int] = None
    alpha26_sialic_count: Optional[int] = None
    warnings: Tuple[str, ...] = field(default_factory=tuple)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "library_id": self.library_id,
            "composition": self.composition,
            "residue_counts": dict(self.residue_counts),
            "glycan_class": self.glycan_class.value,
            "antennary_hypothesis": self.antennary_hypothesis,
            "is_fucosylated": self.is_fucosylated,
            "is_sialylated": self.is_sialylated,
            "contains_neu5gc": self.contains_neu5gc,
            "theoretical_neutral_mass": round(self.theoretical_neutral_mass, 8),
            "elemental_formula": dict(self.elemental_formula),
            "hill_formula": self.hill_formula,
            "adduct_mzs": {k: round(v, 8) for k, v in self.adduct_mzs.items()},
            "labeling": self.labeling,
            "default_charge": self.default_charge,
            "expected_isotope_mzs_na": [round(x, 6) for x in self.expected_isotope_mzs_na],
            "expected_isotope_abundances_na": [round(x, 4) for x in self.expected_isotope_abundances_na],
            "structural_annotation_status": self.structural_annotation_status.value,
            "isomer_group_note": self.isomer_group_note,
            "source_provenance": self.source_provenance,
            "reference_id": self.reference_id,
            "license_info": self.license_info,
            "validation_status": self.validation_status.value,
            "library_version": self.library_version,
            "alpha23_sialic_count": self.alpha23_sialic_count,
            "alpha26_sialic_count": self.alpha26_sialic_count,
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class GlycanLibrary:
    """Complete versioned glycan library with provenance and generation constraints."""

    library_name: str
    library_version: str
    created_at_utc: str
    labeling: str
    constraints_used: Optional[CompositionConstraints]
    entries: Tuple[GlycanLibraryEntry, ...]
    provenance_summary: str
    warnings: Tuple[str, ...] = field(default_factory=tuple)

    @property
    def entry_count(self) -> int:
        return len(self.entries)


def classify_n_glycan_composition(
    counts: Mapping[str, int],
    constraints: Optional[CompositionConstraints] = None,
) -> Tuple[bool, GlycanClass, str, str, List[str]]:
    """Evaluate biosynthetic validity, N-glycan class, antennary hypothesis, and isomer notes.

    Returns:
      (is_valid_under_rules, glycan_class, antennary_hypothesis, isomer_group_note, warnings)
    """
    h = counts.get("Hex", 0)
    n = counts.get("HexNAc", 0)
    f = counts.get("dHex", 0)
    s_ac = counts.get("Neu5Ac", 0)
    s_gc = counts.get("Neu5Gc", 0)
    s_tot = s_ac + s_gc

    enforce = constraints.enforce_biosynthetic_rules if constraints is not None else True
    inc_pauci = constraints.include_paucimannose if constraints is not None else False
    warnings: List[str] = []

    if s_gc > 0:
        warnings.append(
            "Contains Neu5Gc (non-human sialic acid; in human samples may reflect dietary trace, biotherapeutic origin, or oxidation artifact)."
        )

    # 1. Paucimannose / truncated core (Hex 1..3, HexNAc 2, Fuc 0..1, Sia == 0)
    if n == 2 and 1 <= h <= 3 and s_tot == 0 and f <= 1:
        if h < 3 and not inc_pauci and enforce:
            return False, GlycanClass.PAUCIMANNOSE, "Truncated core", "", warnings
        if h < 3:
            return (
                True,
                GlycanClass.PAUCIMANNOSE,
                "Paucimannose / truncated core",
                "Composition-level candidate; regioisomers of core mannose/fucose not resolved by MS1.",
                warnings,
            )
        # Hex3HexNAc2 (or +Fuc1) is the canonical trimannosyl core
        return (
            True,
            GlycanClass.PAUCIMANNOSE if h == 3 and f >= 0 else GlycanClass.HIGH_MANNOSE,
            "Trimannosyl N-glycan core",
            "Canonical N-glycan core composition (Man3GlcNAc2 +/- Fuc).",
            warnings,
        )

    # Enforce core requirement (Hex >= 3, HexNAc >= 2) when biosynthetic rules are active
    if enforce and (h < 3 or n < 2):
        return False, GlycanClass.UNCLASSIFIED, "Invalid core", "", warnings

    # 2. High-mannose N-glycans: HexNAc == 2, Hex >= 4, Sia == 0, Fuc <= 1
    if n == 2:
        if enforce and (s_tot > 0 or f > 1 or h < 4):
            return False, GlycanClass.HIGH_MANNOSE, "Invalid high-mannose", "", warnings
        return (
            True,
            GlycanClass.HIGH_MANNOSE,
            f"Oligomannose (Man{h}GlcNAc2)",
            "Isomeric mannose branching topologies (e.g., D1/D2/D3 arm isomers) share identical MS1 mass.",
            warnings,
        )

    # From here on, HexNAc >= 3 (n_ant = n - 2 extra GlcNAc residues beyond the chitobiose core)
    n_extra_hexnac = n - 2
    n_extra_hex = h - 3

    # 3. Hybrid N-glycans: typically HexNAc == 3 (monoantennary or bisecting on oligomannose) with Hex >= 4,
    # or HexNAc == 4 with high mannose content (Hex >= 6 and Sia <= 1)
    if n == 3:
        if enforce:
            # With 1 extra GlcNAc, at most 1 antenna can be galactosylated and sialylated (Sia <= 1),
            # and Fuc can be core + at most 1 antennary fucose (Fuc <= 2), and Hex >= 3.
            if s_tot > 1 or f > 2 or h < 3:
                return False, GlycanClass.HYBRID, "Invalid hybrid/monoantennary", "", warnings
            # If sialylated (s_tot == 1), requires at least 1 galactose on the antenna (h >= 4)
            if s_tot == 1 and h < 4:
                return False, GlycanClass.HYBRID, "Sialylated without galactose", "", warnings
        g_class = GlycanClass.HYBRID if h >= 5 else GlycanClass.COMPLEX
        hyp = "Monoantennary / Hybrid hypothesis" if h >= 5 else "Monoantennary / Bisected core hypothesis"
        note = (
            "MS1 composition cannot distinguish hybrid structures from monoantennary or bisected isomers."
        )
        return True, g_class, hyp, note, warnings

    # 4. Complex N-glycans (HexNAc >= 4)
    if enforce:
        # Galactose / Hexose bounds:
        # In complex glycans (HexNAc >= 4), extra Hex (h - 3) represents antennary Gal (or hybrid Man / alpha-Gal).
        # Under standard human complex rules:
        # - Every sialic acid requires a Gal residue: s_tot <= max(0, h - 3)
        # - Also s_tot cannot exceed the number of antennary GlcNAc branches: s_tot <= n_extra_hexnac
        # - For HexNAc >= 4, extra Hex (h - 3) normally does not exceed n_extra_hexnac (up to 1 Gal per GlcNAc,
        #   allowing +1 for hybrid/polylactosamine if enabled, but bounded so arbitrary combinations are rejected).
        if s_tot > max(0, n_extra_hex) or s_tot > n_extra_hexnac:
            return False, GlycanClass.COMPLEX, "Sialic acid exceeds antennae/galactose", "", warnings

        max_allowed_extra_hex = n_extra_hexnac + (1 if (constraints and constraints.allow_polylacnac) else 0)
        if n_extra_hex > max_allowed_extra_hex:
            return False, GlycanClass.COMPLEX, "Hexose count exceeds complex antennary capacity", "", warnings

        # Fucose bound: at most 1 core fucose + 1 fucose per antennary GlcNAc (and <= total HexNAc - 1)
        if f > min(n - 1, 1 + n_extra_hexnac):
            return False, GlycanClass.COMPLEX, "Excess deoxyhexose/fucose for antennary count", "", warnings

    # Formulate antennary hypothesis label (explicitly noting bisecting / polylactosamine ambiguity)
    if n == 4:
        hyp = "Biantennary (or monoantennary + bisecting GlcNAc) composition hypothesis"
    elif n == 5:
        hyp = "Triantennary / Biantennary + bisecting GlcNAc isomeric composition group"
    elif n == 6:
        hyp = "Tetraantennary / Triantennary + bisecting / poly-LacNAc isomeric composition group"
    else:
        hyp = f"Higher-antennary / poly-LacNAc composition group (HexNAc={n})"

    note = (
        "Composition-level candidate only; branching topology, bisecting GlcNAc vs additional antennae, "
        "core vs antennary fucosylation, and alpha2,3/alpha2,6 sialic acid linkages (unless derivatized) "
        "are unresolved by MS1."
    )
    return True, GlycanClass.COMPLEX, hyp, note, warnings


def _build_entry_from_composition(
    entry_index: int,
    counts: Mapping[str, int],
    glycan_class: GlycanClass,
    antennary_hypothesis: str,
    isomer_note: str,
    entry_warnings: Sequence[str],
    constraints: CompositionConstraints,
    library_version: str,
    alpha23_count: Optional[int] = None,
    alpha26_count: Optional[int] = None,
) -> GlycanLibraryEntry:
    """Calculate exact masses, adduct m/z table, and isotope envelope for a single library entry."""
    mass_res = calculate_released_glycan_mass(
        composition=counts,
        labeling=constraints.labeling,
        alpha23_count=alpha23_count,
        alpha26_count=alpha26_count,
    )

    adduct_mzs: Dict[str, float] = {}
    for adduct_code in constraints.adduct_codes:
        ion_res = calculate_mz_from_neutral_mass(
            neutral_mass=mass_res.neutral_monoisotopic_mass,
            adduct_code=adduct_code,
            neutral_formula=mass_res.elemental_formula,
        )
        adduct_mzs[adduct_code] = ion_res.theoretical_mz

    # Compute expected [M+Na]+ 3-peak isotope envelope for reference
    na_ion = calculate_mz_from_neutral_mass(
        neutral_mass=mass_res.neutral_monoisotopic_mass,
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

    comp_str = format_glycan_composition(counts)
    if alpha23_count is not None and alpha26_count is not None:
        comp_str = f"{comp_str}[a2,3:{alpha23_count}/a2,6:{alpha26_count}]"
        ann_status = StructuralAnnotationStatus.LINKAGE_SPECIFIC_DERIVATIZED
    else:
        ann_status = StructuralAnnotationStatus.COMPOSITION_CANDIDATE

    lib_id = f"MGID_HNG_{entry_index:04d}"

    return GlycanLibraryEntry(
        library_id=lib_id,
        composition=comp_str,
        residue_counts=dict(counts),
        glycan_class=glycan_class,
        antennary_hypothesis=antennary_hypothesis,
        is_fucosylated=counts.get("dHex", 0) > 0,
        is_sialylated=(counts.get("Neu5Ac", 0) + counts.get("Neu5Gc", 0)) > 0,
        contains_neu5gc=counts.get("Neu5Gc", 0) > 0,
        theoretical_neutral_mass=mass_res.neutral_monoisotopic_mass,
        elemental_formula=mass_res.elemental_formula,
        hill_formula=mass_res.hill_formula,
        adduct_mzs=adduct_mzs,
        labeling=constraints.labeling.upper().strip(),
        default_charge=1,
        expected_isotope_mzs_na=tuple(p.mz for p in iso_env),
        expected_isotope_abundances_na=tuple(p.relative_abundance for p in iso_env),
        structural_annotation_status=ann_status,
        isomer_group_note=isomer_note,
        source_provenance="MALDI-GlycoID Biosynthetic Human N-Glycan Composition Generator (Theoretical)",
        reference_id=None,
        license_info="CC0-1.0 (Theoretical Biosynthetic Composition)",
        validation_status=ValidationStatus.VALIDATED_ENABLED,
        library_version=library_version,
        alpha23_sialic_count=alpha23_count,
        alpha26_sialic_count=alpha26_count,
        warnings=tuple(entry_warnings),
    )


def generate_human_n_glycan_library(
    constraints: Optional[CompositionConstraints] = None,
    library_name: str = "Human N-Glycan Biosynthetic Composition Library",
    library_version: str = "1.0.0",
) -> GlycanLibrary:
    """Generate a reproducible, biosynthetically constrained human N-glycan composition library offline."""
    if constraints is None:
        constraints = CompositionConstraints()

    entries: List[GlycanLibraryEntry] = []
    idx = 1

    max_gc = constraints.max_neu5gc if constraints.include_neu5gc else 0
    label_upper = constraints.labeling.upper().strip()
    is_linkage_sialic = label_upper in ("EEC_SIALIC", "DMA_SIALIC") and constraints.sialic_linkage_expansion

    for n in range(constraints.min_hexnac, constraints.max_hexnac + 1):
        for h in range(constraints.min_hex, constraints.max_hex + 1):
            for f in range(constraints.min_dhex, constraints.max_dhex + 1):
                for s_ac in range(constraints.min_neu5ac, constraints.max_neu5ac + 1):
                    for s_gc in range(0, max_gc + 1):
                        counts: Dict[str, int] = {"Hex": h, "HexNAc": n}
                        if f > 0:
                            counts["dHex"] = f
                        if s_ac > 0:
                            counts["Neu5Ac"] = s_ac
                        if s_gc > 0:
                            counts["Neu5Gc"] = s_gc

                        is_valid, g_class, hyp, note, warns = classify_n_glycan_composition(
                            counts, constraints=constraints
                        )
                        if not is_valid:
                            continue
                        if constraints.allowed_classes and g_class not in constraints.allowed_classes:
                            continue

                        total_sia = s_ac + s_gc
                        if is_linkage_sialic and total_sia > 0:
                            # Expand into alpha2,3 / alpha2,6 linkage combinations
                            for a23 in range(total_sia + 1):
                                a26 = total_sia - a23
                                entries.append(
                                    _build_entry_from_composition(
                                        entry_index=idx,
                                        counts=counts,
                                        glycan_class=g_class,
                                        antennary_hypothesis=hyp,
                                        isomer_note=note,
                                        entry_warnings=warns,
                                        constraints=constraints,
                                        library_version=library_version,
                                        alpha23_count=a23,
                                        alpha26_count=a26,
                                    )
                                )
                                idx += 1
                        else:
                            entries.append(
                                _build_entry_from_composition(
                                    entry_index=idx,
                                    counts=counts,
                                    glycan_class=g_class,
                                    antennary_hypothesis=hyp,
                                    isomer_note=note,
                                    entry_warnings=warns,
                                    constraints=constraints,
                                    library_version=library_version,
                                )
                            )
                            idx += 1

    # Sort deterministically by theoretical neutral mass ascending, then composition string
    entries.sort(key=lambda e: (round(e.theoretical_neutral_mass, 8), e.composition))

    # Re-index sequentially after sorting so IDs are strictly ordered by mass
    reindexed: List[GlycanLibraryEntry] = []
    for new_i, e in enumerate(entries, start=1):
        reindexed.append(
            GlycanLibraryEntry(
                library_id=f"MGID_HNG_{new_i:04d}",
                composition=e.composition,
                residue_counts=e.residue_counts,
                glycan_class=e.glycan_class,
                antennary_hypothesis=e.antennary_hypothesis,
                is_fucosylated=e.is_fucosylated,
                is_sialylated=e.is_sialylated,
                contains_neu5gc=e.contains_neu5gc,
                theoretical_neutral_mass=e.theoretical_neutral_mass,
                elemental_formula=e.elemental_formula,
                hill_formula=e.hill_formula,
                adduct_mzs=e.adduct_mzs,
                labeling=e.labeling,
                default_charge=e.default_charge,
                expected_isotope_mzs_na=e.expected_isotope_mzs_na,
                expected_isotope_abundances_na=e.expected_isotope_abundances_na,
                structural_annotation_status=e.structural_annotation_status,
                isomer_group_note=e.isomer_group_note,
                source_provenance=e.source_provenance,
                reference_id=e.reference_id,
                license_info=e.license_info,
                validation_status=e.validation_status,
                library_version=e.library_version,
                alpha23_sialic_count=e.alpha23_sialic_count,
                alpha26_sialic_count=e.alpha26_sialic_count,
                warnings=e.warnings,
            )
        )

    lib_warnings: List[str] = [
        "Generated theoretical compositions are composition-level hypotheses and must not be represented as experimentally verified reference standards unless matched to validated reference standards."
    ]
    if constraints.include_neu5gc:
        lib_warnings.append(
            "Library includes Neu5Gc-containing compositions; Neu5Gc is non-endogenous in humans."
        )

    return GlycanLibrary(
        library_name=library_name,
        library_version=library_version,
        created_at_utc=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        labeling=label_upper,
        constraints_used=constraints,
        entries=tuple(reindexed),
        provenance_summary=(
            "Generated offline by MALDI-GlycoID biosynthetic human N-glycan composition engine "
            f"(Hex {constraints.min_hex}-{constraints.max_hex}, HexNAc {constraints.min_hexnac}-{constraints.max_hexnac}, "
            f"dHex {constraints.min_dhex}-{constraints.max_dhex}, Neu5Ac {constraints.min_neu5ac}-{constraints.max_neu5ac}, "
            f"labeling={label_upper})."
        ),
        warnings=tuple(lib_warnings),
    )

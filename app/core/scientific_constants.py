"""Authoritative scientific constants for MALDI-GlycoID.

All elemental masses, physical constants, monosaccharide residue masses,
adduct definitions, derivatization reactions, and amino acid residue masses
are defined in this single authoritative module.

No downstream module may define competing mass constants.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Dict, Mapping, Optional, Tuple


# ============================================================================
# 1. FUNDAMENTAL PHYSICAL & ELEMENTAL CONSTANTS (IUPAC / CODATA / Unimod)
# ============================================================================

ELECTRON_MASS_DA: float = 0.000548579909065

# Monoisotopic atomic masses (Da) harmonized with IUPAC/NIST/Unimod 12-decimal definitions
# such that composite monosaccharide formulas match the 10-decimal reference table:
#   Hex    (C6H10O5)    = 162.0528234315 Da
#   HexNAc (C8H13NO5)   = 203.0793725330 Da
#   dHex   (C6H10O4)    = 146.0579088090 Da
#   Neu5Ac (C11H17NO8)  = 291.0954165270 Da
#   Neu5Gc (C11H17NO9)  = 307.0903311490 Da
ELEMENT_MONOISOTOPIC_MASSES: Mapping[str, float] = {
    "C": 12.0000000000,
    "H": 1.00782503215,
    "N": 14.00307400525,
    "O": 15.99491462200,
    "S": 31.9720711744,
    "P": 30.9737619984,
    "Na": 22.9897692820,
    "K": 38.9637064864,
    "Li": 7.0160034366,
    "Cl": 34.968852682,
    "F": 18.9984031627,
    "I": 126.9044719,
}

PROTON_MASS_DA: float = ELEMENT_MONOISOTOPIC_MASSES["H"] - ELECTRON_MASS_DA  # 1.007276452241 Da (H - e-)
NEUTRON_MASS_SHIFT_DA: float = 1.00335483507  # Effective 13C - 12C isotopic spacing


@dataclass(frozen=True)
class IsotopeSpecification:
    """Stable isotope mass shift (relative to monoisotopic) and natural abundance."""

    mass_number: int
    exact_mass: float
    abundance: float
    neutron_offset: int


# Natural isotopic abundances (IUPAC representative isotopic compositions)
ELEMENT_ISOTOPES: Mapping[str, Tuple[IsotopeSpecification, ...]] = {
    "C": (
        IsotopeSpecification(12, 12.0000000000, 0.9893, 0),
        IsotopeSpecification(13, 13.0033548351, 0.0107, 1),
    ),
    "H": (
        IsotopeSpecification(1, 1.00782503219, 0.999885, 0),
        IsotopeSpecification(2, 2.01410177812, 0.000115, 1),
    ),
    "N": (
        IsotopeSpecification(14, 14.00307400509, 0.99636, 0),
        IsotopeSpecification(15, 15.0001088989, 0.00364, 1),
    ),
    "O": (
        IsotopeSpecification(16, 15.99491462210, 0.99757, 0),
        IsotopeSpecification(17, 16.99913175650, 0.00038, 1),
        IsotopeSpecification(18, 17.99915961286, 0.00205, 2),
    ),
    "S": (
        IsotopeSpecification(32, 31.9720711744, 0.9499, 0),
        IsotopeSpecification(33, 32.9714589098, 0.0075, 1),
        IsotopeSpecification(34, 33.967867004, 0.0425, 2),
        IsotopeSpecification(36, 35.96708071, 0.0001, 4),
    ),
    "Na": (
        IsotopeSpecification(23, 22.9897692820, 1.0000, 0),
    ),
    "K": (
        IsotopeSpecification(39, 38.9637064864, 0.932581, 0),
        IsotopeSpecification(40, 39.963998166, 0.000117, 1),
        IsotopeSpecification(41, 40.9618252579, 0.067302, 2),
    ),
    "Cl": (
        IsotopeSpecification(35, 34.968852682, 0.7576, 0),
        IsotopeSpecification(37, 36.965902602, 0.2424, 2),
    ),
}

# Common small-molecule neutral masses
WATER_FORMULA: Mapping[str, int] = {"H": 2, "O": 1}
WATER_MASS_DA: float = 2 * ELEMENT_MONOISOTOPIC_MASSES["H"] + ELEMENT_MONOISOTOPIC_MASSES["O"]  # 18.01056468648 Da

METHYLENE_FORMULA: Mapping[str, int] = {"C": 1, "H": 2}
METHYLENE_MASS_DA: float = ELEMENT_MONOISOTOPIC_MASSES["C"] + 2 * ELEMENT_MONOISOTOPIC_MASSES["H"]  # 14.01565006438 Da

AMMONIA_MASS_DA: float = ELEMENT_MONOISOTOPIC_MASSES["N"] + 3 * ELEMENT_MONOISOTOPIC_MASSES["H"]


# ============================================================================
# 2. MONOSACCHARIDE RESIDUE & FREE SUGAR DEFINITIONS
# ============================================================================


@dataclass(frozen=True)
class MonosaccharideResidue:
    """Authoritative definition of a dehydratively linked monosaccharide residue.

    Note on convention:
    - `residue_mass` is the anhydro-monosaccharide mass (free sugar minus H2O).
    - `free_mass` is the intact unlinked monosaccharide mass (residue_mass + H2O).
    - `permethylation_sites` is the number of free -OH and -NH- hydrogens per internal
      residue replaced by -CH3 (+CH2 = +14.0156500644 Da each) during complete permethylation.
       Note: In a branched glycan tree, every glycosidic linkage consumes one -OH on the
      acceptor residue and frees no new -OH, while the non-reducing terminal residue retains
      one additional free -OH that an internal residue lacks; across the entire intact tree,
      the sum of free hydroxyl/amide sites on the dehydrated residues equals the sum of
      `permethylation_sites` per residue plus 2 sites on the terminal reducing-end H2O
      (or 3 sites if reduced to an alditol).
    """

    name: str
    symbol: str
    aliases: Tuple[str, ...]
    residue_formula: Mapping[str, int]
    residue_mass: float
    free_mass: float
    permethylation_sites: int
    is_sialic_acid: bool = False


MONOSACCHARIDE_RESIDUES: Mapping[str, MonosaccharideResidue] = {
    "Hex": MonosaccharideResidue(
        name="Hexose",
        symbol="Hex",
        aliases=("Hex", "H", "Man", "Gal", "Glc"),
        residue_formula={"C": 6, "H": 10, "O": 5},
        residue_mass=162.0528234315,
        free_mass=162.0528234315 + WATER_MASS_DA,
        permethylation_sites=3,
        is_sialic_acid=False,
    ),
    "HexNAc": MonosaccharideResidue(
        name="N-acetylhexosamine",
        symbol="HexNAc",
        aliases=("HexNAc", "N", "GlcNAc", "GalNAc"),
        residue_formula={"C": 8, "H": 13, "N": 1, "O": 5},
        residue_mass=203.0793725330,
        free_mass=203.0793725330 + WATER_MASS_DA,
        permethylation_sites=3,  # Two free -OH groups + one acetamido -NH- group
        is_sialic_acid=False,
    ),
    "dHex": MonosaccharideResidue(
        name="Deoxyhexose (Fucose)",
        symbol="dHex",
        aliases=("dHex", "Fuc", "F"),
        residue_formula={"C": 6, "H": 10, "O": 4},
        residue_mass=146.0579088090,
        free_mass=146.0579088090 + WATER_MASS_DA,
        permethylation_sites=2,
        is_sialic_acid=False,
    ),
    "Neu5Ac": MonosaccharideResidue(
        name="N-acetylneuraminic acid",
        symbol="Neu5Ac",
        aliases=("Neu5Ac", "NeuAc", "S", "Sia"),
        residue_formula={"C": 11, "H": 17, "N": 1, "O": 8},
        residue_mass=291.0954165270,
        free_mass=291.0954165270 + WATER_MASS_DA,
        permethylation_sites=5,  # Three -OH + one -COOH + one acetamido -NH-
        is_sialic_acid=True,
    ),
    "Neu5Gc": MonosaccharideResidue(
        name="N-glycolylneuraminic acid",
        symbol="Neu5Gc",
        aliases=("Neu5Gc", "NeuGc", "G"),
        residue_formula={"C": 11, "H": 17, "N": 1, "O": 9},
        residue_mass=307.0903311490,
        free_mass=307.0903311490 + WATER_MASS_DA,
        permethylation_sites=6,  # Four -OH (including glycolyl -OH) + one -COOH + one amide -NH-
        is_sialic_acid=True,
    ),
    "Pent": MonosaccharideResidue(
        name="Pentose (Xylose)",
        symbol="Pent",
        aliases=("Pent", "Xyl", "X"),
        residue_formula={"C": 5, "H": 8, "O": 4},
        residue_mass=(
            5 * ELEMENT_MONOISOTOPIC_MASSES["C"]
            + 8 * ELEMENT_MONOISOTOPIC_MASSES["H"]
            + 4 * ELEMENT_MONOISOTOPIC_MASSES["O"]
        ),
        free_mass=(
            5 * ELEMENT_MONOISOTOPIC_MASSES["C"]
            + 8 * ELEMENT_MONOISOTOPIC_MASSES["H"]
            + 4 * ELEMENT_MONOISOTOPIC_MASSES["O"]
            + WATER_MASS_DA
        ),
        permethylation_sites=2,
        is_sialic_acid=False,
    ),
    "HexA": MonosaccharideResidue(
        name="Hexuronic acid",
        symbol="HexA",
        aliases=("HexA", "GlcA", "IdoA"),
        residue_formula={"C": 6, "H": 8, "O": 6},
        residue_mass=(
            6 * ELEMENT_MONOISOTOPIC_MASSES["C"]
            + 8 * ELEMENT_MONOISOTOPIC_MASSES["H"]
            + 6 * ELEMENT_MONOISOTOPIC_MASSES["O"]
        ),
        free_mass=(
            6 * ELEMENT_MONOISOTOPIC_MASSES["C"]
            + 8 * ELEMENT_MONOISOTOPIC_MASSES["H"]
            + 6 * ELEMENT_MONOISOTOPIC_MASSES["O"]
            + WATER_MASS_DA
        ),
        permethylation_sites=3,
        is_sialic_acid=False,
    ),
}

# Alias lookup mapping (e.g., "Fuc" -> "dHex", "NeuAc" -> "Neu5Ac")
RESIDUE_ALIAS_MAP: Dict[str, str] = {}
for _canonical_symbol, _res in MONOSACCHARIDE_RESIDUES.items():
    RESIDUE_ALIAS_MAP[_canonical_symbol] = _canonical_symbol
    RESIDUE_ALIAS_MAP[_canonical_symbol.lower()] = _canonical_symbol
    for _alias in _res.aliases:
        RESIDUE_ALIAS_MAP[_alias] = _canonical_symbol
        RESIDUE_ALIAS_MAP[_alias.lower()] = _canonical_symbol


# ============================================================================
# 3. IONIZATION POLARITY & ADDUCT DEFINITIONS
# ============================================================================


class IonPolarity(str, Enum):
    """Mass spectrometry ionization polarity."""

    POSITIVE = "positive"
    NEGATIVE = "negative"


@dataclass(frozen=True)
class AdductDefinition:
    """Authoritative definition of an ion species / adduct.

    For charge state z = abs(charge):
      m/z = (neutral_mass * mult + net_carrier_mass) / z
    where `net_carrier_mass` already accounts for the loss or gain of z electrons!
    """

    code: str
    name: str
    polarity: IonPolarity
    charge: int  # Signed integer (+1, +2, +3, -1, -2, etc.)
    net_carrier_mass: float  # Exact cationic/anionic mass added to M (already includes electron offset)
    formula_added: Mapping[str, int]
    formula_removed: Mapping[str, int]
    description: str


# Exact ionic masses of individual charge carriers (neutral atom mass minus electron mass for cations,
# neutral atom mass plus electron mass for anions)
CATION_H_MASS: float = PROTON_MASS_DA  # 1.007276466621 Da (H - e-)
CATION_NA_MASS: float = ELEMENT_MONOISOTOPIC_MASSES["Na"] - ELECTRON_MASS_DA  # 22.9892207021 Da (Na - e-)
CATION_K_MASS: float = ELEMENT_MONOISOTOPIC_MASSES["K"] - ELECTRON_MASS_DA  # 38.9631579065 Da (K - e-)
CATION_LI_MASS: float = ELEMENT_MONOISOTOPIC_MASSES["Li"] - ELECTRON_MASS_DA  # 7.0154548567 Da (Li - e-)
CATION_NH4_MASS: float = (
    ELEMENT_MONOISOTOPIC_MASSES["N"] + 4 * ELEMENT_MONOISOTOPIC_MASSES["H"] - ELECTRON_MASS_DA
)  # 18.0338255533 Da (NH4 - e-)

ANION_DEPROTONATED_MASS: float = -(ELEMENT_MONOISOTOPIC_MASSES["H"] - ELECTRON_MASS_DA)  # -1.007276466621 Da (-H+)
ANION_CL_MASS: float = ELEMENT_MONOISOTOPIC_MASSES["Cl"] + ELECTRON_MASS_DA  # 34.9694012619 Da (Cl + e-)
ANION_FORMATE_MASS: float = (
    ELEMENT_MONOISOTOPIC_MASSES["C"]
    + 2 * ELEMENT_MONOISOTOPIC_MASSES["H"]
    + 2 * ELEMENT_MONOISOTOPIC_MASSES["O"]
    + ELECTRON_MASS_DA
)  # HCOO-
ANION_ACETATE_MASS: float = (
    2 * ELEMENT_MONOISOTOPIC_MASSES["C"]
    + 3 * ELEMENT_MONOISOTOPIC_MASSES["H"]
    + 2 * ELEMENT_MONOISOTOPIC_MASSES["O"]
    + ELECTRON_MASS_DA
)  # CH3COO-

ADDUCT_DEFINITIONS: Mapping[str, AdductDefinition] = {
    "[M+H]+": AdductDefinition(
        code="[M+H]+",
        name="Protonated (+1)",
        polarity=IonPolarity.POSITIVE,
        charge=1,
        net_carrier_mass=CATION_H_MASS,
        formula_added={"H": 1},
        formula_removed={},
        description="Singly charged protonated ion: m/z = M + m(H+)",
    ),
    "[M+Na]+": AdductDefinition(
        code="[M+Na]+",
        name="Sodium adduct (+1)",
        polarity=IonPolarity.POSITIVE,
        charge=1,
        net_carrier_mass=CATION_NA_MASS,
        formula_added={"Na": 1},
        formula_removed={},
        description="Singly charged sodium adduct: m/z = M + m(Na+)",
    ),
    "[M+K]+": AdductDefinition(
        code="[M+K]+",
        name="Potassium adduct (+1)",
        polarity=IonPolarity.POSITIVE,
        charge=1,
        net_carrier_mass=CATION_K_MASS,
        formula_added={"K": 1},
        formula_removed={},
        description="Singly charged potassium adduct: m/z = M + m(K+)",
    ),
    "[M+Li]+": AdductDefinition(
        code="[M+Li]+",
        name="Lithium adduct (+1)",
        polarity=IonPolarity.POSITIVE,
        charge=1,
        net_carrier_mass=CATION_LI_MASS,
        formula_added={"Li": 1},
        formula_removed={},
        description="Singly charged lithium adduct: m/z = M + m(Li+)",
    ),
    "[M+NH4]+": AdductDefinition(
        code="[M+NH4]+",
        name="Ammonium adduct (+1)",
        polarity=IonPolarity.POSITIVE,
        charge=1,
        net_carrier_mass=CATION_NH4_MASS,
        formula_added={"N": 1, "H": 4},
        formula_removed={},
        description="Singly charged ammonium adduct: m/z = M + m(NH4+)",
    ),
    "[M+2Na-H]+": AdductDefinition(
        code="[M+2Na-H]+",
        name="Disodium-monodeprotonated adduct (+1)",
        polarity=IonPolarity.POSITIVE,
        charge=1,
        net_carrier_mass=2 * CATION_NA_MASS - CATION_H_MASS,
        formula_added={"Na": 2},
        formula_removed={"H": 1},
        description="Singly charged sodium salt adduct (common in sialylated native glycans): m/z = M + 2*m(Na+) - m(H+)",
    ),
    "[M+2H]2+": AdductDefinition(
        code="[M+2H]2+",
        name="Doubly protonated (+2)",
        polarity=IonPolarity.POSITIVE,
        charge=2,
        net_carrier_mass=2 * CATION_H_MASS,
        formula_added={"H": 2},
        formula_removed={},
        description="Doubly charged protonated ion: m/z = (M + 2*m(H+)) / 2",
    ),
    "[M+2Na]2+": AdductDefinition(
        code="[M+2Na]2+",
        name="Doubly sodiated (+2)",
        polarity=IonPolarity.POSITIVE,
        charge=2,
        net_carrier_mass=2 * CATION_NA_MASS,
        formula_added={"Na": 2},
        formula_removed={},
        description="Doubly charged sodium adduct: m/z = (M + 2*m(Na+)) / 2",
    ),
    "[M+H+Na]2+": AdductDefinition(
        code="[M+H+Na]2+",
        name="Protonated + sodiated (+2)",
        polarity=IonPolarity.POSITIVE,
        charge=2,
        net_carrier_mass=CATION_H_MASS + CATION_NA_MASS,
        formula_added={"H": 1, "Na": 1},
        formula_removed={},
        description="Doubly charged mixed proton/sodium adduct: m/z = (M + m(H+) + m(Na+)) / 2",
    ),
    "[M+3H]3+": AdductDefinition(
        code="[M+3H]3+",
        name="Triply protonated (+3)",
        polarity=IonPolarity.POSITIVE,
        charge=3,
        net_carrier_mass=3 * CATION_H_MASS,
        formula_added={"H": 3},
        formula_removed={},
        description="Triply charged protonated ion: m/z = (M + 3*m(H+)) / 3",
    ),
    "[M+3Na]3+": AdductDefinition(
        code="[M+3Na]3+",
        name="Triply sodiated (+3)",
        polarity=IonPolarity.POSITIVE,
        charge=3,
        net_carrier_mass=3 * CATION_NA_MASS,
        formula_added={"Na": 3},
        formula_removed={},
        description="Triply charged sodium adduct: m/z = (M + 3*m(Na+)) / 3",
    ),
    "[M-H]-": AdductDefinition(
        code="[M-H]-",
        name="Deprotonated (-1)",
        polarity=IonPolarity.NEGATIVE,
        charge=-1,
        net_carrier_mass=ANION_DEPROTONATED_MASS,
        formula_added={},
        formula_removed={"H": 1},
        description="Singly charged deprotonated negative ion: m/z = M - m(H+)",
    ),
    "[M+Cl]-": AdductDefinition(
        code="[M+Cl]-",
        name="Chloride adduct (-1)",
        polarity=IonPolarity.NEGATIVE,
        charge=-1,
        net_carrier_mass=ANION_CL_MASS,
        formula_added={"Cl": 1},
        formula_removed={},
        description="Singly charged chloride adduct: m/z = M + m(Cl-)",
    ),
    "[M+HCOO]-": AdductDefinition(
        code="[M+HCOO]-",
        name="Formate adduct (-1)",
        polarity=IonPolarity.NEGATIVE,
        charge=-1,
        net_carrier_mass=ANION_FORMATE_MASS,
        formula_added={"C": 1, "H": 1, "O": 2},
        formula_removed={},
        description="Singly charged formate adduct: m/z = M + m(HCOO-)",
    ),
    "[M-2H]2-": AdductDefinition(
        code="[M-2H]2-",
        name="Doubly deprotonated (-2)",
        polarity=IonPolarity.NEGATIVE,
        charge=-2,
        net_carrier_mass=2 * ANION_DEPROTONATED_MASS,
        formula_added={},
        formula_removed={"H": 2},
        description="Doubly charged deprotonated negative ion: m/z = (M - 2*m(H+)) / 2",
    ),
}


# ============================================================================
# 4. GLYCAN FORM & DERIVATIZATION / LABELING REGISTRY
# ============================================================================


class GlycanForm(str, Enum):
    """Distinguishes how the glycan exists chemically."""

    RELEASED_FREE = "released_free"  # Released reducing-end glycan (includes reducing-end H2O)
    RELEASED_LABELED = "released_labeled"  # Released glycan modified with a reducing-end tag
    PERNETHYLATED = "permethylated"  # Site-specifically permethylated released glycan
    PEPTIDE_ATTACHED = "peptide_attached"  # Dehydratively linked glycosyl moiety on a peptide (NO terminal H2O)


class ValidationStatus(str, Enum):
    """Validation status of a derivatization chemistry or scientific rule."""

    VALIDATED_ENABLED = "VALIDATED_ENABLED"
    EXPERIMENTAL = "EXPERIMENTAL"
    NOT_IMPLEMENTED = "NOT_IMPLEMENTED"


@dataclass(frozen=True)
class LabelingMethodDefinition:
    """Structured registry entry for glycan reducing-end labeling & derivatization reactions.

    `net_mass_shift_from_free_reducing` is the exact monoisotopic mass delta added to a
    native free reducing glycan (which already includes the reducing-end H2O).
    For reductive amination:
      Glycan-CHO + H2N-R --(NaBH3CN, -H2O, +2H)--> Glycan-CH2-NH-R
      Net elemental delta relative to Free Glycan = Formula(Tag) - H2O + H2 = Formula(Tag) - O.
    """

    code: str
    name: str
    category: str  # "native", "reducing_end_label", "permethylation", "linkage_specific_sialic"
    net_mass_shift_from_free_reducing: Optional[float]
    net_formula_added: Mapping[str, int]
    net_formula_removed: Mapping[str, int]
    is_site_dependent: bool
    applicable_to_peptide_attached: bool
    validation_status: ValidationStatus
    reaction_description: str
    literature_reference: str


# Exact net shifts calculated from elemental composition:
# 1. Reduced Alditol (NaBH4 reduction of free reducing end: +H2):
#    +2 * 1.00782503219 = +2.0156500644 Da
# 2. 2-AB (2-Aminobenzamide, C7H8N2O, 136.06366288 Da):
#    Reductive amination (-H2O + H2 = -O): net formula added = C7H8N2
#    7*12.0 + 8*1.00782503219 + 2*14.00307400509 = +120.0687482677 Da
# 3. 2-AA (2-Aminobenzoic acid / anthranilic acid, C7H7NO2, 137.04767847 Da):
#    Reductive amination (-H2O + H2 = -O): net formula added = C7H7NO
#    7*12.0 + 7*1.00782503219 + 14.00307400509 + 15.99491462210 = +121.0527638525 Da
# 4. Procainamide (4-amino-N-(2-(diethylamino)ethyl)benzamide, C13H21N3O, 235.1684623 Da):
#    Reductive amination (-H2O + H2 = -O): net formula added = C13H21N3
#    13*12.0 + 21*1.00782503219 + 3*14.00307400509 = +219.1735476913 Da
# 5. PA (2-Aminopyridine, C5H6N2, 94.0530982 Da):
#    Reductive amination (-H2O + H2 = -O): net formula added = C5H6N2 - O
#    5*12.0 + 6*1.00782503219 + 2*14.00307400509 - 15.99491462210 = +78.0581835812 Da
# 6. Girard's Reagent P (Hydrazide labeling of reducing end, C7H10N3O+ Cl-, cation C7H10N3O+):
#    Hydrazone formation (-H2O): net formula added for neutral-equivalent mass calculation

LABELING_REGISTRY: Mapping[str, LabelingMethodDefinition] = {
    "NONE": LabelingMethodDefinition(
        code="NONE",
        name="Native (Unlabeled Free Reducing End)",
        category="native",
        net_mass_shift_from_free_reducing=0.0,
        net_formula_added={},
        net_formula_removed={},
        is_site_dependent=False,
        applicable_to_peptide_attached=True,
        validation_status=ValidationStatus.VALIDATED_ENABLED,
        reaction_description="Unmodified native N-glycan (with reducing-end H2O when released; without terminal H2O when peptide-attached).",
        literature_reference="Standard IUPAC/Consortium for Functional Glycomics (CFG) native mass convention.",
    ),
    "REDUCED": LabelingMethodDefinition(
        code="REDUCED",
        name="Reduced Alditol (NaBH4)",
        category="reducing_end_label",
        net_mass_shift_from_free_reducing=2 * ELEMENT_MONOISOTOPIC_MASSES["H"],
        net_formula_added={"H": 2},
        net_formula_removed={},
        is_site_dependent=False,
        applicable_to_peptide_attached=False,
        validation_status=ValidationStatus.VALIDATED_ENABLED,
        reaction_description="Borohydride reduction of the reducing-end hemiacetal to an alditol (+H2 = +2.015650 Da).",
        literature_reference="Harvey DJ, Mass Spectrom Rev. 1999;18(6):349-451.",
    ),
    "2-AB": LabelingMethodDefinition(
        code="2-AB",
        name="2-Aminobenzamide (2-AB)",
        category="reducing_end_label",
        net_mass_shift_from_free_reducing=(
            7 * ELEMENT_MONOISOTOPIC_MASSES["C"]
            + 8 * ELEMENT_MONOISOTOPIC_MASSES["H"]
            + 2 * ELEMENT_MONOISOTOPIC_MASSES["N"]
        ),
        net_formula_added={"C": 7, "H": 8, "N": 2},
        net_formula_removed={},
        is_site_dependent=False,
        applicable_to_peptide_attached=False,
        validation_status=ValidationStatus.VALIDATED_ENABLED,
        reaction_description="Reductive amination with 2-aminobenzamide (C7H8N2O - H2O + H2 = +C7H8N2, +120.068748 Da).",
        literature_reference="Bigge JC et al., Anal Biochem. 1995;230(2):229-238.",
    ),
    "2-AA": LabelingMethodDefinition(
        code="2-AA",
        name="2-Aminobenzoic Acid (2-AA)",
        category="reducing_end_label",
        net_mass_shift_from_free_reducing=(
            7 * ELEMENT_MONOISOTOPIC_MASSES["C"]
            + 7 * ELEMENT_MONOISOTOPIC_MASSES["H"]
            + 1 * ELEMENT_MONOISOTOPIC_MASSES["N"]
            + 1 * ELEMENT_MONOISOTOPIC_MASSES["O"]
        ),
        net_formula_added={"C": 7, "H": 7, "N": 1, "O": 1},
        net_formula_removed={},
        is_site_dependent=False,
        applicable_to_peptide_attached=False,
        validation_status=ValidationStatus.VALIDATED_ENABLED,
        reaction_description="Reductive amination with 2-aminobenzoic acid (C7H7NO2 - H2O + H2 = +C7H7NO, +121.052764 Da).",
        literature_reference="Anumula KR, Anal Biochem. 2006;350(1):1-23.",
    ),
    "PROCAINAMIDE": LabelingMethodDefinition(
        code="PROCAINAMIDE",
        name="Procainamide (ProA)",
        category="reducing_end_label",
        net_mass_shift_from_free_reducing=(
            13 * ELEMENT_MONOISOTOPIC_MASSES["C"]
            + 21 * ELEMENT_MONOISOTOPIC_MASSES["H"]
            + 3 * ELEMENT_MONOISOTOPIC_MASSES["N"]
        ),
        net_formula_added={"C": 13, "H": 21, "N": 3},
        net_formula_removed={},
        is_site_dependent=False,
        applicable_to_peptide_attached=False,
        validation_status=ValidationStatus.VALIDATED_ENABLED,
        reaction_description="Reductive amination with procainamide (C13H21N3O - H2O + H2 = +C13H21N3, +219.173548 Da).",
        literature_reference="Kozak RP et al., Anal Biochem. 2015;486:38-40.",
    ),
    "PA": LabelingMethodDefinition(
        code="PA",
        name="2-Aminopyridine (PA)",
        category="reducing_end_label",
        net_mass_shift_from_free_reducing=(
            5 * ELEMENT_MONOISOTOPIC_MASSES["C"]
            + 6 * ELEMENT_MONOISOTOPIC_MASSES["H"]
            + 2 * ELEMENT_MONOISOTOPIC_MASSES["N"]
            - 1 * ELEMENT_MONOISOTOPIC_MASSES["O"]
        ),
        net_formula_added={"C": 5, "H": 6, "N": 2},
        net_formula_removed={"O": 1},
        is_site_dependent=False,
        applicable_to_peptide_attached=False,
        validation_status=ValidationStatus.VALIDATED_ENABLED,
        reaction_description="Reductive amination with 2-aminopyridine (C5H6N2 - H2O + H2 = +C5H6N2 - O, +78.058184 Da).",
        literature_reference="Hase S et al., J Biochem. 1978;84(5):1325-1327.",
    ),
    "PERMETHYLATED": LabelingMethodDefinition(
        code="PERMETHYLATED",
        name="Permethylated (Native Non-Reduced)",
        category="permethylation",
        net_mass_shift_from_free_reducing=None,  # Calculated dynamically from residue counts
        net_formula_added={},
        net_formula_removed={},
        is_site_dependent=True,
        applicable_to_peptide_attached=False,
        validation_status=ValidationStatus.VALIDATED_ENABLED,
        reaction_description=(
            "Complete O- and N-permethylation of non-reduced free glycan (+CH2 = +14.015650 Da per site: "
            "Hex=3, HexNAc=3, dHex=2, Neu5Ac=5, Neu5Gc=6, Pent=2, HexA=3, plus +2*CH2 = +28.031300 Da for reducing-end H2O)."
        ),
        literature_reference="Ciucanu I, Kerek F, Carbohydr Res. 1984;131(2):209-217.",
    ),
    "PERMETHYLATED_REDUCED": LabelingMethodDefinition(
        code="PERMETHYLATED_REDUCED",
        name="Permethylated (Reduced Alditol)",
        category="permethylation",
        net_mass_shift_from_free_reducing=None,  # Calculated dynamically from residue counts
        net_formula_added={},
        net_formula_removed={},
        is_site_dependent=True,
        applicable_to_peptide_attached=False,
        validation_status=ValidationStatus.VALIDATED_ENABLED,
        reaction_description=(
            "NaBH4 reduction (+H2) followed by complete O- and N-permethylation (+CH2 per site on residues "
            "plus +3*CH2 = +42.046950 Da on the opened alditol reducing end, net +44.062600 Da terminal contribution above native H2O)."
        ),
        literature_reference="Ciucanu I, Costello CE, J Am Chem Soc. 2003;125(52):16213-16219.",
    ),
    "EEC_SIALIC": LabelingMethodDefinition(
        code="EEC_SIALIC",
        name="Linkage-Specific Sialic Acid Ethyl Esterification / Lactonization (ETH-EDC / HOBt)",
        category="linkage_specific_sialic",
        net_mass_shift_from_free_reducing=None,  # Calculated from alpha2,3 vs alpha2,6 sialic acid counts
        net_formula_added={},
        net_formula_removed={},
        is_site_dependent=True,
        applicable_to_peptide_attached=True,
        validation_status=ValidationStatus.VALIDATED_ENABLED,
        reaction_description=(
            "Linkage-specific sialic acid derivatization in ethanol with EDC/HOBt: "
            "alpha2,6-linked sialic acids undergo ethyl esterification (+C2H4 = +28.031300 Da per residue), "
            "whereas alpha2,3-linked sialic acids undergo intramolecular lactonization (-H2O = -18.010565 Da per residue)."
        ),
        literature_reference="Reiding KR et al., Anal Chem. 2014;86(12):5784-5793.",
    ),
    "DMA_SIALIC": LabelingMethodDefinition(
        code="DMA_SIALIC",
        name="Linkage-Specific Sialic Acid Dimethylamidation / Lactonization",
        category="linkage_specific_sialic",
        net_mass_shift_from_free_reducing=None,  # Calculated from alpha2,3 vs alpha2,6 sialic acid counts
        net_formula_added={},
        net_formula_removed={},
        is_site_dependent=True,
        applicable_to_peptide_attached=True,
        validation_status=ValidationStatus.VALIDATED_ENABLED,
        reaction_description=(
            "Linkage-specific dimethylamidation: alpha2,6-linked sialic acids form a dimethylamide "
            "(+C2H5N - O = +27.047285 Da per residue), whereas alpha2,3-linked sialic acids form a lactone "
            "(-H2O = -18.010565 Da per residue) or ammonia/methylamide product depending on second-step nucleophile."
        ),
        literature_reference="de Haan N et al., Anal Chem. 2015;87(16):8284-8291.",
    ),
}

# Linkage-specific sialic acid mass shifts per residue
SIALIC_EEC_ALPHA26_SHIFT_DA: float = (
    2 * ELEMENT_MONOISOTOPIC_MASSES["C"] + 4 * ELEMENT_MONOISOTOPIC_MASSES["H"]
)  # +C2H4 = +28.0313001288 Da
SIALIC_EEC_ALPHA23_SHIFT_DA: float = -WATER_MASS_DA  # -H2O = -18.0105646865 Da

SIALIC_DMA_ALPHA26_SHIFT_DA: float = (
    2 * ELEMENT_MONOISOTOPIC_MASSES["C"]
    + 5 * ELEMENT_MONOISOTOPIC_MASSES["H"]
    + 1 * ELEMENT_MONOISOTOPIC_MASSES["N"]
    - 1 * ELEMENT_MONOISOTOPIC_MASSES["O"]
)  # +C2H5N - O = +27.0472845440 Da
SIALIC_DMA_ALPHA23_SHIFT_DA: float = -WATER_MASS_DA  # -H2O = -18.0105646865 Da


# ============================================================================
# 5. AMINO ACID RESIDUE MASSES & PEPTIDE MODIFICATIONS (FOR STAGE 6)
# ============================================================================

@dataclass(frozen=True)
class AminoAcidResidue:
    """Standard amino acid residue definition (dehydrated peptide bond residue)."""

    code1: str
    code3: str
    name: str
    formula: Mapping[str, int]
    monoisotopic_mass: float


def _calc_formula_mass(formula: Mapping[str, int]) -> float:
    return sum(ELEMENT_MONOISOTOPIC_MASSES[el] * count for el, count in formula.items())


AMINO_ACID_RESIDUES: Mapping[str, AminoAcidResidue] = {
    "A": AminoAcidResidue("A", "Ala", "Alanine", {"C": 3, "H": 5, "N": 1, "O": 1}, _calc_formula_mass({"C": 3, "H": 5, "N": 1, "O": 1})),
    "R": AminoAcidResidue("R", "Arg", "Arginine", {"C": 6, "H": 12, "N": 4, "O": 1}, _calc_formula_mass({"C": 6, "H": 12, "N": 4, "O": 1})),
    "N": AminoAcidResidue("N", "Asn", "Asparagine", {"C": 4, "H": 6, "N": 2, "O": 2}, _calc_formula_mass({"C": 4, "H": 6, "N": 2, "O": 2})),
    "D": AminoAcidResidue("D", "Asp", "Aspartic acid", {"C": 4, "H": 5, "N": 1, "O": 3}, _calc_formula_mass({"C": 4, "H": 5, "N": 1, "O": 3})),
    "C": AminoAcidResidue("C", "Cys", "Cysteine", {"C": 3, "H": 5, "N": 1, "O": 1, "S": 1}, _calc_formula_mass({"C": 3, "H": 5, "N": 1, "O": 1, "S": 1})),
    "E": AminoAcidResidue("E", "Glu", "Glutamic acid", {"C": 5, "H": 7, "N": 1, "O": 3}, _calc_formula_mass({"C": 5, "H": 7, "N": 1, "O": 3})),
    "Q": AminoAcidResidue("Q", "Gln", "Glutamine", {"C": 5, "H": 8, "N": 2, "O": 2}, _calc_formula_mass({"C": 5, "H": 8, "N": 2, "O": 2})),
    "G": AminoAcidResidue("G", "Gly", "Glycine", {"C": 2, "H": 3, "N": 1, "O": 1}, _calc_formula_mass({"C": 2, "H": 3, "N": 1, "O": 1})),
    "H": AminoAcidResidue("H", "His", "Histidine", {"C": 6, "H": 7, "N": 3, "O": 1}, _calc_formula_mass({"C": 6, "H": 7, "N": 3, "O": 1})),
    "I": AminoAcidResidue("I", "Ile", "Isoleucine", {"C": 6, "H": 11, "N": 1, "O": 1}, _calc_formula_mass({"C": 6, "H": 11, "N": 1, "O": 1})),
    "L": AminoAcidResidue("L", "Leu", "Leucine", {"C": 6, "H": 11, "N": 1, "O": 1}, _calc_formula_mass({"C": 6, "H": 11, "N": 1, "O": 1})),
    "K": AminoAcidResidue("K", "Lys", "Lysine", {"C": 6, "H": 12, "N": 2, "O": 1}, _calc_formula_mass({"C": 6, "H": 12, "N": 2, "O": 1})),
    "M": AminoAcidResidue("M", "Met", "Methionine", {"C": 5, "H": 9, "N": 1, "O": 1, "S": 1}, _calc_formula_mass({"C": 5, "H": 9, "N": 1, "O": 1, "S": 1})),
    "F": AminoAcidResidue("F", "Phe", "Phenylalanine", {"C": 9, "H": 9, "N": 1, "O": 1}, _calc_formula_mass({"C": 9, "H": 9, "N": 1, "O": 1})),
    "P": AminoAcidResidue("P", "Pro", "Proline", {"C": 5, "H": 7, "N": 1, "O": 1}, _calc_formula_mass({"C": 5, "H": 7, "N": 1, "O": 1})),
    "S": AminoAcidResidue("S", "Ser", "Serine", {"C": 3, "H": 5, "N": 1, "O": 2}, _calc_formula_mass({"C": 3, "H": 5, "N": 1, "O": 2})),
    "T": AminoAcidResidue("T", "Thr", "Threonine", {"C": 4, "H": 7, "N": 1, "O": 2}, _calc_formula_mass({"C": 4, "H": 7, "N": 1, "O": 2})),
    "W": AminoAcidResidue("W", "Trp", "Tryptophan", {"C": 11, "H": 10, "N": 2, "O": 1}, _calc_formula_mass({"C": 11, "H": 10, "N": 2, "O": 1})),
    "Y": AminoAcidResidue("Y", "Tyr", "Tyrosine", {"C": 9, "H": 9, "N": 1, "O": 2}, _calc_formula_mass({"C": 9, "H": 9, "N": 1, "O": 2})),
    "V": AminoAcidResidue("V", "Val", "Valine", {"C": 5, "H": 9, "N": 1, "O": 1}, _calc_formula_mass({"C": 5, "H": 9, "N": 1, "O": 1})),
}


@dataclass(frozen=True)
class PeptideModificationDefinition:
    """Authoritative definition of a peptide chemical or post-translational modification."""

    code: str
    name: str
    target_residues: Tuple[str, ...]
    formula_delta: Mapping[str, int]
    mass_delta: float
    is_fixed_by_default: bool = False


PEPTIDE_MODIFICATIONS: Mapping[str, PeptideModificationDefinition] = {
    "Carbamidomethyl_C": PeptideModificationDefinition(
        code="Carbamidomethyl_C",
        name="Carbamidomethyl (Cys)",
        target_residues=("C",),
        formula_delta={"C": 2, "H": 3, "N": 1, "O": 1},
        mass_delta=_calc_formula_mass({"C": 2, "H": 3, "N": 1, "O": 1}),  # +57.02146372 Da
        is_fixed_by_default=True,
    ),
    "Oxidation_M": PeptideModificationDefinition(
        code="Oxidation_M",
        name="Oxidation (Met)",
        target_residues=("M",),
        formula_delta={"O": 1},
        mass_delta=ELEMENT_MONOISOTOPIC_MASSES["O"],  # +15.99491462 Da
        is_fixed_by_default=False,
    ),
    "Deamidation_NQ": PeptideModificationDefinition(
        code="Deamidation_NQ",
        name="Deamidation (Asn/Gln)",
        target_residues=("N", "Q"),
        formula_delta={"O": 1, "N": -1, "H": -1},
        mass_delta=(
            ELEMENT_MONOISOTOPIC_MASSES["O"]
            - ELEMENT_MONOISOTOPIC_MASSES["N"]
            - ELEMENT_MONOISOTOPIC_MASSES["H"]
        ),  # +0.98401558 Da
        is_fixed_by_default=False,
    ),
    "Pyro_glu_Q": PeptideModificationDefinition(
        code="Pyro_glu_Q",
        name="Pyro-glu from N-terminal Gln (-NH3)",
        target_residues=("Q",),
        formula_delta={"N": -1, "H": -3},
        mass_delta=-AMMONIA_MASS_DA,  # -17.02654910 Da
        is_fixed_by_default=False,
    ),
}

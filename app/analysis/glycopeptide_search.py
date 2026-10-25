"""Glycopeptide screening, in-silico protease digestion, N-glycosylation sequon detection, and MS1/MS2 matching engine.

CRITICAL SCIENTIFIC RULES:
1. Glycopeptide mass calculation MUST use the peptide-attached glycan moiety mass
   (`calculate_attached_glycan_moiety_mass` / `calculate_glycopeptide_neutral_mass`),
   which omits the released glycan's terminal reducing-end H2O (18.0105646865 Da).
2. Canonical N-glycosylation sequons (`N-X-S/T` where `X != P`) are reported as
   candidate attachment sites (`1-based` peptide and protein coordinates); a sequon
   is a potential site, not proof of site occupancy.
3. Evidence levels (`GlycopeptideConfidenceCategory`) are strictly enforced:
   - `MS1_MASS_MATCH` (custom precursor table / mass-only match)
   - `MS1_COMPOSITION_CANDIDATE` (precursor m/z match to peptide + glycan composition)
   - `MSMS_SUPPORTED` (requires observed MS/MS spectrum with verified B-type oxonium ions
     AND Y-type peptide/core-GlcNAc fragment ions)
   - `SITE_LOCALIZED` (requires MS/MS support AND unambiguous single sequon or site-determining fragment)
   - `REFERENCE_VALIDATED` (matched to an externally validated reference standard)
"""

from __future__ import annotations

import itertools
import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

from app.core.ion_calculator import (
    build_homogeneous_adduct_for_charge,
    calculate_mz_from_neutral_mass,
    resolve_adduct_definition,
)
from app.core.mass_calculator import (
    calculate_attached_glycan_moiety_mass,
    calculate_glycopeptide_neutral_mass,
    calculate_peptide_neutral_mass,
    parse_glycan_composition,
)
from app.core.mass_error import (
    MassToleranceUnit,
    evaluate_mass_error,
    is_within_mass_tolerance,
)
from app.core.scientific_constants import (
    AMINO_ACID_RESIDUES,
    CATION_H_MASS,
    IonPolarity,
    MONOSACCHARIDE_RESIDUES,
    PEPTIDE_MODIFICATIONS,
)
from app.libraries.composition_generator import GlycanLibrary
from app.parsers.base_models import SpectrumData, decode_text_bytes


class ProteaseEnzyme(str, Enum):
    """Supported proteases for in-silico protein digestion."""

    TRYPSIN = "trypsin"  # Cleaves C-terminal to K/R unless followed by P
    TRYPSIN_P = "trypsin_p"  # Cleaves C-terminal to K/R even before P
    CHYMOTRYPSIN = "chymotrypsin"  # Cleaves C-terminal to F/Y/W/L unless followed by P
    GLU_C = "glu_c"  # Cleaves C-terminal to E (and D) unless followed by P
    LYS_C = "lys_c"  # Cleaves C-terminal to K unless followed by P
    NO_CLEAVAGE = "no_cleavage"  # Intact sequence (user-supplied peptides)


class GlycopeptideConfidenceCategory(str, Enum):
    """Explicit identification confidence categories for glycopeptide screening."""

    MS1_MASS_MATCH = "MS1_MASS_MATCH"
    MS1_COMPOSITION_CANDIDATE = "MS1_COMPOSITION_CANDIDATE"
    MSMS_SUPPORTED = "MSMS_SUPPORTED"
    SITE_LOCALIZED = "SITE_LOCALIZED"
    REFERENCE_VALIDATED = "REFERENCE_VALIDATED"


@dataclass(frozen=True)
class FastaProteinRecord:
    """Single parsed protein sequence from a FASTA input."""

    protein_id: str
    description: str
    sequence: str


@dataclass(frozen=True)
class SequonSite:
    """Candidate N-linked glycosylation sequon (`N-X-S/T`, `X != P`)."""

    motif: str  # e.g., "NST", "NKT"
    peptide_position_1based: int  # 1-based index of Asn ('N') within the peptide
    protein_position_1based: Optional[int]  # 1-based index of Asn ('N') within the parent protein


@dataclass(frozen=True)
class DigestedPeptideCandidate:
    """In-silico digested or user-supplied peptide with modification state and sequon positions."""

    protein_id: str
    protein_sequence: str
    peptide_sequence: str
    start_pos_1based: int
    end_pos_1based: int
    missed_cleavages: int
    sequons: Tuple[SequonSite, ...]
    modifications: Mapping[str, int]
    modification_summary: str
    peptide_neutral_mass: float


@dataclass(frozen=True)
class PrecursorObservation:
    """Observed MS1 precursor signal for glycopeptide screening."""

    sample_id: str
    observed_mz: float
    intensity: float = 0.0
    known_charge: Optional[int] = None
    retention_time_sec: Optional[float] = None
    ms2_spectrum: Optional[SpectrumData] = None


@dataclass(frozen=True)
class MS2GlycopeptideEvidence:
    """Tested oxonium and Y-ion fragment evidence from an MS/MS spectrum."""

    evaluated: bool
    msms_supported: bool
    oxonium_ions_matched: Tuple[str, ...]
    peptide_y_ions_matched: Tuple[str, ...]
    status_summary: str


@dataclass(frozen=True)
class GlycopeptideCandidateMatch:
    """Complete result record for a candidate glycopeptide match."""

    sample_id: str
    protein_id: str
    protein_sequence: str
    peptide_sequence: str
    peptide_start_1based: int
    peptide_end_1based: int
    candidate_glycosylation_sites: Tuple[str, ...]  # e.g. ("N4 (Prot N297: NST)",)
    sequon_count_in_peptide: int
    glycan_composition: str
    peptide_neutral_mass: float
    attached_glycan_moiety_mass: float
    theoretical_glycopeptide_neutral_mass: float
    observed_precursor_mz: float
    theoretical_precursor_mz: float
    charge_state: int
    adduct: str
    signed_error_da: float
    abs_error_da: float
    signed_error_ppm: float
    abs_error_ppm: float
    missed_cleavages: int
    modification_info: str
    retention_time_sec: Optional[float]
    msms_evidence_status: str
    msms_oxonium_matched: Tuple[str, ...]
    msms_y_ions_matched: Tuple[str, ...]
    confidence_category: GlycopeptideConfidenceCategory
    candidate_rank: int
    alternative_candidates: Tuple[str, ...]
    warnings: Tuple[str, ...]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "protein_id": self.protein_id,
            "protein_sequence": self.protein_sequence,
            "peptide_sequence": self.peptide_sequence,
            "peptide_start_1based": self.peptide_start_1based,
            "peptide_end_1based": self.peptide_end_1based,
            "candidate_glycosylation_sites": list(self.candidate_glycosylation_sites),
            "sequon_count_in_peptide": self.sequon_count_in_peptide,
            "glycan_composition": self.glycan_composition,
            "peptide_neutral_mass": round(self.peptide_neutral_mass, 6),
            "attached_glycan_moiety_mass": round(self.attached_glycan_moiety_mass, 6),
            "theoretical_glycopeptide_neutral_mass": round(self.theoretical_glycopeptide_neutral_mass, 6),
            "observed_precursor_mz": round(self.observed_precursor_mz, 6),
            "theoretical_precursor_mz": round(self.theoretical_precursor_mz, 6),
            "charge_state": self.charge_state,
            "adduct": self.adduct,
            "signed_error_da": round(self.signed_error_da, 6),
            "abs_error_da": round(self.abs_error_da, 6),
            "signed_error_ppm": round(self.signed_error_ppm, 4),
            "abs_error_ppm": round(self.abs_error_ppm, 4),
            "missed_cleavages": self.missed_cleavages,
            "modification_info": self.modification_info,
            "retention_time_sec": self.retention_time_sec,
            "msms_evidence_status": self.msms_evidence_status,
            "msms_oxonium_matched": list(self.msms_oxonium_matched),
            "msms_y_ions_matched": list(self.msms_y_ions_matched),
            "confidence_category": self.confidence_category.value,
            "candidate_rank": self.candidate_rank,
            "alternative_candidates": list(self.alternative_candidates),
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class GlycopeptideSearchSettings:
    """User-configurable parameters for Stage 6 glycopeptide screening."""

    enzyme: ProteaseEnzyme = ProteaseEnzyme.TRYPSIN
    max_missed_cleavages: int = 2
    min_peptide_length: int = 4
    max_peptide_length: int = 50
    require_n_glycan_sequon: bool = True
    fixed_carbamidomethyl_cys: bool = True
    variable_oxidation_met: bool = False
    variable_deamidation_nq: bool = False
    tolerance: float = 20.0
    tolerance_unit: MassToleranceUnit = MassToleranceUnit.PPM
    adduct_codes: Tuple[str, ...] = ("[M+H]+", "[M+Na]+")
    charge_states: Tuple[int, ...] = (1, 2, 3)
    sialic_derivatization: str = "NONE"
    ms2_fragment_tolerance_da: float = 0.05
    candidate_limit_per_precursor: int = 10


@dataclass(frozen=True)
class GlycopeptideScreeningResult:
    """Complete result of a glycopeptide screening run."""

    sample_id: str
    proteins_screened: int
    peptides_generated: int
    sequon_peptides_count: int
    glycan_compositions_screened: int
    precursors_evaluated: int
    matched_precursors_count: int
    total_candidate_matches: int
    matches: Tuple[GlycopeptideCandidateMatch, ...]
    warnings: Tuple[str, ...]

    def to_summary_dict(self) -> Dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "proteins_screened": self.proteins_screened,
            "peptides_generated": self.peptides_generated,
            "sequon_peptides_count": self.sequon_peptides_count,
            "glycan_compositions_screened": self.glycan_compositions_screened,
            "precursors_evaluated": self.precursors_evaluated,
            "matched_precursors_count": self.matched_precursors_count,
            "total_candidate_matches": self.total_candidate_matches,
            "matches": [m.to_dict() for m in self.matches],
            "warnings": list(self.warnings),
        }


# Canonical N-linked glycosylation sequon regex: N followed by any AA except P, followed by S or T
# Uses lookahead `(?=...)` so overlapping sequons (e.g. NNST or NSTNKT) are all detected!
_N_SEQUON_RE = re.compile(r"(?=(N[^P][ST]))")

# Diagnostic B-type carbohydrate oxonium ions (protonated [M+H]+ m/z)
DIAGNOSTIC_OXONIUM_IONS: Mapping[str, float] = {
    "HexNAc+ (204.0866)": MONOSACCHARIDE_RESIDUES["HexNAc"].residue_mass + CATION_H_MASS,
    "HexNAc-H2O+ (186.0761)": MONOSACCHARIDE_RESIDUES["HexNAc"].residue_mass - 18.0105646865 + CATION_H_MASS,
    "Hex+ (163.0601)": MONOSACCHARIDE_RESIDUES["Hex"].residue_mass + CATION_H_MASS,
    "HexHexNAc+ (366.1395)": (
        MONOSACCHARIDE_RESIDUES["Hex"].residue_mass
        + MONOSACCHARIDE_RESIDUES["HexNAc"].residue_mass
        + CATION_H_MASS
    ),
    "Neu5Ac+ (292.1027)": MONOSACCHARIDE_RESIDUES["Neu5Ac"].residue_mass + CATION_H_MASS,
    "Neu5Ac-H2O+ (274.0921)": MONOSACCHARIDE_RESIDUES["Neu5Ac"].residue_mass - 18.0105646865 + CATION_H_MASS,
}


def parse_fasta_sequences(fasta_input: Union[str, bytes, Path]) -> List[FastaProteinRecord]:
    """Parse one or more protein sequences from FASTA text or raw sequence lines."""
    if isinstance(fasta_input, Path):
        text, _ = decode_text_bytes(fasta_input.read_bytes())
    elif isinstance(fasta_input, bytes):
        text, _ = decode_text_bytes(fasta_input)
    else:
        text = fasta_input

    cleaned = text.strip()
    if not cleaned:
        raise ValueError("FASTA / protein sequence input cannot be empty.")

    records: List[FastaProteinRecord] = []
    if not cleaned.startswith(">"):
        # Treat as a single raw amino acid sequence
        seq = re.sub(r"[\s0-9]+", "", cleaned).upper()
        for ch in seq:
            if ch not in AMINO_ACID_RESIDUES:
                raise ValueError(f"Invalid amino acid character '{ch}' in protein sequence.")
        return [FastaProteinRecord(protein_id="USER_PROTEIN_1", description="User-supplied sequence", sequence=seq)]

    curr_id = ""
    curr_desc = ""
    curr_seq_parts: List[str] = []

    for raw_ln in cleaned.splitlines():
        ln = raw_ln.strip()
        if not ln:
            continue
        if ln.startswith(">"):
            if curr_id and curr_seq_parts:
                full_seq = "".join(curr_seq_parts).upper()
                records.append(FastaProteinRecord(protein_id=curr_id, description=curr_desc, sequence=full_seq))
            header = ln[1:].strip()
            parts = header.split(None, 1)
            curr_id = parts[0] if parts else f"PROTEIN_{len(records) + 1}"
            curr_desc = parts[1] if len(parts) > 1 else ""
            curr_seq_parts = []
        else:
            seq_clean = re.sub(r"[\s0-9*]+", "", ln).upper()
            for ch in seq_clean:
                if ch not in AMINO_ACID_RESIDUES:
                    raise ValueError(f"Invalid amino acid character '{ch}' in protein '{curr_id}'.")
            curr_seq_parts.append(seq_clean)

    if curr_id and curr_seq_parts:
        full_seq = "".join(curr_seq_parts).upper()
        records.append(FastaProteinRecord(protein_id=curr_id, description=curr_desc, sequence=full_seq))

    if not records:
        raise ValueError("No valid protein sequences found in FASTA input.")

    return records


def find_n_glycosylation_sequons(
    peptide_sequence: str,
    peptide_start_in_protein_1based: Optional[int] = 1,
    protein_sequence: Optional[str] = None,
) -> Tuple[SequonSite, ...]:
    """Identify all canonical `N-X-S/T` (`X != P`) sequons in a peptide sequence.

    Also checks the C-terminal boundary against `protein_sequence` if an Asn at the last
    or second-to-last residue of the peptide completes an `N-X-S/T` motif across the cleavage site.
    """
    pep = peptide_sequence.strip().upper()
    sites: List[SequonSite] = []

    # Check within peptide sequence
    for m in _N_SEQUON_RE.finditer(pep):
        idx0 = m.start(1)
        motif = m.group(1)
        pep_pos1 = idx0 + 1
        prot_pos1 = (
            (peptide_start_in_protein_1based + idx0)
            if peptide_start_in_protein_1based is not None
            else None
        )
        sites.append(
            SequonSite(
                motif=motif,
                peptide_position_1based=pep_pos1,
                protein_position_1based=prot_pos1,
            )
        )

    # Check if Asn at len(pep)-2 or len(pep)-1 forms a sequon in the intact parent protein
    if protein_sequence and peptide_start_in_protein_1based is not None:
        prot = protein_sequence.strip().upper()
        for offset in (len(pep) - 2, len(pep) - 1):
            if offset < 0:
                continue
            if pep[offset] == "N":
                prot_idx0 = (peptide_start_in_protein_1based - 1) + offset
                if prot_idx0 + 2 < len(prot):
                    tri = prot[prot_idx0 : prot_idx0 + 3]
                    if tri[1] != "P" and tri[2] in ("S", "T"):
                        pep_pos1 = offset + 1
                        if not any(s.peptide_position_1based == pep_pos1 for s in sites):
                            sites.append(
                                SequonSite(
                                    motif=tri,
                                    peptide_position_1based=pep_pos1,
                                    protein_position_1based=prot_idx0 + 1,
                                )
                            )

    return tuple(sites)


def _get_cleavage_indices(sequence: str, enzyme: ProteaseEnzyme) -> List[int]:
    """Return 0-based cleavage positions (index after the cleaved residue)."""
    n = len(sequence)
    if enzyme == ProteaseEnzyme.NO_CLEAVAGE:
        return [0, n]

    cuts = [0]
    for i in range(n - 1):
        aa = sequence[i]
        next_aa = sequence[i + 1]
        if enzyme == ProteaseEnzyme.TRYPSIN:
            if aa in ("K", "R") and next_aa != "P":
                cuts.append(i + 1)
        elif enzyme == ProteaseEnzyme.TRYPSIN_P:
            if aa in ("K", "R"):
                cuts.append(i + 1)
        elif enzyme == ProteaseEnzyme.CHYMOTRYPSIN:
            if aa in ("F", "Y", "W", "L") and next_aa != "P":
                cuts.append(i + 1)
        elif enzyme == ProteaseEnzyme.GLU_C:
            if aa in ("E", "D") and next_aa != "P":
                cuts.append(i + 1)
        elif enzyme == ProteaseEnzyme.LYS_C:
            if aa == "K" and next_aa != "P":
                cuts.append(i + 1)
    cuts.append(n)
    return cuts


def digest_protein_in_silico(
    protein: FastaProteinRecord,
    settings: Optional[GlycopeptideSearchSettings] = None,
) -> List[DigestedPeptideCandidate]:
    """Perform in-silico protease digestion and generate modified/unmodified peptide candidates."""
    if settings is None:
        settings = GlycopeptideSearchSettings()

    seq = protein.sequence
    cuts = _get_cleavage_indices(seq, settings.enzyme)
    num_segments = len(cuts) - 1

    candidates: List[DigestedPeptideCandidate] = []
    max_mc = 0 if settings.enzyme == ProteaseEnzyme.NO_CLEAVAGE else settings.max_missed_cleavages

    for i in range(num_segments):
        for mc in range(0, max_mc + 1):
            j = i + 1 + mc
            if j >= len(cuts):
                break
            start0 = cuts[i]
            end0 = cuts[j]
            pep_seq = seq[start0:end0]
            if not (settings.min_peptide_length <= len(pep_seq) <= settings.max_peptide_length):
                continue

            sequons = find_n_glycosylation_sequons(
                peptide_sequence=pep_seq,
                peptide_start_in_protein_1based=start0 + 1,
                protein_sequence=seq,
            )
            if settings.require_n_glycan_sequon and not sequons:
                continue

            # Enumerate variable modification states (0 or all-oxidized M, 0 or 1 deamidation)
            ox_options = [0]
            if settings.variable_oxidation_met and "M" in pep_seq:
                ox_options.extend(range(1, pep_seq.count("M") + 1))

            deam_options = [0]
            if settings.variable_deamidation_nq and ("N" in pep_seq or "Q" in pep_seq):
                deam_options.append(1)

            for ox_cnt, deam_cnt in itertools.product(ox_options, deam_options):
                mods: Dict[str, int] = {}
                if ox_cnt > 0:
                    mods["Oxidation_M"] = ox_cnt
                if deam_cnt > 0:
                    mods["Deamidation_NQ"] = deam_cnt

                pep_mass_res = calculate_peptide_neutral_mass(
                    sequence=pep_seq,
                    modifications=mods,
                    fixed_carbamidomethyl_cys=settings.fixed_carbamidomethyl_cys,
                )
                mod_summary = (
                    ", ".join(pep_mass_res.modifications_applied)
                    if pep_mass_res.modifications_applied
                    else "Unmodified"
                )
                candidates.append(
                    DigestedPeptideCandidate(
                        protein_id=protein.protein_id,
                        protein_sequence=seq,
                        peptide_sequence=pep_seq,
                        start_pos_1based=start0 + 1,
                        end_pos_1based=end0,
                        missed_cleavages=mc,
                        sequons=sequons,
                        modifications=mods,
                        modification_summary=mod_summary,
                        peptide_neutral_mass=pep_mass_res.neutral_monoisotopic_mass,
                    )
                )

    return candidates


def evaluate_ms2_glycopeptide_fragments(
    ms2_spectrum: Optional[SpectrumData],
    peptide_neutral_mass: float,
    tolerance_da: float = 0.05,
) -> MS2GlycopeptideEvidence:
    """Evaluate MS/MS spectrum for diagnostic B-type oxonium ions and Y0/Y1 peptide-core ions.

    Only reports `msms_supported=True` when BOTH at least one carbohydrate oxonium ion
    (e.g. m/z 204.0866 HexNAc+, 366.1395 HexHexNAc+) AND at least one Y-type peptide/core
    ion (Y0 intact peptide [M+H]+ or Y1 [Peptide+HexNAc+H]+ in z=1, 2, or 3) are observed
    within `tolerance_da`.
    """
    if ms2_spectrum is None or not ms2_spectrum.peaks:
        return MS2GlycopeptideEvidence(
            evaluated=False,
            msms_supported=False,
            oxonium_ions_matched=(),
            peptide_y_ions_matched=(),
            status_summary="Precursor-only MS1 screening (no MS/MS spectrum provided)",
        )

    obs_mzs = [p.mz for p in ms2_spectrum.peaks if p.intensity > 0]
    if not obs_mzs:
        return MS2GlycopeptideEvidence(
            evaluated=False,
            msms_supported=False,
            oxonium_ions_matched=(),
            peptide_y_ions_matched=(),
            status_summary="MS/MS spectrum contained no positive-intensity fragment peaks",
        )

    matched_oxonium: List[str] = []
    for ox_name, ox_mz in DIAGNOSTIC_OXONIUM_IONS.items():
        if any(abs(m - ox_mz) <= tolerance_da for m in obs_mzs):
            matched_oxonium.append(ox_name)

    # Check Y0 (peptide), Y1 (peptide + HexNAc), Y2 (peptide + 2 HexNAc), Y1F (peptide + HexNAc + dHex)
    hexnac_mass = MONOSACCHARIDE_RESIDUES["HexNAc"].residue_mass
    dhex_mass = MONOSACCHARIDE_RESIDUES["dHex"].residue_mass
    y_neutral_masses = {
        "Y0 (Peptide)": peptide_neutral_mass,
        "Y1 (Peptide+HexNAc)": peptide_neutral_mass + hexnac_mass,
        "Y1F (Peptide+HexNAc+Fuc)": peptide_neutral_mass + hexnac_mass + dhex_mass,
        "Y2 (Peptide+HexNAc2)": peptide_neutral_mass + 2 * hexnac_mass,
    }

    matched_y: List[str] = []
    for y_label, y_neutral in y_neutral_masses.items():
        for z in (1, 2, 3):
            theo_y_mz = (y_neutral + z * CATION_H_MASS) / z
            if any(abs(m - theo_y_mz) <= tolerance_da for m in obs_mzs):
                matched_y.append(f"{y_label} [{z}+] ({theo_y_mz:.4f})")

    is_supported = len(matched_oxonium) >= 1 and len(matched_y) >= 1
    if is_supported:
        summary = (
            f"MS/MS supported ({len(matched_oxonium)} oxonium ion(s), {len(matched_y)} Y-type peptide-core ion(s) matched)"
        )
    else:
        summary = (
            f"Insufficient MS/MS fragment evidence (oxonium={len(matched_oxonium)}, Y-ions={len(matched_y)}); "
            "retained as MS1 precursor match only"
        )

    return MS2GlycopeptideEvidence(
        evaluated=True,
        msms_supported=is_supported,
        oxonium_ions_matched=tuple(matched_oxonium),
        peptide_y_ions_matched=tuple(matched_y),
        status_summary=summary,
    )


def screen_glycopeptides(
    precursors: Sequence[PrecursorObservation],
    fasta_or_proteins: Optional[Union[str, bytes, Path, Sequence[FastaProteinRecord]]] = None,
    peptide_sequences: Optional[Sequence[str]] = None,
    glycan_compositions: Optional[Union[Sequence[str], GlycanLibrary]] = None,
    settings: Optional[GlycopeptideSearchSettings] = None,
    user_candidate_table: Optional[Sequence[Mapping[str, Any]]] = None,
) -> GlycopeptideScreeningResult:
    """Screen observed precursor ions against in-silico digested proteins/peptides + N-glycan compositions."""
    if settings is None:
        settings = GlycopeptideSearchSettings()

    if not precursors:
        raise ValueError("At least one observed precursor m/z must be provided for glycopeptide screening.")

    sample_id = precursors[0].sample_id or "glycopeptide_sample"
    warnings: List[str] = [
        "Scientific note: An N-X-S/T sequon represents a potential N-glycosylation site, not proof of site occupancy.",
        "Attached glycan moiety masses are calculated as anhydro-residue sums WITHOUT the released glycan's terminal reducing-end H2O.",
    ]

    # 1. Build list of DigestedPeptideCandidates from FASTA and/or user-supplied peptides
    proteins: List[FastaProteinRecord] = []
    if fasta_or_proteins is not None:
        if isinstance(fasta_or_proteins, (str, bytes, Path)):
            proteins.extend(parse_fasta_sequences(fasta_or_proteins))
        else:
            proteins.extend(fasta_or_proteins)

    if peptide_sequences:
        for idx, p_seq in enumerate(peptide_sequences, start=1):
            clean_p = p_seq.strip().upper()
            if clean_p:
                proteins.append(
                    FastaProteinRecord(
                        protein_id=f"USER_PEPTIDE_{idx}",
                        description="User-supplied peptide sequence",
                        sequence=clean_p,
                    )
                )

    digested_peptides: List[DigestedPeptideCandidate] = []
    for prot in proteins:
        if prot.protein_id.startswith("USER_PEPTIDE_"):
            # Do not cleave user-supplied discrete peptide sequences unless requested
            no_cleave_settings = GlycopeptideSearchSettings(
                **{**settings.__dict__, "enzyme": ProteaseEnzyme.NO_CLEAVAGE}
            )
            digested_peptides.extend(digest_protein_in_silico(prot, no_cleave_settings))
        else:
            digested_peptides.extend(digest_protein_in_silico(prot, settings))

    # 2. Build list of attached glycan moiety masses
    if isinstance(glycan_compositions, GlycanLibrary):
        comp_strings = [e.composition.split("[")[0] for e in glycan_compositions.entries]
    elif glycan_compositions:
        comp_strings = list(glycan_compositions)
    else:
        # Default core set of common human N-glycans if none specified
        comp_strings = [
            "Hex5HexNAc2",
            "Hex3HexNAc4dHex1",
            "Hex4HexNAc4dHex1",
            "Hex5HexNAc4dHex1",
            "Hex5HexNAc4Neu5Ac1",
            "Hex5HexNAc4dHex1Neu5Ac1",
            "Hex5HexNAc4Neu5Ac2",
            "Hex5HexNAc4dHex1Neu5Ac2",
        ]

    # Deduplicate compositions and precompute attached glycan moiety masses
    attached_glycans: List[Tuple[str, float]] = []
    seen_comps = set()
    for c_str in comp_strings:
        counts = parse_glycan_composition(c_str)
        att_res = calculate_attached_glycan_moiety_mass(
            counts,
            sialic_derivatization=settings.sialic_derivatization,
        )
        if att_res.canonical_composition not in seen_comps:
            seen_comps.add(att_res.canonical_composition)
            attached_glycans.append((att_res.canonical_composition, att_res.neutral_monoisotopic_mass))

    # 3. Expand adduct and charge state combinations
    adduct_charge_defs: List[Tuple[str, int]] = []
    seen_ac = set()
    for add_code in settings.adduct_codes:
        base_def = resolve_adduct_definition(add_code)
        for z in settings.charge_states:
            abs_z = abs(z)
            if abs_z == abs(base_def.charge):
                key = (base_def.code, base_def.charge)
            else:
                try:
                    scaled = build_homogeneous_adduct_for_charge(add_code, abs_z)
                    key = (scaled.code, scaled.charge)
                except ValueError:
                    continue
            if key not in seen_ac:
                seen_ac.add(key)
                adduct_charge_defs.append(key)

    all_matches: List[GlycopeptideCandidateMatch] = []
    matched_precursors_count = 0

    for prec in precursors:
        prec_hits: List[Tuple[float, GlycopeptideCandidateMatch]] = []

        for pep in digested_peptides:
            site_labels = (
                tuple(
                    f"N{s.peptide_position_1based} (Prot N{s.protein_position_1based}: {s.motif})"
                    if s.protein_position_1based is not None
                    else f"N{s.peptide_position_1based} ({s.motif})"
                    for s in pep.sequons
                )
                if pep.sequons
                else ("No canonical N-X-S/T sequon",)
            )

            for gly_comp, att_gly_mass in attached_glycans:
                gp_neutral = pep.peptide_neutral_mass + att_gly_mass

                for adduct_code, signed_z in adduct_charge_defs:
                    if prec.known_charge is not None and abs(prec.known_charge) != abs(signed_z):
                        continue

                    ion_res = calculate_mz_from_neutral_mass(
                        neutral_mass=gp_neutral,
                        adduct_code=adduct_code,
                    )
                    if not is_within_mass_tolerance(
                        observed=prec.observed_mz,
                        theoretical=ion_res.theoretical_mz,
                        tolerance=settings.tolerance,
                        unit=settings.tolerance_unit,
                    ):
                        continue

                    err = evaluate_mass_error(
                        observed=prec.observed_mz,
                        theoretical=ion_res.theoretical_mz,
                    )
                    ms2_ev = evaluate_ms2_glycopeptide_fragments(
                        ms2_spectrum=prec.ms2_spectrum,
                        peptide_neutral_mass=pep.peptide_neutral_mass,
                        tolerance_da=settings.ms2_fragment_tolerance_da,
                    )

                    # Assign strict evidence confidence category
                    if ms2_ev.msms_supported:
                        if len(pep.sequons) == 1:
                            conf_cat = GlycopeptideConfidenceCategory.SITE_LOCALIZED
                        else:
                            conf_cat = GlycopeptideConfidenceCategory.MSMS_SUPPORTED
                    elif pep.sequons:
                        conf_cat = GlycopeptideConfidenceCategory.MS1_COMPOSITION_CANDIDATE
                    else:
                        conf_cat = GlycopeptideConfidenceCategory.MS1_MASS_MATCH

                    cand_warns: List[str] = []
                    if len(pep.sequons) > 1:
                        cand_warns.append(
                            f"Peptide contains {len(pep.sequons)} candidate sequons ({', '.join(site_labels)}); "
                            "site occupancy and site-specific localization require site-determining ETD/EThcD/HCD fragments."
                        )
                    elif len(pep.sequons) == 1 and not ms2_ev.msms_supported:
                        cand_warns.append(
                            "Precursor m/z match only; peptide sequence and glycan composition require MS/MS confirmation."
                        )

                    match_draft = GlycopeptideCandidateMatch(
                        sample_id=prec.sample_id or sample_id,
                        protein_id=pep.protein_id,
                        protein_sequence=pep.protein_sequence,
                        peptide_sequence=pep.peptide_sequence,
                        peptide_start_1based=pep.start_pos_1based,
                        peptide_end_1based=pep.end_pos_1based,
                        candidate_glycosylation_sites=site_labels,
                        sequon_count_in_peptide=len(pep.sequons),
                        glycan_composition=gly_comp,
                        peptide_neutral_mass=pep.peptide_neutral_mass,
                        attached_glycan_moiety_mass=att_gly_mass,
                        theoretical_glycopeptide_neutral_mass=gp_neutral,
                        observed_precursor_mz=prec.observed_mz,
                        theoretical_precursor_mz=ion_res.theoretical_mz,
                        charge_state=signed_z,
                        adduct=adduct_code,
                        signed_error_da=err.signed_error_da,
                        abs_error_da=err.abs_error_da,
                        signed_error_ppm=err.signed_error_ppm,
                        abs_error_ppm=err.abs_error_ppm,
                        missed_cleavages=pep.missed_cleavages,
                        modification_info=pep.modification_summary,
                        retention_time_sec=prec.retention_time_sec,
                        msms_evidence_status=ms2_ev.status_summary,
                        msms_oxonium_matched=ms2_ev.oxonium_ions_matched,
                        msms_y_ions_matched=ms2_ev.peptide_y_ions_matched,
                        confidence_category=conf_cat,
                        candidate_rank=1,
                        alternative_candidates=(),
                        warnings=tuple(cand_warns),
                    )
                    prec_hits.append((err.abs_error_ppm, match_draft))

        # Also evaluate optional user-provided glycopeptide candidate table rows if supplied
        if user_candidate_table:
            for row_i, row in enumerate(user_candidate_table, start=1):
                row_pep = str(row.get("peptide_sequence", "")).strip().upper()
                row_gly = str(row.get("glycan_composition", "")).strip()
                if not row_pep or not row_gly:
                    continue
                is_ref_val = bool(row.get("reference_validated", False))
                gp_neutral, pep_res, gly_res, _ = calculate_glycopeptide_neutral_mass(
                    peptide_sequence=row_pep,
                    glycan_composition=row_gly,
                    fixed_carbamidomethyl_cys=settings.fixed_carbamidomethyl_cys,
                )
                sequons = find_n_glycosylation_sequons(row_pep, 1, row_pep)
                site_labels = tuple(f"N{s.peptide_position_1based} ({s.motif})" for s in sequons)
                for adduct_code, signed_z in adduct_charge_defs:
                    if prec.known_charge is not None and abs(prec.known_charge) != abs(signed_z):
                        continue
                    ion_res = calculate_mz_from_neutral_mass(gp_neutral, adduct_code)
                    if not is_within_mass_tolerance(
                        prec.observed_mz, ion_res.theoretical_mz, settings.tolerance, settings.tolerance_unit
                    ):
                        continue
                    err = evaluate_mass_error(prec.observed_mz, ion_res.theoretical_mz)
                    conf_cat = (
                        GlycopeptideConfidenceCategory.REFERENCE_VALIDATED
                        if is_ref_val
                        else GlycopeptideConfidenceCategory.MS1_COMPOSITION_CANDIDATE
                    )
                    match_draft = GlycopeptideCandidateMatch(
                        sample_id=prec.sample_id or sample_id,
                        protein_id=str(row.get("protein_id", f"USER_TABLE_{row_i}")),
                        protein_sequence=row_pep,
                        peptide_sequence=row_pep,
                        peptide_start_1based=1,
                        peptide_end_1based=len(row_pep),
                        candidate_glycosylation_sites=site_labels,
                        sequon_count_in_peptide=len(sequons),
                        glycan_composition=gly_res.canonical_composition,
                        peptide_neutral_mass=pep_res.neutral_monoisotopic_mass,
                        attached_glycan_moiety_mass=gly_res.neutral_monoisotopic_mass,
                        theoretical_glycopeptide_neutral_mass=gp_neutral,
                        observed_precursor_mz=prec.observed_mz,
                        theoretical_precursor_mz=ion_res.theoretical_mz,
                        charge_state=signed_z,
                        adduct=adduct_code,
                        signed_error_da=err.signed_error_da,
                        abs_error_da=err.abs_error_da,
                        signed_error_ppm=err.signed_error_ppm,
                        abs_error_ppm=err.abs_error_ppm,
                        missed_cleavages=0,
                        modification_info=", ".join(pep_res.modifications_applied) or "Unmodified",
                        retention_time_sec=prec.retention_time_sec,
                        msms_evidence_status="Matched against user-supplied candidate table",
                        msms_oxonium_matched=(),
                        msms_y_ions_matched=(),
                        confidence_category=conf_cat,
                        candidate_rank=1,
                        alternative_candidates=(),
                        warnings=(),
                    )
                    prec_hits.append((err.abs_error_ppm, match_draft))

        if not prec_hits:
            continue

        matched_precursors_count += 1
        # Sort by MS/MS support first, then abs_error_ppm ascending
        prec_hits.sort(
            key=lambda x: (
                0 if x[1].confidence_category in (
                    GlycopeptideConfidenceCategory.SITE_LOCALIZED,
                    GlycopeptideConfidenceCategory.MSMS_SUPPORTED,
                    GlycopeptideConfidenceCategory.REFERENCE_VALIDATED,
                ) else 1,
                x[0],
                x[1].peptide_sequence,
                x[1].glycan_composition,
            )
        )

        retained = prec_hits[: max(1, settings.candidate_limit_per_precursor)]
        for rank_i, (_, m_obj) in enumerate(retained, start=1):
            alt_list = tuple(
                f"{other.peptide_sequence}+{other.glycan_composition} ({other.adduct}, {other.signed_error_ppm:+.2f} ppm)"
                for _, other in prec_hits
                if (other.peptide_sequence, other.glycan_composition, other.adduct)
                != (m_obj.peptide_sequence, m_obj.glycan_composition, m_obj.adduct)
            )
            extra_warns = list(m_obj.warnings)
            if len(prec_hits) > 1:
                extra_warns.append(
                    f"Ambiguous precursor match: {len(prec_hits)} peptide+glycan combinations match m/z {prec.observed_mz:.4f} within tolerance."
                )
            all_matches.append(
                GlycopeptideCandidateMatch(
                    **{
                        **m_obj.__dict__,
                        "candidate_rank": rank_i,
                        "alternative_candidates": alt_list,
                        "warnings": tuple(extra_warns),
                    }
                )
            )

    sequon_peptides = sum(1 for p in digested_peptides if p.sequons)

    return GlycopeptideScreeningResult(
        sample_id=sample_id,
        proteins_screened=len(proteins),
        peptides_generated=len(digested_peptides),
        sequon_peptides_count=sequon_peptides,
        glycan_compositions_screened=len(attached_glycans),
        precursors_evaluated=len(precursors),
        matched_precursors_count=matched_precursors_count,
        total_candidate_matches=len(all_matches),
        matches=tuple(all_matches),
        warnings=tuple(warnings),
    )

"""GlycoQuest-like MS1 glycan search, candidate ranking, and ambiguity reporting engine.

Implements:
- Mode A: Targeted glycan search (by compositions, library IDs, theoretical neutral masses, or expected m/z values).
- Mode B: Library-based glycan search against a reference/composition library.
- Multi-adduct and multi-charge matching using Stage 1 authoritative calculators.
- Transparent heuristic candidate ranking vs. optional isotope dot-product (`idotp`) scoring.
- Explicit one-to-many ambiguity reporting, adduct-vs-composition overlap warnings,
  and strict evidence level classification (`MS1_EXACT_MASS_ONLY`, `MS1_COMPOSITION_MATCH`,
  `MS1_ISOTOPE_SUPPORTED`, `REFERENCE_LIBRARY_COMPOSITION_MATCH`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np

from app.analysis.peak_processing import (
    PeakPreprocessingConfig,
    ProcessedSpectrumResult,
    preprocess_spectrum,
)
from app.analysis.spectrum_corrections import find_isotope_satellite_indices
from app.core.ion_calculator import (
    build_homogeneous_adduct_for_charge,
    calculate_mz_from_neutral_mass,
    resolve_adduct_definition,
)
from app.core.isotope_scoring import (
    generate_theoretical_isotope_envelope,
    score_envelope_against_peaks,
)
from app.core.mass_calculator import (
    calculate_released_glycan_mass,
    format_glycan_composition,
    parse_glycan_composition,
)
from app.core.mass_error import (
    MassToleranceUnit,
    evaluate_mass_error,
    is_within_mass_tolerance,
)
from app.core.scientific_constants import (
    ADDUCT_DEFINITIONS,
    IonPolarity,
    ValidationStatus,
)
from app.core.validation import validate_search_parameters
from app.libraries.composition_generator import (
    GlycanClass,
    GlycanLibrary,
    GlycanLibraryEntry,
    StructuralAnnotationStatus,
    classify_n_glycan_composition,
)
from app.parsers.base_models import SpectrumData, SpectrumPeak


class SearchMode(str, Enum):
    """MS1 glycan search mode."""

    TARGETED = "targeted"  # Mode A
    LIBRARY = "library"    # Mode B


class GlycanEvidenceLevel(str, Enum):
    """Explicit identification evidence level for MS1 glycan search candidates.

    Strictly distinguishes a theoretical mass/composition match from isotope-supported
    or reference-curated matches. Never implies full topological or linkage resolution
    unless linkage-specific derivatization was used.
    """

    MS1_MASS_MATCH_ONLY = "MS1_MASS_MATCH_ONLY"
    MS1_COMPOSITION_CANDIDATE = "MS1_COMPOSITION_CANDIDATE"
    MS1_ISOTOPE_SUPPORTED_CANDIDATE = "MS1_ISOTOPE_SUPPORTED_CANDIDATE"
    MS1_REFERENCE_COMPOSITION_CANDIDATE = "MS1_REFERENCE_COMPOSITION_CANDIDATE"


@dataclass(frozen=True)
class TargetedQueryItem:
    """Single target specification for Mode A (Targeted MS1 Search)."""

    query_label: str
    target_type: str  # "composition", "library_id", "neutral_mass", "expected_mz"
    composition: Optional[str] = None
    library_id: Optional[str] = None
    neutral_mass: Optional[float] = None
    expected_mz: Optional[float] = None
    adduct_code: Optional[str] = None
    charge: Optional[int] = None


@dataclass(frozen=True)
class MS1SearchSettings:
    """User-configurable parameters for Mode A and Mode B MS1 glycan searches."""

    search_mode: SearchMode = SearchMode.LIBRARY
    tolerance: float = 20.0
    tolerance_unit: MassToleranceUnit = MassToleranceUnit.PPM
    polarity: IonPolarity = IonPolarity.POSITIVE
    adduct_codes: Tuple[str, ...] = ("[M+Na]+",)
    charge_states: Tuple[int, ...] = (1,)
    labeling: str = "NONE"
    min_intensity: float = 0.0
    min_relative_intensity_pct: float = 0.0
    candidate_limit_per_peak: int = 10
    enable_isotope_scoring: bool = True
    preprocessing: PeakPreprocessingConfig = field(default_factory=PeakPreprocessingConfig)
    # --- match-quality filters (all optional; 0 / False = disabled) ---
    min_snr: float = 0.0  # drop peaks below this S/N (instrument S/N, else local-median estimate)
    min_score: float = 0.0  # drop candidates whose composite heuristic score is below this (0-100)
    min_idotp: float = 0.0  # drop candidates whose AVAILABLE idotp is below this (0-1)
    require_isotope_support: bool = False  # drop candidates without a computable idotp
    exclude_isotope_satellites: bool = False  # do not annotate M+1/M+2 satellites as separate species

    def to_dict(self) -> Dict[str, Any]:
        return {
            "search_mode": self.search_mode.value,
            "tolerance": self.tolerance,
            "tolerance_unit": self.tolerance_unit.value,
            "polarity": self.polarity.value,
            "adduct_codes": list(self.adduct_codes),
            "charge_states": list(self.charge_states),
            "labeling": self.labeling,
            "min_intensity": self.min_intensity,
            "min_relative_intensity_pct": self.min_relative_intensity_pct,
            "candidate_limit_per_peak": self.candidate_limit_per_peak,
            "enable_isotope_scoring": self.enable_isotope_scoring,
            "min_snr": self.min_snr,
            "min_score": self.min_score,
            "min_idotp": self.min_idotp,
            "require_isotope_support": self.require_isotope_support,
            "exclude_isotope_satellites": self.exclude_isotope_satellites,
        }


def build_preprocessing_config(settings: "MS1SearchSettings") -> PeakPreprocessingConfig:
    """Merge search-level thresholds with the nested `PeakPreprocessingConfig`."""
    pre = settings.preprocessing
    snr_candidates = [v for v in (settings.min_snr, pre.min_snr) if v is not None and v > 0]
    return PeakPreprocessingConfig(
        min_intensity=max(settings.min_intensity, pre.min_intensity),
        min_relative_intensity_pct=max(settings.min_relative_intensity_pct, pre.min_relative_intensity_pct),
        min_mz=pre.min_mz,
        max_mz=pre.max_mz,
        min_snr=max(snr_candidates) if snr_candidates else None,
        apply_profile_baseline_subtraction=pre.apply_profile_baseline_subtraction,
        apply_profile_peak_detection=pre.apply_profile_peak_detection,
        integrate_profile_peak_area=pre.integrate_profile_peak_area,
        baseline_window_points=pre.baseline_window_points,
        baseline_percentile=pre.baseline_percentile,
    )


@dataclass(frozen=True)
class ScoringBreakdown:
    """Transparent heuristic ranking score breakdown (0 to 100 scale, NOT a probability)."""

    mass_accuracy_score: float  # [0..100]: 100 * (1 - abs_error / tolerance)
    intensity_score: float  # [0..100]: log/linear relative intensity contribution
    isotope_similarity_score: Optional[float]  # [0..100]: 100 * idotp if available, else None
    composite_heuristic_score: float  # Weighted heuristic ranking score in [0..100]
    scoring_formula_note: str


@dataclass(frozen=True)
class MS1CandidateMatch:
    """Single candidate match for an observed MS1 peak (row in the search results table)."""

    sample_id: str
    peak_index: int
    observed_mz: float
    observed_intensity: float
    relative_intensity_pct: float
    observed_area: Optional[float]
    observed_snr: Optional[float]
    library_id: str
    glycan_composition: str
    candidate_class: str
    antennary_hypothesis: str
    theoretical_neutral_mass: float
    theoretical_mz: float
    assigned_adduct: str
    charge_state: int
    signed_error_da: float
    abs_error_da: float
    signed_error_ppm: float
    abs_error_ppm: float
    candidate_rank: int
    mass_error_rank: int
    scoring_components: ScoringBreakdown
    isotope_score_idotp: Optional[float]
    isotope_score_status: str
    evidence_level: GlycanEvidenceLevel
    structural_annotation_status: str
    total_candidates_for_peak: int
    alternative_candidates: Tuple[str, ...]
    is_ambiguous: bool
    warnings: Tuple[str, ...]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "peak_index": self.peak_index,
            "observed_mz": round(self.observed_mz, 6),
            "observed_intensity": round(self.observed_intensity, 4),
            "relative_intensity_pct": round(self.relative_intensity_pct, 4),
            "observed_area": round(self.observed_area, 4) if self.observed_area is not None else None,
            "observed_snr": round(self.observed_snr, 2) if self.observed_snr is not None else None,
            "library_id": self.library_id,
            "glycan_composition": self.glycan_composition,
            "candidate_class": self.candidate_class,
            "antennary_hypothesis": self.antennary_hypothesis,
            "theoretical_neutral_mass": round(self.theoretical_neutral_mass, 6),
            "theoretical_mz": round(self.theoretical_mz, 6),
            "assigned_adduct": self.assigned_adduct,
            "charge_state": self.charge_state,
            "signed_error_da": round(self.signed_error_da, 6),
            "abs_error_da": round(self.abs_error_da, 6),
            "signed_error_ppm": round(self.signed_error_ppm, 4),
            "abs_error_ppm": round(self.abs_error_ppm, 4),
            "candidate_rank": self.candidate_rank,
            "mass_error_rank": self.mass_error_rank,
            "composite_heuristic_score": round(self.scoring_components.composite_heuristic_score, 2),
            "scoring_components": {
                "mass_accuracy_score": round(self.scoring_components.mass_accuracy_score, 2),
                "intensity_score": round(self.scoring_components.intensity_score, 2),
                "isotope_similarity_score": (
                    round(self.scoring_components.isotope_similarity_score, 2)
                    if self.scoring_components.isotope_similarity_score is not None
                    else None
                ),
                "composite_heuristic_score": round(self.scoring_components.composite_heuristic_score, 2),
                "scoring_formula_note": self.scoring_components.scoring_formula_note,
            },
            "isotope_score_idotp": self.isotope_score_idotp,
            "isotope_score_status": self.isotope_score_status,
            "evidence_level": self.evidence_level.value,
            "structural_annotation_status": self.structural_annotation_status,
            "total_candidates_for_peak": self.total_candidates_for_peak,
            "alternative_candidates": list(self.alternative_candidates),
            "is_ambiguous": self.is_ambiguous,
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class MS1SearchResult:
    """Complete result of an MS1 glycan search (Mode A or Mode B)."""

    sample_id: str
    filename: str
    search_mode: SearchMode
    settings: MS1SearchSettings
    library_name: str
    library_version: str
    total_imported_peaks: int
    eligible_peaks_searched: int
    matched_peaks_count: int
    unmatched_peaks_count: int
    total_candidate_matches: int
    ambiguous_peaks_count: int
    matches: Tuple[MS1CandidateMatch, ...]  # All ranked candidate rows across matched peaks
    top_matches_per_peak: Tuple[MS1CandidateMatch, ...]  # Rank-1 candidate per matched peak
    unmatched_peaks: Tuple[SpectrumPeak, ...]
    warnings: Tuple[str, ...]

    def to_summary_dict(self) -> Dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "filename": self.filename,
            "search_mode": self.search_mode.value,
            "settings": self.settings.to_dict(),
            "library_name": self.library_name,
            "library_version": self.library_version,
            "total_imported_peaks": self.total_imported_peaks,
            "eligible_peaks_searched": self.eligible_peaks_searched,
            "matched_peaks_count": self.matched_peaks_count,
            "unmatched_peaks_count": self.unmatched_peaks_count,
            "total_candidate_matches": self.total_candidate_matches,
            "ambiguous_peaks_count": self.ambiguous_peaks_count,
            "matches": [m.to_dict() for m in self.matches],
            "top_matches_per_peak": [m.to_dict() for m in self.top_matches_per_peak],
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class _TheoreticalIonCandidate:
    """Precomputed theoretical ion candidate for rapid MS1 comparison."""

    library_id: str
    composition: str
    glycan_class: str
    antennary_hypothesis: str
    neutral_mass: float
    theoretical_mz: float
    adduct_code: str
    charge: int
    ion_formula: Optional[Mapping[str, int]]
    structural_annotation_status: str
    isomer_group_note: str
    entry_warnings: Tuple[str, ...]


def _expand_adduct_charge_pairs(
    adduct_codes: Sequence[str],
    charge_states: Sequence[int],
    polarity: IonPolarity,
) -> List[Tuple[str, int]]:
    """Expand user-selected adduct codes and charge states into valid (adduct_code, signed_charge) pairs."""
    pairs: List[Tuple[str, int]] = []
    seen = set()

    for code in adduct_codes:
        base_def = resolve_adduct_definition(code)
        if base_def.polarity != polarity:
            continue
        base_abs_z = abs(base_def.charge)
        if base_abs_z in charge_states or not charge_states:
            key = (base_def.code, base_def.charge)
            if key not in seen:
                pairs.append(key)
                seen.add(key)

        # Also scale homogeneous single-charged carriers if user requested higher charge states
        for req_z in charge_states:
            if req_z != base_abs_z:
                try:
                    scaled_def = build_homogeneous_adduct_for_charge(code, abs(req_z))
                    key = (scaled_def.code, scaled_def.charge)
                    if key not in seen:
                        pairs.append(key)
                        seen.add(key)
                except ValueError:
                    continue
    return pairs


_ION_CACHE_MAX = 8
_ION_CACHE: "Dict[Tuple[Any, ...], Tuple[GlycanLibrary, List[_TheoreticalIonCandidate]]]" = {}


def _build_library_theoretical_ions(
    library: GlycanLibrary,
    settings: MS1SearchSettings,
) -> List[_TheoreticalIonCandidate]:
    """Cached wrapper: theoretical ions only depend on the library and adduct/charge/label settings.

    Re-running a search (e.g. after tolerance / S/N / calibration changes) therefore skips the
    expensive per-entry mass + formula computation.
    """
    if library.library_name == "Targeted Subset":
        return _build_library_theoretical_ions_uncached(library, settings)
    key = (
        id(library),
        library.entry_count,
        tuple(settings.adduct_codes),
        tuple(settings.charge_states),
        settings.polarity.value,
        settings.labeling.upper().strip(),
    )
    hit = _ION_CACHE.get(key)
    if hit is not None and hit[0] is library:
        return hit[1]
    ions = _build_library_theoretical_ions_uncached(library, settings)
    if len(_ION_CACHE) >= _ION_CACHE_MAX:
        _ION_CACHE.pop(next(iter(_ION_CACHE)))
    _ION_CACHE[key] = (library, ions)
    return ions


def _build_library_theoretical_ions_uncached(
    library: GlycanLibrary,
    settings: MS1SearchSettings,
) -> List[_TheoreticalIonCandidate]:
    """Generate all eligible theoretical ion candidates from a `GlycanLibrary`."""
    adduct_pairs = _expand_adduct_charge_pairs(
        adduct_codes=settings.adduct_codes,
        charge_states=settings.charge_states,
        polarity=settings.polarity,
    )
    target_label = settings.labeling.upper().strip()
    candidates: List[_TheoreticalIonCandidate] = []

    for entry in library.entries:
        # If the user selected a different labeling chemistry than the library's precomputed default,
        # recalculate neutral mass and formula using Stage 1 `calculate_released_glycan_mass`.
        if entry.labeling.upper() == target_label:
            neutral_mass = entry.theoretical_neutral_mass
            neutral_formula = entry.elemental_formula
        else:
            mass_res = calculate_released_glycan_mass(
                composition=entry.residue_counts,
                labeling=target_label,
                alpha23_count=entry.alpha23_sialic_count,
                alpha26_count=entry.alpha26_sialic_count,
            )
            neutral_mass = mass_res.neutral_monoisotopic_mass
            neutral_formula = mass_res.elemental_formula

        for adduct_code, _ in adduct_pairs:
            ion_res = calculate_mz_from_neutral_mass(
                neutral_mass=neutral_mass,
                adduct_code=adduct_code,
                neutral_formula=neutral_formula,
            )
            candidates.append(
                _TheoreticalIonCandidate(
                    library_id=entry.library_id,
                    composition=entry.composition,
                    glycan_class=entry.glycan_class.value,
                    antennary_hypothesis=entry.antennary_hypothesis,
                    neutral_mass=neutral_mass,
                    theoretical_mz=ion_res.theoretical_mz,
                    adduct_code=ion_res.adduct_code,
                    charge=ion_res.charge,
                    ion_formula=ion_res.ion_formula,
                    structural_annotation_status=entry.structural_annotation_status.value,
                    isomer_group_note=entry.isomer_group_note,
                    entry_warnings=entry.warnings,
                )
            )
    return candidates


def _build_targeted_theoretical_ions(
    targets: Sequence[TargetedQueryItem],
    settings: MS1SearchSettings,
    reference_library: Optional[GlycanLibrary] = None,
) -> List[_TheoreticalIonCandidate]:
    """Generate theoretical ion candidates from Mode A user-supplied targets."""
    adduct_pairs = _expand_adduct_charge_pairs(
        adduct_codes=settings.adduct_codes,
        charge_states=settings.charge_states,
        polarity=settings.polarity,
    )
    target_label = settings.labeling.upper().strip()
    lib_by_id: Dict[str, GlycanLibraryEntry] = {}
    if reference_library is not None:
        lib_by_id = {e.library_id.upper(): e for e in reference_library.entries}

    candidates: List[_TheoreticalIonCandidate] = []

    for idx, t in enumerate(targets, start=1):
        t_type = t.target_type.lower().strip()

        if t_type == "library_id":
            lid = (t.library_id or t.query_label).strip().upper()
            if lid not in lib_by_id:
                raise ValueError(f"Target library_id '{lid}' not found in selected library.")
            entry = lib_by_id[lid]
            sub_lib = GlycanLibrary(
                library_name="Targeted Subset",
                library_version="1.0.0",
                created_at_utc="",
                labeling=entry.labeling,
                constraints_used=None,
                entries=(entry,),
                provenance_summary="Targeted lookup",
            )
            candidates.extend(_build_library_theoretical_ions(sub_lib, settings))

        elif t_type == "composition":
            comp_str = t.composition or t.query_label
            counts = parse_glycan_composition(comp_str)
            mass_res = calculate_released_glycan_mass(counts, labeling=target_label)
            _, g_class, hyp, note, warns = classify_n_glycan_composition(counts)
            for adduct_code, _ in adduct_pairs:
                ion_res = calculate_mz_from_neutral_mass(
                    neutral_mass=mass_res.neutral_monoisotopic_mass,
                    adduct_code=adduct_code,
                    neutral_formula=mass_res.elemental_formula,
                )
                candidates.append(
                    _TheoreticalIonCandidate(
                        library_id=f"TARGET_COMP_{idx:03d}",
                        composition=mass_res.canonical_composition,
                        glycan_class=g_class.value,
                        antennary_hypothesis=hyp,
                        neutral_mass=mass_res.neutral_monoisotopic_mass,
                        theoretical_mz=ion_res.theoretical_mz,
                        adduct_code=ion_res.adduct_code,
                        charge=ion_res.charge,
                        ion_formula=ion_res.ion_formula,
                        structural_annotation_status=StructuralAnnotationStatus.COMPOSITION_CANDIDATE.value,
                        isomer_group_note=note,
                        entry_warnings=tuple(warns),
                    )
                )

        elif t_type == "neutral_mass":
            if t.neutral_mass is None or t.neutral_mass <= 0:
                raise ValueError(f"Targeted neutral_mass must be positive for target '{t.query_label}'")
            for adduct_code, _ in adduct_pairs:
                ion_res = calculate_mz_from_neutral_mass(
                    neutral_mass=t.neutral_mass,
                    adduct_code=adduct_code,
                )
                candidates.append(
                    _TheoreticalIonCandidate(
                        library_id=f"TARGET_MASS_{idx:03d}",
                        composition=t.query_label or f"NeutralMass_{t.neutral_mass:.4f}",
                        glycan_class=GlycanClass.UNCLASSIFIED.value,
                        antennary_hypothesis="Targeted neutral mass query (no composition specified)",
                        neutral_mass=t.neutral_mass,
                        theoretical_mz=ion_res.theoretical_mz,
                        adduct_code=ion_res.adduct_code,
                        charge=ion_res.charge,
                        ion_formula=None,
                        structural_annotation_status=StructuralAnnotationStatus.COMPOSITION_CANDIDATE.value,
                        isomer_group_note="Target neutral mass match only.",
                        entry_warnings=("Target specified by neutral mass without elemental composition.",),
                    )
                )

        elif t_type == "expected_mz":
            if t.expected_mz is None or t.expected_mz <= 0:
                raise ValueError(f"Targeted expected_mz must be positive for target '{t.query_label}'")
            adduct_used = t.adduct_code or settings.adduct_codes[0]
            charge_used = t.charge or settings.charge_states[0]
            candidates.append(
                _TheoreticalIonCandidate(
                    library_id=f"TARGET_MZ_{idx:03d}",
                    composition=t.query_label or f"TargetMz_{t.expected_mz:.4f}",
                    glycan_class=GlycanClass.UNCLASSIFIED.value,
                    antennary_hypothesis="Targeted m/z query (no composition specified)",
                    neutral_mass=0.0,
                    theoretical_mz=t.expected_mz,
                    adduct_code=adduct_used,
                    charge=charge_used,
                    ion_formula=None,
                    structural_annotation_status=StructuralAnnotationStatus.COMPOSITION_CANDIDATE.value,
                    isomer_group_note="Target m/z match only.",
                    entry_warnings=("Target specified by expected m/z only.",),
                )
            )
        else:
            raise ValueError(f"Unsupported targeted search item type: '{t.target_type}'")

    return candidates


def _compute_heuristic_score(
    abs_error: float,
    tolerance: float,
    relative_intensity_pct: float,
    idotp: Optional[float],
) -> ScoringBreakdown:
    """Calculate a transparent heuristic ranking score (0..100).

    Explicitly documented as a heuristic ranking metric, NOT a calibrated identification probability.
    Components:
      - `mass_accuracy_score` = max(0, 100 * (1 - abs_error / tolerance))
      - `intensity_score` = min(100, 20 * log10(1 + relative_intensity_pct * 999))  [maps 0.1%..100% into ~40..100]
      - `isotope_similarity_score` = 100 * idotp (when available)
    Weights:
      - When isotope score is available: 0.65 * mass_accuracy + 0.25 * isotope_similarity + 0.10 * intensity
      - When isotope score is unavailable: 0.85 * mass_accuracy + 0.15 * intensity
    """
    import math

    tol_safe = max(tolerance, 1e-9)
    mass_score = max(0.0, min(100.0, 100.0 * (1.0 - (abs_error / tol_safe))))
    int_score = max(0.0, min(100.0, 25.0 * math.log10(1.0 + max(0.0, relative_intensity_pct) * 99.0)))

    if idotp is not None:
        iso_score = max(0.0, min(100.0, 100.0 * idotp))
        composite = 0.65 * mass_score + 0.25 * iso_score + 0.10 * int_score
        note = "Heuristic rank score = 0.65*MassAccuracy + 0.25*IsotopeSimilarity + 0.10*Intensity (not a probability)"
    else:
        iso_score = None
        composite = 0.85 * mass_score + 0.15 * int_score
        note = "Heuristic rank score = 0.85*MassAccuracy + 0.15*Intensity; isotope similarity unavailable (not a probability)"

    return ScoringBreakdown(
        mass_accuracy_score=mass_score,
        intensity_score=int_score,
        isotope_similarity_score=iso_score,
        composite_heuristic_score=composite,
        scoring_formula_note=note,
    )


def _execute_matching_core(
    spectrum: SpectrumData,
    theoretical_ions: Sequence[_TheoreticalIonCandidate],
    settings: MS1SearchSettings,
    library_name: str,
    library_version: str,
) -> MS1SearchResult:
    """Core matching loop shared by Mode A (targeted) and Mode B (library) MS1 searches."""
    val_rep = validate_search_parameters(
        tolerance=settings.tolerance,
        tolerance_unit=settings.tolerance_unit,
        polarity=settings.polarity,
        adduct_codes=list(settings.adduct_codes),
        charge_states=list(settings.charge_states),
        labeling=settings.labeling,
    )
    if not val_rep.is_valid:
        raise ValueError(f"Invalid MS1 search parameters: {'; '.join(val_rep.errors)}")

    # Sync PeakPreprocessingConfig thresholds (incl. S/N) with search settings
    prep_cfg = build_preprocessing_config(settings)
    proc: ProcessedSpectrumResult = preprocess_spectrum(spectrum, config=prep_cfg)
    warnings: List[str] = list(val_rep.warnings) + list(proc.warnings)

    # Isotope scoring pool: all peaks that passed m/z + S/N filters (sorted by m/z), so that weak
    # M+1/M+2 satellites are available even when the relative-intensity cut removed them from annotation.
    pool_peaks = sorted(proc.isotope_pool_peaks or proc.eligible_peaks, key=lambda pk: pk.mz)
    pool_mz = np.array([p.mz for p in pool_peaks], dtype=float)
    pool_int = np.array([p.intensity for p in pool_peaks], dtype=float)

    # Optionally skip peaks that are M+1/M+2 satellites of a lower peak for annotation purposes.
    satellite_peak_ids: set = set()
    if settings.exclude_isotope_satellites and len(proc.eligible_peaks) > 1:
        e_mz = np.array([p.mz for p in proc.eligible_peaks], dtype=float)
        e_int = np.array([p.intensity for p in proc.eligible_peaks], dtype=float)
        tol_ppm = settings.tolerance if settings.tolerance_unit == MassToleranceUnit.PPM else 25.0
        charges = tuple(sorted({abs(int(z)) for z in settings.charge_states} | {1}))
        # Satellites are identified among eligible peaks; they remain in the isotope pool for idotp
        sat_idx = find_isotope_satellite_indices(e_mz, e_int, tolerance_ppm=tol_ppm, charges=charges)
        satellite_peak_ids = {id(proc.eligible_peaks[i]) for i in sat_idx}
        if satellite_peak_ids:
            warnings.append(
                f"{len(satellite_peak_ids)} peaks were treated as isotope satellites (M+1/M+2) and not annotated separately."
            )

    # Theoretical ions sorted by m/z for O(log n) window lookup per observed peak
    ion_order = sorted(range(len(theoretical_ions)), key=lambda k: theoretical_ions[k].theoretical_mz)
    ion_mz_sorted = np.array([theoretical_ions[k].theoretical_mz for k in ion_order], dtype=float)
    envelope_cache: Dict[int, Tuple[Any, ...]] = {}
    use_ppm = settings.tolerance_unit == MassToleranceUnit.PPM
    n_isotopes = 4

    all_matches: List[MS1CandidateMatch] = []
    top_matches: List[MS1CandidateMatch] = []
    unmatched_peaks: List[SpectrumPeak] = []
    ambiguous_peaks_count = 0

    for p_idx, (peak, rel_int_pct) in enumerate(
        zip(proc.eligible_peaks, proc.relative_intensities_pct), start=1
    ):
        if id(peak) in satellite_peak_ids:
            continue

        raw_hits: List[Tuple[_TheoreticalIonCandidate, Any, Any, ScoringBreakdown]] = []

        # Binary-search window around the observed m/z (padded 1.5x; exact test follows)
        if use_ppm:
            half_w = peak.mz * settings.tolerance * 1e-6 * 1.5 + 1e-9
        else:
            half_w = settings.tolerance * 1.5 + 1e-9
        lo = int(np.searchsorted(ion_mz_sorted, peak.mz - half_w, side="left"))
        hi = int(np.searchsorted(ion_mz_sorted, peak.mz + half_w, side="right"))
        cand_indices = sorted(ion_order[lo:hi])  # restore original (deterministic) ordering

        for cand_idx in cand_indices:
            cand = theoretical_ions[cand_idx]
            if not is_within_mass_tolerance(
                observed=peak.mz,
                theoretical=cand.theoretical_mz,
                tolerance=settings.tolerance,
                unit=settings.tolerance_unit,
            ):
                continue

            err = evaluate_mass_error(observed=peak.mz, theoretical=cand.theoretical_mz)

            # Evaluate optional isotope envelope similarity (idotp); envelopes are cached per candidate
            iso_res = None
            if settings.enable_isotope_scoring and cand.ion_formula is not None:
                env = envelope_cache.get(cand_idx)
                if env is None:
                    env = generate_theoretical_isotope_envelope(
                        ion_formula=cand.ion_formula,
                        monoisotopic_mz=cand.theoretical_mz,
                        charge=cand.charge,
                        num_peaks=n_isotopes,
                    )
                    envelope_cache[cand_idx] = env
                iso_res = score_envelope_against_peaks(
                    theo_envelope=env,
                    obs_mz_sorted=pool_mz,
                    obs_int_sorted=pool_int,
                    tolerance=settings.tolerance,
                    tolerance_unit=settings.tolerance_unit,
                )

            idotp_val = iso_res.idotp if (iso_res and iso_res.available) else None

            # --- optional isotope-based quality filters ---
            if settings.require_isotope_support and idotp_val is None:
                continue
            if settings.min_idotp > 0 and idotp_val is not None and idotp_val < settings.min_idotp:
                continue

            abs_err_metric = (
                err.abs_error_ppm
                if settings.tolerance_unit == MassToleranceUnit.PPM
                else err.abs_error_da
            )
            breakdown = _compute_heuristic_score(
                abs_error=abs_err_metric,
                tolerance=settings.tolerance,
                relative_intensity_pct=rel_int_pct,
                idotp=idotp_val,
            )
            if settings.min_score > 0 and breakdown.composite_heuristic_score < settings.min_score:
                continue
            raw_hits.append((cand, err, iso_res, breakdown))


        if not raw_hits:
            unmatched_peaks.append(peak)
            continue

        # Compute mass_error_rank (1 = smallest abs_error_ppm)
        sorted_by_mass_err = sorted(raw_hits, key=lambda x: (x[1].abs_error_ppm, x[0].library_id))
        mass_rank_map: Dict[Tuple[str, str], int] = {
            (item[0].library_id, item[0].adduct_code): r_i
            for r_i, item in enumerate(sorted_by_mass_err, start=1)
        }

        # Sort by composite_heuristic_score descending (within 1e-4), then prefer canonical
        # human (non-Neu5Gc) compositions when two candidates are exact elemental isomers,
        # then abs_error_ppm ascending
        raw_hits.sort(
            key=lambda x: (
                -round(x[3].composite_heuristic_score, 4),
                "Neu5Gc" in x[0].composition,
                round(x[1].abs_error_ppm, 5),
                x[0].library_id,
            )
        )

        total_cands = len(raw_hits)
        is_ambig = total_cands > 1
        if is_ambig:
            ambiguous_peaks_count += 1

        # Detect adduct confusion vs composition ambiguity among candidates matching the same peak
        distinct_adducts = {h[0].adduct_code for h in raw_hits}
        distinct_comps = {h[0].composition for h in raw_hits}

        retained_hits = raw_hits[: max(1, settings.candidate_limit_per_peak)]

        for rank_idx, (cand, err, iso_res, breakdown) in enumerate(retained_hits, start=1):
            alt_labels = tuple(
                f"{other[0].composition} ({other[0].adduct_code}, {other[1].signed_error_ppm:+.2f} ppm)"
                for other in raw_hits
                if (other[0].library_id, other[0].adduct_code) != (cand.library_id, cand.adduct_code)
            )

            cand_warnings: List[str] = list(cand.entry_warnings)
            cand_warnings.append(cand.isomer_group_note)

            if is_ambig:
                if len(distinct_adducts) > 1 and len(distinct_comps) > 1:
                    cand_warnings.append(
                        f"Ambiguous peak assignment ({total_cands} candidates): near-isobaric overlap between different glycan compositions and adduct species ({', '.join(sorted(distinct_adducts))})."
                    )
                elif len(distinct_comps) > 1:
                    cand_warnings.append(
                        f"Multiple compositions ({total_cands}) match within tolerance; MS1 mass cannot resolve structural or isobaric ambiguity."
                    )
                else:
                    cand_warnings.append(
                        f"Multiple library records ({total_cands}) share this composition/mass."
                    )

            if iso_res is not None:
                iso_status = iso_res.status_message
                idotp_out = iso_res.idotp
            else:
                iso_status = (
                    "Insufficient isotope information for isotope scoring"
                    if settings.enable_isotope_scoring
                    else "Isotope scoring disabled by user"
                )
                idotp_out = None

            # Determine explicit evidence level
            if cand.ion_formula is None and cand.neutral_mass == 0.0:
                ev_level = GlycanEvidenceLevel.MS1_MASS_MATCH_ONLY
            elif cand.structural_annotation_status == StructuralAnnotationStatus.STRUCTURALLY_VALIDATED_REFERENCE.value:
                ev_level = GlycanEvidenceLevel.MS1_REFERENCE_COMPOSITION_CANDIDATE
            elif idotp_out is not None and idotp_out >= 0.85:
                ev_level = GlycanEvidenceLevel.MS1_ISOTOPE_SUPPORTED_CANDIDATE
            else:
                ev_level = GlycanEvidenceLevel.MS1_COMPOSITION_CANDIDATE

            match_obj = MS1CandidateMatch(
                sample_id=spectrum.sample_id,
                peak_index=p_idx,
                observed_mz=peak.mz,
                observed_intensity=peak.intensity,
                relative_intensity_pct=rel_int_pct,
                observed_area=peak.area,
                observed_snr=peak.snr,
                library_id=cand.library_id,
                glycan_composition=cand.composition,
                candidate_class=cand.glycan_class,
                antennary_hypothesis=cand.antennary_hypothesis,
                theoretical_neutral_mass=cand.neutral_mass,
                theoretical_mz=cand.theoretical_mz,
                assigned_adduct=cand.adduct_code,
                charge_state=cand.charge,
                signed_error_da=err.signed_error_da,
                abs_error_da=err.abs_error_da,
                signed_error_ppm=err.signed_error_ppm,
                abs_error_ppm=err.abs_error_ppm,
                candidate_rank=rank_idx,
                mass_error_rank=mass_rank_map[(cand.library_id, cand.adduct_code)],
                scoring_components=breakdown,
                isotope_score_idotp=idotp_out,
                isotope_score_status=iso_status,
                evidence_level=ev_level,
                structural_annotation_status=cand.structural_annotation_status,
                total_candidates_for_peak=total_cands,
                alternative_candidates=alt_labels,
                is_ambiguous=is_ambig,
                warnings=tuple(cand_warnings),
            )
            all_matches.append(match_obj)
            if rank_idx == 1:
                top_matches.append(match_obj)

    if not all_matches:
        warnings.append("No candidates found within the selected tolerance")

    warnings.append(
        "Scientific note: MS1 mass matches represent candidate monosaccharide compositions, not definitive structural identifications or calibrated probabilities."
    )

    return MS1SearchResult(
        sample_id=spectrum.sample_id,
        filename=spectrum.filename,
        search_mode=settings.search_mode,
        settings=settings,
        library_name=library_name,
        library_version=library_version,
        total_imported_peaks=spectrum.valid_points_count,
        eligible_peaks_searched=len(proc.eligible_peaks),
        matched_peaks_count=len(top_matches),
        unmatched_peaks_count=len(unmatched_peaks),
        total_candidate_matches=len(all_matches),
        ambiguous_peaks_count=ambiguous_peaks_count,
        matches=tuple(all_matches),
        top_matches_per_peak=tuple(top_matches),
        unmatched_peaks=tuple(unmatched_peaks),
        warnings=tuple(warnings),
    )


def search_ms1_library(
    spectrum: SpectrumData,
    library: GlycanLibrary,
    settings: Optional[MS1SearchSettings] = None,
) -> MS1SearchResult:
    """Mode B — Library-based MS1 glycan search against all eligible candidates in `library`."""
    if settings is None:
        settings = MS1SearchSettings(search_mode=SearchMode.LIBRARY)
    elif settings.search_mode != SearchMode.LIBRARY:
        settings = MS1SearchSettings(
            **{**settings.__dict__, "search_mode": SearchMode.LIBRARY}
        )

    theo_ions = _build_library_theoretical_ions(library=library, settings=settings)
    return _execute_matching_core(
        spectrum=spectrum,
        theoretical_ions=theo_ions,
        settings=settings,
        library_name=library.library_name,
        library_version=library.library_version,
    )


def search_ms1_targeted(
    spectrum: SpectrumData,
    targets: Sequence[TargetedQueryItem],
    settings: Optional[MS1SearchSettings] = None,
    reference_library: Optional[GlycanLibrary] = None,
) -> MS1SearchResult:
    """Mode A — Targeted MS1 glycan search for user-specified compositions, IDs, neutral masses, or m/z values."""
    if not targets:
        raise ValueError("Targeted MS1 search requires at least one target item.")

    if settings is None:
        settings = MS1SearchSettings(search_mode=SearchMode.TARGETED)
    elif settings.search_mode != SearchMode.TARGETED:
        settings = MS1SearchSettings(
            **{**settings.__dict__, "search_mode": SearchMode.TARGETED}
        )

    theo_ions = _build_targeted_theoretical_ions(
        targets=targets,
        settings=settings,
        reference_library=reference_library,
    )
    lib_name = reference_library.library_name if reference_library else "User Targeted Query List"
    lib_ver = reference_library.library_version if reference_library else "1.0.0"

    return _execute_matching_core(
        spectrum=spectrum,
        theoretical_ions=theo_ions,
        settings=settings,
        library_name=lib_name,
        library_version=lib_ver,
    )

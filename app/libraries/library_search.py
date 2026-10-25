"""Offline search and filtering engine for N-glycan reference and composition libraries.

Supports:
- Composition search (exact or minimum/maximum residue counts).
- Mass-range search (neutral mass or adduct m/z).
- Exact mass search (with tolerance in ppm or Da).
- Glycan class filters (`high_mannose`, `hybrid`, `complex`, `paucimannose`).
- Fucosylation and sialylation filters.
- Labeling and adduct filters.
- Grouping of isomeric / isobaric entries sharing the same composition or mass.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

from app.core.ion_calculator import calculate_mz_from_neutral_mass
from app.core.mass_calculator import format_glycan_composition, parse_glycan_composition
from app.core.mass_error import MassToleranceUnit, evaluate_mass_error, is_within_mass_tolerance
from app.libraries.composition_generator import (
    GlycanClass,
    GlycanLibrary,
    GlycanLibraryEntry,
)


@dataclass(frozen=True)
class LibrarySearchFilter:
    """Configurable query parameters for searching a `GlycanLibrary`."""

    composition_query: Optional[str] = None
    min_neutral_mass: Optional[float] = None
    max_neutral_mass: Optional[float] = None
    target_mass: Optional[float] = None
    target_mass_is_mz: bool = False
    adduct_code: str = "[M+Na]+"
    tolerance: float = 20.0
    tolerance_unit: MassToleranceUnit = MassToleranceUnit.PPM
    glycan_classes: Optional[Sequence[GlycanClass]] = None
    labeling: Optional[str] = None
    fucosylated_only: Optional[bool] = None
    sialylated_only: Optional[bool] = None
    exclude_neu5gc: bool = False


@dataclass(frozen=True)
class LibraryQueryMatch:
    """Single library entry matched by a search query, including optional mass error."""

    entry: GlycanLibraryEntry
    queried_adduct: str
    theoretical_mz: float
    signed_error_da: Optional[float] = None
    abs_error_da: Optional[float] = None
    signed_error_ppm: Optional[float] = None
    abs_error_ppm: Optional[float] = None


def get_entry_adduct_mz(entry: GlycanLibraryEntry, adduct_code: str) -> float:
    """Return precomputed adduct m/z from entry or compute it on the fly using Stage 1 `ion_calculator`."""
    if adduct_code in entry.adduct_mzs:
        return entry.adduct_mzs[adduct_code]
    res = calculate_mz_from_neutral_mass(
        neutral_mass=entry.theoretical_neutral_mass,
        adduct_code=adduct_code,
        neutral_formula=entry.elemental_formula,
    )
    return res.theoretical_mz


def query_glycan_library(
    library: GlycanLibrary,
    search_filter: LibrarySearchFilter,
) -> List[LibraryQueryMatch]:
    """Search and filter a `GlycanLibrary` offline according to `search_filter`."""
    target_canon_comp: Optional[str] = None
    target_counts: Optional[Dict[str, int]] = None
    if search_filter.composition_query and search_filter.composition_query.strip():
        q_str = search_filter.composition_query.strip()
        try:
            target_counts = parse_glycan_composition(q_str)
            target_canon_comp = format_glycan_composition(target_counts)
        except ValueError:
            # Allow substring match on composition or library_id if not a full composition string
            target_canon_comp = q_str

    allowed_classes_set = set(search_filter.glycan_classes) if search_filter.glycan_classes else None
    label_filter = search_filter.labeling.upper().strip() if search_filter.labeling else None

    matches: List[LibraryQueryMatch] = []

    for entry in library.entries:
        # 1. Labeling filter
        if label_filter and entry.labeling.upper() != label_filter:
            continue

        # 2. Class filter
        if allowed_classes_set is not None and entry.glycan_class not in allowed_classes_set:
            continue

        # 3. Fucosylation / sialylation / Neu5Gc filters
        if search_filter.fucosylated_only is not None and entry.is_fucosylated != search_filter.fucosylated_only:
            continue
        if search_filter.sialylated_only is not None and entry.is_sialylated != search_filter.sialylated_only:
            continue
        if search_filter.exclude_neu5gc and entry.contains_neu5gc:
            continue

        # 4. Composition filter
        if target_counts is not None:
            if dict(entry.residue_counts) != target_counts:
                continue
        elif target_canon_comp is not None:
            if (
                target_canon_comp.lower() not in entry.composition.lower()
                and target_canon_comp.lower() not in entry.library_id.lower()
            ):
                continue

        # 5. Neutral mass range filter
        if search_filter.min_neutral_mass is not None and entry.theoretical_neutral_mass < search_filter.min_neutral_mass:
            continue
        if search_filter.max_neutral_mass is not None and entry.theoretical_neutral_mass > search_filter.max_neutral_mass:
            continue

        theo_mz = get_entry_adduct_mz(entry, search_filter.adduct_code)

        # 6. Exact target mass / m/z filter
        if search_filter.target_mass is not None:
            ref_mass = theo_mz if search_filter.target_mass_is_mz else entry.theoretical_neutral_mass
            if not is_within_mass_tolerance(
                observed=search_filter.target_mass,
                theoretical=ref_mass,
                tolerance=search_filter.tolerance,
                unit=search_filter.tolerance_unit,
            ):
                continue
            err = evaluate_mass_error(observed=search_filter.target_mass, theoretical=ref_mass)
            matches.append(
                LibraryQueryMatch(
                    entry=entry,
                    queried_adduct=search_filter.adduct_code,
                    theoretical_mz=theo_mz,
                    signed_error_da=err.signed_error_da,
                    abs_error_da=err.abs_error_da,
                    signed_error_ppm=err.signed_error_ppm,
                    abs_error_ppm=err.abs_error_ppm,
                )
            )
        else:
            matches.append(
                LibraryQueryMatch(
                    entry=entry,
                    queried_adduct=search_filter.adduct_code,
                    theoretical_mz=theo_mz,
                )
            )

    if search_filter.target_mass is not None:
        matches.sort(key=lambda m: (m.abs_error_ppm if m.abs_error_ppm is not None else float("inf"), m.entry.library_id))

    return matches

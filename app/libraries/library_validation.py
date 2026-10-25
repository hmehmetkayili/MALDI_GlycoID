"""Validation and duplicate-detection utilities for N-glycan libraries.

Verifies:
- Schema completeness and unique internal IDs.
- Theoretical neutral mass and adduct m/z consistency against Stage 1 calculators.
- Detection of exact duplicate compositions (retaining structural isomers as alternative
  candidates while flagging exact redundant records).
- Retention of source and license provenance on imported entries.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

from app.core.ion_calculator import calculate_mz_from_neutral_mass
from app.core.mass_calculator import calculate_released_glycan_mass
from app.libraries.composition_generator import (
    GlycanLibrary,
    GlycanLibraryEntry,
    StructuralAnnotationStatus,
)


@dataclass(frozen=True)
class DuplicateGroup:
    """Group of library entries sharing identical composition, labeling, and linkage counts."""

    composition_key: str
    library_ids: Tuple[str, ...]
    reference_ids: Tuple[str, ...]
    are_structural_isomers: bool
    description: str


@dataclass(frozen=True)
class LibraryValidationResult:
    """Comprehensive validation report for a `GlycanLibrary`."""

    is_valid: bool
    entry_count: int
    errors: Tuple[str, ...]
    warnings: Tuple[str, ...]
    duplicate_groups: Tuple[DuplicateGroup, ...]


def detect_library_duplicates(entries: Sequence[GlycanLibraryEntry]) -> Tuple[DuplicateGroup, ...]:
    """Detect entries sharing identical composition, labeling, and sialic linkage state.

    Distinguishes:
    - True redundant duplicates (identical composition and no distinct reference structure ID).
    - Validated structural isomers (identical composition and mass, but distinct `reference_id`
      or topology hypotheses retained intentionally as an isomeric composition group).
    """
    by_key: Dict[str, List[GlycanLibraryEntry]] = {}
    for entry in entries:
        key = f"{entry.composition}|{entry.labeling}|{entry.alpha23_sialic_count}|{entry.alpha26_sialic_count}"
        by_key.setdefault(key, []).append(entry)

    groups: List[DuplicateGroup] = []
    for key, group in by_key.items():
        if len(group) <= 1:
            continue
        lib_ids = tuple(e.library_id for e in group)
        ref_ids = tuple(e.reference_id or "NONE" for e in group)
        distinct_refs = len(set(r for r in ref_ids if r != "NONE")) == len(group)
        if distinct_refs:
            desc = (
                f"Isomeric composition group ({group[0].composition}): {len(group)} distinct reference structures "
                f"share the same composition and theoretical mass ({group[0].theoretical_neutral_mass:.4f} Da)."
            )
        else:
            desc = (
                f"Duplicate composition records detected for {group[0].composition} ({', '.join(lib_ids)}) "
                "without distinct structural reference IDs."
            )
        groups.append(
            DuplicateGroup(
                composition_key=key,
                library_ids=lib_ids,
                reference_ids=ref_ids,
                are_structural_isomers=distinct_refs,
                description=desc,
            )
        )
    return tuple(groups)


def validate_glycan_library(
    library: GlycanLibrary,
    mass_tolerance_da: float = 1e-4,
) -> LibraryValidationResult:
    """Validate all entries in a `GlycanLibrary` against Stage 1 scientific mass rules."""
    errors: List[str] = []
    warnings: List[str] = list(library.warnings)

    if not library.library_version.strip():
        errors.append("Library is missing a version string (`library_version`).")

    if not library.entries:
        errors.append("Library contains zero entries.")
        return LibraryValidationResult(False, 0, tuple(errors), tuple(warnings), ())

    seen_ids = set()
    for idx, entry in enumerate(library.entries, start=1):
        if not entry.library_id:
            errors.append(f"Entry #{idx} has an empty `library_id`.")
        elif entry.library_id in seen_ids:
            errors.append(f"Duplicate `library_id` '{entry.library_id}' found at entry #{idx}.")
        seen_ids.add(entry.library_id)

        if not entry.source_provenance or not entry.source_provenance.strip():
            errors.append(f"Entry '{entry.library_id}' is missing `source_provenance`.")

        if (
            entry.structural_annotation_status == StructuralAnnotationStatus.STRUCTURALLY_VALIDATED_REFERENCE
            and not entry.reference_id
        ):
            errors.append(
                f"Entry '{entry.library_id}' is marked STRUCTURALLY_VALIDATED_REFERENCE but lacks a `reference_id`."
            )

        # Verify neutral mass against Stage 1 `calculate_released_glycan_mass`
        try:
            calc_res = calculate_released_glycan_mass(
                composition=entry.residue_counts,
                labeling=entry.labeling,
                alpha23_count=entry.alpha23_sialic_count,
                alpha26_count=entry.alpha26_sialic_count,
            )
            diff = abs(calc_res.neutral_monoisotopic_mass - entry.theoretical_neutral_mass)
            if diff > mass_tolerance_da:
                errors.append(
                    f"Entry '{entry.library_id}' ({entry.composition}) neutral mass ({entry.theoretical_neutral_mass:.6f} Da) "
                    f"deviates from Stage 1 authoritative mass ({calc_res.neutral_monoisotopic_mass:.6f} Da) by {diff:.6f} Da."
                )

            # Verify stored adduct m/z values
            for adduct_code, stored_mz in entry.adduct_mzs.items():
                ion_res = calculate_mz_from_neutral_mass(
                    neutral_mass=calc_res.neutral_monoisotopic_mass,
                    adduct_code=adduct_code,
                )
                mz_diff = abs(ion_res.theoretical_mz - stored_mz)
                if mz_diff > mass_tolerance_da:
                    errors.append(
                        f"Entry '{entry.library_id}' adduct {adduct_code} m/z ({stored_mz:.6f}) "
                        f"deviates from Stage 1 calculation ({ion_res.theoretical_mz:.6f}) by {mz_diff:.6f} Da."
                    )
        except Exception as exc:
            errors.append(f"Entry '{entry.library_id}' failed Stage 1 mass verification: {exc}")

    dup_groups = detect_library_duplicates(library.entries)
    for dg in dup_groups:
        if dg.are_structural_isomers:
            warnings.append(dg.description)
        else:
            warnings.append(f"Redundant composition warning: {dg.description}")

    return LibraryValidationResult(
        is_valid=len(errors) == 0,
        entry_count=len(library.entries),
        errors=tuple(errors),
        warnings=tuple(warnings),
        duplicate_groups=dup_groups,
    )

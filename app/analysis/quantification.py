"""Glycan intensity, relative abundance, peak-area handling, and replicate quantification engine.

Strictly distinguishes and documents:
1. Raw observed peak intensity (`observed_intensity`, instrument counts / a.u.).
2. Integrated peak area (`observed_area`, available ONLY from profile spectra or pre-integrated peak tables; never fabricated on centroid-only lists).
3. Base-peak-normalized relative intensity (`relative_intensity_pct` = observed_intensity / max_denominator_intensity * 100), with an explicit `DenominatorSetPolicy`.
4. Total-annotated-glycan relative abundance (`relative_abundance_pct` = glycan_signal / sum_all_annotated_glycan_signals * 100).
5. Absolute concentration: NEVER claimed from uncalibrated MALDI MS1 intensities.

Also implements configurable composition-level aggregation rules (`GlycanAggregationPolicy`)
that never blindly sum different adducts, charge states, or ambiguous overlapping peaks
without explicit rule selection and full peak-level provenance tracking.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from app.analysis.glycan_search import MS1CandidateMatch, MS1SearchResult


class DenominatorSetPolicy(str, Enum):
    """Defines which peak set is used as the denominator for base-peak relative intensity."""

    ALL_ELIGIBLE_SPECTRUM_PEAKS = "all_eligible_spectrum_peaks"
    ANNOTATED_GLYCAN_PEAKS_ONLY = "annotated_glycan_peaks_only"


class GlycanAggregationPolicy(str, Enum):
    """Rule for summarizing peak-level signals to the glycan composition level."""

    NO_AGGREGATION_PEAK_LEVEL = "no_aggregation_peak_level"
    PRIMARY_ADDUCT_ONLY = "primary_adduct_only"
    HIGHEST_INTENSITY_PEAK_PER_COMPOSITION = "highest_intensity_peak_per_composition"
    EXPLICIT_SUM_DISTINCT_ADDUCTS = "explicit_sum_distinct_adducts"


@dataclass(frozen=True)
class PeakQuantificationRecord:
    """Peak-level quantitative record preserving original intensity, area status, and normalization."""

    sample_id: str
    peak_index: int
    observed_mz: float
    library_id: str
    glycan_composition: str
    candidate_class: str
    assigned_adduct: str
    charge_state: int
    is_ambiguous: bool
    observed_intensity: float
    observed_area: Optional[float]
    area_status_note: str
    base_peak_relative_intensity_pct: float
    denominator_policy: DenominatorSetPolicy
    denominator_max_intensity: float
    fraction_of_annotated_intensity_pct: float
    fraction_of_annotated_area_pct: Optional[float]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "peak_index": self.peak_index,
            "observed_mz": round(self.observed_mz, 6),
            "library_id": self.library_id,
            "glycan_composition": self.glycan_composition,
            "candidate_class": self.candidate_class,
            "assigned_adduct": self.assigned_adduct,
            "charge_state": self.charge_state,
            "is_ambiguous": self.is_ambiguous,
            "observed_intensity": round(self.observed_intensity, 4),
            "observed_area": round(self.observed_area, 4) if self.observed_area is not None else None,
            "area_status_note": self.area_status_note,
            "base_peak_relative_intensity_pct": round(self.base_peak_relative_intensity_pct, 4),
            "denominator_policy": self.denominator_policy.value,
            "denominator_max_intensity": round(self.denominator_max_intensity, 4),
            "fraction_of_annotated_intensity_pct": round(self.fraction_of_annotated_intensity_pct, 4),
            "fraction_of_annotated_area_pct": (
                round(self.fraction_of_annotated_area_pct, 4)
                if self.fraction_of_annotated_area_pct is not None
                else None
            ),
        }


@dataclass(frozen=True)
class GlycanAbundanceRecord:
    """Composition-level summarized abundance record with full peak-level provenance."""

    sample_id: str
    glycan_composition: str
    candidate_class: str
    aggregation_policy_used: GlycanAggregationPolicy
    primary_adduct_filter: Optional[str]
    aggregated_intensity: float
    aggregated_area: Optional[float]
    base_peak_relative_intensity_pct: float
    relative_abundance_intensity_pct: float  # % of total quantified glycan intensity in sample
    relative_abundance_area_pct: Optional[float]  # % of total quantified glycan area in sample (if available)
    contributing_peaks_count: int
    contributing_adducts: Tuple[str, ...]
    contributing_mzs: Tuple[float, ...]
    peak_provenance_history: Tuple[str, ...]
    contains_ambiguous_peak: bool
    warnings: Tuple[str, ...]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "glycan_composition": self.glycan_composition,
            "candidate_class": self.candidate_class,
            "aggregation_policy_used": self.aggregation_policy_used.value,
            "primary_adduct_filter": self.primary_adduct_filter,
            "aggregated_intensity": round(self.aggregated_intensity, 4),
            "aggregated_area": round(self.aggregated_area, 4) if self.aggregated_area is not None else None,
            "base_peak_relative_intensity_pct": round(self.base_peak_relative_intensity_pct, 4),
            "relative_abundance_intensity_pct": round(self.relative_abundance_intensity_pct, 4),
            "relative_abundance_area_pct": (
                round(self.relative_abundance_area_pct, 4)
                if self.relative_abundance_area_pct is not None
                else None
            ),
            "contributing_peaks_count": self.contributing_peaks_count,
            "contributing_adducts": list(self.contributing_adducts),
            "contributing_mzs": [round(m, 6) for m in self.contributing_mzs],
            "peak_provenance_history": list(self.peak_provenance_history),
            "contains_ambiguous_peak": self.contains_ambiguous_peak,
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class SampleQuantificationResult:
    """Complete quantitative output for a single sample."""

    sample_id: str
    denominator_policy: DenominatorSetPolicy
    aggregation_policy: GlycanAggregationPolicy
    primary_adduct: str
    exclude_ambiguous_from_aggregation: bool
    denominator_max_intensity: float
    total_annotated_intensity: float
    total_annotated_area: Optional[float]
    has_peak_area_data: bool
    peak_records: Tuple[PeakQuantificationRecord, ...]
    glycan_abundances: Tuple[GlycanAbundanceRecord, ...]
    warnings: Tuple[str, ...]

    def to_summary_dict(self) -> Dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "denominator_policy": self.denominator_policy.value,
            "aggregation_policy": self.aggregation_policy.value,
            "primary_adduct": self.primary_adduct,
            "exclude_ambiguous_from_aggregation": self.exclude_ambiguous_from_aggregation,
            "denominator_max_intensity": round(self.denominator_max_intensity, 4),
            "total_annotated_intensity": round(self.total_annotated_intensity, 4),
            "total_annotated_area": (
                round(self.total_annotated_area, 4) if self.total_annotated_area is not None else None
            ),
            "has_peak_area_data": self.has_peak_area_data,
            "peak_records": [r.to_dict() for r in self.peak_records],
            "glycan_abundances": [g.to_dict() for g in self.glycan_abundances],
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class ReplicateGlycanSummary:
    """Replicate statistics for a single glycan composition across multiple samples (no silent imputation)."""

    glycan_composition: str
    candidate_class: str
    total_replicates: int
    observed_replicates: int
    missing_replicates: int
    missingness_pct: float
    mean_relative_abundance_pct: Optional[float]  # Computed over observed replicates
    sd_relative_abundance_pct: Optional[float]  # Sample SD (ddof=1 when n>=2) over observed replicates
    cv_pct: Optional[float]  # (SD / Mean) * 100 when n>=2 and mean > 0
    mean_base_peak_intensity_pct: Optional[float]
    sample_values_relative_abundance_pct: Mapping[str, Optional[float]]  # None for missing samples!
    imputation_applied: str = "NONE (Missing values preserved as None/NaN)"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "glycan_composition": self.glycan_composition,
            "candidate_class": self.candidate_class,
            "total_replicates": self.total_replicates,
            "observed_replicates": self.observed_replicates,
            "missing_replicates": self.missing_replicates,
            "missingness_pct": round(self.missingness_pct, 2),
            "mean_relative_abundance_pct": (
                round(self.mean_relative_abundance_pct, 4)
                if self.mean_relative_abundance_pct is not None
                else None
            ),
            "sd_relative_abundance_pct": (
                round(self.sd_relative_abundance_pct, 4)
                if self.sd_relative_abundance_pct is not None
                else None
            ),
            "cv_pct": round(self.cv_pct, 2) if self.cv_pct is not None else None,
            "mean_base_peak_intensity_pct": (
                round(self.mean_base_peak_intensity_pct, 4)
                if self.mean_base_peak_intensity_pct is not None
                else None
            ),
            "sample_values_relative_abundance_pct": {
                k: (round(v, 4) if v is not None else None)
                for k, v in self.sample_values_relative_abundance_pct.items()
            },
            "imputation_applied": self.imputation_applied,
        }


@dataclass(frozen=True)
class MultiSampleQuantificationMatrix:
    """Multi-sample quantification matrix and replicate statistical summary."""

    sample_ids: Tuple[str, ...]
    glycan_compositions: Tuple[str, ...]
    abundance_matrix: Mapping[str, Mapping[str, Optional[float]]]  # sample_id -> {composition: rel_abund_pct or None}
    intensity_matrix: Mapping[str, Mapping[str, Optional[float]]]  # sample_id -> {composition: raw_intensity or None}
    replicate_summaries: Tuple[ReplicateGlycanSummary, ...]
    aggregation_policy_used: GlycanAggregationPolicy
    denominator_policy_used: DenominatorSetPolicy
    warnings: Tuple[str, ...]

    def to_summary_dict(self) -> Dict[str, Any]:
        return {
            "sample_ids": list(self.sample_ids),
            "glycan_compositions": list(self.glycan_compositions),
            "abundance_matrix": {
                sid: {
                    comp: (round(val, 4) if val is not None else None)
                    for comp, val in row.items()
                }
                for sid, row in self.abundance_matrix.items()
            },
            "intensity_matrix": {
                sid: {
                    comp: (round(val, 4) if val is not None else None)
                    for comp, val in row.items()
                }
                for sid, row in self.intensity_matrix.items()
            },
            "replicate_summaries": [s.to_dict() for s in self.replicate_summaries],
            "aggregation_policy_used": self.aggregation_policy_used.value,
            "denominator_policy_used": self.denominator_policy_used.value,
            "warnings": list(self.warnings),
        }


def quantify_ms1_search_result(
    search_result: MS1SearchResult,
    denominator_policy: DenominatorSetPolicy = DenominatorSetPolicy.ALL_ELIGIBLE_SPECTRUM_PEAKS,
    aggregation_policy: GlycanAggregationPolicy = GlycanAggregationPolicy.PRIMARY_ADDUCT_ONLY,
    primary_adduct: str = "[M+Na]+",
    exclude_ambiguous_from_aggregation: bool = False,
) -> SampleQuantificationResult:
    """Compute peak-level and composition-level relative quantification from an `MS1SearchResult`."""
    warnings: List[str] = [
        "Normalized peak intensities and relative abundances are semi-quantitative compositional metrics and do NOT represent absolute molar glycan concentrations."
    ]

    top_matches: Sequence[MS1CandidateMatch] = search_result.top_matches_per_peak
    if not top_matches:
        warnings.append("No matched glycan peaks available for quantification.")
        return SampleQuantificationResult(
            sample_id=search_result.sample_id,
            denominator_policy=denominator_policy,
            aggregation_policy=aggregation_policy,
            primary_adduct=primary_adduct,
            exclude_ambiguous_from_aggregation=exclude_ambiguous_from_aggregation,
            denominator_max_intensity=0.0,
            total_annotated_intensity=0.0,
            total_annotated_area=None,
            has_peak_area_data=False,
            peak_records=(),
            glycan_abundances=(),
            warnings=tuple(warnings),
        )

    # Determine denominator maximum peak intensity
    max_annotated_int = max((m.observed_intensity for m in top_matches), default=0.0)
    max_unmatched_int = max((p.intensity for p in search_result.unmatched_peaks), default=0.0)

    if denominator_policy == DenominatorSetPolicy.ALL_ELIGIBLE_SPECTRUM_PEAKS:
        denom_max = max(max_annotated_int, max_unmatched_int)
    else:
        denom_max = max_annotated_int

    total_annot_int = sum(m.observed_intensity for m in top_matches)
    all_have_area = all(m.observed_area is not None and m.observed_area > 0.0 for m in top_matches)
    total_annot_area: Optional[float] = (
        sum(m.observed_area for m in top_matches if m.observed_area is not None) if all_have_area else None
    )

    if not all_have_area:
        warnings.append(
            "Insufficient profile data for peak-area calculation on one or more peaks; area-based relative abundance is reported as None (not fabricated)."
        )

    peak_records: List[PeakQuantificationRecord] = []
    for m in top_matches:
        bp_rel = (m.observed_intensity / denom_max) * 100.0 if denom_max > 0 else 0.0
        frac_int = (m.observed_intensity / total_annot_int) * 100.0 if total_annot_int > 0 else 0.0
        frac_area: Optional[float] = None
        if all_have_area and total_annot_area and total_annot_area > 0 and m.observed_area is not None:
            frac_area = (m.observed_area / total_annot_area) * 100.0

        area_note = (
            "Integrated profile/table peak area"
            if m.observed_area is not None
            else "Unavailable (centroid-only peak list)"
        )
        peak_records.append(
            PeakQuantificationRecord(
                sample_id=search_result.sample_id,
                peak_index=m.peak_index,
                observed_mz=m.observed_mz,
                library_id=m.library_id,
                glycan_composition=m.glycan_composition,
                candidate_class=m.candidate_class,
                assigned_adduct=m.assigned_adduct,
                charge_state=m.charge_state,
                is_ambiguous=m.is_ambiguous,
                observed_intensity=m.observed_intensity,
                observed_area=m.observed_area,
                area_status_note=area_note,
                base_peak_relative_intensity_pct=bp_rel,
                denominator_policy=denominator_policy,
                denominator_max_intensity=denom_max,
                fraction_of_annotated_intensity_pct=frac_int,
                fraction_of_annotated_area_pct=frac_area,
            )
        )

    # Group eligible peak matches by glycan composition according to `aggregation_policy`
    eligible_for_agg: List[MS1CandidateMatch] = []
    for m in top_matches:
        if exclude_ambiguous_from_aggregation and m.is_ambiguous:
            warnings.append(
                f"Excluded ambiguous peak at m/z {m.observed_mz:.4f} ({m.glycan_composition}) from composition-level aggregation."
            )
            continue
        eligible_for_agg.append(m)

    raw_groups: Dict[str, List[MS1CandidateMatch]] = {}
    if aggregation_policy == GlycanAggregationPolicy.NO_AGGREGATION_PEAK_LEVEL:
        # Each peak remains its own entry keyed by composition + adduct + peak_index
        for m in eligible_for_agg:
            key = f"{m.glycan_composition}@{m.assigned_adduct}#p{m.peak_index}"
            raw_groups[key] = [m]
    elif aggregation_policy == GlycanAggregationPolicy.PRIMARY_ADDUCT_ONLY:
        for m in eligible_for_agg:
            if m.assigned_adduct != primary_adduct:
                continue
            raw_groups.setdefault(m.glycan_composition, []).append(m)
        # If more than 1 peak matched the same composition and same primary adduct, retain only the highest intensity peak
        for comp_key, grp in list(raw_groups.items()):
            if len(grp) > 1:
                best = max(grp, key=lambda x: x.observed_intensity)
                raw_groups[comp_key] = [best]
                warnings.append(
                    f"Multiple peaks matched {comp_key} with primary adduct {primary_adduct}; retained highest-intensity peak (m/z {best.observed_mz:.4f}) to avoid summing unrelated peaks."
                )
    elif aggregation_policy == GlycanAggregationPolicy.HIGHEST_INTENSITY_PEAK_PER_COMPOSITION:
        for m in eligible_for_agg:
            raw_groups.setdefault(m.glycan_composition, []).append(m)
        for comp_key, grp in list(raw_groups.items()):
            best = max(grp, key=lambda x: x.observed_intensity)
            raw_groups[comp_key] = [best]
    elif aggregation_policy == GlycanAggregationPolicy.EXPLICIT_SUM_DISTINCT_ADDUCTS:
        for m in eligible_for_agg:
            raw_groups.setdefault(m.glycan_composition, []).append(m)
        # Deduplicate per adduct so we only sum across distinct adducts, never duplicate peaks of the same adduct
        for comp_key, grp in list(raw_groups.items()):
            by_adduct: Dict[str, MS1CandidateMatch] = {}
            for item in grp:
                if item.assigned_adduct not in by_adduct or item.observed_intensity > by_adduct[item.assigned_adduct].observed_intensity:
                    by_adduct[item.assigned_adduct] = item
            raw_groups[comp_key] = list(by_adduct.values())

    # Compute total aggregated intensity and area across all composition groups
    total_agg_int = sum(sum(m.observed_intensity for m in grp) for grp in raw_groups.values())
    all_groups_have_area = bool(raw_groups) and all(
        all(m.observed_area is not None and m.observed_area > 0.0 for m in grp)
        for grp in raw_groups.values()
    )
    total_agg_area: Optional[float] = (
        sum(sum(m.observed_area for m in grp if m.observed_area is not None) for grp in raw_groups.values())
        if all_groups_have_area
        else None
    )

    max_agg_int = max(
        (sum(m.observed_intensity for m in grp) for grp in raw_groups.values()),
        default=0.0,
    )
    bp_denom_for_agg = denom_max if denominator_policy == DenominatorSetPolicy.ALL_ELIGIBLE_SPECTRUM_PEAKS else max_agg_int

    abundance_records: List[GlycanAbundanceRecord] = []
    for key, grp in raw_groups.items():
        comp_name = grp[0].glycan_composition
        cand_class = grp[0].candidate_class
        agg_int = sum(m.observed_intensity for m in grp)
        agg_area = (
            sum(m.observed_area for m in grp if m.observed_area is not None)
            if all(m.observed_area is not None for m in grp)
            else None
        )
        bp_pct = (agg_int / bp_denom_for_agg) * 100.0 if bp_denom_for_agg > 0 else 0.0
        rel_int_pct = (agg_int / total_agg_int) * 100.0 if total_agg_int > 0 else 0.0
        rel_area_pct = (
            (agg_area / total_agg_area) * 100.0
            if (all_groups_have_area and total_agg_area and total_agg_area > 0 and agg_area is not None)
            else None
        )

        provenance = tuple(
            f"peak#{m.peak_index}: m/z={m.observed_mz:.4f} ({m.assigned_adduct}, z={m.charge_state:+d}, int={m.observed_intensity:.2f})"
            for m in grp
        )
        rec_warnings: List[str] = []
        if len(grp) > 1:
            rec_warnings.append(
                f"Aggregated {len(grp)} distinct adduct signals ({', '.join(m.assigned_adduct for m in grp)}) under policy '{aggregation_policy.value}'."
            )
        if any(m.is_ambiguous for m in grp):
            rec_warnings.append("Includes at least one peak with multiple alternative candidate matches within tolerance.")

        abundance_records.append(
            GlycanAbundanceRecord(
                sample_id=search_result.sample_id,
                glycan_composition=comp_name if aggregation_policy != GlycanAggregationPolicy.NO_AGGREGATION_PEAK_LEVEL else key,
                candidate_class=cand_class,
                aggregation_policy_used=aggregation_policy,
                primary_adduct_filter=primary_adduct if aggregation_policy == GlycanAggregationPolicy.PRIMARY_ADDUCT_ONLY else None,
                aggregated_intensity=agg_int,
                aggregated_area=agg_area,
                base_peak_relative_intensity_pct=bp_pct,
                relative_abundance_intensity_pct=rel_int_pct,
                relative_abundance_area_pct=rel_area_pct,
                contributing_peaks_count=len(grp),
                contributing_adducts=tuple(m.assigned_adduct for m in grp),
                contributing_mzs=tuple(m.observed_mz for m in grp),
                peak_provenance_history=provenance,
                contains_ambiguous_peak=any(m.is_ambiguous for m in grp),
                warnings=tuple(rec_warnings),
            )
        )

    abundance_records.sort(key=lambda r: (-r.relative_abundance_intensity_pct, r.glycan_composition))

    return SampleQuantificationResult(
        sample_id=search_result.sample_id,
        denominator_policy=denominator_policy,
        aggregation_policy=aggregation_policy,
        primary_adduct=primary_adduct,
        exclude_ambiguous_from_aggregation=exclude_ambiguous_from_aggregation,
        denominator_max_intensity=denom_max,
        total_annotated_intensity=total_annot_int,
        total_annotated_area=total_annot_area,
        has_peak_area_data=all_have_area,
        peak_records=tuple(peak_records),
        glycan_abundances=tuple(abundance_records),
        warnings=tuple(warnings),
    )


def summarize_replicates_and_build_matrix(
    sample_results: Sequence[SampleQuantificationResult],
) -> MultiSampleQuantificationMatrix:
    """Compute replicate statistics (Mean, SD, CV, Missingness %) and build sample-by-glycan matrices.

    CRITICAL SCIENTIFIC RULE:
    Missing glycan signals in a replicate are NEVER silently imputed as zero.
    They are stored as `None` in the sample-by-glycan matrix, counted in `missing_replicates`
    and `missingness_pct`, and Mean/SD/CV are computed across the observed replicates.
    """
    if not sample_results:
        raise ValueError("Cannot build replicate summary from an empty list of sample quantification results.")

    sample_ids = tuple(r.sample_id for r in sample_results)
    n_samples = len(sample_ids)

    # Collect union of all observed glycan compositions and their classes
    comp_to_class: Dict[str, str] = {}
    for res in sample_results:
        for g in res.glycan_abundances:
            if g.glycan_composition not in comp_to_class:
                comp_to_class[g.glycan_composition] = g.candidate_class

    sorted_comps = tuple(sorted(comp_to_class.keys()))

    # Build lookup maps per sample
    abundance_matrix: Dict[str, Dict[str, Optional[float]]] = {sid: {} for sid in sample_ids}
    intensity_matrix: Dict[str, Dict[str, Optional[float]]] = {sid: {} for sid in sample_ids}
    bp_matrix: Dict[str, Dict[str, Optional[float]]] = {sid: {} for sid in sample_ids}

    for res in sample_results:
        sid = res.sample_id
        by_comp = {g.glycan_composition: g for g in res.glycan_abundances}
        for comp in sorted_comps:
            if comp in by_comp:
                abundance_matrix[sid][comp] = by_comp[comp].relative_abundance_intensity_pct
                intensity_matrix[sid][comp] = by_comp[comp].aggregated_intensity
                bp_matrix[sid][comp] = by_comp[comp].base_peak_relative_intensity_pct
            else:
                abundance_matrix[sid][comp] = None
                intensity_matrix[sid][comp] = None
                bp_matrix[sid][comp] = None

    summaries: List[ReplicateGlycanSummary] = []
    for comp in sorted_comps:
        sample_vals: Dict[str, Optional[float]] = {sid: abundance_matrix[sid][comp] for sid in sample_ids}
        obs_vals = [v for v in sample_vals.values() if v is not None]
        obs_bp = [bp_matrix[sid][comp] for sid in sample_ids if bp_matrix[sid][comp] is not None]

        n_obs = len(obs_vals)
        n_miss = n_samples - n_obs
        miss_pct = (n_miss / n_samples) * 100.0 if n_samples > 0 else 0.0

        mean_rel: Optional[float] = None
        sd_rel: Optional[float] = None
        cv_pct: Optional[float] = None
        mean_bp: Optional[float] = None

        if n_obs >= 1:
            arr = np.array(obs_vals, dtype=np.float64)
            mean_rel = float(np.mean(arr))
            mean_bp = float(np.mean(np.array(obs_bp, dtype=np.float64)))
        if n_obs >= 2 and mean_rel is not None:
            sd_rel = float(np.std(np.array(obs_vals, dtype=np.float64), ddof=1))
            if mean_rel > 1e-12:
                cv_pct = (sd_rel / mean_rel) * 100.0

        summaries.append(
            ReplicateGlycanSummary(
                glycan_composition=comp,
                candidate_class=comp_to_class[comp],
                total_replicates=n_samples,
                observed_replicates=n_obs,
                missing_replicates=n_miss,
                missingness_pct=miss_pct,
                mean_relative_abundance_pct=mean_rel,
                sd_relative_abundance_pct=sd_rel,
                cv_pct=cv_pct,
                mean_base_peak_intensity_pct=mean_bp,
                sample_values_relative_abundance_pct=sample_vals,
            )
        )

    warnings = (
        "Missing glycan values across replicates are preserved as None (not silently imputed as zero).",
        "Relative abundance (%) is normalized to the sum of quantified glycan intensities within each sample.",
    )

    return MultiSampleQuantificationMatrix(
        sample_ids=sample_ids,
        glycan_compositions=sorted_comps,
        abundance_matrix=abundance_matrix,
        intensity_matrix=intensity_matrix,
        replicate_summaries=tuple(summaries),
        aggregation_policy_used=sample_results[0].aggregation_policy,
        denominator_policy_used=sample_results[0].denominator_policy,
        warnings=warnings,
    )

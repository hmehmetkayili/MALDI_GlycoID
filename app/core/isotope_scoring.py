"""Theoretical isotopic distribution generator and Isotope Dot Product (idotp) scoring.

Implements exact discrete convolution of elemental isotope distributions to generate
theoretical [M], [M+1], [M+2], ... envelopes, and computes normalized dot-product
similarity (`idotp`) ONLY when observed isotopic peaks are genuinely present.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from app.core.mass_error import MassToleranceUnit, is_within_mass_tolerance
from app.core.scientific_constants import (
    ELEMENT_ISOTOPES,
    NEUTRON_MASS_SHIFT_DA,
)


@dataclass(frozen=True)
class IsotopicPeak:
    """Single peak in a theoretical or observed isotopic envelope."""

    neutron_offset: int  # 0 for monoisotopic M0, 1 for M+1, 2 for M+2, etc.
    mz: float
    relative_abundance: float  # Normalized so sum(abundances) == 1.0
    base_peak_ratio: float  # Normalized so max(abundances) == 100.0


@dataclass(frozen=True)
class IsotopeScoreResult:
    """Result of comparing an observed spectrum against a theoretical isotopic envelope."""

    available: bool
    idotp: Optional[float]  # Normalized dot product in [0.0, 1.0], or None if unavailable
    matched_isotope_peaks: int
    expected_isotope_peaks: int
    status_message: str
    theoretical_envelope: Tuple[IsotopicPeak, ...]
    observed_intensities_matched: Tuple[float, ...]
    observed_mzs_matched: Tuple[Optional[float], ...] = ()

    def to_comparison_dict(self, ion_formula_str: str = "") -> Dict[str, object]:
        """Serialize theoretical vs experimental isotope envelope for interactive UI plotting."""
        max_obs = max(self.observed_intensities_matched, default=0.0)
        tot_obs = sum(self.observed_intensities_matched)
        peaks_list = []
        for idx, theo_pk in enumerate(self.theoretical_envelope):
            obs_int = self.observed_intensities_matched[idx] if idx < len(self.observed_intensities_matched) else 0.0
            obs_mz = (
                self.observed_mzs_matched[idx]
                if idx < len(self.observed_mzs_matched)
                else None
            )
            obs_base_ratio = (obs_int / max_obs * 100.0) if max_obs > 0 and obs_int > 0 else 0.0
            obs_rel_abund = (obs_int / tot_obs) if tot_obs > 0 and obs_int > 0 else 0.0
            err_ppm = (
                round(((obs_mz - theo_pk.mz) / theo_pk.mz) * 1e6, 2)
                if (obs_mz is not None and theo_pk.mz > 0)
                else None
            )
            peaks_list.append(
                {
                    "neutron_offset": theo_pk.neutron_offset,
                    "label": f"M+{theo_pk.neutron_offset}",
                    "theoretical_mz": round(theo_pk.mz, 5),
                    "theoretical_abundance_pct": round(theo_pk.relative_abundance * 100.0, 2),
                    "theoretical_base_ratio_pct": round(theo_pk.base_peak_ratio, 2),
                    "observed_mz": round(obs_mz, 5) if obs_mz is not None else None,
                    "observed_intensity": round(obs_int, 2),
                    "observed_abundance_pct": round(obs_rel_abund * 100.0, 2),
                    "observed_base_ratio_pct": round(obs_base_ratio, 2),
                    "delta_base_ratio_pct": round(obs_base_ratio - theo_pk.base_peak_ratio, 2)
                    if obs_int > 0
                    else None,
                    "error_ppm": err_ppm,
                }
            )
        return {
            "available": self.available,
            "idotp": self.idotp,
            "matched_isotope_peaks": self.matched_isotope_peaks,
            "expected_isotope_peaks": self.expected_isotope_peaks,
            "status_message": self.status_message,
            "ion_formula": ion_formula_str,
            "isotopologues": peaks_list,
        }


def _convolve_distributions(
    dist_a: List[Tuple[float, float]],
    dist_b: List[Tuple[float, float]],
    max_peaks: int,
) -> List[Tuple[float, float]]:
    """Convolve two discrete (mass_shift, probability) distributions grouped by nominal neutron offset."""
    len_a = len(dist_a)
    len_b = len(dist_b)
    out_len = min(max_peaks, len_a + len_b - 1)
    probs = [0.0] * out_len
    weighted_masses = [0.0] * out_len

    for i, (mass_a, prob_a) in enumerate(dist_a):
        if prob_a == 0.0:
            continue
        for j, (mass_b, prob_b) in enumerate(dist_b):
            k = i + j
            if k >= out_len or prob_b == 0.0:
                continue
            p = prob_a * prob_b
            probs[k] += p
            weighted_masses[k] += p * (mass_a + mass_b)

    result: List[Tuple[float, float]] = []
    for k in range(out_len):
        if probs[k] > 0.0:
            result.append((weighted_masses[k] / probs[k], probs[k]))
        else:
            result.append((k * NEUTRON_MASS_SHIFT_DA, 0.0))
    return result


def _power_convolve(
    single_atom_dist: List[Tuple[float, float]],
    atom_count: int,
    max_peaks: int,
) -> List[Tuple[float, float]]:
    """Compute the n-fold convolution of a single atom's isotopic distribution via binary exponentiation."""
    result: List[Tuple[float, float]] = [(0.0, 1.0)]
    base = list(single_atom_dist)
    n = atom_count
    while n > 0:
        if n & 1:
            result = _convolve_distributions(result, base, max_peaks)
        n >>= 1
        if n > 0:
            base = _convolve_distributions(base, base, max_peaks)
    return result


def generate_theoretical_isotope_envelope(
    ion_formula: Mapping[str, int],
    monoisotopic_mz: float,
    charge: int = 1,
    num_peaks: int = 5,
) -> Tuple[IsotopicPeak, ...]:
    """Generate the theoretical isotopic envelope for an ion given its elemental formula and monoisotopic m/z.

    Uses exact polynomial convolution of IUPAC natural isotopic abundances and isotope mass defects.
    """
    abs_z = abs(charge)
    if abs_z < 1:
        raise ValueError(f"Charge state must be non-zero, got {charge}")
    if num_peaks < 1:
        raise ValueError(f"num_peaks must be >= 1, got {num_peaks}")

    overall_dist: List[Tuple[float, float]] = [(0.0, 1.0)]

    for el, count in ion_formula.items():
        if count <= 0:
            continue
        if el not in ELEMENT_ISOTOPES:
            continue
        iso_specs = ELEMENT_ISOTOPES[el]
        if len(iso_specs) <= 1:
            continue

        mono_mass = iso_specs[0].exact_mass
        max_offset = max(spec.neutron_offset for spec in iso_specs)
        single_atom: List[Tuple[float, float]] = [
            (i * NEUTRON_MASS_SHIFT_DA, 0.0) for i in range(max_offset + 1)
        ]
        for spec in iso_specs:
            mass_defect_shift = spec.exact_mass - mono_mass
            single_atom[spec.neutron_offset] = (mass_defect_shift, spec.abundance)

        el_dist = _power_convolve(single_atom, count, num_peaks)
        overall_dist = _convolve_distributions(overall_dist, el_dist, num_peaks)

    total_prob = sum(p for _, p in overall_dist)
    if total_prob <= 0:
        return (IsotopicPeak(0, monoisotopic_mz, 1.0, 100.0),)

    max_prob = max(p for _, p in overall_dist)
    peaks: List[IsotopicPeak] = []
    for offset, (mass_shift, prob) in enumerate(overall_dist):
        if prob <= 1e-9:
            continue
        peak_mz = monoisotopic_mz + (mass_shift / abs_z)
        peaks.append(
            IsotopicPeak(
                neutron_offset=offset,
                mz=peak_mz,
                relative_abundance=prob / total_prob,
                base_peak_ratio=(prob / max_prob) * 100.0 if max_prob > 0 else 0.0,
            )
        )
    return tuple(peaks)


def score_envelope_against_peaks(
    theo_envelope: Sequence[IsotopicPeak],
    obs_mz_sorted: np.ndarray,
    obs_int_sorted: np.ndarray,
    tolerance: float = 20.0,
    tolerance_unit: MassToleranceUnit = MassToleranceUnit.PPM,
    min_matched_isotopes: int = 2,
) -> IsotopeScoreResult:
    """Score a precomputed theoretical envelope against m/z-sorted observed peaks (fast path).

    Nearest-peak lookup uses binary search (O(log n) per isotopologue) instead of a full
    array scan, so scoring thousands of candidates against a large peak list stays cheap.
    """
    theo_envelope = tuple(theo_envelope)
    n_obs = len(obs_mz_sorted)

    def _unavailable(matched_count: int, ints: Sequence[float], mzs: Sequence[Optional[float]]) -> IsotopeScoreResult:
        return IsotopeScoreResult(
            available=False,
            idotp=None,
            matched_isotope_peaks=matched_count,
            expected_isotope_peaks=len(theo_envelope),
            status_message="Insufficient isotope information for isotope scoring",
            theoretical_envelope=theo_envelope,
            observed_intensities_matched=tuple(ints),
            observed_mzs_matched=tuple(mzs),
        )

    # Fewer observed peaks than required isotopologues: score is meaningless
    if n_obs < min_matched_isotopes or len(theo_envelope) < min_matched_isotopes:
        return _unavailable(min(n_obs, 1), (), ())

    matched_ints: List[float] = []
    matched_mzs: List[Optional[float]] = []
    matched_count = 0

    for iso_peak in theo_envelope:
        pos = int(np.searchsorted(obs_mz_sorted, iso_peak.mz))
        best_idx = -1
        best_diff = math.inf
        for j in (pos - 1, pos):
            if 0 <= j < n_obs:
                d = abs(float(obs_mz_sorted[j]) - iso_peak.mz)
                if d < best_diff:
                    best_diff, best_idx = d, j
        if best_idx >= 0:
            cand_mz = float(obs_mz_sorted[best_idx])
            cand_int = float(obs_int_sorted[best_idx])
            if cand_int > 0 and is_within_mass_tolerance(
                observed=cand_mz,
                theoretical=iso_peak.mz,
                tolerance=tolerance,
                unit=tolerance_unit,
            ):
                matched_ints.append(cand_int)
                matched_mzs.append(cand_mz)
                matched_count += 1
                continue
        matched_ints.append(0.0)
        matched_mzs.append(None)

    if matched_count < min_matched_isotopes:
        return _unavailable(matched_count, matched_ints, matched_mzs)

    theo_vec = np.array([p.relative_abundance for p in theo_envelope], dtype=float)
    obs_vec = np.array(matched_ints, dtype=float)
    norm_theo = float(np.linalg.norm(theo_vec))
    norm_obs = float(np.linalg.norm(obs_vec))
    if norm_theo <= 0.0 or norm_obs <= 0.0:
        return _unavailable(matched_count, matched_ints, matched_mzs)

    dot = float(np.dot(theo_vec, obs_vec) / (norm_theo * norm_obs))
    dot = max(0.0, min(1.0, dot))
    return IsotopeScoreResult(
        available=True,
        idotp=round(dot, 6),
        matched_isotope_peaks=matched_count,
        expected_isotope_peaks=len(theo_envelope),
        status_message=f"Isotope envelope scored across {matched_count}/{len(theo_envelope)} peaks",
        theoretical_envelope=theo_envelope,
        observed_intensities_matched=tuple(matched_ints),
        observed_mzs_matched=tuple(matched_mzs),
    )


def calculate_isotope_dot_product(
    observed_mzs: Sequence[float],
    observed_intensities: Sequence[float],
    ion_formula: Mapping[str, int],
    monoisotopic_mz: float,
    charge: int = 1,
    tolerance: float = 20.0,
    tolerance_unit: MassToleranceUnit = MassToleranceUnit.PPM,
    num_isotopes: int = 4,
    min_matched_isotopes: int = 2,
) -> IsotopeScoreResult:
    """Compare observed spectrum peaks against the theoretical isotope distribution (`idotp`).

    CRITICAL SCIENTIFIC RULE:
    If only a single centroid peak (or fewer than `min_matched_isotopes` peaks) is present
    in the observed spectrum for the candidate's isotopic envelope, `idotp` cannot be
    meaningfully calculated and MUST be returned as `None` with `available=False` and
    message `'Insufficient isotope information for isotope scoring'`.
    """
    theo_envelope = generate_theoretical_isotope_envelope(
        ion_formula=ion_formula,
        monoisotopic_mz=monoisotopic_mz,
        charge=charge,
        num_peaks=num_isotopes,
    )

    if len(observed_mzs) != len(observed_intensities):
        raise ValueError("observed_mzs and observed_intensities must have equal length")

    obs_mz_arr = np.asarray(observed_mzs, dtype=float)
    obs_int_arr = np.asarray(observed_intensities, dtype=float)
    if len(obs_mz_arr) > 1 and np.any(np.diff(obs_mz_arr) < 0):
        order = np.argsort(obs_mz_arr, kind="stable")
        obs_mz_arr = obs_mz_arr[order]
        obs_int_arr = obs_int_arr[order]

    return score_envelope_against_peaks(
        theo_envelope=theo_envelope,
        obs_mz_sorted=obs_mz_arr,
        obs_int_sorted=obs_int_arr,
        tolerance=tolerance,
        tolerance_unit=tolerance_unit,
        min_matched_isotopes=min_matched_isotopes,
    )

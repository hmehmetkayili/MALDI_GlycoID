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

    # If the input only provided a single peak (or empty), isotope scoring is unavailable
    if len(observed_mzs) < min_matched_isotopes or len(theo_envelope) < min_matched_isotopes:
        return IsotopeScoreResult(
            available=False,
            idotp=None,
            matched_isotope_peaks=min(len(observed_mzs), 1),
            expected_isotope_peaks=len(theo_envelope),
            status_message="Insufficient isotope information for isotope scoring",
            theoretical_envelope=theo_envelope,
            observed_intensities_matched=(),
        )

    matched_obs_intensities: List[float] = []
    matched_count = 0

    obs_mz_arr = np.asarray(observed_mzs, dtype=float)
    obs_int_arr = np.asarray(observed_intensities, dtype=float)

    for iso_peak in theo_envelope:
        # Find closest observed peak within tolerance
        diffs = np.abs(obs_mz_arr - iso_peak.mz)
        best_idx = int(np.argmin(diffs))
        candidate_mz = float(obs_mz_arr[best_idx])
        candidate_int = float(obs_int_arr[best_idx])

        if candidate_int > 0 and is_within_mass_tolerance(
            observed=candidate_mz,
            theoretical=iso_peak.mz,
            tolerance=tolerance,
            unit=tolerance_unit,
        ):
            matched_obs_intensities.append(candidate_int)
            matched_count += 1
        else:
            matched_obs_intensities.append(0.0)

    if matched_count < min_matched_isotopes:
        return IsotopeScoreResult(
            available=False,
            idotp=None,
            matched_isotope_peaks=matched_count,
            expected_isotope_peaks=len(theo_envelope),
            status_message="Insufficient isotope information for isotope scoring",
            theoretical_envelope=theo_envelope,
            observed_intensities_matched=tuple(matched_obs_intensities),
        )

    # Compute normalized cosine similarity (dot product) between theoretical and observed vectors
    theo_vec = np.array([p.relative_abundance for p in theo_envelope], dtype=float)
    obs_vec = np.array(matched_obs_intensities, dtype=float)

    norm_theo = float(np.linalg.norm(theo_vec))
    norm_obs = float(np.linalg.norm(obs_vec))
    if norm_theo <= 0.0 or norm_obs <= 0.0:
        return IsotopeScoreResult(
            available=False,
            idotp=None,
            matched_isotope_peaks=matched_count,
            expected_isotope_peaks=len(theo_envelope),
            status_message="Insufficient isotope information for isotope scoring",
            theoretical_envelope=theo_envelope,
            observed_intensities_matched=tuple(matched_obs_intensities),
        )

    dot = float(np.dot(theo_vec, obs_vec) / (norm_theo * norm_obs))
    dot = max(0.0, min(1.0, dot))

    return IsotopeScoreResult(
        available=True,
        idotp=round(dot, 6),
        matched_isotope_peaks=matched_count,
        expected_isotope_peaks=len(theo_envelope),
        status_message=f"Isotope envelope scored across {matched_count}/{len(theo_envelope)} peaks",
        theoretical_envelope=theo_envelope,
        observed_intensities_matched=tuple(matched_obs_intensities),
    )

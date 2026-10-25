"""Deterministic synthetic MALDI-TOF N-glycan spectra with realistic complexity (for tests and demos).

The generated *profile* spectra deliberately contain everything that makes real MALDI data hard:

- ~45 N-glycans spanning high-mannose, hybrid, complex, fucosylated and sialylated species
  with a realistic abundance distribution (spanning ~3 orders of magnitude)
- exact isotopic envelopes (convolved IUPAC abundances) for every ion
- [M+Na]+ with [M+K]+ satellites in positive mode, [M-H]- in negative mode
- TOF peak shape (Gaussian, FWHM = m/z / resolution) on a non-uniform TOF-like m/z grid
- exponentially decaying chemical-noise baseline + Gaussian + shot noise
- low-mass matrix ions (DHB / CHCA clusters)
- a systematic, m/z-dependent mass-calibration error (offset + slope, in ppm)
- near-isobaric interferers and spurious noise spikes

The ground truth is returned together with the spectrum so tests can verify identification.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, replace
from typing import Dict, List, Optional, Tuple

import numpy as np

from app.analysis.spectrum_corrections import matrix_ion_list
from app.core.ion_calculator import calculate_mz_from_neutral_mass
from app.core.isotope_scoring import generate_theoretical_isotope_envelope
from app.core.mass_calculator import calculate_released_glycan_mass
from app.parsers.base_models import SpectrumData, SpectrumDataType, SpectrumPeak

# (composition, relative abundance %) — canonical, biologically plausible N-glycan panel.
# The first eight are the default calibrants of `get_default_nglycan_calibrants`.
_CORE_PANEL: List[Tuple[str, float]] = [
    ("Hex5HexNAc2", 45.0),
    ("Hex3HexNAc4dHex1", 35.0),
    ("Hex4HexNAc4dHex1", 40.0),
    ("Hex5HexNAc4", 38.0),
    ("Hex5HexNAc4dHex1", 100.0),
    ("Hex5HexNAc4Neu5Ac1", 55.0),
    ("Hex5HexNAc4dHex1Neu5Ac1", 70.0),
    ("Hex5HexNAc4Neu5Ac2", 48.0),
]

_EXTRA_POOL: List[str] = [
    "Hex6HexNAc2", "Hex7HexNAc2", "Hex8HexNAc2", "Hex9HexNAc2",
    "Hex3HexNAc3", "Hex4HexNAc3", "Hex5HexNAc3", "Hex3HexNAc4", "Hex4HexNAc4",
    "Hex3HexNAc5", "Hex4HexNAc5", "Hex5HexNAc5", "Hex3HexNAc3dHex1", "Hex4HexNAc3dHex1",
    "Hex5HexNAc3dHex1", "Hex3HexNAc5dHex1", "Hex4HexNAc5dHex1", "Hex5HexNAc5dHex1",
    "Hex6HexNAc3", "Hex6HexNAc3dHex1", "Hex5HexNAc4Neu5Ac1dHex2", "Hex6HexNAc5",
    "Hex6HexNAc5dHex1", "Hex6HexNAc5Neu5Ac1", "Hex6HexNAc5dHex1Neu5Ac1",
    "Hex6HexNAc5Neu5Ac2", "Hex6HexNAc5dHex1Neu5Ac2", "Hex6HexNAc5Neu5Ac3",
    "Hex7HexNAc6", "Hex7HexNAc6dHex1", "Hex7HexNAc6Neu5Ac2", "Hex7HexNAc6Neu5Ac3",
    "Hex7HexNAc6dHex1Neu5Ac3", "Hex7HexNAc6Neu5Ac4", "Hex6HexNAc4Neu5Ac1",
    "Hex6HexNAc4dHex1Neu5Ac1", "Hex5HexNAc5Neu5Ac1", "Hex4HexNAc4Neu5Ac1",
]


@dataclass(frozen=True)
class TruthEntry:
    composition: str
    adduct: str
    theoretical_mz: float
    observed_mz_nominal: float  # after the injected calibration error (apex of the mono peak)
    abundance: float
    is_calibrant: bool


@dataclass(frozen=True)
class SyntheticSpectrum:
    spectrum: SpectrumData
    truth: Tuple[TruthEntry, ...]
    ppm_offset: float
    ppm_slope_per_1000da: float
    polarity: str

    def inject_error_ppm(self, mz: float) -> float:
        return self.ppm_offset + self.ppm_slope_per_1000da * (mz - 2000.0) / 1000.0


def _ion(comp: str, adduct: str):
    mass = calculate_released_glycan_mass(comp, labeling="NONE")
    ion = calculate_mz_from_neutral_mass(mass.neutral_monoisotopic_mass, adduct, neutral_formula=mass.elemental_formula)
    return ion


def generate_complex_nglycan_spectrum(
    polarity: str = "positive",
    seed: int = 7,
    n_extra_glycans: int = 37,
    mz_range: Tuple[float, float] = (180.0, 3900.0),
    resolution: float = 12000.0,
    ppm_offset: float = 55.0,
    ppm_slope_per_1000da: float = 6.0,
    noise_sigma: float = 35.0,
    baseline_amplitude: float = 900.0,
    include_matrix: bool = True,
    points_per_fwhm: float = 5.0,
    sample_id: Optional[str] = None,
) -> SyntheticSpectrum:
    """Build a complex profile MALDI-TOF N-glycan spectrum (see module docstring)."""
    rng = random.Random(seed)
    nrng = np.random.default_rng(seed)
    neg = polarity.lower().startswith("neg")
    main_adduct = "[M-H]-" if neg else "[M+Na]+"
    sec_adduct = None if neg else "[M+K]+"

    # --- choose glycans & abundances -------------------------------------------------
    panel = list(_CORE_PANEL)
    extras = rng.sample(_EXTRA_POOL, k=min(n_extra_glycans, len(_EXTRA_POOL)))
    for comp in extras:
        # log-uniform 2.0 .. 35 % so every spiked glycan is above the 5x S/N detection floor
        ab = float(np.exp(rng.uniform(np.log(2.0), np.log(35.0))))
        panel.append((comp, ab))

    core_set = {c for c, _ in _CORE_PANEL}
    ions = []  # (comp, adduct, ion_result, abundance, is_cal)
    for comp, ab in panel:
        ions.append((comp, main_adduct, _ion(comp, main_adduct), ab, comp in core_set))
        if sec_adduct and ab > 12.0:
            ions.append((comp, sec_adduct, _ion(comp, sec_adduct), ab * 0.12, False))

    # --- grid: TOF-like non-uniform sampling (step proportional to local FWHM) ---------
    lo, hi = mz_range
    grid: List[float] = []
    x = lo
    while x < hi:
        grid.append(x)
        x += (x / resolution) / points_per_fwhm
    mz_grid = np.array(grid, dtype=np.float64)
    y = np.zeros_like(mz_grid)

    def apply_error(mz: float) -> float:
        return mz * (1.0 + (ppm_offset + ppm_slope_per_1000da * (mz - 2000.0) / 1000.0) * 1e-6)

    amp_scale = 9000.0  # intensity of a 100 % glycan monoisotopic-ish apex (counts)

    def add_peak(mz_obs: float, height: float) -> None:
        fwhm = mz_obs / resolution
        sigma = fwhm / 2.3548
        i0 = int(np.searchsorted(mz_grid, mz_obs - 5 * sigma))
        i1 = int(np.searchsorted(mz_grid, mz_obs + 5 * sigma))
        if i1 > i0:
            y[i0:i1] += height * np.exp(-0.5 * ((mz_grid[i0:i1] - mz_obs) / sigma) ** 2)

    truth: List[TruthEntry] = []
    for comp, adduct, ion, ab, is_cal in ions:
        env = generate_theoretical_isotope_envelope(
            ion_formula=ion.ion_formula,
            monoisotopic_mz=ion.theoretical_mz,
            charge=ion.charge,
            num_peaks=7,
        )
        base = max(p.base_peak_ratio for p in env) or 1.0
        # Peak-height scaling: envelope's tallest isotopologue gets amp_scale * ab/100
        for pk in env:
            h = amp_scale * (ab / 100.0) * (pk.base_peak_ratio / base)
            if h > 1.0:
                add_peak(apply_error(pk.mz), h)
        truth.append(
            TruthEntry(
                composition=comp,
                adduct=adduct,
                theoretical_mz=ion.theoretical_mz,
                observed_mz_nominal=apply_error(ion.theoretical_mz),
                abundance=ab,
                is_calibrant=is_cal,
            )
        )

    # --- matrix ions (strong, low mass) -----------------------------------------------
    if include_matrix:
        for label, mz in matrix_ion_list("negative" if neg else "positive"):
            if lo < mz < hi and "2M" not in label[:6]:
                add_peak(apply_error(mz), rng.uniform(1500, 6000))

    # --- near-isobaric interferer: +0.036 Da (Na vs. H2-like fine mass difference) ------
    for comp, _ in panel[:3]:
        ion = _ion(comp, main_adduct)
        add_peak(apply_error(ion.theoretical_mz + 0.0364 * (-1 if neg else 1)), amp_scale * 0.02)

    # --- baseline + noise --------------------------------------------------------------
    baseline = baseline_amplitude * np.exp(-(mz_grid - lo) / 700.0) + 25.0
    shot = nrng.normal(0.0, 1.0, size=len(mz_grid)) * np.sqrt(np.maximum(y, 0.0)) * 0.6
    gauss = nrng.normal(0.0, noise_sigma, size=len(mz_grid))
    total = np.maximum(0.0, y + baseline + shot + gauss)

    # --- spurious noise spikes ----------------------------------------------------------
    n_spikes = 60
    for _ in range(n_spikes):
        add_idx = int(rng.uniform(0, len(mz_grid) - 1))
        total[add_idx] += rng.uniform(4, 9) * noise_sigma

    peaks = tuple(
        SpectrumPeak(mz=float(m), intensity=float(i), row_index=k + 1)
        for k, (m, i) in enumerate(zip(mz_grid, total))
    )
    spectrum = SpectrumData(
        sample_id=sample_id or f"SYNTH_COMPLEX_{'NEG' if neg else 'POS'}",
        filename=f"synthetic_complex_nglycans_{'neg' if neg else 'pos'}.csv",
        source_format="CSV",
        data_type=SpectrumDataType.PROFILE,
        user_specified_data_type=SpectrumDataType.PROFILE,
        peaks=peaks,
        total_rows_read=len(peaks),
        valid_points_count=len(peaks),
        min_mz=float(mz_grid[0]),
        max_mz=float(mz_grid[-1]),
        max_intensity=float(total.max()),
        metadata={"synthetic": True, "polarity": "negative" if neg else "positive"},
    )
    return SyntheticSpectrum(
        spectrum=spectrum,
        truth=tuple(truth),
        ppm_offset=ppm_offset,
        ppm_slope_per_1000da=ppm_slope_per_1000da,
        polarity="negative" if neg else "positive",
    )


def centroid_from_profile(spec: SyntheticSpectrum, min_snr: float = 3.0) -> SpectrumData:
    """Convenience: centroid a synthetic profile with the library's own peak picker."""
    from app.analysis.spectrum_corrections import CorrectionConfig, apply_corrections

    return apply_corrections(spec.spectrum, CorrectionConfig(pick_peaks=True, pick_min_snr=min_snr)).spectrum

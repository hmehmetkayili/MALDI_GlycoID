"""Mass error calculation and tolerance-window filtering engine.

Conventions:
- Signed mass error (Da)  = observed_mass - theoretical_mass
- Absolute mass error (Da) = |observed_mass - theoretical_mass|
- Signed mass error (ppm) = ((observed_mass - theoretical_mass) / theoretical_mass) * 1e6
- Absolute mass error (ppm) = |Signed mass error (ppm)|
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Tuple


class MassToleranceUnit(str, Enum):
    """Supported units for mass accuracy tolerance windows."""

    PPM = "ppm"
    DA = "Da"


@dataclass(frozen=True)
class MassErrorResult:
    """Comprehensive mass accuracy metrics comparing an observed and theoretical mass or m/z."""

    observed: float
    theoretical: float
    signed_error_da: float
    abs_error_da: float
    signed_error_ppm: float
    abs_error_ppm: float


def calculate_signed_error_da(observed: float, theoretical: float) -> float:
    """Calculate signed mass error in Daltons (Da): observed - theoretical."""
    if theoretical <= 0 or observed <= 0:
        raise ValueError(f"Observed ({observed}) and theoretical ({theoretical}) masses must be > 0")
    return observed - theoretical


def calculate_abs_error_da(observed: float, theoretical: float) -> float:
    """Calculate absolute mass error in Daltons (Da): |observed - theoretical|."""
    return abs(calculate_signed_error_da(observed, theoretical))


def calculate_signed_error_ppm(observed: float, theoretical: float) -> float:
    """Calculate signed mass error in parts-per-million (ppm): ((observed - theoretical) / theoretical) * 1e6."""
    if theoretical <= 0 or observed <= 0:
        raise ValueError(f"Observed ({observed}) and theoretical ({theoretical}) masses must be > 0")
    return ((observed - theoretical) / theoretical) * 1_000_000.0


def calculate_abs_error_ppm(observed: float, theoretical: float) -> float:
    """Calculate absolute mass error in parts-per-million (ppm)."""
    return abs(calculate_signed_error_ppm(observed, theoretical))


def evaluate_mass_error(observed: float, theoretical: float) -> MassErrorResult:
    """Compute signed and absolute mass errors in both Da and ppm."""
    signed_da = calculate_signed_error_da(observed, theoretical)
    signed_ppm = (signed_da / theoretical) * 1_000_000.0
    return MassErrorResult(
        observed=observed,
        theoretical=theoretical,
        signed_error_da=signed_da,
        abs_error_da=abs(signed_da),
        signed_error_ppm=signed_ppm,
        abs_error_ppm=abs(signed_ppm),
    )


def get_tolerance_window_bounds(
    target_mass: float,
    tolerance: float,
    unit: MassToleranceUnit = MassToleranceUnit.PPM,
) -> Tuple[float, float]:
    """Return (lower_bound, upper_bound) around `target_mass` for a given tolerance."""
    if target_mass <= 0:
        raise ValueError(f"Target mass must be positive, got {target_mass}")
    if tolerance < 0:
        raise ValueError(f"Mass tolerance must be non-negative, got {tolerance}")

    if unit == MassToleranceUnit.PPM:
        delta = target_mass * (tolerance / 1_000_000.0)
    elif unit == MassToleranceUnit.DA:
        delta = tolerance
    else:
        raise ValueError(f"Unsupported tolerance unit: {unit}")

    return (target_mass - delta, target_mass + delta)


def is_within_mass_tolerance(
    observed: float,
    theoretical: float,
    tolerance: float,
    unit: MassToleranceUnit = MassToleranceUnit.PPM,
) -> bool:
    """Return True if `observed` is within `tolerance` (ppm or Da) of `theoretical`."""
    if tolerance < 0:
        raise ValueError(f"Mass tolerance must be non-negative, got {tolerance}")
    err = evaluate_mass_error(observed, theoretical)
    if unit == MassToleranceUnit.PPM:
        return err.abs_error_ppm <= tolerance + 1e-12
    if unit == MassToleranceUnit.DA:
        return err.abs_error_da <= tolerance + 1e-12
    raise ValueError(f"Unsupported tolerance unit: {unit}")

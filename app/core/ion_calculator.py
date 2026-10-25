"""Ion and adduct m/z calculation engine.

Implements exact charge-state and adduct transformations without mixing
neutral hydrogen atom mass (1.007825 Da) and proton mass (1.007276 Da).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence

from app.core.formula_calculator import add_formulas, format_hill_formula
from app.core.scientific_constants import (
    ADDUCT_DEFINITIONS,
    ANION_CL_MASS,
    ANION_DEPROTONATED_MASS,
    ANION_FORMATE_MASS,
    AdductDefinition,
    CATION_H_MASS,
    CATION_K_MASS,
    CATION_LI_MASS,
    CATION_NA_MASS,
    CATION_NH4_MASS,
    IonPolarity,
)

# Maps generic carrier symbols to their single-charge ionic mass and elemental formula delta
CARRIER_ION_SPECS: Mapping[str, Tuple[float, Mapping[str, int], Mapping[str, int], IonPolarity]] = {
    "H+": (CATION_H_MASS, {"H": 1}, {}, IonPolarity.POSITIVE),
    "Na+": (CATION_NA_MASS, {"Na": 1}, {}, IonPolarity.POSITIVE),
    "K+": (CATION_K_MASS, {"K": 1}, {}, IonPolarity.POSITIVE),
    "Li+": (CATION_LI_MASS, {"Li": 1}, {}, IonPolarity.POSITIVE),
    "NH4+": (CATION_NH4_MASS, {"N": 1, "H": 4}, {}, IonPolarity.POSITIVE),
    "-H+": (ANION_DEPROTONATED_MASS, {}, {"H": 1}, IonPolarity.NEGATIVE),
    "Cl-": (ANION_CL_MASS, {"Cl": 1}, {}, IonPolarity.NEGATIVE),
    "HCOO-": (ANION_FORMATE_MASS, {"C": 1, "H": 1, "O": 2}, {}, IonPolarity.NEGATIVE),
}


@dataclass(frozen=True)
class IonCalculationResult:
    """Result of an ion m/z calculation for a given neutral mass and adduct."""

    neutral_mass: float
    adduct_code: str
    polarity: IonPolarity
    charge: int  # Signed charge (+1, +2, -1, etc.)
    abs_charge: int
    theoretical_mz: float
    net_carrier_mass: float
    ion_formula: Optional[Mapping[str, int]] = None
    ion_hill_formula: Optional[str] = None


def resolve_adduct_definition(
    adduct_code: str,
    custom_charge: Optional[int] = None,
) -> AdductDefinition:
    """Resolve an adduct code from ADDUCT_DEFINITIONS or build a multi-charged homogeneous adduct."""
    code = adduct_code.strip()
    if code in ADDUCT_DEFINITIONS:
        adduct = ADDUCT_DEFINITIONS[code]
        if custom_charge is not None and abs(custom_charge) != abs(adduct.charge):
            # If user requested e.g. "[M+H]+" with custom_charge=2 -> upgrade to [M+2H]2+
            return build_homogeneous_adduct_for_charge(code, abs(custom_charge))
        return adduct

    # Handle shorthand carrier names ("H+", "Na+", "K+", "-H+", etc.)
    if code in CARRIER_ION_SPECS:
        z = abs(custom_charge) if custom_charge is not None else 1
        return build_homogeneous_adduct_from_carrier(code, z)

    raise ValueError(
        f"Unsupported adduct code: '{adduct_code}'. Supported codes: {', '.join(ADDUCT_DEFINITIONS.keys())}"
    )


def build_homogeneous_adduct_from_carrier(carrier: str, abs_charge: int) -> AdductDefinition:
    """Construct a z-charged homogeneous adduct from a single carrier ion (e.g., H+, Na+, K+, -H+)."""
    if abs_charge < 1 or abs_charge > 10:
        raise ValueError(f"Invalid charge state magnitude: {abs_charge} (must be between 1 and 10)")
    if carrier not in CARRIER_ION_SPECS:
        raise ValueError(f"Unsupported charge carrier '{carrier}'")

    unit_mass, added, removed, polarity = CARRIER_ION_SPECS[carrier]
    signed_charge = abs_charge if polarity == IonPolarity.POSITIVE else -abs_charge
    z_str = f"{abs_charge}" if abs_charge > 1 else ""
    sign_str = "+" if polarity == IonPolarity.POSITIVE else "-"

    if carrier == "-H+":
        code = f"[M-{z_str}H]{z_str}-"
    else:
        base_sym = carrier.rstrip("+-")
        code = f"[M+{z_str}{base_sym}]{z_str}{sign_str}"

    if code in ADDUCT_DEFINITIONS:
        return ADDUCT_DEFINITIONS[code]

    scaled_added = {k: v * abs_charge for k, v in added.items()}
    scaled_removed = {k: v * abs_charge for k, v in removed.items()}
    return AdductDefinition(
        code=code,
        name=f"{carrier} adduct ({signed_charge:+d})",
        polarity=polarity,
        charge=signed_charge,
        net_carrier_mass=unit_mass * abs_charge,
        formula_added=scaled_added,
        formula_removed=scaled_removed,
        description=f"Charge {signed_charge:+d} ion with {abs_charge} x {carrier}",
    )


def build_homogeneous_adduct_for_charge(adduct_code: str, abs_charge: int) -> AdductDefinition:
    """Map a singly charged adduct code (e.g. '[M+Na]+') to a target charge state z."""
    mapping = {
        "[M+H]+": "H+",
        "[M+Na]+": "Na+",
        "[M+K]+": "K+",
        "[M+Li]+": "Li+",
        "[M+NH4]+": "NH4+",
        "[M-H]-": "-H+",
        "[M+Cl]-": "Cl-",
        "[M+HCOO]-": "HCOO-",
    }
    if adduct_code in mapping:
        return build_homogeneous_adduct_from_carrier(mapping[adduct_code], abs_charge)
    if adduct_code in ADDUCT_DEFINITIONS and abs(ADDUCT_DEFINITIONS[adduct_code].charge) == abs_charge:
        return ADDUCT_DEFINITIONS[adduct_code]
    raise ValueError(
        f"Cannot automatically scale heterogeneous adduct '{adduct_code}' to charge state {abs_charge}"
    )


def calculate_mz_from_neutral_mass(
    neutral_mass: float,
    adduct_code: str = "[M+Na]+",
    custom_charge: Optional[int] = None,
    neutral_formula: Optional[Mapping[str, int]] = None,
) -> IonCalculationResult:
    """Calculate the exact theoretical m/z of an ion from its neutral monoisotopic mass.

    Formula:
      m/z = (neutral_mass + net_carrier_mass) / |z|
    where `net_carrier_mass` uses true ionic carrier masses (accounting for electron mass).
    """
    if neutral_mass <= 0:
        raise ValueError(f"Neutral mass must be positive, got {neutral_mass}")

    adduct = resolve_adduct_definition(adduct_code, custom_charge=custom_charge)
    abs_z = abs(adduct.charge)
    mz = (neutral_mass + adduct.net_carrier_mass) / abs_z

    ion_formula: Optional[Dict[str, int]] = None
    ion_hill: Optional[str] = None
    if neutral_formula is not None:
        temp = dict(neutral_formula)
        if adduct.formula_added:
            temp = add_formulas(temp, adduct.formula_added, scale=1, allow_negative=False)
        if adduct.formula_removed:
            temp = add_formulas(temp, adduct.formula_removed, scale=-1, allow_negative=False)
        ion_formula = temp
        ion_hill = format_hill_formula(temp)

    return IonCalculationResult(
        neutral_mass=neutral_mass,
        adduct_code=adduct.code,
        polarity=adduct.polarity,
        charge=adduct.charge,
        abs_charge=abs_z,
        theoretical_mz=mz,
        net_carrier_mass=adduct.net_carrier_mass,
        ion_formula=ion_formula,
        ion_hill_formula=ion_hill,
    )


def calculate_neutral_mass_from_mz(
    observed_mz: float,
    adduct_code: str = "[M+Na]+",
    custom_charge: Optional[int] = None,
) -> float:
    """Back-calculate neutral monoisotopic mass from an observed m/z and assumed adduct.

    Formula:
      M = observed_mz * |z| - net_carrier_mass
    """
    if observed_mz <= 0:
        raise ValueError(f"Observed m/z must be positive, got {observed_mz}")
    adduct = resolve_adduct_definition(adduct_code, custom_charge=custom_charge)
    abs_z = abs(adduct.charge)
    return observed_mz * abs_z - adduct.net_carrier_mass


def calculate_all_adduct_mzs(
    neutral_mass: float,
    adduct_codes: Optional[Sequence[str]] = None,
    polarity: Optional[IonPolarity] = None,
    neutral_formula: Optional[Mapping[str, int]] = None,
) -> Dict[str, IonCalculationResult]:
    """Calculate theoretical m/z values across a list of adducts (or all adducts for a polarity)."""
    if adduct_codes is None:
        codes = [
            code
            for code, defn in ADDUCT_DEFINITIONS.items()
            if polarity is None or defn.polarity == polarity
        ]
    else:
        codes = list(adduct_codes)

    results: Dict[str, IonCalculationResult] = {}
    for code in codes:
        res = calculate_mz_from_neutral_mass(
            neutral_mass=neutral_mass,
            adduct_code=code,
            neutral_formula=neutral_formula,
        )
        results[res.adduct_code] = res
    return results

"""Elemental formula parsing, stoichiometry arithmetic, and monoisotopic mass calculation."""

from __future__ import annotations

import re
from typing import Dict, Mapping, Union

from app.core.scientific_constants import ELEMENT_MONOISOTOPIC_MASSES

_FORMULA_TOKEN_RE = re.compile(r"([A-Z][a-z]?)([+-]?\d*)")


def parse_elemental_formula(formula: Union[str, Mapping[str, int]]) -> Dict[str, int]:
    """Parse an elemental formula string (e.g. 'C6H10O5') or mapping into a normalized dict.

    Raises ValueError if unknown element symbols or malformed characters are present.
    """
    if isinstance(formula, Mapping):
        result: Dict[str, int] = {}
        for el, count in formula.items():
            if el not in ELEMENT_MONOISOTOPIC_MASSES:
                raise ValueError(f"Unknown chemical element symbol: '{el}'")
            if not isinstance(count, int):
                raise ValueError(f"Element count for '{el}' must be an integer, got {count!r}")
            if count != 0:
                result[el] = result.get(el, 0) + count
        return result

    if not isinstance(formula, str):
        raise ValueError(f"Formula must be a string or mapping, got {type(formula).__name__}")

    cleaned = formula.strip().replace(" ", "")
    if not cleaned:
        return {}

    pos = 0
    counts: Dict[str, int] = {}
    for match in _FORMULA_TOKEN_RE.finditer(cleaned):
        if match.start() != pos:
            raise ValueError(f"Invalid syntax in chemical formula '{formula}' near position {pos}")
        pos = match.end()
        element = match.group(1)
        count_str = match.group(2)
        if element not in ELEMENT_MONOISOTOPIC_MASSES:
            raise ValueError(f"Unknown chemical element symbol '{element}' in formula '{formula}'")
        if count_str in ("", "+"):
            count = 1
        elif count_str == "-":
            count = -1
        else:
            count = int(count_str)
        if count != 0:
            counts[element] = counts.get(element, 0) + count
            if counts[element] == 0:
                del counts[element]

    if pos != len(cleaned):
        raise ValueError(f"Invalid trailing characters in chemical formula '{formula}'")

    return counts


def add_formulas(
    base: Mapping[str, int],
    delta: Mapping[str, int],
    scale: int = 1,
    allow_negative: bool = False,
) -> Dict[str, int]:
    """Combine two elemental formula dictionaries: base + scale * delta."""
    combined: Dict[str, int] = dict(base)
    for el, count in delta.items():
        if el not in ELEMENT_MONOISOTOPIC_MASSES:
            raise ValueError(f"Unknown chemical element symbol: '{el}'")
        new_val = combined.get(el, 0) + count * scale
        if new_val < 0 and not allow_negative:
            raise ValueError(
                f"Elemental stoichiometry cannot have negative count for element '{el}' ({new_val})"
            )
        if new_val == 0:
            combined.pop(el, None)
        else:
            combined[el] = new_val
    return combined


def calculate_formula_monoisotopic_mass(
    formula: Union[str, Mapping[str, int]],
    allow_negative: bool = False,
) -> float:
    """Calculate the exact monoisotopic neutral mass (Da) from an elemental formula."""
    parsed = parse_elemental_formula(formula)
    total_mass = 0.0
    for el, count in parsed.items():
        if count < 0 and not allow_negative:
            raise ValueError(f"Negative atom count for '{el}' ({count}) is not allowed for a neutral molecule")
        total_mass += ELEMENT_MONOISOTOPIC_MASSES[el] * count
    return total_mass


def format_hill_formula(formula: Mapping[str, int]) -> str:
    """Format an elemental composition dictionary using the standard Hill system (C, H, then alphabetical)."""
    if not formula:
        return ""
    parts = []
    if "C" in formula and formula["C"] != 0:
        c_cnt = formula["C"]
        parts.append(f"C{c_cnt if c_cnt != 1 else ''}")
        if "H" in formula and formula["H"] != 0:
            h_cnt = formula["H"]
            parts.append(f"H{h_cnt if h_cnt != 1 else ''}")
        remaining = sorted(k for k in formula.keys() if k not in ("C", "H") and formula[k] != 0)
    else:
        remaining = sorted(k for k in formula.keys() if formula[k] != 0)

    for el in remaining:
        cnt = formula[el]
        parts.append(f"{el}{cnt if cnt != 1 else ''}")
    return "".join(parts)

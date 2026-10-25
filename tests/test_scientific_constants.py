"""Unit tests for Stage 1: Authoritative scientific constants and chemical definitions."""

from __future__ import annotations

import unittest

from app.core.formula_calculator import calculate_formula_monoisotopic_mass
from app.core.scientific_constants import (
    ADDUCT_DEFINITIONS,
    AMINO_ACID_RESIDUES,
    CATION_H_MASS,
    CATION_K_MASS,
    CATION_NA_MASS,
    ELECTRON_MASS_DA,
    ELEMENT_MONOISOTOPIC_MASSES,
    IonPolarity,
    LABELING_REGISTRY,
    METHYLENE_MASS_DA,
    MONOSACCHARIDE_RESIDUES,
    PROTON_MASS_DA,
    ValidationStatus,
    WATER_MASS_DA,
)


class TestScientificConstants(unittest.TestCase):
    """Verify authoritative residue, atomic, and ionic mass definitions."""

    def test_required_monosaccharide_residue_masses(self) -> None:
        """Verify the 5 mandatory monosaccharide residue masses to 10 decimal places."""
        self.assertAlmostEqual(MONOSACCHARIDE_RESIDUES["Hex"].residue_mass, 162.0528234315, places=10)
        self.assertAlmostEqual(MONOSACCHARIDE_RESIDUES["HexNAc"].residue_mass, 203.0793725330, places=10)
        self.assertAlmostEqual(MONOSACCHARIDE_RESIDUES["dHex"].residue_mass, 146.0579088090, places=10)
        self.assertAlmostEqual(MONOSACCHARIDE_RESIDUES["Neu5Ac"].residue_mass, 291.0954165270, places=10)
        self.assertAlmostEqual(MONOSACCHARIDE_RESIDUES["Neu5Gc"].residue_mass, 307.0903311490, places=10)

    def test_residue_vs_free_monosaccharide_mass_distinction(self) -> None:
        """Ensure free monosaccharide mass equals residue mass + H2O (18.0105646865 Da)."""
        for sym, res in MONOSACCHARIDE_RESIDUES.items():
            expected_free = res.residue_mass + WATER_MASS_DA
            self.assertAlmostEqual(
                res.free_mass,
                expected_free,
                places=9,
                msg=f"Free mass mismatch for {sym}",
            )
            self.assertGreater(res.free_mass - res.residue_mass, 18.01)

    def test_formula_mass_matches_residue_mass(self) -> None:
        """Verify elemental formulas for each monosaccharide compute to the exact residue mass."""
        for sym in ("Hex", "HexNAc", "dHex", "Neu5Ac", "Neu5Gc"):
            res = MONOSACCHARIDE_RESIDUES[sym]
            calc_mass = calculate_formula_monoisotopic_mass(res.residue_formula)
            self.assertAlmostEqual(
                calc_mass,
                res.residue_mass,
                places=8,
                msg=f"Elemental formula mass does not match residue mass for {sym}",
            )

    def test_proton_vs_neutral_hydrogen_distinction(self) -> None:
        """Ensure proton mass and neutral hydrogen atom mass are strictly distinguished by electron mass."""
        h_atom = ELEMENT_MONOISOTOPIC_MASSES["H"]
        self.assertAlmostEqual(h_atom - PROTON_MASS_DA, ELECTRON_MASS_DA, places=9)
        self.assertAlmostEqual(CATION_H_MASS, PROTON_MASS_DA, places=12)
        self.assertNotAlmostEqual(CATION_H_MASS, h_atom, places=3)

    def test_alkali_metal_cation_masses(self) -> None:
        """Verify Na+ and K+ cation masses subtract the electron mass from neutral atomic masses."""
        self.assertAlmostEqual(
            CATION_NA_MASS,
            ELEMENT_MONOISOTOPIC_MASSES["Na"] - ELECTRON_MASS_DA,
            places=10,
        )
        self.assertAlmostEqual(
            CATION_K_MASS,
            ELEMENT_MONOISOTOPIC_MASSES["K"] - ELECTRON_MASS_DA,
            places=10,
        )

    def test_adduct_registry_completeness(self) -> None:
        """Verify positive and negative adducts are properly registered."""
        self.assertIn("[M+H]+", ADDUCT_DEFINITIONS)
        self.assertIn("[M+Na]+", ADDUCT_DEFINITIONS)
        self.assertIn("[M+K]+", ADDUCT_DEFINITIONS)
        self.assertIn("[M-H]-", ADDUCT_DEFINITIONS)
        self.assertEqual(ADDUCT_DEFINITIONS["[M+Na]+"].polarity, IonPolarity.POSITIVE)
        self.assertEqual(ADDUCT_DEFINITIONS["[M-H]-"].polarity, IonPolarity.NEGATIVE)

    def test_labeling_registry_validation_status(self) -> None:
        """Ensure all registered labeling chemistries declare a valid ValidationStatus."""
        required_codes = (
            "NONE",
            "REDUCED",
            "2-AB",
            "2-AA",
            "PROCAINAMIDE",
            "PA",
            "PERMETHYLATED",
            "PERMETHYLATED_REDUCED",
            "EEC_SIALIC",
            "DMA_SIALIC",
        )
        for code in required_codes:
            self.assertIn(code, LABELING_REGISTRY)
            self.assertEqual(LABELING_REGISTRY[code].validation_status, ValidationStatus.VALIDATED_ENABLED)

    def test_methylene_and_water_masses(self) -> None:
        """Check water and methylene (permethylation unit) exact masses."""
        self.assertAlmostEqual(WATER_MASS_DA, 18.0105646865, places=9)
        self.assertAlmostEqual(METHYLENE_MASS_DA, 14.0156500644, places=9)

    def test_amino_acid_table(self) -> None:
        """Verify canonical 20 amino acids are present and match their elemental formulas."""
        self.assertEqual(len(AMINO_ACID_RESIDUES), 20)
        self.assertAlmostEqual(AMINO_ACID_RESIDUES["N"].monoisotopic_mass, 114.042927447, places=7)


if __name__ == "__main__":
    unittest.main()

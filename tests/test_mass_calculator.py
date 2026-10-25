"""Unit tests for Stage 1: Mass calculator, ion calculator, derivatization, ppm errors, and isotope scoring."""

from __future__ import annotations

import unittest

from app.core.formula_calculator import (
    calculate_formula_monoisotopic_mass,
    format_hill_formula,
    parse_elemental_formula,
)
from app.core.ion_calculator import (
    calculate_all_adduct_mzs,
    calculate_mz_from_neutral_mass,
    calculate_neutral_mass_from_mz,
)
from app.core.isotope_scoring import (
    calculate_isotope_dot_product,
    generate_theoretical_isotope_envelope,
)
from app.core.mass_calculator import (
    calculate_attached_glycan_moiety_mass,
    calculate_glycopeptide_neutral_mass,
    calculate_peptide_neutral_mass,
    calculate_released_glycan_mass,
    format_glycan_composition,
    parse_glycan_composition,
)
from app.core.mass_error import (
    MassToleranceUnit,
    calculate_abs_error_da,
    calculate_abs_error_ppm,
    calculate_signed_error_da,
    calculate_signed_error_ppm,
    evaluate_mass_error,
    get_tolerance_window_bounds,
    is_within_mass_tolerance,
)
from app.core.scientific_constants import (
    CATION_H_MASS,
    CATION_K_MASS,
    CATION_NA_MASS,
    GlycanForm,
    IonPolarity,
    METHYLENE_MASS_DA,
    WATER_MASS_DA,
)
from app.core.validation import (
    validate_human_n_glycan_composition,
    validate_search_parameters,
)


class TestMassAndIonCalculators(unittest.TestCase):
    """Comprehensive verification of Stage 1 mass, ion, derivatization, error, and isotope calculations."""

    def test_composition_parsing_formats(self) -> None:
        """Verify inline, parenthesized, shorthand, and dictionary composition parsing."""
        c1 = parse_glycan_composition("Hex5HexNAc4dHex1Neu5Ac2")
        c2 = parse_glycan_composition("Hex(5)HexNAc(4)Fuc(1)NeuAc(2)")
        c3 = parse_glycan_composition("H5N4F1S2")
        c4 = parse_glycan_composition({"Hex": 5, "HexNAc": 4, "Fuc": 1, "Neu5Ac": 2})
        self.assertEqual(c1, {"Hex": 5, "HexNAc": 4, "dHex": 1, "Neu5Ac": 2})
        self.assertEqual(c1, c2)
        self.assertEqual(c1, c3)
        self.assertEqual(c1, c4)
        self.assertEqual(format_glycan_composition(c1), "Hex5HexNAc4dHex1Neu5Ac2")
        self.assertEqual(format_glycan_composition(c1, use_fuc_alias=True), "Hex5HexNAc4Fuc1Neu5Ac2")

    def test_native_released_glycan_masses_independent_reference(self) -> None:
        """Compare native released N-glycan masses against independent reference calculations."""
        # Man5 (Hex5HexNAc2): 5*162.0528234315 + 2*203.0793725330 + 18.0105646865 = 1234.4334269100 Da
        man5 = calculate_released_glycan_mass("Hex5HexNAc2", labeling="NONE")
        self.assertEqual(man5.glycan_form, GlycanForm.RELEASED_FREE)
        self.assertAlmostEqual(man5.neutral_monoisotopic_mass, 1234.4334269100, places=8)
        self.assertEqual(man5.hill_formula, "C46H78N2O36")

        # G2F / A2G2F1 (Hex5HexNAc4dHex1):
        # 5*162.0528234315 + 4*203.0793725330 + 146.0579088090 + 18.0105646865 = 1786.6500807850 Da
        g2f = calculate_released_glycan_mass("Hex5HexNAc4Fuc1", labeling="NONE")
        self.assertAlmostEqual(g2f.neutral_monoisotopic_mass, 1786.6500807850, places=8)
        self.assertEqual(g2f.hill_formula, "C68H114N4O50")

        # A2G2S2 (Hex5HexNAc4Neu5Ac2):
        # 5*162.0528234315 + 4*203.0793725330 + 2*291.0954165270 + 18.0105646865 = 2222.7830050300 Da
        a2g2s2 = calculate_released_glycan_mass("Hex5HexNAc4Neu5Ac2", labeling="NONE")
        self.assertAlmostEqual(a2g2s2.neutral_monoisotopic_mass, 2222.7830050300, places=8)
        self.assertEqual(a2g2s2.hill_formula, "C84H138N6O62")

    def test_released_vs_peptide_attached_mass_difference(self) -> None:
        """Ensure peptide-attached glycan moiety is strictly 18.0105646865 Da (H2O) lighter than released glycan."""
        released = calculate_released_glycan_mass("Hex5HexNAc4dHex1Neu5Ac2", labeling="NONE")
        attached = calculate_attached_glycan_moiety_mass("Hex5HexNAc4dHex1Neu5Ac2")
        self.assertEqual(attached.glycan_form, GlycanForm.PEPTIDE_ATTACHED)
        self.assertAlmostEqual(
            released.neutral_monoisotopic_mass - attached.neutral_monoisotopic_mass,
            WATER_MASS_DA,
            places=10,
        )
        self.assertEqual(attached.terminal_mass_contribution, 0.0)

    def test_reducing_end_labeling_calculations(self) -> None:
        """Verify 2-AB, 2-AA, and Procainamide reductive amination mass shifts and formulas."""
        native = calculate_released_glycan_mass("Hex5HexNAc2", labeling="NONE")
        ab2 = calculate_released_glycan_mass("Hex5HexNAc2", labeling="2-AB")
        aa2 = calculate_released_glycan_mass("Hex5HexNAc2", labeling="2-AA")
        proa = calculate_released_glycan_mass("Hex5HexNAc2", labeling="PROCAINAMIDE")

        self.assertAlmostEqual(
            ab2.neutral_monoisotopic_mass - native.neutral_monoisotopic_mass,
            120.0687482677,
            places=8,
        )
        self.assertAlmostEqual(
            aa2.neutral_monoisotopic_mass - native.neutral_monoisotopic_mass,
            121.0527638525,
            places=8,
        )
        self.assertAlmostEqual(
            proa.neutral_monoisotopic_mass - native.neutral_monoisotopic_mass,
            219.1735476913,
            places=8,
        )
        # Verify formula mass matches neutral mass for each labeled glycan
        for res in (ab2, aa2, proa):
            self.assertAlmostEqual(
                calculate_formula_monoisotopic_mass(res.elemental_formula),
                res.neutral_monoisotopic_mass,
                places=8,
            )

    def test_permethylation_site_dependent_calculation(self) -> None:
        """Verify stoichiometric permethylation of Man5 (Hex5HexNAc2) non-reduced and reduced."""
        # Man5: 5 Hex (5*3=15 sites) + 2 HexNAc (2*3=6 sites) = 21 residue sites.
        # Non-reduced adds +2 terminal sites = 23 CH2 groups (+23 * 14.0156500644 = +322.35995148 Da).
        # Total neutral mass = 1234.43342691 + 322.35995148 = 1556.79337839 Da.
        # [M+Na]+ m/z = 1556.79337839 + 22.98922070 = 1579.78259909 (classic Man5 permethylated [M+Na]+ @ 1579.78).
        perm = calculate_released_glycan_mass("Hex5HexNAc2", labeling="PERMETHYLATED")
        self.assertEqual(perm.permethylation_sites_total, 23)
        self.assertAlmostEqual(
            perm.derivatization_mass_shift,
            23 * METHYLENE_MASS_DA,
            places=8,
        )
        self.assertAlmostEqual(perm.neutral_monoisotopic_mass, 1556.79337839, places=7)

        mz_na = calculate_mz_from_neutral_mass(perm.neutral_monoisotopic_mass, "[M+Na]+")
        self.assertAlmostEqual(mz_na.theoretical_mz, 1579.78259909, places=7)

        # Reduced permethylated Man5: +H2 (+2.015650 Da) and 21 + 3 = 24 CH2 groups
        perm_red = calculate_released_glycan_mass("Hex5HexNAc2", labeling="PERMETHYLATED_REDUCED")
        self.assertEqual(perm_red.permethylation_sites_total, 24)
        self.assertAlmostEqual(
            perm_red.neutral_monoisotopic_mass - perm.neutral_monoisotopic_mass,
            2 * 1.00782503219 + METHYLENE_MASS_DA,
            places=8,
        )

    def test_linkage_specific_sialic_acid_derivatization(self) -> None:
        """Verify ETH-EDC (EEC_SIALIC) distinguishes alpha2,3 (-18.0106 Da) from alpha2,6 (+28.0313 Da)."""
        # A2G2S2 (Hex5HexNAc4Neu5Ac2) with:
        # Case 1: two alpha2,6 sialic acids
        s2_26_26 = calculate_released_glycan_mass(
            "Hex5HexNAc4Neu5Ac2",
            labeling="EEC_SIALIC",
            alpha23_count=0,
            alpha26_count=2,
        )
        # Case 2: one alpha2,3 and one alpha2,6 sialic acid
        s2_23_26 = calculate_released_glycan_mass(
            "Hex5HexNAc4Neu5Ac2",
            labeling="EEC_SIALIC",
            alpha23_count=1,
            alpha26_count=1,
        )
        # Case 3: two alpha2,3 sialic acids
        s2_23_23 = calculate_released_glycan_mass(
            "Hex5HexNAc4Neu5Ac2",
            labeling="EEC_SIALIC",
            alpha23_count=2,
            alpha26_count=0,
        )

        # Mass difference between alpha2,6 ethyl ester (+C2H4) and alpha2,3 lactone (-H2O) is +46.041865 Da
        expected_step = 28.0313001288 - (-18.0105646865)
        self.assertAlmostEqual(
            s2_26_26.neutral_monoisotopic_mass - s2_23_26.neutral_monoisotopic_mass,
            expected_step,
            places=7,
        )
        self.assertAlmostEqual(
            s2_23_26.neutral_monoisotopic_mass - s2_23_23.neutral_monoisotopic_mass,
            expected_step,
            places=7,
        )

    def test_ion_and_adduct_calculations(self) -> None:
        """Verify singly and multiply charged positive and negative ion m/z calculations."""
        neutral = 1786.6500807850  # G2F neutral mass

        ion_h = calculate_mz_from_neutral_mass(neutral, "[M+H]+")
        ion_na = calculate_mz_from_neutral_mass(neutral, "[M+Na]+")
        ion_k = calculate_mz_from_neutral_mass(neutral, "[M+K]+")
        ion_2h = calculate_mz_from_neutral_mass(neutral, "[M+2H]2+")
        ion_mh_neg = calculate_mz_from_neutral_mass(neutral, "[M-H]-")

        self.assertAlmostEqual(ion_h.theoretical_mz, neutral + CATION_H_MASS, places=9)
        self.assertAlmostEqual(ion_na.theoretical_mz, neutral + CATION_NA_MASS, places=9)
        self.assertAlmostEqual(ion_k.theoretical_mz, neutral + CATION_K_MASS, places=9)
        self.assertAlmostEqual(ion_2h.theoretical_mz, (neutral + 2 * CATION_H_MASS) / 2.0, places=9)
        self.assertAlmostEqual(ion_mh_neg.theoretical_mz, neutral - CATION_H_MASS, places=9)

        # Round-trip neutral mass inversion
        back_neutral = calculate_neutral_mass_from_mz(ion_2h.theoretical_mz, "[M+2H]2+")
        self.assertAlmostEqual(back_neutral, neutral, places=9)

        all_pos = calculate_all_adduct_mzs(neutral, polarity=IonPolarity.POSITIVE)
        self.assertIn("[M+Na]+", all_pos)
        self.assertNotIn("[M-H]-", all_pos)

    def test_mass_error_and_tolerance_filtering(self) -> None:
        """Verify Da and ppm signed/absolute mass error and tolerance boundary checks."""
        theo = 2000.0000
        obs_high = 2000.0200  # +0.0200 Da = +10.0 ppm
        obs_low = 1999.9700   # -0.0300 Da = -15.0 ppm

        self.assertAlmostEqual(calculate_signed_error_da(obs_high, theo), 0.0200, places=8)
        self.assertAlmostEqual(calculate_abs_error_da(obs_low, theo), 0.0300, places=8)
        self.assertAlmostEqual(calculate_signed_error_ppm(obs_high, theo), 10.0, places=6)
        self.assertAlmostEqual(calculate_signed_error_ppm(obs_low, theo), -15.0, places=6)
        self.assertAlmostEqual(calculate_abs_error_ppm(obs_low, theo), 15.0, places=6)

        err_obj = evaluate_mass_error(obs_low, theo)
        self.assertAlmostEqual(err_obj.signed_error_ppm, -15.0, places=6)
        self.assertAlmostEqual(err_obj.abs_error_ppm, 15.0, places=6)

        self.assertTrue(is_within_mass_tolerance(obs_high, theo, 10.0, MassToleranceUnit.PPM))
        self.assertFalse(is_within_mass_tolerance(obs_low, theo, 10.0, MassToleranceUnit.PPM))
        self.assertTrue(is_within_mass_tolerance(obs_low, theo, 15.0, MassToleranceUnit.PPM))
        self.assertTrue(is_within_mass_tolerance(obs_low, theo, 0.03, MassToleranceUnit.DA))

        low_b, high_b = get_tolerance_window_bounds(1000.0, 20.0, MassToleranceUnit.PPM)
        self.assertAlmostEqual(low_b, 999.98, places=8)
        self.assertAlmostEqual(high_b, 1000.02, places=8)

    def test_isotope_envelope_and_idotp_scoring(self) -> None:
        """Verify theoretical isotope envelope and strict unavailability on single centroid peaks."""
        man5 = calculate_released_glycan_mass("Hex5HexNAc2", labeling="NONE")
        ion = calculate_mz_from_neutral_mass(
            man5.neutral_monoisotopic_mass,
            adduct_code="[M+Na]+",
            neutral_formula=man5.elemental_formula,
        )
        assert ion.ion_formula is not None
        envelope = generate_theoretical_isotope_envelope(
            ion_formula=ion.ion_formula,
            monoisotopic_mz=ion.theoretical_mz,
            charge=1,
            num_peaks=4,
        )
        self.assertEqual(len(envelope), 4)
        self.assertAlmostEqual(envelope[0].mz, ion.theoretical_mz, places=8)
        # Spacing between M0 and M+1 for z=1 is ~1.00335 Da
        self.assertAlmostEqual(envelope[1].mz - envelope[0].mz, 1.00335, places=3)
        self.assertAlmostEqual(sum(p.relative_abundance for p in envelope), 1.0, places=6)

        # Case A: Single centroid peak -> idotp MUST be None and available=False
        single_peak_score = calculate_isotope_dot_product(
            observed_mzs=[ion.theoretical_mz],
            observed_intensities=[5000.0],
            ion_formula=ion.ion_formula,
            monoisotopic_mz=ion.theoretical_mz,
            charge=1,
        )
        self.assertFalse(single_peak_score.available)
        self.assertIsNone(single_peak_score.idotp)
        self.assertEqual(
            single_peak_score.status_message,
            "Insufficient isotope information for isotope scoring",
        )

        # Case B: Observed multiplet matching theoretical ratios -> idotp ~ 1.0
        obs_mzs = [p.mz for p in envelope]
        obs_ints = [p.relative_abundance * 10000.0 for p in envelope]
        multiplet_score = calculate_isotope_dot_product(
            observed_mzs=obs_mzs,
            observed_intensities=obs_ints,
            ion_formula=ion.ion_formula,
            monoisotopic_mz=ion.theoretical_mz,
            charge=1,
            tolerance=10.0,
        )
        self.assertTrue(multiplet_score.available)
        assert multiplet_score.idotp is not None
        self.assertGreater(multiplet_score.idotp, 0.999)

    def test_glycopeptide_neutral_mass_calculation(self) -> None:
        """Verify IgG1 EEQYNSTYR glycopeptide neutral mass uses attached glycan moiety convention."""
        # Peptide EEQYNSTYR = 1188.504731 Da
        # Attached G0F (Hex3HexNAc4dHex1) = 3*162.0528234 + 4*203.0793725 + 146.0579088 = 1444.533869 Da
        # Total neutral mass = 2633.038600 Da (NOT +18.010565 Da higher!)
        total_mass, pep_res, gly_res, comb_formula = calculate_glycopeptide_neutral_mass(
            peptide_sequence="EEQYNSTYR",
            glycan_composition="Hex3HexNAc4dHex1",
        )
        self.assertAlmostEqual(pep_res.neutral_monoisotopic_mass, 1188.504731, places=5)
        self.assertAlmostEqual(gly_res.neutral_monoisotopic_mass, 1444.533869, places=5)
        self.assertAlmostEqual(total_mass, 2633.038600, places=5)
        self.assertAlmostEqual(
            calculate_formula_monoisotopic_mass(comb_formula),
            total_mass,
            places=7,
        )

    def test_formula_and_validation_utilities(self) -> None:
        """Test formula parsing and biosynthetic validation rules."""
        f = parse_elemental_formula("C6H12O6")
        self.assertEqual(format_hill_formula(f), "C6H12O6")

        rep_valid = validate_human_n_glycan_composition({"Hex": 5, "HexNAc": 4, "dHex": 1, "Neu5Ac": 2})
        self.assertTrue(rep_valid.is_valid)
        self.assertEqual(len(rep_valid.warnings), 0)

        rep_neu5gc = validate_human_n_glycan_composition({"Hex": 5, "HexNAc": 4, "Neu5Gc": 1})
        self.assertTrue(rep_neu5gc.is_valid)
        self.assertTrue(any("Neu5Gc" in w for w in rep_neu5gc.warnings))

        param_rep = validate_search_parameters(
            tolerance=20.0,
            tolerance_unit=MassToleranceUnit.PPM,
            polarity=IonPolarity.POSITIVE,
            adduct_codes=["[M+Na]+"],
            charge_states=[1],
            labeling="NONE",
        )
        self.assertTrue(param_rep.is_valid)


if __name__ == "__main__":
    unittest.main()

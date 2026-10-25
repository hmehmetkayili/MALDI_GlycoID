"""Unit tests for Stage 6: Glycopeptide screening engine (digestion, sequon detection, attached glycan moiety mass, MS1/MS2 matching)."""

from __future__ import annotations

import unittest

from app.analysis import (
    GlycopeptideConfidenceCategory,
    GlycopeptideSearchSettings,
    PrecursorObservation,
    ProteaseEnzyme,
    digest_protein_in_silico,
    find_n_glycosylation_sequons,
    parse_fasta_sequences,
    screen_glycopeptides,
)
from app.config import EXAMPLES_DIR
from app.core.ion_calculator import calculate_mz_from_neutral_mass
from app.core.mass_calculator import (
    calculate_glycopeptide_neutral_mass,
    calculate_released_glycan_mass,
)
from app.core.mass_error import MassToleranceUnit
from app.core.scientific_constants import WATER_MASS_DA
from app.parsers import parse_txt_spectrum


class TestGlycopeptideScreening(unittest.TestCase):
    """Comprehensive verification of Stage 6 glycopeptide digestion, sequon rules, and precursor matching."""

    def test_canonical_sequon_detection_including_proline_exclusion(self) -> None:
        """Verify N-X-S/T (X != P) detection, overlapping sequons, and proline rejection."""
        # 'EEQYNSTYR' -> N at position 5 ('NST')
        s1 = find_n_glycosylation_sequons("EEQYNSTYR", peptide_start_in_protein_1based=57)
        self.assertEqual(len(s1), 1)
        self.assertEqual(s1[0].motif, "NST")
        self.assertEqual(s1[0].peptide_position_1based, 5)
        self.assertEqual(s1[0].protein_position_1based, 61)

        # 'EEQYNPTYR' -> X is Proline ('P'), MUST NOT be detected as an N-glycosylation sequon!
        s_pro = find_n_glycosylation_sequons("EEQYNPTYR")
        self.assertEqual(len(s_pro), 0)

        # Multiple sequons: 'ANSTGANKT' -> 'NST' at pos 2 and 'NKT' at pos 7
        s_multi = find_n_glycosylation_sequons("ANSTGANKT", peptide_start_in_protein_1based=10)
        self.assertEqual(len(s_multi), 2)
        self.assertEqual(s_multi[0].peptide_position_1based, 2)
        self.assertEqual(s_multi[1].peptide_position_1based, 7)

    def test_in_silico_tryptic_digestion_and_modifications(self) -> None:
        """Verify tryptic digestion (K/R not before P), missed cleavages, and Cys/Met modifications."""
        proteins = parse_fasta_sequences(EXAMPLES_DIR / "synthetic_igg1_fc.fasta")
        self.assertEqual(len(proteins), 1)

        # Note in synthetic_igg1_fc.fasta: '...HNAKTKPREEQYNSTYRVVSVL...'
        # Notice 'KPR': K followed by P is NOT cleaved by Trypsin! Cleavage happens after R -> 'EEQYNSTYR'
        peps_0mc = digest_protein_in_silico(
            proteins[0],
            GlycopeptideSearchSettings(
                enzyme=ProteaseEnzyme.TRYPSIN,
                max_missed_cleavages=0,
                require_n_glycan_sequon=True,
            ),
        )
        seqs_0mc = [p.peptide_sequence for p in peps_0mc]
        self.assertIn("EEQYNSTYR", seqs_0mc)

        # With max_missed_cleavages=1, we also get 'TKPREEQYNSTYR'
        peps_1mc = digest_protein_in_silico(
            proteins[0],
            GlycopeptideSearchSettings(
                enzyme=ProteaseEnzyme.TRYPSIN,
                max_missed_cleavages=1,
                require_n_glycan_sequon=True,
            ),
        )
        seqs_1mc = [p.peptide_sequence for p in peps_1mc]
        self.assertIn("TKPREEQYNSTYR", seqs_1mc)

        # Test Cys Carbamidomethyl (+57.021464 Da) and Met Oxidation (+15.994915 Da)
        mod_prot = parse_fasta_sequences(">MOD_TEST\nACMNSTYR\n")[0]
        mod_peps = digest_protein_in_silico(
            mod_prot,
            GlycopeptideSearchSettings(
                fixed_carbamidomethyl_cys=True,
                variable_oxidation_met=True,
            ),
        )
        # Should generate 2 forms: Carbamidomethyl_C x1 (unoxidized M) and Carbamidomethyl_C x1 + Oxidation_M x1
        self.assertEqual(len(mod_peps), 2)
        self.assertAlmostEqual(
            mod_peps[1].peptide_neutral_mass - mod_peps[0].peptide_neutral_mass,
            15.99491462,
            places=6,
        )

    def test_glycopeptide_attached_moiety_vs_released_glycan_water_guard(self) -> None:
        """Verify that screening uses the attached glycan moiety mass and rejects released-glycan +H2O mass."""
        # Peptide EEQYNSTYR (1188.504731 Da) + attached G0F (Hex3HexNAc4dHex1, 1444.533869 Da)
        # True glycopeptide neutral mass = 2633.038600 Da
        # Singly protonated [M+H]+ m/z = 2633.038600 + 1.00727645 = 2634.045876 Da
        gp_neutral, pep_res, gly_res, _ = calculate_glycopeptide_neutral_mass("EEQYNSTYR", "Hex3HexNAc4dHex1")
        true_mh_mz = calculate_mz_from_neutral_mass(gp_neutral, "[M+H]+").theoretical_mz

        # If someone mistakenly added released G0F (which has +18.010565 Da H2O), m/z would be 2652.056441 Da
        rel_gly = calculate_released_glycan_mass("Hex3HexNAc4dHex1", labeling="NONE")
        wrong_mh_mz = calculate_mz_from_neutral_mass(
            pep_res.neutral_monoisotopic_mass + rel_gly.neutral_monoisotopic_mass,
            "[M+H]+",
        ).theoretical_mz
        self.assertAlmostEqual(wrong_mh_mz - true_mh_mz, WATER_MASS_DA, places=7)

        res = screen_glycopeptides(
            precursors=[
                PrecursorObservation(sample_id="IgG1_Sample", observed_mz=true_mh_mz, known_charge=1),
                PrecursorObservation(sample_id="IgG1_Sample", observed_mz=wrong_mh_mz, known_charge=1),
            ],
            fasta_or_proteins=EXAMPLES_DIR / "synthetic_igg1_fc.fasta",
            glycan_compositions=["Hex3HexNAc4dHex1"],
            settings=GlycopeptideSearchSettings(
                tolerance=10.0,
                tolerance_unit=MassToleranceUnit.PPM,
                adduct_codes=("[M+H]+",),
                charge_states=(1,),
            ),
        )
        # Only the true attached-moiety precursor (2634.0459) matches; the +18 Da wrong precursor is rejected!
        self.assertEqual(res.matched_precursors_count, 1)
        self.assertEqual(len(res.matches), 1)
        m = res.matches[0]
        self.assertEqual(m.peptide_sequence, "EEQYNSTYR")
        self.assertEqual(m.glycan_composition, "Hex3HexNAc4dHex1")
        self.assertAlmostEqual(m.observed_precursor_mz, true_mh_mz, places=6)
        self.assertLess(m.abs_error_ppm, 0.5)
        self.assertEqual(m.confidence_category, GlycopeptideConfidenceCategory.MS1_COMPOSITION_CANDIDATE)

    def test_ms2_oxonium_and_y_ion_evidence_and_confidence_categories(self) -> None:
        """Verify MSMS_SUPPORTED and SITE_LOCALIZED categories when oxonium + Y1 ions are present."""
        gp_neutral, pep_res, _, _ = calculate_glycopeptide_neutral_mass("EEQYNSTYR", "Hex3HexNAc4dHex1")
        prec_3h_mz = calculate_mz_from_neutral_mass(gp_neutral, "[M+3H]3+").theoretical_mz

        # Build synthetic MS2 spectrum containing:
        # - HexNAc+ oxonium ion (204.0866)
        # - HexHexNAc+ oxonium ion (366.1395)
        # - Y0 intact peptide [M+H]+ (1188.5047 + 1.0073 = 1189.5120)
        # - Y1 [Peptide+HexNAc+H]+ (1188.5047 + 203.0794 + 1.0073 = 1392.5914)
        ms2_txt = """m/z intensity
204.0866 15000.0
366.1395 12000.0
1189.5120 4500.0
1392.5914 8200.0
"""
        ms2_spec = parse_txt_spectrum(ms2_txt, filename="ms2_hcd.txt")

        res_ms2 = screen_glycopeptides(
            precursors=[
                PrecursorObservation(
                    sample_id="IgG1_LCMS",
                    observed_mz=prec_3h_mz,
                    known_charge=3,
                    retention_time_sec=845.0,
                    ms2_spectrum=ms2_spec,
                )
            ],
            peptide_sequences=["EEQYNSTYR"],
            glycan_compositions=["Hex3HexNAc4dHex1"],
            settings=GlycopeptideSearchSettings(
                tolerance=10.0,
                adduct_codes=("[M+H]+",),
                charge_states=(3,),
            ),
        )
        self.assertEqual(len(res_ms2.matches), 1)
        hit = res_ms2.matches[0]
        # Single sequon + verified MS/MS oxonium & Y0/Y1 ions -> SITE_LOCALIZED
        self.assertEqual(hit.confidence_category, GlycopeptideConfidenceCategory.SITE_LOCALIZED)
        self.assertGreaterEqual(len(hit.msms_oxonium_matched), 2)
        self.assertGreaterEqual(len(hit.msms_y_ions_matched), 2)

        # Now test a peptide with 2 sequons ('ANSTGANKT') with the same MS/MS oxonium + Y0 ions:
        # Because there are 2 sequons, MS/MS oxonium + Y0 proves glycopeptide identity (MSMS_SUPPORTED)
        # but cannot localize which of the 2 sequons carries the glycan without site-specific c/z ions!
        gp2_neutral, pep2_res, _, _ = calculate_glycopeptide_neutral_mass("ANSTGANKT", "Hex5HexNAc2")
        prec2_mz = calculate_mz_from_neutral_mass(gp2_neutral, "[M+2H]2+").theoretical_mz
        y0_pep2 = pep2_res.neutral_monoisotopic_mass + 1.007276
        ms2_2seq = parse_txt_spectrum(f"m/z intensity\n204.0866 9000.0\n{y0_pep2:.4f} 5000.0\n", filename="ms2_2seq.txt")

        res_2seq = screen_glycopeptides(
            precursors=[
                PrecursorObservation(
                    sample_id="MultiSequon",
                    observed_mz=prec2_mz,
                    known_charge=2,
                    ms2_spectrum=ms2_2seq,
                )
            ],
            peptide_sequences=["ANSTGANKT"],
            glycan_compositions=["Hex5HexNAc2"],
            settings=GlycopeptideSearchSettings(adduct_codes=("[M+H]+",), charge_states=(2,)),
        )
        self.assertEqual(len(res_2seq.matches), 1)
        self.assertEqual(res_2seq.matches[0].sequon_count_in_peptide, 2)
        self.assertEqual(res_2seq.matches[0].confidence_category, GlycopeptideConfidenceCategory.MSMS_SUPPORTED)
        self.assertTrue(any("2 candidate sequons" in w for w in res_2seq.matches[0].warnings))


if __name__ == "__main__":
    unittest.main()

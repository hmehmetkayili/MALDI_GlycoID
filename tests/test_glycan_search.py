"""Unit tests for Stage 4: GlycoQuest-like MS1 search and matching engine."""

from __future__ import annotations

import unittest

from app.analysis import (
    GlycanEvidenceLevel,
    MS1SearchSettings,
    PeakPreprocessingConfig,
    SearchMode,
    TargetedQueryItem,
    preprocess_spectrum,
    search_ms1_library,
    search_ms1_targeted,
)
from app.config import EXAMPLES_DIR
from app.core.ion_calculator import calculate_mz_from_neutral_mass
from app.core.mass_calculator import calculate_released_glycan_mass
from app.core.mass_error import MassToleranceUnit
from app.core.scientific_constants import IonPolarity
from app.libraries import CompositionConstraints, generate_human_n_glycan_library
from app.parsers import SpectrumDataType, parse_spectrum_file, parse_txt_spectrum


class TestMS1GlycanSearch(unittest.TestCase):
    """Comprehensive tests for Mode A (targeted), Mode B (library), tolerances, and ambiguity reporting."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.library = generate_human_n_glycan_library(
            constraints=CompositionConstraints(
                min_hex=3,
                max_hex=9,
                min_hexnac=2,
                max_hexnac=6,
                min_dhex=0,
                max_dhex=2,
                min_neu5ac=0,
                max_neu5ac=3,
                include_neu5gc=True,
                max_neu5gc=1,
                labeling="NONE",
            ),
            library_name="Stage 4 Test Library",
            library_version="1.0.0",
        )

    def test_mode_b_library_search_synthetic_spectrum(self) -> None:
        """Test Mode B library search on the synthetic native N-glycan spectrum (with M+1 isotopes)."""
        spec = parse_spectrum_file(EXAMPLES_DIR / "synthetic_native_nglycans_centroid.txt")
        settings = MS1SearchSettings(
            search_mode=SearchMode.LIBRARY,
            tolerance=15.0,
            tolerance_unit=MassToleranceUnit.PPM,
            polarity=IonPolarity.POSITIVE,
            adduct_codes=("[M+Na]+",),
            charge_states=(1,),
            labeling="NONE",
            enable_isotope_scoring=True,
        )
        res = search_ms1_library(spectrum=spec, library=self.library, settings=settings)

        # The synthetic file has 6 monoisotopic peaks (1257.4226, 1485.5336, 1647.5865, 1809.6393, 1954.6768, 2245.7722)
        # plus 4 M+1 isotopic peaks (1258.4260, 1486.5370, 1648.5899, 1810.6427).
        self.assertEqual(res.matched_peaks_count, 6)
        self.assertEqual(res.unmatched_peaks_count, 4)

        matched_comps = {m.glycan_composition for m in res.top_matches_per_peak}
        self.assertEqual(
            matched_comps,
            {
                "Hex5HexNAc2",
                "Hex3HexNAc4dHex1",
                "Hex4HexNAc4dHex1",
                "Hex5HexNAc4dHex1",
                "Hex5HexNAc4Neu5Ac1",
                "Hex5HexNAc4Neu5Ac2",
            },
        )

        # Verify isotope scoring worked on Man5 (which has M0=1257.4226 and M+1=1258.4260 in the file)
        man5_match = next(m for m in res.top_matches_per_peak if m.glycan_composition == "Hex5HexNAc2")
        self.assertIsNotNone(man5_match.isotope_score_idotp)
        assert man5_match.isotope_score_idotp is not None
        self.assertGreater(man5_match.isotope_score_idotp, 0.90)
        self.assertEqual(man5_match.evidence_level, GlycanEvidenceLevel.MS1_ISOTOPE_SUPPORTED_CANDIDATE)

        # Verify A2G2S2 (2245.7722) has NO M+1 peak in the file -> idotp is None (unavailable)
        s2_match = next(m for m in res.top_matches_per_peak if m.glycan_composition == "Hex5HexNAc4Neu5Ac2")
        self.assertIsNone(s2_match.isotope_score_idotp)
        self.assertEqual(s2_match.isotope_score_status, "Insufficient isotope information for isotope scoring")
        self.assertEqual(s2_match.evidence_level, GlycanEvidenceLevel.MS1_COMPOSITION_CANDIDATE)

    def test_mode_a_targeted_search_all_query_types(self) -> None:
        """Test Mode A targeted search using composition, library_id, neutral_mass, and expected_mz."""
        spec = parse_spectrum_file(EXAMPLES_DIR / "synthetic_native_nglycans.xy")
        targets = [
            TargetedQueryItem(
                query_label="Man5 by composition",
                target_type="composition",
                composition="Hex5HexNAc2",
            ),
            TargetedQueryItem(
                query_label="G2F by neutral mass",
                target_type="neutral_mass",
                neutral_mass=1786.65008079,
            ),
            TargetedQueryItem(
                query_label="A2G2S2 by expected m/z",
                target_type="expected_mz",
                expected_mz=2245.7722,
                adduct_code="[M+Na]+",
                charge=1,
            ),
        ]
        settings = MS1SearchSettings(
            search_mode=SearchMode.TARGETED,
            tolerance=10.0,
            tolerance_unit=MassToleranceUnit.PPM,
            adduct_codes=("[M+Na]+",),
            charge_states=(1,),
            labeling="NONE",
        )
        res = search_ms1_targeted(spectrum=spec, targets=targets, settings=settings)
        self.assertEqual(res.matched_peaks_count, 3)

    def test_tolerance_boundaries_and_no_matches_warning(self) -> None:
        """Verify strict tolerance rejection outside boundary and explicit 'No candidates found' warning."""
        # G2F [M+Na]+ exact m/z = 1809.63930149
        # At +25 ppm -> 1809.63930149 * (1 + 25e-6) = 1809.684542
        shifted_txt = "m/z intensity\n1809.684542 10000.0\n"
        spec = parse_txt_spectrum(shifted_txt, filename="shifted.txt")

        # 1. Search at 20 ppm -> should reject
        res_20 = search_ms1_library(
            spectrum=spec,
            library=self.library,
            settings=MS1SearchSettings(tolerance=20.0, tolerance_unit=MassToleranceUnit.PPM),
        )
        self.assertEqual(res_20.matched_peaks_count, 0)
        self.assertIn("No candidates found within the selected tolerance", res_20.warnings)

        # 2. Search at 30 ppm -> should match G2F with ~ +25.0 ppm error
        res_30 = search_ms1_library(
            spectrum=spec,
            library=self.library,
            settings=MS1SearchSettings(tolerance=30.0, tolerance_unit=MassToleranceUnit.PPM),
        )
        self.assertEqual(res_30.matched_peaks_count, 1)
        self.assertEqual(res_30.top_matches_per_peak[0].glycan_composition, "Hex5HexNAc4dHex1")
        self.assertAlmostEqual(res_30.top_matches_per_peak[0].signed_error_ppm, 25.0, places=1)

    def test_adduct_and_composition_ambiguity_reporting(self) -> None:
        """Verify near-isobaric overlap between [M+K]+ of fucosylated and [M+Na]+ of non-fucosylated glycans."""
        # Compare:
        # Candidate 1: Hex4HexNAc4dHex1 (G1F) + [M+K]+
        #   Neutral = 1624.59725735 Da, [M+K]+ = 1624.59725735 + 38.96315791 = 1663.56041526
        # Candidate 2: Hex5HexNAc4 (G2, non-fucosylated) + [M+Na]+
        #   Neutral = 1640.59217198 Da, [M+Na]+ = 1640.59217198 + 22.98922070 = 1663.58139268
        # Difference is 0.02098 Da (~12.6 ppm)!
        # Also Candidate 3: Hex4HexNAc4Neu5Gc1 vs Hex5HexNAc4Neu5Ac1 (identical elemental formula C73H121N5O54!)
        #   Hex5 (5*C6H10O5 = C30H50O25) + Neu5Ac1 (C11H17NO8) = C41H67NO33
        #   Hex4dHex1 (C24H40O20 + C6H10O4 = C30H50O24) + Neu5Gc1 (C11H17NO9) = C41H67NO33 -> EXACT ISOBARS (0.00 ppm)!
        midpoint_mz = 1663.5709
        ambig_txt = f"m/z intensity\n{midpoint_mz:.6f} 12000.0\n1954.6768 15000.0\n"
        spec = parse_txt_spectrum(ambig_txt, filename="ambig.txt")

        settings = MS1SearchSettings(
            tolerance=15.0,
            tolerance_unit=MassToleranceUnit.PPM,
            adduct_codes=("[M+Na]+", "[M+K]+"),
            charge_states=(1,),
            labeling="NONE",
        )
        res = search_ms1_library(spectrum=spec, library=self.library, settings=settings)

        # Both peaks should be flagged as ambiguous (one due to [M+Na]+ vs [M+K]+, one due to Neu5Ac+Hex vs Neu5Gc+dHex)
        self.assertEqual(res.matched_peaks_count, 2)
        self.assertEqual(res.ambiguous_peaks_count, 2)

        peak1_matches = [m for m in res.matches if m.peak_index == 1]
        self.assertGreaterEqual(len(peak1_matches), 2)
        self.assertTrue(peak1_matches[0].is_ambiguous)
        self.assertGreaterEqual(len(peak1_matches[0].alternative_candidates), 1)
        self.assertTrue(any("Ambiguous peak assignment" in w for w in peak1_matches[0].warnings))

        # Peak 2 (1954.6768) matches both Hex5HexNAc4Neu5Ac1 and Hex4HexNAc4dHex1Neu5Gc1
        peak2_comps = {m.glycan_composition for m in res.matches if m.peak_index == 2}
        self.assertIn("Hex5HexNAc4Neu5Ac1", peak2_comps)
        self.assertIn("Hex4HexNAc4dHex1Neu5Gc1", peak2_comps)

    def test_multiply_charged_and_labeled_glycan_search(self) -> None:
        """Verify 2-AB labeled and doubly charged [M+2H]2+ search calculations."""
        # Calculate exact m/z for 2-AB labeled Man5 (Hex5HexNAc2) as [M+2H]2+
        man5_2ab = calculate_released_glycan_mass("Hex5HexNAc2", labeling="2-AB")
        ion_2h = calculate_mz_from_neutral_mass(man5_2ab.neutral_monoisotopic_mass, "[M+2H]2+")

        spec = parse_txt_spectrum(f"m/z intensity\n{ion_2h.theoretical_mz:.6f} 9500.0\n", filename="man5_2ab.txt")
        settings = MS1SearchSettings(
            tolerance=5.0,
            tolerance_unit=MassToleranceUnit.PPM,
            polarity=IonPolarity.POSITIVE,
            adduct_codes=("[M+H]+",),
            charge_states=(2,),
            labeling="2-AB",
        )
        res = search_ms1_library(spectrum=spec, library=self.library, settings=settings)
        self.assertEqual(res.matched_peaks_count, 1)
        top = res.top_matches_per_peak[0]
        self.assertEqual(top.glycan_composition, "Hex5HexNAc2")
        self.assertEqual(top.assigned_adduct, "[M+2H]2+")
        self.assertEqual(top.charge_state, 2)
        self.assertLess(top.abs_error_ppm, 0.5)

    def test_centroid_vs_profile_preprocessing_guard(self) -> None:
        """Verify requesting profile peak integration on a centroid list records the required warning."""
        spec = parse_spectrum_file(EXAMPLES_DIR / "synthetic_native_nglycans.xy")
        self.assertEqual(spec.data_type, SpectrumDataType.CENTROID)
        proc = preprocess_spectrum(
            spec,
            PeakPreprocessingConfig(apply_profile_baseline_subtraction=True, apply_profile_peak_detection=True),
        )
        self.assertFalse(proc.was_profile_processed)
        self.assertTrue(
            any("Insufficient profile data for peak-area calculation" in w for w in proc.warnings)
        )


if __name__ == "__main__":
    unittest.main()

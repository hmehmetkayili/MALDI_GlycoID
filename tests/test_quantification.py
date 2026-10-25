"""Unit tests for Stage 5: Glycan intensity, relative abundance, peak area, aggregation, and replicate quantification."""

from __future__ import annotations

import unittest

from app.analysis import (
    DenominatorSetPolicy,
    GlycanAggregationPolicy,
    MS1SearchSettings,
    quantify_ms1_search_result,
    search_ms1_library,
    summarize_replicates_and_build_matrix,
)
from app.config import EXAMPLES_DIR
from app.core.mass_error import MassToleranceUnit
from app.libraries import CompositionConstraints, generate_human_n_glycan_library
from app.parsers import parse_spectrum_file, parse_txt_spectrum


class TestGlycanQuantification(unittest.TestCase):
    """Comprehensive verification of Stage 5 quantification, aggregation rules, and replicate statistics."""

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
                include_neu5gc=False,
                labeling="NONE",
            ),
            library_name="Stage 5 Quantification Test Library",
            library_version="1.0.0",
        )

    def test_base_peak_normalization_and_denominator_policies(self) -> None:
        """Verify base-peak relative intensity under ALL_ELIGIBLE vs ANNOTATED_GLYCAN_PEAKS_ONLY."""
        # Create a spectrum with:
        # - Man5 [M+Na]+ (1257.4226) at intensity 5000
        # - G2F [M+Na]+ (1809.6393) at intensity 10000 (max annotated peak)
        # - An unmatched matrix/contaminant peak (999.9999) at intensity 20000 (overall spectrum base peak)
        raw_txt = "m/z intensity\n999.9999 20000.0\n1257.4226 5000.0\n1809.6393 10000.0\n"
        spec = parse_txt_spectrum(raw_txt, filename="denom_test.txt")
        search_res = search_ms1_library(
            spectrum=spec,
            library=self.library,
            settings=MS1SearchSettings(tolerance=15.0, adduct_codes=("[M+Na]+",)),
        )

        # Policy 1: ALL_ELIGIBLE_SPECTRUM_PEAKS (denominator = 20000.0)
        q_all = quantify_ms1_search_result(
            search_res,
            denominator_policy=DenominatorSetPolicy.ALL_ELIGIBLE_SPECTRUM_PEAKS,
        )
        self.assertAlmostEqual(q_all.denominator_max_intensity, 20000.0, places=2)
        g2f_rec = next(g for g in q_all.glycan_abundances if g.glycan_composition == "Hex5HexNAc4dHex1")
        man5_rec = next(g for g in q_all.glycan_abundances if g.glycan_composition == "Hex5HexNAc2")
        self.assertAlmostEqual(g2f_rec.base_peak_relative_intensity_pct, 50.0, places=4)
        self.assertAlmostEqual(man5_rec.base_peak_relative_intensity_pct, 25.0, places=4)
        # Relative abundance % (sum of annotated = 15000): G2F = 66.6667%, Man5 = 33.3333%
        self.assertAlmostEqual(g2f_rec.relative_abundance_intensity_pct, 200.0 / 3.0, places=3)
        self.assertAlmostEqual(man5_rec.relative_abundance_intensity_pct, 100.0 / 3.0, places=3)

        # Policy 2: ANNOTATED_GLYCAN_PEAKS_ONLY (denominator = 10000.0)
        q_annot = quantify_ms1_search_result(
            search_res,
            denominator_policy=DenominatorSetPolicy.ANNOTATED_GLYCAN_PEAKS_ONLY,
        )
        self.assertAlmostEqual(q_annot.denominator_max_intensity, 10000.0, places=2)
        g2f_annot = next(g for g in q_annot.glycan_abundances if g.glycan_composition == "Hex5HexNAc4dHex1")
        self.assertAlmostEqual(g2f_annot.base_peak_relative_intensity_pct, 100.0, places=4)

    def test_peak_area_vs_centroid_intensity_distinction(self) -> None:
        """Ensure peak area is reported when present (CSV example) and strictly None on centroid-only XY."""
        # 1. Centroid-only XY file -> has_peak_area_data is False, area is None
        xy_spec = parse_spectrum_file(EXAMPLES_DIR / "synthetic_native_nglycans.xy")
        xy_search = search_ms1_library(xy_spec, self.library)
        xy_quant = quantify_ms1_search_result(xy_search)
        self.assertFalse(xy_quant.has_peak_area_data)
        self.assertIsNone(xy_quant.total_annotated_area)
        self.assertTrue(all(g.aggregated_area is None for g in xy_quant.glycan_abundances))
        self.assertTrue(all(g.relative_abundance_area_pct is None for g in xy_quant.glycan_abundances))

        # 2. CSV file with integrated peak areas -> has_peak_area_data is True, area % sums to 100%
        csv_spec = parse_spectrum_file(EXAMPLES_DIR / "synthetic_native_nglycans.csv")
        csv_search = search_ms1_library(csv_spec, self.library)
        csv_quant = quantify_ms1_search_result(csv_search)
        self.assertTrue(csv_quant.has_peak_area_data)
        self.assertIsNotNone(csv_quant.total_annotated_area)
        area_pct_sum = sum(g.relative_abundance_area_pct or 0.0 for g in csv_quant.glycan_abundances)
        self.assertAlmostEqual(area_pct_sum, 100.0, places=4)

    def test_adduct_aggregation_policies_and_provenance(self) -> None:
        """Verify PRIMARY_ADDUCT_ONLY vs EXPLICIT_SUM_DISTINCT_ADDUCTS and provenance history."""
        # Spectrum containing Man5 as both [M+Na]+ (1257.4226, int=8000) and [M+K]+ (1273.3966, int=2000)
        # plus G0F as [M+Na]+ (1485.5336, int=10000)
        txt = "m/z intensity\n1257.4226 8000.0\n1273.3966 2000.0\n1485.5336 10000.0\n"
        spec = parse_txt_spectrum(txt, filename="multi_adduct.txt")
        search_res = search_ms1_library(
            spec,
            self.library,
            settings=MS1SearchSettings(
                tolerance=15.0,
                tolerance_unit=MassToleranceUnit.PPM,
                adduct_codes=("[M+Na]+", "[M+K]+"),
            ),
        )

        # Policy A: PRIMARY_ADDUCT_ONLY ([M+Na]+) -> does NOT sum [M+K]+! Man5 = 8000, G0F = 10000
        q_prim = quantify_ms1_search_result(
            search_res,
            aggregation_policy=GlycanAggregationPolicy.PRIMARY_ADDUCT_ONLY,
            primary_adduct="[M+Na]+",
        )
        man5_prim = next(g for g in q_prim.glycan_abundances if g.glycan_composition == "Hex5HexNAc2")
        self.assertAlmostEqual(man5_prim.aggregated_intensity, 8000.0, places=2)
        self.assertEqual(man5_prim.contributing_peaks_count, 1)
        self.assertEqual(man5_prim.contributing_adducts, ("[M+Na]+",))

        # Policy B: EXPLICIT_SUM_DISTINCT_ADDUCTS -> sums [M+Na]+ (8000) + [M+K]+ (2000) = 10000 with full provenance
        q_sum = quantify_ms1_search_result(
            search_res,
            aggregation_policy=GlycanAggregationPolicy.EXPLICIT_SUM_DISTINCT_ADDUCTS,
        )
        man5_sum = next(g for g in q_sum.glycan_abundances if g.glycan_composition == "Hex5HexNAc2")
        self.assertAlmostEqual(man5_sum.aggregated_intensity, 10000.0, places=2)
        self.assertEqual(man5_sum.contributing_peaks_count, 2)
        self.assertEqual(set(man5_sum.contributing_adducts), {"[M+Na]+", "[M+K]+"})
        self.assertEqual(len(man5_sum.peak_provenance_history), 2)
        self.assertAlmostEqual(man5_sum.relative_abundance_intensity_pct, 50.0, places=4)

    def test_replicates_and_missing_data_without_silent_imputation(self) -> None:
        """Verify replicate Mean, SD, CV, and missingness % without silent zero-imputation."""
        # Replicate 1: Man5 (5000), G0F (5000) -> 50% each
        s1 = parse_txt_spectrum("m/z intensity\n1257.4226 5000\n1485.5336 5000\n", sample_id="Rep1")
        # Replicate 2: Man5 (6000), G0F (4000) -> 60% and 40%
        s2 = parse_txt_spectrum("m/z intensity\n1257.4226 6000\n1485.5336 4000\n", sample_id="Rep2")
        # Replicate 3: Man5 (5500), G0F (4500), plus G2F (1809.6393, 10000) -> G2F present ONLY in Rep3!
        s3 = parse_txt_spectrum(
            "m/z intensity\n1257.4226 5500\n1485.5336 4500\n1809.6393 10000\n",
            sample_id="Rep3",
        )

        q1 = quantify_ms1_search_result(search_ms1_library(s1, self.library))
        q2 = quantify_ms1_search_result(search_ms1_library(s2, self.library))
        q3 = quantify_ms1_search_result(search_ms1_library(s3, self.library))

        matrix = summarize_replicates_and_build_matrix([q1, q2, q3])
        self.assertEqual(matrix.sample_ids, ("Rep1", "Rep2", "Rep3"))

        # Check G2F (Hex5HexNAc4dHex1): observed in 1/3 replicates, missing in 2/3 (66.67% missingness)
        g2f_stat = next(s for s in matrix.replicate_summaries if s.glycan_composition == "Hex5HexNAc4dHex1")
        self.assertEqual(g2f_stat.total_replicates, 3)
        self.assertEqual(g2f_stat.observed_replicates, 1)
        self.assertEqual(g2f_stat.missing_replicates, 2)
        self.assertAlmostEqual(g2f_stat.missingness_pct, 200.0 / 3.0, places=2)
        # Crucial check: Rep1 and Rep2 values in matrix are None (NOT imputed as 0.0!)
        self.assertIsNone(matrix.abundance_matrix["Rep1"]["Hex5HexNAc4dHex1"])
        self.assertIsNone(matrix.abundance_matrix["Rep2"]["Hex5HexNAc4dHex1"])
        self.assertAlmostEqual(matrix.abundance_matrix["Rep3"]["Hex5HexNAc4dHex1"] or 0.0, 50.0, places=4)
        # Mean over observed replicates is 50.0% (not 50/3 = 16.67%), and SD/CV are None for n_obs=1
        self.assertAlmostEqual(g2f_stat.mean_relative_abundance_pct or 0.0, 50.0, places=4)
        self.assertIsNone(g2f_stat.sd_relative_abundance_pct)
        self.assertIsNone(g2f_stat.cv_pct)

        # Check Man5 (Hex5HexNAc2): observed in all 3 replicates (50%, 60%, 27.5%)
        man5_stat = next(s for s in matrix.replicate_summaries if s.glycan_composition == "Hex5HexNAc2")
        self.assertEqual(man5_stat.observed_replicates, 3)
        self.assertEqual(man5_stat.missing_replicates, 0)
        self.assertAlmostEqual(man5_stat.missingness_pct, 0.0, places=4)
        self.assertAlmostEqual(man5_stat.mean_relative_abundance_pct or 0.0, (50.0 + 60.0 + 27.5) / 3.0, places=4)
        self.assertIsNotNone(man5_stat.sd_relative_abundance_pct)
        self.assertIsNotNone(man5_stat.cv_pct)


if __name__ == "__main__":
    unittest.main()

"""Unit tests for Stage 9 Statistics, Visualization, Derived Traits, and Leakage-Free ML (`app/analysis/statistics.py`)."""

from __future__ import annotations

from pathlib import Path
import unittest

from app.analysis.glycan_search import MS1SearchSettings, search_ms1_library
from app.analysis.quantification import (
    DenominatorSetPolicy,
    GlycanAbundanceRecord,
    GlycanAggregationPolicy,
    SampleQuantificationResult,
    quantify_ms1_search_result,
    summarize_replicates_and_build_matrix,
)
from app.analysis.statistics import (
    MissingValueImputationMethod,
    build_mixomics_diablo_export_tables,
    compute_correlation_matrix,
    compute_derived_glycan_traits,
    compute_descriptive_statistics,
    compute_univariate_roc,
    render_annotated_spectrum_png,
    render_pca_plot_png,
    run_group_comparisons,
    run_hierarchical_clustering,
    run_leakage_free_classification,
    run_pca_analysis,
)
from app.libraries.library_io import import_library_from_json
from app.parsers import parse_spectrum_file


def _make_synthetic_sample_quant(
    sample_id: str,
    abundances: dict[str, tuple[str, float]],
) -> SampleQuantificationResult:
    """Helper to construct a `SampleQuantificationResult` with specified relative abundances."""
    records = []
    for comp, (g_class, rel_pct) in abundances.items():
        records.append(
            GlycanAbundanceRecord(
                sample_id=sample_id,
                glycan_composition=comp,
                candidate_class=g_class,
                aggregation_policy_used=GlycanAggregationPolicy.PRIMARY_ADDUCT_ONLY,
                primary_adduct_filter="[M+Na]+",
                aggregated_intensity=rel_pct * 100.0,
                aggregated_area=None,
                base_peak_relative_intensity_pct=rel_pct,
                relative_abundance_intensity_pct=rel_pct,
                relative_abundance_area_pct=None,
                contributing_peaks_count=1,
                contributing_adducts=("[M+Na]+",),
                contributing_mzs=(1500.0,),
                peak_provenance_history=("peak_1",),
                contains_ambiguous_peak=False,
                warnings=(),
            )
        )
    return SampleQuantificationResult(
        sample_id=sample_id,
        denominator_policy=DenominatorSetPolicy.ALL_ELIGIBLE_SPECTRUM_PEAKS,
        aggregation_policy=GlycanAggregationPolicy.PRIMARY_ADDUCT_ONLY,
        primary_adduct="[M+Na]+",
        exclude_ambiguous_from_aggregation=False,
        denominator_max_intensity=5000.0,
        total_annotated_intensity=10000.0,
        total_annotated_area=None,
        has_peak_area_data=False,
        peak_records=(),
        glycan_abundances=tuple(records),
        warnings=(),
    )


class TestStatisticsAndML(unittest.TestCase):
    """Test suite for Stage 9 multivariate statistics, group tests, derived traits, plots, and ML."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.library = import_library_from_json(Path("data/libraries/human_nglycan_native_v1.json"))

        # Build a 6-sample synthetic cohort (3 Control, 3 Disease) with distinct glycan shifts
        cls.quants = [
            _make_synthetic_sample_quant(
                "Ctrl_1",
                {
                    "Hex5HexNAc2": ("high_mannose", 40.0),
                    "Hex3HexNAc4dHex1": ("complex", 30.0),
                    "Hex5HexNAc4Neu5Ac1": ("complex", 20.0),
                    "Hex5HexNAc4dHex1Neu5Ac2": ("complex", 10.0),
                },
            ),
            _make_synthetic_sample_quant(
                "Ctrl_2",
                {
                    "Hex5HexNAc2": ("high_mannose", 42.0),
                    "Hex3HexNAc4dHex1": ("complex", 28.0),
                    "Hex5HexNAc4Neu5Ac1": ("complex", 21.0),
                    "Hex5HexNAc4dHex1Neu5Ac2": ("complex", 9.0),
                },
            ),
            _make_synthetic_sample_quant(
                "Ctrl_3",
                {
                    "Hex5HexNAc2": ("high_mannose", 38.0),
                    "Hex3HexNAc4dHex1": ("complex", 32.0),
                    "Hex5HexNAc4Neu5Ac1": ("complex", 19.0),
                    # Missing Hex5HexNAc4dHex1Neu5Ac2 in Ctrl_3 to test missing-value handling
                },
            ),
            _make_synthetic_sample_quant(
                "Dis_1",
                {
                    "Hex5HexNAc2": ("high_mannose", 15.0),
                    "Hex3HexNAc4dHex1": ("complex", 15.0),
                    "Hex5HexNAc4Neu5Ac1": ("complex", 35.0),
                    "Hex5HexNAc4dHex1Neu5Ac2": ("complex", 35.0),
                },
            ),
            _make_synthetic_sample_quant(
                "Dis_2",
                {
                    "Hex5HexNAc2": ("high_mannose", 14.0),
                    "Hex3HexNAc4dHex1": ("complex", 16.0),
                    "Hex5HexNAc4Neu5Ac1": ("complex", 36.0),
                    "Hex5HexNAc4dHex1Neu5Ac2": ("complex", 34.0),
                },
            ),
            _make_synthetic_sample_quant(
                "Dis_3",
                {
                    "Hex5HexNAc2": ("high_mannose", 16.0),
                    "Hex3HexNAc4dHex1": ("complex", 14.0),
                    "Hex5HexNAc4Neu5Ac1": ("complex", 34.0),
                    "Hex5HexNAc4dHex1Neu5Ac2": ("complex", 36.0),
                },
            ),
        ]
        cls.matrix = summarize_replicates_and_build_matrix(cls.quants)
        cls.groups = {
            "Ctrl_1": "Control",
            "Ctrl_2": "Control",
            "Ctrl_3": "Control",
            "Dis_1": "Disease",
            "Dis_2": "Disease",
            "Dis_3": "Disease",
        }

    def test_descriptive_stats_pca_clustering_and_correlations(self) -> None:
        """Verify descriptive statistics, PCA, hierarchical clustering, and correlation matrices."""
        desc = compute_descriptive_statistics(self.matrix, self.groups)
        self.assertGreaterEqual(len(desc), 12)
        # Check missingness recorded for Hex5HexNAc4dHex1Neu5Ac2 in Control group
        ctrl_s2 = [
            d for d in desc
            if d.glycan_composition == "Hex5HexNAc4dHex1Neu5Ac2" and d.group_label == "Control"
        ][0]
        self.assertEqual(ctrl_s2.n_missing, 1)
        self.assertAlmostEqual(ctrl_s2.missingness_pct, 33.33, places=2)

        # PCA
        pca_res = run_pca_analysis(
            self.matrix,
            sample_groups=self.groups,
            imputation_method=MissingValueImputationMethod.HALF_MIN_OBSERVED,
        )
        self.assertEqual(pca_res.n_components, 2)
        self.assertGreater(pca_res.explained_variance_ratio[0], 0.80)
        self.assertTrue(any("Imputed 1 missing matrix entries" in w for w in pca_res.warnings))

        # Hierarchical clustering
        hc_res = run_hierarchical_clustering(self.matrix)
        self.assertEqual(len(hc_res.ordered_sample_ids), 6)
        self.assertEqual(len(hc_res.ordered_feature_names), 4)

        # Correlation analysis
        corr_res = compute_correlation_matrix(self.matrix, method="spearman")
        self.assertGreaterEqual(len(corr_res), 6)
        self.assertTrue(all(0.0 <= c.fdr_q_value <= 1.0 for c in corr_res))

    def test_group_comparisons_bh_fdr_and_univariate_roc(self) -> None:
        """Verify two-group Welch t-test / Mann-Whitney U, 3-group ANOVA / Kruskal-Wallis, and ROC."""
        # 1. Two-group comparison
        comp_2 = run_group_comparisons(self.matrix, self.groups)
        self.assertEqual(len(comp_2), 4)
        man5_row = [r for r in comp_2 if r.glycan_composition == "Hex5HexNAc2"][0]
        self.assertEqual(man5_row.comparison_type, "two_group")
        self.assertLess(man5_row.log2_fold_change, -1.0)
        self.assertLess(man5_row.parametric_fdr_q_value, 0.01)
        self.assertTrue(man5_row.significant_fdr_0_05)

        # 2. Multi-group (3 groups) comparison
        groups_3 = {
            "Ctrl_1": "G1",
            "Ctrl_2": "G1",
            "Ctrl_3": "G2",
            "Dis_1": "G2",
            "Dis_2": "G3",
            "Dis_3": "G3",
        }
        comp_3 = run_group_comparisons(self.matrix, groups_3)
        self.assertEqual(comp_3[0].comparison_type, "multi_group")
        self.assertEqual(comp_3[0].parametric_test_name, "One-way ANOVA")

        # 3. Univariate ROC
        roc_list = compute_univariate_roc(
            self.matrix,
            self.groups,
            positive_group="Disease",
            negative_group="Control",
        )
        self.assertEqual(len(roc_list), 4)
        self.assertEqual(roc_list[0].auc, 1.0)
        self.assertEqual(roc_list[0].youden_j_index, 1.0)

    def test_leakage_free_cross_validation_and_small_sample_guard(self) -> None:
        """Verify leakage-free scikit-learn pipeline and minimum class size safeguard."""
        cv_res = run_leakage_free_classification(
            self.matrix,
            self.groups,
            classifier_type="logistic_regression",
            n_splits=3,
            max_selected_features=3,
        )
        self.assertEqual(cv_res.oof_accuracy, 1.0)
        self.assertEqual(cv_res.oof_roc_auc, 1.0)
        self.assertIn("fitted strictly inside each training fold", cv_res.leakage_prevention_statement)
        self.assertTrue(any("Small cohort warning" in w for w in cv_res.warnings))

        # Verify guard when a class has < 3 samples
        tiny_groups = {"Ctrl_1": "A", "Ctrl_2": "A", "Dis_1": "B", "Dis_2": "B"}
        with self.assertRaises(ValueError) as ctx:
            run_leakage_free_classification(self.matrix, tiny_groups)
        self.assertIn("Insufficient samples per class", str(ctx.exception))

    def test_hypothesis_labeled_derived_traits_plots_and_mixomics_export(self) -> None:
        """Verify derived glycan traits, static PNG plot generation, and mixOmics/DIABLO export tables."""
        traits_ctrl1 = compute_derived_glycan_traits(self.quants[0])
        self.assertAlmostEqual(traits_ctrl1.trait_dict["TM_pct"], 40.0, places=2)
        self.assertAlmostEqual(traits_ctrl1.trait_dict["TC_pct"], 60.0, places=2)
        self.assertAlmostEqual(traits_ctrl1.trait_dict["TFuc_pct"], 40.0, places=2)
        self.assertAlmostEqual(traits_ctrl1.trait_dict["TSia_pct"], 30.0, places=2)
        self.assertAlmostEqual(traits_ctrl1.trait_dict["A2_hypothesis_pct"], 60.0, places=2)

        # Verify branching/class traits carry explicit composition-inferred hypothesis flag
        a2_trait = [t for t in traits_ctrl1.traits if t.trait_code == "A2_hypothesis_pct"][0]
        self.assertTrue(a2_trait.is_composition_inferred_hypothesis)

        # Verify PNG plot generators produce valid PNG byte streams
        spec = parse_spectrum_file(Path("data/examples/synthetic_native_nglycans.csv"), sample_id="Plot_Spec")
        ms1_res = search_ms1_library(spec, self.library, MS1SearchSettings(tolerance=25.0))
        spec_png = render_annotated_spectrum_png(spec, ms1_res)
        self.assertTrue(spec_png.startswith(b"\x89PNG"))

        pca_res = run_pca_analysis(self.matrix, sample_groups=self.groups)
        pca_png = render_pca_plot_png(pca_res)
        self.assertTrue(pca_png.startswith(b"\x89PNG"))

        # Verify R mixOmics / DIABLO export tables
        diablo_tables = build_mixomics_diablo_export_tables(
            self.matrix,
            self.groups,
            self.quants,
        )
        self.assertEqual(len(diablo_tables["X_glycans"]), 6)
        self.assertEqual(len(diablo_tables["X_derived_traits"]), 6)
        self.assertEqual(len(diablo_tables["Y_metadata"]), 6)


if __name__ == "__main__":
    unittest.main()

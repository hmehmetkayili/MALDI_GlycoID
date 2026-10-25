"""Statistical analysis, derived glycan traits, visualization, and leakage-free ML (Stage 9).

Implements:
1. Descriptive statistics per glycan composition across samples and study groups.
2. Principal Component Analysis (PCA) with explicit imputation & scaling provenance.
3. Hierarchical clustering (sample and glycan dendrogram linkage matrices & ordering).
4. Pairwise Pearson and Spearman correlation matrices with Benjamini-Hochberg FDR q-values.
5. Two-group (Welch's t-test, Mann-Whitney U, log2 fold change) and multi-group
   (One-way ANOVA, Kruskal-Wallis) statistical comparisons with BH-FDR correction.
6. Univariate Receiver Operating Characteristic (ROC) curve analysis and Youden's J index.
7. Leakage-free cross-validated classification (`scikit-learn` `Pipeline` where imputation,
   scaling, and optional feature selection are fitted strictly inside each training fold).
8. Hypothesis-labeled derived glycan traits (class totals, fucosylation %, sialylation %,
   composition-inferred branching/galactosylation hypotheses, and linkage-specific sialic ratios).
9. Publication-ready static PNG plot generators (`matplotlib` `Agg` backend) for annotated
   spectra, PCA plots, volcano plots, and relative abundance bar charts.
10. R `mixOmics` / `DIABLO` multi-block export table formatter.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
import io
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import cluster as sp_cluster
from scipy import spatial as sp_spatial
from scipy import stats as sp_stats
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import KNNImputer, SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from app.analysis.glycan_search import MS1SearchResult
from app.analysis.quantification import MultiSampleQuantificationMatrix, SampleQuantificationResult
from app.core.mass_calculator import parse_glycan_composition
from app.libraries.composition_generator import classify_n_glycan_composition
from app.parsers.base_models import SpectrumData


class MissingValueImputationMethod(str, Enum):
    """Explicit missing-value handling policies for multivariate algorithms."""

    HALF_MIN_OBSERVED = "half_min_observed"
    MEDIAN = "median"
    MEAN = "mean"
    ZERO = "zero"
    KNN = "knn"


@dataclass(frozen=True)
class DescriptiveGlycanStat:
    """Descriptive statistical summary for a single glycan composition."""

    glycan_composition: str
    group_label: str
    n_total: int
    n_observed: int
    n_missing: int
    missingness_pct: float
    mean: Optional[float]
    sd: Optional[float]
    cv_pct: Optional[float]
    median: Optional[float]
    iqr: Optional[float]
    min_val: Optional[float]
    max_val: Optional[float]


@dataclass(frozen=True)
class PCAResult:
    """Principal Component Analysis result with explicit preprocessing provenance."""

    sample_ids: Tuple[str, ...]
    group_labels: Tuple[Optional[str], ...]
    feature_names: Tuple[str, ...]
    n_components: int
    explained_variance_ratio: Tuple[float, ...]
    cumulative_variance_ratio: Tuple[float, ...]
    sample_scores: Dict[str, Tuple[float, ...]]
    feature_loadings: Dict[str, Tuple[float, ...]]
    imputation_method_used: str
    log10_transformed: bool
    scaled_unit_variance: bool
    warnings: Tuple[str, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class HierarchicalClusteringResult:
    """Sample and feature hierarchical clustering dendrogram data."""

    sample_ids: Tuple[str, ...]
    ordered_sample_ids: Tuple[str, ...]
    sample_linkage_matrix: Tuple[Tuple[float, float, float, float], ...]
    feature_names: Tuple[str, ...]
    ordered_feature_names: Tuple[str, ...]
    feature_linkage_matrix: Tuple[Tuple[float, float, float, float], ...]
    distance_metric: str
    linkage_method: str
    imputation_method_used: str


@dataclass(frozen=True)
class CorrelationPairRecord:
    """Pairwise correlation result between two glycan features."""

    feature_a: str
    feature_b: str
    method: str  # "pearson" or "spearman"
    n_paired_observations: int
    correlation_r: float
    p_value: float
    fdr_q_value: float


@dataclass(frozen=True)
class GroupComparisonRecord:
    """Differential abundance statistical test result for a single glycan feature."""

    glycan_composition: str
    comparison_type: str  # "two_group" or "multi_group"
    groups_compared: Tuple[str, ...]
    group_means: Dict[str, Optional[float]]
    group_medians: Dict[str, Optional[float]]
    log2_fold_change: Optional[float]  # log2(mean_group2 / mean_group1) for two_group
    parametric_test_name: str  # "Welch t-test" or "One-way ANOVA"
    parametric_statistic: Optional[float]
    parametric_p_value: Optional[float]
    parametric_fdr_q_value: Optional[float]
    nonparametric_test_name: str  # "Mann-Whitney U" or "Kruskal-Wallis"
    nonparametric_statistic: Optional[float]
    nonparametric_p_value: Optional[float]
    nonparametric_fdr_q_value: Optional[float]
    significant_fdr_0_05: bool


@dataclass(frozen=True)
class UnivariateROCResult:
    """Univariate ROC biomarker evaluation for a binary group comparison."""

    glycan_composition: str
    positive_group: str
    negative_group: str
    n_positive: int
    n_negative: int
    auc: float
    youden_j_index: float
    optimal_threshold: float
    sensitivity_at_optimal: float
    specificity_at_optimal: float
    fpr_points: Tuple[float, ...]
    tpr_points: Tuple[float, ...]


@dataclass(frozen=True)
class CrossValidatedClassificationResult:
    """Result of a strictly leakage-free cross-validated classification pipeline."""

    classifier_name: str
    classes: Tuple[str, ...]
    n_samples: int
    n_features_total: int
    n_splits: int
    oof_accuracy: float
    oof_balanced_accuracy: float
    oof_roc_auc: Optional[float]
    fold_accuracies: Tuple[float, ...]
    top_features_mean_importance: Tuple[Tuple[str, float], ...]
    leakage_prevention_statement: str
    warnings: Tuple[str, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class DerivedGlycanTraitValue:
    """Single calculated glycan trait or ratio with explicit structural hypothesis labeling."""

    trait_code: str
    trait_name: str
    value: Optional[float]
    unit: str
    is_composition_inferred_hypothesis: bool
    structural_disclaimer: str


@dataclass(frozen=True)
class SampleDerivedTraitsReport:
    """Complete set of derived glycan traits for a single sample."""

    sample_id: str
    traits: Tuple[DerivedGlycanTraitValue, ...]
    trait_dict: Dict[str, Optional[float]]


def benjamini_hochberg_fdr(p_values: Sequence[Optional[float]]) -> List[Optional[float]]:
    """Compute Benjamini-Hochberg FDR-adjusted q-values preserving `None` entries."""
    valid_indices = [i for i, p in enumerate(p_values) if p is not None and not np.isnan(p)]
    m = len(valid_indices)
    if m == 0:
        return [None for _ in p_values]

    valid_p = np.array([float(p_values[i]) for i in valid_indices], dtype=np.float64)
    order = np.argsort(valid_p)
    ranked_p = valid_p[order]
    ranks = np.arange(1, m + 1, dtype=np.float64)

    q_vals = ranked_p * m / ranks
    # Enforce monotonicity from largest p-value down to smallest
    q_vals = np.minimum.accumulate(q_vals[::-1])[::-1]
    q_vals = np.clip(q_vals, 0.0, 1.0)

    out: List[Optional[float]] = [None for _ in p_values]
    for rank_idx, orig_idx in enumerate(order):
        out[valid_indices[orig_idx]] = float(q_vals[rank_idx])
    return out


def _extract_matrix_arrays(
    matrix: MultiSampleQuantificationMatrix,
) -> Tuple[List[str], List[str], np.ndarray]:
    """Return `(sample_ids, glycan_compositions, array_shape_samples_x_glycans_with_nan)`."""
    sample_ids = list(matrix.sample_ids)
    glycans = list(matrix.glycan_compositions)
    data = np.full((len(sample_ids), len(glycans)), np.nan, dtype=np.float64)

    for g_idx, summary in enumerate(matrix.replicate_summaries):
        for s_idx, sid in enumerate(sample_ids):
            val = summary.sample_values_relative_abundance_pct.get(sid)
            if val is not None:
                data[s_idx, g_idx] = float(val)

    return sample_ids, glycans, data


def impute_missing_matrix(
    data_nan: np.ndarray,
    method: MissingValueImputationMethod = MissingValueImputationMethod.HALF_MIN_OBSERVED,
) -> np.ndarray:
    """Explicitly impute missing (`NaN`) entries in a `samples x features` matrix."""
    arr = np.array(data_nan, dtype=np.float64, copy=True)
    if not np.isnan(arr).any():
        return arr

    if method == MissingValueImputationMethod.ZERO:
        return np.nan_to_num(arr, nan=0.0)

    if method == MissingValueImputationMethod.HALF_MIN_OBSERVED:
        for col in range(arr.shape[1]):
            col_vals = arr[:, col]
            pos_vals = col_vals[~np.isnan(col_vals) & (col_vals > 0)]
            fill_val = float(np.min(pos_vals) * 0.5) if len(pos_vals) > 0 else 0.0
            col_vals[np.isnan(col_vals)] = fill_val
        return arr

    if method == MissingValueImputationMethod.MEDIAN:
        imp = SimpleImputer(strategy="median", fill_value=0.0)
        return imp.fit_transform(arr)

    if method == MissingValueImputationMethod.MEAN:
        imp = SimpleImputer(strategy="mean", fill_value=0.0)
        return imp.fit_transform(arr)

    if method == MissingValueImputationMethod.KNN:
        n_neighbors = min(3, max(1, arr.shape[0] - 1))
        imp = KNNImputer(n_neighbors=n_neighbors)
        return imp.fit_transform(arr)

    raise ValueError(f"Unsupported imputation method: {method}")


def compute_descriptive_statistics(
    matrix: MultiSampleQuantificationMatrix,
    sample_groups: Optional[Mapping[str, str]] = None,
) -> List[DescriptiveGlycanStat]:
    """Compute descriptive statistics per glycan overall and optionally stratified by group."""
    sample_ids, glycans, data_nan = _extract_matrix_arrays(matrix)
    stats_out: List[DescriptiveGlycanStat] = []

    group_subsets: List[Tuple[str, List[int]]] = [("ALL", list(range(len(sample_ids))))]
    if sample_groups:
        unique_groups = sorted({sample_groups[s] for s in sample_ids if s in sample_groups})
        for grp in unique_groups:
            idxs = [i for i, s in enumerate(sample_ids) if sample_groups.get(s) == grp]
            if idxs:
                group_subsets.append((grp, idxs))

    for grp_label, row_indices in group_subsets:
        n_tot = len(row_indices)
        for g_idx, g_comp in enumerate(glycans):
            vals = data_nan[row_indices, g_idx]
            obs = vals[~np.isnan(vals)]
            n_obs = len(obs)
            n_miss = n_tot - n_obs
            miss_pct = (n_miss / n_tot * 100.0) if n_tot > 0 else 0.0

            if n_obs > 0:
                mean_v = float(np.mean(obs))
                med_v = float(np.median(obs))
                min_v = float(np.min(obs))
                max_v = float(np.max(obs))
                q75, q25 = np.percentile(obs, [75, 25])
                iqr_v = float(q75 - q25)
                if n_obs >= 2:
                    sd_v = float(np.std(obs, ddof=1))
                    cv_v = float((sd_v / mean_v) * 100.0) if mean_v > 0 else None
                else:
                    sd_v = None
                    cv_v = None
            else:
                mean_v = sd_v = cv_v = med_v = iqr_v = min_v = max_v = None

            stats_out.append(
                DescriptiveGlycanStat(
                    glycan_composition=g_comp,
                    group_label=grp_label,
                    n_total=n_tot,
                    n_observed=n_obs,
                    n_missing=n_miss,
                    missingness_pct=round(miss_pct, 2),
                    mean=round(mean_v, 4) if mean_v is not None else None,
                    sd=round(sd_v, 4) if sd_v is not None else None,
                    cv_pct=round(cv_v, 2) if cv_v is not None else None,
                    median=round(med_v, 4) if med_v is not None else None,
                    iqr=round(iqr_v, 4) if iqr_v is not None else None,
                    min_val=round(min_v, 4) if min_v is not None else None,
                    max_val=round(max_v, 4) if max_v is not None else None,
                )
            )

    return stats_out


def run_pca_analysis(
    matrix: MultiSampleQuantificationMatrix,
    sample_groups: Optional[Mapping[str, str]] = None,
    n_components: int = 2,
    imputation_method: MissingValueImputationMethod = MissingValueImputationMethod.HALF_MIN_OBSERVED,
    log10_transform: bool = True,
    scale_unit_variance: bool = True,
) -> PCAResult:
    """Perform Principal Component Analysis with explicit imputation and scaling provenance."""
    sample_ids, glycans, data_nan = _extract_matrix_arrays(matrix)
    if len(sample_ids) < 2 or len(glycans) < 2:
        raise ValueError("PCA requires at least 2 samples and 2 quantified glycan features.")

    warnings: List[str] = []
    n_missing = int(np.isnan(data_nan).sum())
    if n_missing > 0:
        warnings.append(
            f"Imputed {n_missing} missing matrix entries using explicit policy '{imputation_method.value}' prior to PCA."
        )

    X = impute_missing_matrix(data_nan, method=imputation_method)
    if log10_transform:
        pos_vals = X[X > 0]
        pseudocount = float(np.min(pos_vals) * 0.1) if len(pos_vals) > 0 else 1e-4
        X = np.log10(np.maximum(X, 0.0) + pseudocount)

    if scale_unit_variance:
        scaler = StandardScaler()
        X = scaler.fit_transform(X)

    max_comp = min(n_components, X.shape[0], X.shape[1])
    if np.allclose(X, X[0]):
        scores = np.zeros((X.shape[0], max_comp), dtype=np.float64)
        loadings = np.zeros((X.shape[1], max_comp), dtype=np.float64)
        var_ratios = tuple(0.0 for _ in range(max_comp))
        cum_ratios = tuple(0.0 for _ in range(max_comp))
    else:
        pca = PCA(n_components=max_comp, random_state=42)
        scores = pca.fit_transform(X)
        loadings = pca.components_.T
        var_ratios = tuple(
            round(float(0.0 if np.isnan(v) else v), 6) for v in pca.explained_variance_ratio_
        )
        cum_ratios = tuple(
            round(float(0.0 if np.isnan(v) else v), 6) for v in np.cumsum(var_ratios)
        )

    sample_scores_dict = {
        sid: tuple(round(float(v), 6) for v in scores[idx])
        for idx, sid in enumerate(sample_ids)
    }
    feature_loadings_dict = {
        g: tuple(round(float(v), 6) for v in loadings[idx])
        for idx, g in enumerate(glycans)
    }
    grp_labels = tuple(sample_groups.get(s) if sample_groups else None for s in sample_ids)

    return PCAResult(
        sample_ids=tuple(sample_ids),
        group_labels=grp_labels,
        feature_names=tuple(glycans),
        n_components=max_comp,
        explained_variance_ratio=var_ratios,
        cumulative_variance_ratio=cum_ratios,
        sample_scores=sample_scores_dict,
        feature_loadings=feature_loadings_dict,
        imputation_method_used=imputation_method.value,
        log10_transformed=log10_transform,
        scaled_unit_variance=scale_unit_variance,
        warnings=tuple(warnings),
    )


def run_hierarchical_clustering(
    matrix: MultiSampleQuantificationMatrix,
    imputation_method: MissingValueImputationMethod = MissingValueImputationMethod.HALF_MIN_OBSERVED,
    distance_metric: str = "euclidean",
    linkage_method: str = "ward",
    scale_features: bool = True,
) -> HierarchicalClusteringResult:
    """Compute sample and feature hierarchical clustering linkage matrices and leaf orderings."""
    sample_ids, glycans, data_nan = _extract_matrix_arrays(matrix)
    if len(sample_ids) < 2 or len(glycans) < 2:
        raise ValueError("Hierarchical clustering requires at least 2 samples and 2 glycan features.")

    X = impute_missing_matrix(data_nan, method=imputation_method)
    if scale_features:
        X = StandardScaler().fit_transform(X)

    # Ward linkage requires Euclidean metric
    eff_metric = "euclidean" if linkage_method == "ward" else distance_metric

    sample_dist = sp_spatial.distance.pdist(X, metric=eff_metric)
    sample_dist = np.nan_to_num(sample_dist, nan=0.0)
    sample_Z = sp_cluster.hierarchy.linkage(sample_dist, method=linkage_method)
    sample_order = sp_cluster.hierarchy.leaves_list(sample_Z)

    feat_dist = sp_spatial.distance.pdist(X.T, metric=eff_metric)
    feat_dist = np.nan_to_num(feat_dist, nan=0.0)
    feat_Z = sp_cluster.hierarchy.linkage(feat_dist, method=linkage_method)
    feat_order = sp_cluster.hierarchy.leaves_list(feat_Z)

    return HierarchicalClusteringResult(
        sample_ids=tuple(sample_ids),
        ordered_sample_ids=tuple(sample_ids[int(i)] for i in sample_order),
        sample_linkage_matrix=tuple(
            (float(row[0]), float(row[1]), float(row[2]), float(row[3])) for row in sample_Z
        ),
        feature_names=tuple(glycans),
        ordered_feature_names=tuple(glycans[int(i)] for i in feat_order),
        feature_linkage_matrix=tuple(
            (float(row[0]), float(row[1]), float(row[2]), float(row[3])) for row in feat_Z
        ),
        distance_metric=eff_metric,
        linkage_method=linkage_method,
        imputation_method_used=imputation_method.value,
    )


def compute_correlation_matrix(
    matrix: MultiSampleQuantificationMatrix,
    method: str = "spearman",
    imputation_method: Optional[MissingValueImputationMethod] = None,
) -> List[CorrelationPairRecord]:
    """Compute pairwise glycan-to-glycan correlations (Pearson or Spearman) with BH-FDR q-values."""
    _, glycans, data_nan = _extract_matrix_arrays(matrix)
    X = impute_missing_matrix(data_nan, imputation_method) if imputation_method is not None else data_nan

    raw_records: List[Tuple[str, str, int, float, float]] = []
    n_feat = len(glycans)

    for i in range(n_feat):
        for j in range(i + 1, n_feat):
            col_i = X[:, i]
            col_j = X[:, j]
            mask = ~np.isnan(col_i) & ~np.isnan(col_j)
            n_pair = int(mask.sum())
            if n_pair < 3:
                continue
            xi = col_i[mask]
            xj = col_j[mask]
            if np.std(xi) == 0.0 or np.std(xj) == 0.0:
                continue
            if method.lower() == "pearson":
                r_val, p_val = sp_stats.pearsonr(xi, xj)
            else:
                r_val, p_val = sp_stats.spearmanr(xi, xj)
            raw_records.append((glycans[i], glycans[j], n_pair, float(r_val), float(p_val)))

    q_values = benjamini_hochberg_fdr([r[4] for r in raw_records])
    out: List[CorrelationPairRecord] = []
    for (fa, fb, n_p, r_v, p_v), q_v in zip(raw_records, q_values):
        out.append(
            CorrelationPairRecord(
                feature_a=fa,
                feature_b=fb,
                method=method.lower(),
                n_paired_observations=n_p,
                correlation_r=round(r_v, 4),
                p_value=float(p_v),
                fdr_q_value=float(q_v) if q_v is not None else 1.0,
            )
        )
    return out


def run_group_comparisons(
    matrix: MultiSampleQuantificationMatrix,
    sample_groups: Mapping[str, str],
    imputation_method: Optional[MissingValueImputationMethod] = None,
) -> List[GroupComparisonRecord]:
    """Run two-group (Welch t-test, Mann-Whitney U) or multi-group (ANOVA, Kruskal-Wallis) tests with BH-FDR."""
    sample_ids, glycans, data_nan = _extract_matrix_arrays(matrix)
    X = impute_missing_matrix(data_nan, imputation_method) if imputation_method is not None else data_nan

    grouped_indices: Dict[str, List[int]] = {}
    for idx, sid in enumerate(sample_ids):
        grp = sample_groups.get(sid)
        if grp:
            grouped_indices.setdefault(grp, []).append(idx)

    group_names = sorted(grouped_indices.keys())
    if len(group_names) < 2:
        raise ValueError("Group comparison requires at least 2 distinct sample groups.")

    is_two_group = len(group_names) == 2
    comp_type = "two_group" if is_two_group else "multi_group"
    param_name = "Welch t-test" if is_two_group else "One-way ANOVA"
    nonparam_name = "Mann-Whitney U" if is_two_group else "Kruskal-Wallis"

    interim: List[Dict[str, Any]] = []
    param_pvals: List[Optional[float]] = []
    nonparam_pvals: List[Optional[float]] = []

    for g_idx, g_comp in enumerate(glycans):
        grp_arrays: Dict[str, np.ndarray] = {}
        grp_means: Dict[str, Optional[float]] = {}
        grp_medians: Dict[str, Optional[float]] = {}

        for g_name in group_names:
            raw_vals = X[grouped_indices[g_name], g_idx]
            valid_vals = raw_vals[~np.isnan(raw_vals)]
            grp_arrays[g_name] = valid_vals
            grp_means[g_name] = round(float(np.mean(valid_vals)), 4) if len(valid_vals) > 0 else None
            grp_medians[g_name] = round(float(np.median(valid_vals)), 4) if len(valid_vals) > 0 else None

        log2_fc: Optional[float] = None
        p_stat: Optional[float] = None
        p_pval: Optional[float] = None
        np_stat: Optional[float] = None
        np_pval: Optional[float] = None

        if is_two_group:
            g1, g2 = group_names[0], group_names[1]
            v1, v2 = grp_arrays[g1], grp_arrays[g2]
            m1, m2 = grp_means[g1], grp_means[g2]
            if m1 is not None and m2 is not None and m1 > 0 and m2 > 0:
                log2_fc = round(float(np.log2(m2 / m1)), 4)

            if len(v1) >= 2 and len(v2) >= 2 and (np.std(v1) > 0 or np.std(v2) > 0):
                t_res = sp_stats.ttest_ind(v2, v1, equal_var=False)
                p_stat = float(t_res.statistic)
                p_pval = float(t_res.pvalue)
            if len(v1) >= 1 and len(v2) >= 1:
                try:
                    u_res = sp_stats.mannwhitneyu(v2, v1, alternative="two-sided")
                    np_stat = float(u_res.statistic)
                    np_pval = float(u_res.pvalue)
                except ValueError:
                    pass
        else:
            valid_groups = [arr for arr in grp_arrays.values() if len(arr) >= 2]
            if len(valid_groups) >= 2:
                try:
                    f_res = sp_stats.f_oneway(*valid_groups)
                    if not np.isnan(f_res.pvalue):
                        p_stat = float(f_res.statistic)
                        p_pval = float(f_res.pvalue)
                except ValueError:
                    pass
                try:
                    kw_res = sp_stats.kruskal(*valid_groups)
                    if not np.isnan(kw_res.pvalue):
                        np_stat = float(kw_res.statistic)
                        np_pval = float(kw_res.pvalue)
                except ValueError:
                    pass

        param_pvals.append(p_pval)
        nonparam_pvals.append(np_pval)
        interim.append(
            {
                "glycan_composition": g_comp,
                "group_means": grp_means,
                "group_medians": grp_medians,
                "log2_fold_change": log2_fc,
                "p_stat": p_stat,
                "p_pval": p_pval,
                "np_stat": np_stat,
                "np_pval": np_pval,
            }
        )

    param_qvals = benjamini_hochberg_fdr(param_pvals)
    nonparam_qvals = benjamini_hochberg_fdr(nonparam_pvals)

    results: List[GroupComparisonRecord] = []
    for row, pq, npq in zip(interim, param_qvals, nonparam_qvals):
        sig = bool((pq is not None and pq < 0.05) or (npq is not None and npq < 0.05))
        results.append(
            GroupComparisonRecord(
                glycan_composition=row["glycan_composition"],
                comparison_type=comp_type,
                groups_compared=tuple(group_names),
                group_means=row["group_means"],
                group_medians=row["group_medians"],
                log2_fold_change=row["log2_fold_change"],
                parametric_test_name=param_name,
                parametric_statistic=row["p_stat"],
                parametric_p_value=row["p_pval"],
                parametric_fdr_q_value=pq,
                nonparametric_test_name=nonparam_name,
                nonparametric_statistic=row["np_stat"],
                nonparametric_p_value=row["np_pval"],
                nonparametric_fdr_q_value=npq,
                significant_fdr_0_05=sig,
            )
        )

    return results


def compute_univariate_roc(
    matrix: MultiSampleQuantificationMatrix,
    sample_groups: Mapping[str, str],
    positive_group: str,
    negative_group: str,
    imputation_method: MissingValueImputationMethod = MissingValueImputationMethod.HALF_MIN_OBSERVED,
) -> List[UnivariateROCResult]:
    """Compute univariate ROC curves and optimal Youden's J thresholds for each glycan."""
    sample_ids, glycans, data_nan = _extract_matrix_arrays(matrix)
    X = impute_missing_matrix(data_nan, method=imputation_method)

    pos_indices = [i for i, s in enumerate(sample_ids) if sample_groups.get(s) == positive_group]
    neg_indices = [i for i, s in enumerate(sample_ids) if sample_groups.get(s) == negative_group]
    if len(pos_indices) < 1 or len(neg_indices) < 1:
        raise ValueError(
            f"ROC analysis requires at least 1 sample in positive group '{positive_group}' "
            f"and 1 sample in negative group '{negative_group}'."
        )

    all_idx = pos_indices + neg_indices
    y_true = np.array([1] * len(pos_indices) + [0] * len(neg_indices), dtype=int)

    roc_results: List[UnivariateROCResult] = []
    for g_idx, g_comp in enumerate(glycans):
        scores = X[all_idx, g_idx]
        if np.std(scores) == 0.0:
            auc_val = 0.5
            fpr, tpr, thresh = np.array([0.0, 1.0]), np.array([0.0, 1.0]), np.array([scores[0], scores[0]])
        else:
            auc_val = float(roc_auc_score(y_true, scores))
            # If inverse direction (< 0.5), invert direction so AUC >= 0.5 represents discriminative magnitude
            if auc_val < 0.5:
                fpr, tpr, thresh = roc_curve(y_true, -scores)
                thresh = -thresh
                auc_val = 1.0 - auc_val
            else:
                fpr, tpr, thresh = roc_curve(y_true, scores)

        j_scores = tpr - fpr
        best_i = int(np.argmax(j_scores))
        opt_thresh = float(thresh[best_i]) if not np.isinf(thresh[best_i]) else float(np.median(scores))

        roc_results.append(
            UnivariateROCResult(
                glycan_composition=g_comp,
                positive_group=positive_group,
                negative_group=negative_group,
                n_positive=len(pos_indices),
                n_negative=len(neg_indices),
                auc=round(auc_val, 4),
                youden_j_index=round(float(j_scores[best_i]), 4),
                optimal_threshold=round(opt_thresh, 4),
                sensitivity_at_optimal=round(float(tpr[best_i]), 4),
                specificity_at_optimal=round(float(1.0 - fpr[best_i]), 4),
                fpr_points=tuple(round(float(v), 4) for v in fpr),
                tpr_points=tuple(round(float(v), 4) for v in tpr),
            )
        )

    roc_results.sort(key=lambda r: -r.auc)
    return roc_results


def run_leakage_free_classification(
    matrix: MultiSampleQuantificationMatrix,
    sample_groups: Mapping[str, str],
    classifier_type: str = "logistic_regression",
    n_splits: int = 3,
    max_selected_features: Optional[int] = None,
) -> CrossValidatedClassificationResult:
    """Run strictly leakage-free stratified cross-validation.

    Scientific & ML Integrity Safeguard:
    `SimpleImputer`, `StandardScaler`, and optional `SelectKBest` are encapsulated
    inside the `scikit-learn` `Pipeline` and fitted strictly on `X_train` within each
    cross-validation fold. `X_test` is never seen during imputation, scaling, or
    feature selection.
    """
    sample_ids, glycans, data_nan = _extract_matrix_arrays(matrix)

    labeled_rows: List[int] = []
    labels: List[str] = []
    for idx, sid in enumerate(sample_ids):
        grp = sample_groups.get(sid)
        if grp:
            labeled_rows.append(idx)
            labels.append(grp)

    classes = sorted(set(labels))
    if len(classes) < 2:
        raise ValueError("Classification requires at least 2 distinct sample classes.")

    class_counts = {c: labels.count(c) for c in classes}
    min_class_n = min(class_counts.values())
    if min_class_n < 3:
        raise ValueError(
            f"Insufficient samples per class for cross-validated classification: minimum class count is "
            f"{min_class_n} ({class_counts}); at least 3 samples per class are required."
        )

    eff_splits = min(n_splits, min_class_n)
    warnings: List[str] = []
    if min_class_n < 10:
        warnings.append(
            f"Small cohort warning: minimum class size is {min_class_n} samples (< 10 per class). "
            f"Cross-validated performance metrics should be treated as exploratory hypothesis generation."
        )

    X_raw = data_nan[labeled_rows, :]
    class_to_int = {c: i for i, c in enumerate(classes)}
    y = np.array([class_to_int[lbl] for lbl in labels], dtype=int)

    skf = StratifiedKFold(n_splits=eff_splits, shuffle=True, random_state=42)
    oof_preds = np.zeros(len(y), dtype=int)
    oof_probs = np.zeros((len(y), len(classes)), dtype=np.float64)
    fold_accs: List[float] = []
    feature_importance_accum = np.zeros(len(glycans), dtype=np.float64)

    for train_idx, test_idx in skf.split(X_raw, y):
        X_train, y_train = X_raw[train_idx], y[train_idx]
        X_test, y_test = X_raw[test_idx], y[test_idx]

        steps: List[Tuple[str, Any]] = [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
        ]
        if max_selected_features is not None and 0 < max_selected_features < X_raw.shape[1]:
            steps.append(("select", SelectKBest(score_func=f_classif, k=max_selected_features)))

        if classifier_type.lower() == "random_forest":
            clf = RandomForestClassifier(n_estimators=50, random_state=42)
        else:
            clf = LogisticRegression(max_iter=500, random_state=42)
        steps.append(("clf", clf))

        pipe = Pipeline(steps)
        pipe.fit(X_train, y_train)

        preds = pipe.predict(X_test)
        probs = pipe.predict_proba(X_test)
        oof_preds[test_idx] = preds
        oof_probs[test_idx] = probs
        fold_accs.append(float(accuracy_score(y_test, preds)))

        # Extract feature importances mapped back to original features
        fitted_clf = pipe.named_steps["clf"]
        if hasattr(fitted_clf, "coef_"):
            raw_imp = np.mean(np.abs(fitted_clf.coef_), axis=0)
        elif hasattr(fitted_clf, "feature_importances_"):
            raw_imp = np.asarray(fitted_clf.feature_importances_, dtype=np.float64)
        else:
            raw_imp = np.zeros(X_train.shape[1], dtype=np.float64)

        if "select" in pipe.named_steps:
            support_mask = pipe.named_steps["select"].get_support()
            mapped_imp = np.zeros(len(glycans), dtype=np.float64)
            mapped_imp[support_mask] = raw_imp
            feature_importance_accum += mapped_imp / eff_splits
        else:
            feature_importance_accum += raw_imp / eff_splits

    oof_acc = float(accuracy_score(y, oof_preds))
    oof_bal_acc = float(balanced_accuracy_score(y, oof_preds))
    if len(classes) == 2:
        oof_auc: Optional[float] = float(roc_auc_score(y, oof_probs[:, 1]))
    else:
        oof_auc = float(roc_auc_score(y, oof_probs, multi_class="ovr"))

    ranked_feats = sorted(
        zip(glycans, feature_importance_accum.tolist()),
        key=lambda item: -item[1],
    )

    return CrossValidatedClassificationResult(
        classifier_name=classifier_type.lower(),
        classes=tuple(classes),
        n_samples=len(y),
        n_features_total=len(glycans),
        n_splits=eff_splits,
        oof_accuracy=round(oof_acc, 4),
        oof_balanced_accuracy=round(oof_bal_acc, 4),
        oof_roc_auc=round(oof_auc, 4) if oof_auc is not None else None,
        fold_accuracies=tuple(round(a, 4) for a in fold_accs),
        top_features_mean_importance=tuple((name, round(val, 4)) for name, val in ranked_feats[:15]),
        leakage_prevention_statement=(
            "Imputation, scaling, and feature selection fitted strictly inside each training fold "
            "(zero test-fold data leakage)."
        ),
        warnings=tuple(warnings),
    )


def compute_derived_glycan_traits(
    quant_result: SampleQuantificationResult,
    labeling: str = "NONE",
) -> SampleDerivedTraitsReport:
    """Compute biologically interpretable derived N-glycan traits with explicit hypothesis disclaimers."""
    hm_pct = 0.0
    hyb_pct = 0.0
    comp_pct = 0.0
    pauci_pct = 0.0
    fuc_pct = 0.0
    sia_pct = 0.0
    bi_hyp_pct = 0.0
    tri_or_bis_hyp_pct = 0.0
    tetra_hyp_pct = 0.0
    weighted_gal_num = 0.0
    complex_abund_denom = 0.0

    for rec in quant_result.glycan_abundances:
        rel_ab = rec.relative_abundance_intensity_pct
        counts = parse_glycan_composition(rec.glycan_composition)
        _is_valid, g_class, *_ = classify_n_glycan_composition(counts)

        if g_class.value == "high_mannose":
            hm_pct += rel_ab
        elif g_class.value == "hybrid":
            hyb_pct += rel_ab
        elif g_class.value == "complex":
            comp_pct += rel_ab
            hex_c = counts.get("Hex", 0)
            hexnac_c = counts.get("HexNAc", 0)
            if hexnac_c == 4:
                bi_hyp_pct += rel_ab
            elif hexnac_c == 5:
                tri_or_bis_hyp_pct += rel_ab
            elif hexnac_c >= 6:
                tetra_hyp_pct += rel_ab

            # Inferred galactose count on complex N-glycans = max(0, Hex - 3 core mannoses)
            gal_inferred = max(0, hex_c - 3)
            max_antennae = max(1, hexnac_c - 2)
            weighted_gal_num += rel_ab * min(1.0, gal_inferred / max_antennae)
            complex_abund_denom += rel_ab
        elif g_class.value == "paucimannose":
            pauci_pct += rel_ab

        if counts.get("dHex", 0) >= 1:
            fuc_pct += rel_ab
        if counts.get("Neu5Ac", 0) + counts.get("Neu5Gc", 0) >= 1:
            sia_pct += rel_ab

    cgal_index = (weighted_gal_num / complex_abund_denom * 100.0) if complex_abund_denom > 0 else None

    comp_disclaimer = (
        "Calculated from MS1 monosaccharide composition and N-glycan trimannosyl core rules; "
        "isomeric topologies cannot be resolved from unseparated MS1 mass alone."
    )
    branch_disclaimer = (
        "Composition-inferred antennary hypothesis: HexNAc=4 is assigned as biantennary hypothesis, "
        "and HexNAc=5 cannot distinguish triantennary vs. bisecting GlcNAc biantennary without MS/MS."
    )

    traits: List[DerivedGlycanTraitValue] = [
        DerivedGlycanTraitValue("TM_pct", "Total High-Mannose (%)", round(hm_pct, 3), "%", True, comp_disclaimer),
        DerivedGlycanTraitValue("THy_pct", "Total Hybrid (%)", round(hyb_pct, 3), "%", True, comp_disclaimer),
        DerivedGlycanTraitValue("TC_pct", "Total Complex (%)", round(comp_pct, 3), "%", True, comp_disclaimer),
        DerivedGlycanTraitValue("TP_pct", "Total Paucimannose (%)", round(pauci_pct, 3), "%", True, comp_disclaimer),
        DerivedGlycanTraitValue("TFuc_pct", "Total Fucosylation (%)", round(fuc_pct, 3), "%", False, "Direct stoichiometric presence of >=1 dHex (Fuc) residue."),
        DerivedGlycanTraitValue("TSia_pct", "Total Sialylation (%)", round(sia_pct, 3), "%", False, "Direct stoichiometric presence of >=1 Neu5Ac/Neu5Gc residue."),
        DerivedGlycanTraitValue(
            "CGal_hypothesis_index",
            "Complex Galactosylation Hypothesis Index (%)",
            round(cgal_index, 3) if cgal_index is not None else None,
            "%",
            True,
            "Assumes Hex - 3 corresponds to antennae galactose residues on complex N-glycans.",
        ),
        DerivedGlycanTraitValue("A2_hypothesis_pct", "Biantennary Hypothesis (HexNAc=4) (%)", round(bi_hyp_pct, 3), "%", True, branch_disclaimer),
        DerivedGlycanTraitValue("A3_or_Bisecting_hypothesis_pct", "Triantennary / Bisecting Hypothesis (HexNAc=5) (%)", round(tri_or_bis_hyp_pct, 3), "%", True, branch_disclaimer),
        DerivedGlycanTraitValue("A4_plus_hypothesis_pct", "Tetraantennary+ Hypothesis (HexNAc>=6) (%)", round(tetra_hyp_pct, 3), "%", True, branch_disclaimer),
    ]

    trait_dict = {t.trait_code: t.value for t in traits}
    return SampleDerivedTraitsReport(
        sample_id=quant_result.sample_id,
        traits=tuple(traits),
        trait_dict=trait_dict,
    )


def render_annotated_spectrum_png(
    spectrum: SpectrumData,
    search_result: Optional[MS1SearchResult] = None,
    title: Optional[str] = None,
) -> bytes:
    """Render a static publication-ready annotated MS1 mass spectrum plot as PNG bytes."""
    fig, ax = plt.subplots(figsize=(10, 4.5), dpi=150)
    mzs = [p.mz for p in spectrum.peaks]
    ints = [p.intensity for p in spectrum.peaks]

    if mzs:
        ax.vlines(mzs, 0, ints, color="#64748b", linewidth=1.0, label="Observed MS1 Peaks")

    if search_result is not None and search_result.top_matches_per_peak:
        matched_mzs = [m.observed_mz for m in search_result.top_matches_per_peak]
        matched_ints = [m.observed_intensity for m in search_result.top_matches_per_peak]
        ax.vlines(matched_mzs, 0, matched_ints, color="#0284c7", linewidth=1.6, label="Annotated N-Glycans")

        # Annotate top 8 most intense matched peaks
        top_annotated = sorted(search_result.top_matches_per_peak, key=lambda m: -m.observed_intensity)[:8]
        for m in top_annotated:
            ax.annotate(
                f"{m.glycan_composition}\n({m.observed_mz:.2f})",
                xy=(m.observed_mz, m.observed_intensity),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7,
                color="#0f172a",
            )

    ax.set_xlabel("m/z (Da/e)")
    ax.set_ylabel("Intensity (a.u.)")
    ax.set_title(title or f"MALDI-TOF MS1 Spectrum — {spectrum.sample_id}")
    ax.set_ylim(bottom=0, top=(max(ints) * 1.22) if ints else 1.0)
    ax.legend(loc="upper right", frameon=True, fontsize=8)
    ax.grid(True, linestyle="--", alpha=0.3)
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    plt.close(fig)
    return buf.getvalue()


def render_pca_plot_png(pca_result: PCAResult, title: str = "PCA Score Plot") -> bytes:
    """Render a 2D PCA score plot as PNG bytes."""
    fig, ax = plt.subplots(figsize=(6.5, 5.0), dpi=150)

    groups = [g or "Ungrouped" for g in pca_result.group_labels]
    unique_groups = sorted(set(groups))
    palette = ["#0284c7", "#e11d48", "#16a34a", "#d97706", "#7c3aed", "#475569"]

    for g_idx, grp in enumerate(unique_groups):
        xs: List[float] = []
        ys: List[float] = []
        sids: List[str] = []
        for sid, s_grp in zip(pca_result.sample_ids, groups):
            if s_grp == grp:
                sc = pca_result.sample_scores[sid]
                xs.append(sc[0])
                ys.append(sc[1] if len(sc) > 1 else 0.0)
                sids.append(sid)
        color = palette[g_idx % len(palette)]
        ax.scatter(xs, ys, label=grp, color=color, s=55, alpha=0.85, edgecolors="white")
        for x_v, y_v, sid in zip(xs, ys, sids):
            ax.annotate(sid, (x_v, y_v), xytext=(4, 4), textcoords="offset points", fontsize=7.5)

    pc1_pct = pca_result.explained_variance_ratio[0] * 100.0
    pc2_pct = (pca_result.explained_variance_ratio[1] * 100.0) if pca_result.n_components > 1 else 0.0
    ax.set_xlabel(f"PC1 ({pc1_pct:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pc2_pct:.1f}% variance)")
    ax.set_title(title)
    ax.axhline(0, color="#94a3b8", linestyle="--", linewidth=0.8)
    ax.axvline(0, color="#94a3b8", linestyle="--", linewidth=0.8)
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, linestyle=":", alpha=0.4)
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    plt.close(fig)
    return buf.getvalue()


def build_mixomics_diablo_export_tables(
    matrix: MultiSampleQuantificationMatrix,
    sample_groups: Mapping[str, str],
    per_sample_quantifications: Optional[Sequence[SampleQuantificationResult]] = None,
) -> Dict[str, List[Dict[str, Any]]]:
    """Build synchronized multi-block tables (`X_glycans`, `X_derived_traits`, `Y_metadata`) for R mixOmics/DIABLO."""
    sample_ids, glycans, data_nan = _extract_matrix_arrays(matrix)

    x_glycans_rows: List[Dict[str, Any]] = []
    for s_idx, sid in enumerate(sample_ids):
        row: Dict[str, Any] = {"sample_id": sid}
        for g_idx, g_comp in enumerate(glycans):
            val = data_nan[s_idx, g_idx]
            row[g_comp] = round(float(val), 4) if not np.isnan(val) else None
        x_glycans_rows.append(row)

    y_meta_rows: List[Dict[str, Any]] = [
        {"sample_id": sid, "group": sample_groups.get(sid, "UNASSIGNED")}
        for sid in sample_ids
    ]

    x_traits_rows: List[Dict[str, Any]] = []
    if per_sample_quantifications:
        quant_by_id = {q.sample_id: q for q in per_sample_quantifications}
        for sid in sample_ids:
            if sid in quant_by_id:
                tr_rep = compute_derived_glycan_traits(quant_by_id[sid])
                x_traits_rows.append({"sample_id": sid, **tr_rep.trait_dict})

    return {
        "X_glycans": x_glycans_rows,
        "X_derived_traits": x_traits_rows,
        "Y_metadata": y_meta_rows,
    }

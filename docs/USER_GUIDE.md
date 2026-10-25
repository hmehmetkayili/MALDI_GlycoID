# User Guide — MALDI-GlycoID

MALDI-GlycoID provides a 12-page web interface (`http://127.0.0.1:8000`) covering the complete analytical workflow:

1. **Dashboard:** Overview of active sample, active N-glycan library, matched peaks, QC status, and one-click synthetic validation workflow launch.
2. **Import Data:** Upload `.txt`, `.xy`, `.csv`, `.tsv`, `.mzML`, `.mzXML`, or Bruker 1D `.fid`+`acqus` spectra. Explicit controls for delimiter, decimal separator, column mapping, duplicate `m/z` handling, and `centroid` vs. `profile` data type classification.
3. **Glycan Library:** Browse and filter the 333-entry default human N-glycan library (`data/libraries/human_nglycan_native_v1.json`), generate custom biosynthetically constrained libraries (`Hex`, `HexNAc`, `dHex`, `Neu5Ac`, `Neu5Gc`, reducing-end labels, permethylation, linkage-specific sialic acid derivatization), and import/export JSON/CSV/TSV libraries.
4. **MS1 Glycan Search:** Run **Mode A (Targeted)** or **Mode B (Library)** MS1 search with configurable ppm/Da tolerance, polarity, adducts (`[M+Na]+`, `[M+H]+`, `[M+K]+`, `[M+2H]2+`, etc.), charge states, and candidate ambiguity inspection.
5. **Glycopeptide Search:** Paste FASTA protein sequences, choose protease (`Trypsin`, `Trypsin/P`, `Chymotrypsin`, `Glu-C`, `Lys-C`, `No Cleavage`), detect canonical `N-X-S/T` ($X \neq P$) sequons, and match precursor and optional MS/MS oxonium/Y-ion fragments.
6. **Batch Processing:** Process multiple spectrum files or a `.zip` archive with per-sample fault isolation and automatic sample-by-glycan matrix generation.
7. **Quantification:** Compute base-peak relative intensity, total-annotated-glycan relative abundance, configurable adduct aggregation (`primary_adduct_only`, `highest_intensity_peak_per_composition`, `explicit_sum_distinct_adducts`, `no_aggregation_peak_level`), and hypothesis-labeled derived glycan traits.
8. **Quality Control:** Apply linear or quadratic internal/external mass calibration with minimum calibrant safeguards and inspect MassyTools-inspired QC checks (`PASS` / `WARN` / `FAIL`) with explicit threshold provenance.
9. **Statistics and Visualization:** Run PCA, hierarchical clustering, Welch's $t$-test / Mann-Whitney $U$ / ANOVA / Kruskal-Wallis with Benjamini-Hochberg FDR correction, univariate ROC analysis, and leakage-free cross-validated classification.
10. **Export Results:** Download any of the 12 standardized export categories as CSV, TSV, multi-sheet `.xlsx` Excel workbook, structured `.json` manifest, multi-page `.pdf` report, or complete `.zip` archive.
11. **Settings:** Configure default MS1 ppm tolerance and session parameters.
12. **Help and Documentation:** Built-in scientific reference tables, Bruker conversion guidance, and ambiguity policies.

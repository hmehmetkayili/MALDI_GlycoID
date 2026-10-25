# Validation & Test Suite Reference — MALDI-GlycoID

MALDI-GlycoID includes a comprehensive automated unit, scientific, parser, UI, API, and end-to-end validation suite in `tests/`:

- `tests/test_scientific_constants.py` — Validates atomic, isotopic, monosaccharide residue, proton vs. hydrogen atom, cation, and derivatization registries.
- `tests/test_mass_calculator.py` — Validates released native, 2-AB, 2-AA, Procainamide, permethylated, linkage-specific sialic acid, peptide, and attached glycan moiety masses, multi-charge ion `m/z`, ppm/Da error, and `idotp` isotope envelope scoring.
- `tests/test_parsers.py` — Validates TXT, XY, CSV, TSV, decimal comma, ambiguous multi-column refusal, centroid vs. profile classification, `mzML`, `mzXML`, and Bruker 1D `.fid`+`acqus` & conversion guidance.
- `tests/test_library_generation.py` — Validates biosynthetic human N-glycan generation, duplicate vs. structural isomer detection, offline search filters, and JSON/CSV/TSV round-trip I/O.
- `tests/test_glycan_search.py` — Validates Mode A (targeted) and Mode B (library) MS1 search, tolerance boundaries, isobaric adduct/composition ambiguity retention, and profile/centroid area guards.
- `tests/test_quantification.py` — Validates base-peak normalization, adduct aggregation policies, and replicate CVs/missingness without silent zero-imputation.
- `tests/test_glycopeptide_search.py` — Validates in-silico protease digestion, overlapping `N-X-S/T` ($X \neq P$) sequon detection, attached glycan moiety mass accuracy, and MS/MS oxonium/Y-ion confidence categories.
- `tests/test_batch_processing.py` — Validates fault-isolated batch execution, failed-sample retry, safe in-memory ZIP extraction (path traversal & nested archive rejection), and reproducibility metadata.
- `tests/test_calibration_qc.py` — Validates linear/quadratic mass calibration, minimum calibrant safeguards, and MassyTools-inspired single-spectrum & batch QC checks.
- `tests/test_statistics.py` — Validates PCA, hierarchical clustering, BH-FDR group comparisons, univariate ROC, leakage-free cross-validation, derived glycan traits, and PNG plot generation.
- `tests/test_exports.py` — Validates all 12 standardized export categories across CSV, TSV, multi-sheet `.xlsx`, JSON, multi-page `.pdf`, and `.zip` bundles.
- `tests/test_ui.py` — Validates the 12-page HTML5/CSS3/JS web interface and endpoint bindings.
- `tests/test_api.py` — Validates all REST API endpoints, security headers, upload limits, and path-traversal protections.

## Running All Tests
```bash
python -m unittest discover -s tests -p "test_*.py" -v
```

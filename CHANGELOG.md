# Changelog

All notable changes, architectural decisions, and stage completions for **MALDI-GlycoID** are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.0.0] — 2026-10-25 (Stage 13: End-to-End Validation and Release)

### Added
- `tests/test_e2e_validation.py`: End-to-end scientific and software validation suite verifying the full analytical pipeline (default 333-entry library validation, released vs. attached glycan moiety `WATER_MASS_DA` invariant, spectrum import, linear mass calibration, Mode A & Mode B MS1 search, quantification, derived glycosylation traits, IgG1 Fc `EEQYNSTYR` glycopeptide screening with MS/MS support, 3-sample batch processing, MassyTools-inspired QC, cohort PCA & leakage-free cross-validation, 12-category ZIP/XLSX/PDF/CSV/TSV exports, raw multipart HTTP wire round-trip, and AST verification across all Python modules).
- `scripts/package_release.py` & `MALDI-GlycoID-v1.0.0.zip`: Packaged complete integrated v1.0.0 release archive.
- Verified all 73 unit and end-to-end tests pass (`Ran 73 tests in 4.568s — OK`).

---

## [0.13.0-rc1] — 2026-10-25 (Stage 12: API, Security, and Deployment)

### Added
- `app/api/compat.py`, `app/api/state.py`, `app/api/health.py`, `app/api/libraries.py`, `app/api/spectra.py`, `app/api/search.py`, `app/api/glycopeptides.py`, `app/api/batch.py`, `app/api/qc.py`, `app/api/exports.py`, `app/main.py`, and `scripts/run_local.py`: Complete REST API backend, static asset & template server, upload size & path-traversal security middleware, and local CLI launcher.
- `requirements.txt`, `.env.example`, `.gitignore`, `.dockerignore`, `Dockerfile`, `docker-compose.yml`, `.github/workflows/tests.yml`, `LICENSE`, `README.md`, and `docs/{INSTALLATION,USER_GUIDE,SCIENTIFIC_METHODS,LIBRARY_FORMAT,VALIDATION,DEPLOYMENT}.md`: Complete deployment, CI, and documentation suite.
- `tests/test_api.py`: 4 comprehensive unit and integration test suites verifying all REST API endpoints, security headers, upload/path-traversal guards, deployment files, and CLI `--check-only` mode (70 total tests passing).

---

## [0.12.0-alpha] — 2026-10-25 (Stage 11: User Interface and Workflow Integration)

### Added
- `app/templates/index.html`: Complete 12-page English scientific web application interface (`Dashboard`, `Import Data`, `Glycan Library`, `MS1 Glycan Search`, `Glycopeptide Search`, `Batch Processing`, `Quantification`, `Quality Control`, `Statistics and Visualization`, `Export Results`, `Settings`, `Help and Documentation`).
- `app/static/css/styles.css`: Responsive, accessible stylesheet with sidebar navigation, KPI cards, drag-and-drop upload zone, data tables, status badges, and mobile/tablet media queries.
- `app/static/js/app.js`: Frontend workflow controller binding all 12 views, interactive SVG stick-spectrum rendering, drag-and-drop spectrum import, library generation/filtering, Mode A & Mode B MS1 search, glycopeptide screening, batch execution, quantification & derived traits, polynomial mass calibration & QC, cohort PCA/statistics/ML, and 12-category export downloads.
- `tests/test_ui.py`: 3 automated unit tests verifying all 12 navigation views, all 12 export category cards, scientific disclosures, and JS API endpoint bindings.

---

## [0.11.0-alpha] — 2026-10-25 (Stage 10: Complete Export and Reporting Center)

### Added
- `app/exports/csv_export.py`: Standardized UTF-8 / UTF-8-BOM CSV and TSV exporters (`export_rows_to_csv_bytes`, `export_rows_to_tsv_bytes`, `export_rows_to_delimited_bytes`).
- `app/exports/excel_export.py`: Multi-sheet `.xlsx` workbook exporter and parser (`export_tables_to_xlsx_bytes`, `read_xlsx_bytes_sheets`) supporting `openpyxl` and a built-in standards-compliant ECMA-376 OpenXML `.xlsx` generator/reader.
- `app/exports/json_export.py`: Structured JSON analytical manifest exporter (`export_analysis_to_json_bytes`) with reproducibility metadata and mandatory scientific disclaimers.
- `app/exports/pdf_export.py`: Multi-page vector PDF analytical summary report generator (`export_summary_report_to_pdf_bytes`) using `matplotlib.backends.backend_pdf.PdfPages`.
- `app/exports/report_builder.py`: Unified 12-category export builder (`build_standardized_export_tables`) and complete downloadable ZIP bundler (`generate_complete_export_zip_bundle`).
- `tests/test_exports.py`: 4 comprehensive unit test suites verifying all 12 export categories, CSV/TSV round-trip, multi-sheet `.xlsx` round-trip, JSON schema, 2-page PDF via `PyPDF2`, and complete ZIP bundle contents.

---

## [0.10.0-alpha] — 2026-10-25 (Stage 9: Statistics, Visualization, and Optional Machine Learning)

### Added
- `app/analysis/statistics.py`: Comprehensive statistical analysis, visualization, derived glycan traits, and leakage-free machine learning engine.
- Descriptive statistics per glycan across samples and groups (`compute_descriptive_statistics`), Principal Component Analysis (`run_pca_analysis`) with explicit imputation/log/scaling provenance, sample/feature hierarchical clustering (`run_hierarchical_clustering`), and Pearson/Spearman correlation matrices (`compute_correlation_matrix`) with Benjamini-Hochberg FDR q-values.
- Two-group (`Welch t-test`, `Mann-Whitney U`, $\log_2$ fold change) and multi-group (`One-way ANOVA`, `Kruskal-Wallis`) statistical comparisons (`run_group_comparisons`), plus univariate ROC biomarker analysis (`compute_univariate_roc`).
- Leakage-free cross-validated classification (`run_leakage_free_classification`) encapsulating `SimpleImputer`, `StandardScaler`, and optional `SelectKBest` strictly inside each `StratifiedKFold` training fold.
- Hypothesis-labeled derived glycan traits (`compute_derived_glycan_traits`), static `matplotlib` PNG plot rendering (`render_annotated_spectrum_png`, `render_pca_plot_png`), and R `mixOmics` / `DIABLO` multi-block export table formatting (`build_mixomics_diablo_export_tables`).
- `tests/test_statistics.py`: 4 comprehensive unit test suites verifying all Stage 9 statistical, ML, plotting, and derived-trait functions.

---

## [0.9.0-alpha] — 2026-10-25 (Stage 8: Mass Calibration and Quality Control)

### Added
- `app/core/calibration.py`: Internal and external mass calibration engine (`calibrate_spectrum`, `match_calibrants_to_spectrum`, `apply_calibration_coefficients`, `get_default_nglycan_calibrants`, `get_default_peptide_calibrants`) supporting linear (`mz_cal = a * mz_obs + b`) and quadratic (`mz_cal = a * mz_obs^2 + b * mz_obs + c`) polynomial least-squares models.
- Minimum calibrant count safeguards (`>= 2` linear, `>= 3` quadratic), degrees-of-freedom warnings, narrow m/z span extrapolation warnings, before/after RMS and median ppm error diagnostics, and strict `never silently calibrate` provenance tracking.
- `app/analysis/quality_control.py`: MassyTools-inspired single-spectrum (`evaluate_spectrum_qc`) and multi-sample batch (`evaluate_batch_qc`) quality control engine with explicit `QCThresholdSource` provenance (`USER_DEFINED`, `INSTRUMENT_SPECIFIC`, `LITERATURE_DERIVED`, `PROVISIONAL_DEFAULT`), profile S/N availability guards, batch missingness %, replicate CV %, and outlier sample flagging.
- `tests/test_calibration_qc.py`: 4 comprehensive unit test suites verifying linear/quadratic calibration recovery, insufficient-calibrant safeguards, QC threshold provenance, profile S/N availability guards, and batch outlier detection.

---

## [0.8.0-alpha] — 2026-10-25 (Stage 7: Batch Processing)

### Added
- `app/analysis/batch_processing.py`: Multi-sample batch execution engine (`execute_batch_analysis`, `extract_batch_zip_safely`, `retry_failed_batch_items`, `sanitize_filename`, `BatchCancellationToken`, `BatchReproducibilityRecord`).
- Safe in-memory ZIP extraction enforcing file count limits, uncompressed byte limits, strict path-traversal rejection (`..`, `/`, `\`, `:`), and nested archive rejection (`.zip`, `.tar`, `.gz`, `.7z`, `.rar`).
- Per-sample fault isolation (`Batch item failed; other samples completed successfully`) ensuring corrupted or unparseable spectrum files are logged in `failed_files_report` while valid files complete MS1 search, quantification, and optional glycopeptide screening.
- Combined batch output tables (`combined_annotated_peaks`, `combined_glycan_candidates`, `abundance_matrix`, `sample_metadata_table`, `qc_summary_table`, `failed_files_report`, `processing_log`, `reproducibility`).
- `tests/test_batch_processing.py`: 3 comprehensive unit test suites verifying fault isolation, failed-item retry/merge, safe ZIP extraction, cancellation tokens, and reproducibility manifests.

---

## [0.7.0-alpha] — 2026-10-25 (Stage 6: Glycopeptide Screening Engine)

### Added
- `app/analysis/glycopeptide_search.py`: Dedicated glycopeptide screening engine (`parse_fasta_sequences`, `find_n_glycosylation_sequons`, `digest_protein_in_silico`, `evaluate_ms2_glycopeptide_fragments`, `screen_glycopeptides`).
- Configurable protease digestion (`Trypsin`, `Trypsin/P`, `Chymotrypsin`, `Glu-C`, `Lys-C`, `No Cleavage`) with missed cleavages and peptide modifications (`Carbamidomethyl_C`, `Oxidation_M`, `Deamidation_NQ`).
- Canonical N-linked glycosylation sequon (`N-X-S/T`, `X != P`) detection with 1-based peptide and protein sequence coordinates.
- Strict use of peptide-attached glycan moiety mass (`calculate_attached_glycan_moiety_mass`, omitting released glycan's terminal $\text{H}_2\text{O}$).
- Extensible MS/MS fragment validator checking B-type carbohydrate oxonium ions (`204.0866`, `366.1395`, `292.1027`, etc.) and Y-type peptide/core ions (`Y0`, `Y1`, `Y1F`, `Y2`) to assign explicit `GlycopeptideConfidenceCategory` levels (`MS1_MASS_MATCH`, `MS1_COMPOSITION_CANDIDATE`, `MSMS_SUPPORTED`, `SITE_LOCALIZED`, `REFERENCE_VALIDATED`).
- `data/examples/synthetic_igg1_fc.fasta` & `tests/test_glycopeptide_search.py`: Synthetic IgG1 Fc region FASTA example and 4 unit test suites.

---

## [0.6.0-alpha] — 2026-10-25 (Stage 5: Glycan Intensity, Relative Abundance, and Quantification)

### Added
- `app/analysis/quantification.py`: Dedicated quantification engine (`quantify_ms1_search_result`, `summarize_replicates_and_build_matrix`) distinguishing raw peak intensity, integrated peak area (available only from profile spectra or area-annotated tables), base-peak relative intensity (`DenominatorSetPolicy.ALL_ELIGIBLE_SPECTRUM_PEAKS` vs `ANNOTATED_GLYCAN_PEAKS_ONLY`), and total-annotated-glycan relative abundance.
- Configurable composition-level aggregation rules (`GlycanAggregationPolicy.NO_AGGREGATION_PEAK_LEVEL`, `PRIMARY_ADDUCT_ONLY`, `HIGHEST_INTENSITY_PEAK_PER_COMPOSITION`, `EXPLICIT_SUM_DISTINCT_ADDUCTS`) with full peak-level provenance tracking.
- Multi-sample replicate statistics (Mean, SD, CV %, observed replicates, missingness %) preserving missing glycan signals as `None` without silent zero-imputation.
- `tests/test_quantification.py`: 4 unit test suites verifying denominator policies, peak area vs. centroid intensity distinction, adduct aggregation provenance, and replicate statistics with missing values.

---

## [0.5.0-alpha] — 2026-10-25 (Stage 4: GlycoQuest-Like MS1 Search and Matching Engine)

### Added
- `app/analysis/peak_processing.py`: Spectrum preprocessing, intensity/m/z/S/N filtering, and profile-only baseline subtraction, local-maxima parabolic centroiding, FWHM, and trapezoidal peak area integration (`detect_and_integrate_profile_peaks`, `preprocess_spectrum`). Enforces `Insufficient profile data for peak-area calculation` when profile operations are requested on centroid lists.
- `app/analysis/glycan_search.py`: Mode A (`search_ms1_targeted`) and Mode B (`search_ms1_library`) MS1 glycan search engine with configurable ppm/Da tolerances, polarity, multi-adduct and multi-charge expansion, derivatization chemistry recalculation, transparent heuristic ranking (`ScoringBreakdown`), optional `idotp` isotope envelope scoring, one-to-many candidate retention, and explicit adduct/composition ambiguity reporting.
- `tests/test_glycan_search.py`: 6 unit test suites testing Mode A, Mode B, tolerance boundaries, `[M+Na]+`/`[M+K]+` & Neu5Gc/Neu5Ac isobaric ambiguity, 2-AB `[M+2H]2+` search, and centroid vs. profile guards.

---

## [0.4.0-alpha] — 2026-10-25 (Stage 3: Automatic Human N-Glycan Library Generation)

### Added
- `app/libraries/composition_generator.py`: Biosynthetically constrained human N-glycan composition generator (`CompositionConstraints`, `GlycanLibrary`, `GlycanLibraryEntry`, `GlycanClass`, `StructuralAnnotationStatus`) supporting high-mannose, hybrid, complex, paucimannose, fucosylated, sialylated, optional Neu5Gc, and linkage-specific sialic-acid ($\alpha 2,3/\alpha 2,6$) compositions.
- `app/libraries/library_search.py`: Offline search engine supporting composition queries, neutral mass ranges, exact target mass/`m/z` queries with ppm/Da tolerances, and class/labeling/fucosylation/sialylation filters.
- `app/libraries/library_validation.py`: Schema and Stage 1 mass consistency validator plus duplicate vs. structural isomer group detector (`detect_library_duplicates`).
- `app/libraries/library_io.py`: Round-trip JSON, CSV, and TSV library exporters and importers preserving provenance, license, and structural annotation status.
- `scripts/validate_library.py` & `data/libraries/human_nglycan_native_v1.{json,csv}`: Generated and validated the default 333-entry human N-glycan composition library.
- `tests/test_library_generation.py`: 5 unit tests verifying biosynthetic rules, linkage expansion, duplicate/isomer detection, offline search, and JSON/CSV/TSV round-trip fidelity.

---

## [0.3.0-alpha] — 2026-10-25 (Stage 2: Input Files and Peak-List Parsing)

### Added
- `app/parsers/base_models.py`: Standardized `SpectrumData`, `SpectrumPeak`, `SpectrumDataType` (`CENTROID`, `PROFILE`, `UNKNOWN`), `DuplicateMzPolicy`, and `AmbiguousColumnError`.
- `app/parsers/text_parser.py`, `app/parsers/xy_parser.py`, `app/parsers/csv_parser.py`: Multi-encoding delimited text, XY, TSV, and CSV parsers supporting comma/tab/semicolon/whitespace delimiters, `.` and `,` decimal separators, comment/header detection, negative/zero/duplicate m/z handling, and explicit refusal to silently reinterpret ambiguous multi-column layouts.
- `app/parsers/mzml_parser.py` & `app/parsers/mzxml_parser.py`: Verified XML/binary parsers for PSI-MS `mzML` and `mzXML` supporting 32-bit and 64-bit IEEE-754 floats, `zlib` decompression, centroid/profile CV metadata, retention time, and MS2 precursor metadata.
- `app/parsers/bruker_adapter.py`: Bruker MALDI-TOF adapter supporting flexAnalysis text/XY/CSV exports, unencrypted 1D `.fid` + JCAMP-DX `acqus` TOF calibration (`$ML1`, `$ML2`, `$ML3`, `$DELAY`, `$DW`), and explicit `BrukerFormatConversionRequiredError` conversion instructions for proprietary `.baf` / `.tdf` containers.
- `data/examples/`: Created clearly labeled synthetic validation spectra (`synthetic_native_nglycans_centroid.txt`, `synthetic_native_nglycans.xy`, `synthetic_native_nglycans.csv`).
- `tests/test_parsers.py`: 9 comprehensive unit tests covering all supported formats and malformed/ambiguous edge cases.

---

## [0.2.0-alpha] — 2026-10-25 (Stage 1: Scientific Constants and Mass Calculation Engine)

### Added
- `app/config.py`: Centralized application configuration, upload size limits, and default mass tolerance presets.
- `app/core/scientific_constants.py`: Single authoritative source for elemental monoisotopic masses, natural isotope abundances, monosaccharide residue vs. free masses (`Hex`, `HexNAc`, `dHex`, `Neu5Ac`, `Neu5Gc`, `Pent`, `HexA`), cationic/anionic adduct definitions (`[M+H]+`, `[M+Na]+`, `[M+K]+`, `[M+Li]+`, `[M+NH4]+`, `[M+2Na-H]+`, `[M+2H]2+`, `[M+2Na]2+`, `[M+H+Na]2+`, `[M+3H]3+`, `[M+3Na]3+`, `[M-H]-`, `[M+Cl]-`, `[M+HCOO]-`, `[M-2H]2-`), reducing-end and sialic-acid derivatization registry (`NONE`, `REDUCED`, `2-AB`, `2-AA`, `PROCAINAMIDE`, `PA`, `PERMETHYLATED`, `PERMETHYLATED_REDUCED`, `EEC_SIALIC`, `DMA_SIALIC`), and 20 canonical amino acid residue definitions.
- `app/core/formula_calculator.py`: Elemental formula parser, stoichiometry addition/scaling, Hill notation formatter, and monoisotopic formula mass calculator.
- `app/core/mass_calculator.py`: Strictly separated neutral mass calculation functions for released native/labeled/permethylated glycans (`calculate_released_glycan_mass`), peptide-attached glycan moieties (`calculate_attached_glycan_moiety_mass`), peptides (`calculate_peptide_neutral_mass`), and intact glycopeptides (`calculate_glycopeptide_neutral_mass`).
- `app/core/ion_calculator.py`: Charge-dependent `m/z` and inverse neutral mass calculators using exact cationic/anionic carrier masses.
- `app/core/mass_error.py`: Signed and absolute mass error calculations in Da and ppm, plus tolerance window filtering.
- `app/core/isotope_scoring.py`: Polynomial convolution generator for theoretical isotopic envelopes and `idotp` cosine similarity scoring (returning `None` / `available=False` when fewer than 2 isotopic peaks are observed).
- `app/core/validation.py`: Search parameter and human N-glycan biosynthetic constraint validators.
- `tests/test_scientific_constants.py` and `tests/test_mass_calculator.py`: 20 automated unit tests verifying all Stage 1 calculations.

---

## [0.1.0-alpha] — 2026-10-25 (Stage 0: Project Planning and Requirements)

### Added
- Created `PROJECT_MANIFEST.md` to track completed stages, file inventory, test metrics, known issues, and remaining milestones.
- Created `IMPLEMENTATION_STATUS.md` categorizing all features by readiness, external dependencies, and MS1 scientific boundaries.
- Created `CHANGELOG.md` to maintain a verifiable history of incremental development.
- Defined the 13-stage architectural blueprint, module dependency DAG, scientific risk register, test strategy, and deployment strategy.

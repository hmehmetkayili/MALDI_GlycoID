# MALDI-GlycoID — Project Manifest

**Project Name:** MALDI-GlycoID  
**Version:** 1.0.0 (Stage 13 — Final Release)  
**Target Release Version:** 1.0.0  
**Current Stage:** Stage 13 — End-to-End Validation and Release (Completed)  
**Next Development Stage:** None (All Stages 0–13 Completed & Validated)  
**Integration Status:** Complete Release Ready (`MALDI-GlycoID-v1.0.0.zip` Packaged & All 73 Unit/E2E Tests Passing)  

---

## 1. Completed Stages

- [x] **Stage 0 — Project Planning and Requirements**
  - Defined system architecture, module boundaries, and dependency order.
  - Established scientific constants reference plan, derivatization rules, and ambiguity reporting policies.
  - Created `PROJECT_MANIFEST.md`, `IMPLEMENTATION_STATUS.md`, and `CHANGELOG.md`.
  - Documented the Scientific Risk Register, Test Strategy, Deployment Strategy, and Technical Assumptions.
- [x] **Stage 1 — Scientific Constants and Mass Calculation Engine**
  - Implemented authoritative elemental, isotopic, monosaccharide residue, adduct, labeling, and amino acid constants (`app/core/scientific_constants.py`).
  - Implemented elemental formula parser and exact monoisotopic mass calculator (`app/core/formula_calculator.py`).
  - Implemented neutral glycan mass calculator strictly separating released native, reducing-end labeled, permethylated (native & reduced alditol), linkage-specific sialic-acid derivatized, and peptide-attached glycan moieties (`app/core/mass_calculator.py`).
  - Implemented multi-charge positive/negative adduct `m/z` calculator (`app/core/ion_calculator.py`).
  - Implemented Da/ppm signed and absolute mass error and tolerance filtering (`app/core/mass_error.py`).
  - Implemented theoretical isotopic distribution convolution and `idotp` scoring (`app/core/isotope_scoring.py`).
  - Implemented search parameter and biosynthetic N-glycan composition validators (`app/core/validation.py`).
  - Verified all 20 unit tests across `tests/test_scientific_constants.py` and `tests/test_mass_calculator.py`.
- [x] **Stage 2 — Input Files and Peak-List Parsing**
  - Implemented standardized `SpectrumData`, `SpectrumPeak`, `SpectrumDataType` (`CENTROID`, `PROFILE`, `UNKNOWN`), `DuplicateMzPolicy`, and `AmbiguousColumnError` (`app/parsers/base_models.py`).
  - Implemented delimited text/TSV (`app/parsers/text_parser.py`), two-column XY (`app/parsers/xy_parser.py`), and CSV (`app/parsers/csv_parser.py`) parsers with multi-encoding support, decimal comma/dot detection, malformed row skipping, duplicate m/z policies, and strict refusal to guess ambiguous multi-column layouts.
  - Implemented verified PSI-MS `mzML` (`app/parsers/mzml_parser.py`) and `mzXML` (`app/parsers/mzxml_parser.py`) parsers supporting 32-bit/64-bit IEEE-754 floats and `zlib` decompression.
  - Implemented Bruker MALDI-TOF adapter (`app/parsers/bruker_adapter.py`) supporting flexAnalysis exports, unencrypted 1D `.fid`+`acqus` TOF calibration, and explicit `BrukerFormatConversionRequiredError` guidance for proprietary `.baf`/`.tdf` containers.
  - Created clearly labeled synthetic example spectra in `data/examples/` and verified 9 comprehensive parser unit tests (`tests/test_parsers.py`).
- [x] **Stage 3 — Automatic Human N-Glycan Library Generation**
  - Implemented biosynthetically constrained human N-glycan composition generator (`app/libraries/composition_generator.py`) with `GlycanClass` (`paucimannose`, `high_mannose`, `hybrid`, `complex`), `StructuralAnnotationStatus` (`COMPOSITION_CANDIDATE`, `TOPOLOGY_HYPOTHESIS`, `LINKAGE_SPECIFIC_DERIVATIZED`, `STRUCTURALLY_VALIDATED_REFERENCE`), optional Neu5Gc flagging, and linkage-specific $\alpha 2,3/\alpha 2,6$ sialic acid expansion.
  - Implemented offline library querying and filtering (`app/libraries/library_search.py`).
  - Implemented library schema, Stage 1 mass consistency, and duplicate vs. structural isomer detection (`app/libraries/library_validation.py`).
  - Implemented JSON, CSV, and TSV library import/export preserving provenance and license metadata (`app/libraries/library_io.py`).
  - Generated default 333-entry human N-glycan library in `data/libraries/` via `scripts/validate_library.py` and verified 5 unit test suites in `tests/test_library_generation.py`.
- [x] **Stage 4 — GlycoQuest-Like MS1 Search and Matching Engine**
  - Implemented spectrum peak preprocessing (`app/analysis/peak_processing.py`) with strict centroid vs. profile separation (`Insufficient profile data for peak-area calculation` guard on centroid lists; parabolic centroiding, baseline subtraction, local MAD S/N, FWHM, and trapezoidal area integration on profile spectra).
  - Implemented Mode A (`search_ms1_targeted`) and Mode B (`search_ms1_library`) MS1 glycan search engine (`app/analysis/glycan_search.py`) supporting configurable ppm/Da tolerances, polarity, multiple adducts, multiple charge states, derivatization chemistries, and candidate reporting limits.
  - Implemented transparent heuristic candidate ranking (`ScoringBreakdown`) explicitly separated from calibrated probability, optional `idotp` isotope envelope scoring, one-to-many alternative candidate retention, and explicit isobaric adduct/composition overlap warnings.
  - Verified 6 comprehensive MS1 search unit tests in `tests/test_glycan_search.py`.
- [x] **Stage 5 — Glycan Intensity, Relative Abundance, and Quantification**
  - Implemented dedicated quantification engine (`app/analysis/quantification.py`) strictly separating raw peak intensity, integrated peak area (never fabricated on centroid-only data), base-peak-normalized relative intensity (`DenominatorSetPolicy`), and total-annotated-glycan relative abundance.
  - Implemented configurable composition-level aggregation rules (`GlycanAggregationPolicy`: `NO_AGGREGATION_PEAK_LEVEL`, `PRIMARY_ADDUCT_ONLY`, `HIGHEST_INTENSITY_PEAK_PER_COMPOSITION`, `EXPLICIT_SUM_DISTINCT_ADDUCTS`) with complete peak-level provenance history.
  - Implemented multi-sample replicate statistical summary (`summarize_replicates_and_build_matrix`: Mean, SD, CV %, observed replicates, missingness %) preserving missing glycan signals as `None` without silent zero-imputation.
  - Verified 4 comprehensive quantification unit tests in `tests/test_quantification.py`.
- [x] **Stage 6 — Glycopeptide Screening Engine**
  - Implemented FASTA/peptide sequence parser, configurable protease digestion (`Trypsin`, `Trypsin/P`, `Chymotrypsin`, `Glu-C`, `Lys-C`, `No Cleavage`), and `N-X-S/T` (`X != P`) sequon detector with 1-based peptide and protein coordinates (`app/analysis/glycopeptide_search.py`).
  - Implemented glycopeptide precursor matching strictly using the peptide-attached glycan moiety mass (`calculate_attached_glycan_moiety_mass` without released glycan's terminal $\text{H}_2\text{O}$).
  - Implemented MS/MS oxonium (`HexNAc+`, `HexHexNAc+`, `Neu5Ac+`, etc.) and Y-ion (`Y0`, `Y1`, `Y1F`, `Y2`) fragment evaluation and strict confidence categories (`MS1_MASS_MATCH`, `MS1_COMPOSITION_CANDIDATE`, `MSMS_SUPPORTED`, `SITE_LOCALIZED`, `REFERENCE_VALIDATED`).
  - Created `data/examples/synthetic_igg1_fc.fasta` and verified 4 comprehensive unit tests in `tests/test_glycopeptide_search.py` (48 total tests passing).
- [x] **Stage 7 — Batch Processing**
  - Implemented multi-sample batch execution, safe in-memory ZIP extraction with path-traversal and nested-archive rejection, per-sample fault isolation (`Batch item failed; other samples completed successfully`), failed-sample retry/merge, cancellation token handling, combined sample-by-glycan matrices, and full reproducibility metadata logging (`app/analysis/batch_processing.py`).
  - Verified 3 comprehensive batch processing test suites in `tests/test_batch_processing.py` (51 total tests passing).
- [x] **Stage 8 — Mass Calibration and Quality Control**
  - Implemented linear and quadratic least-squares mass calibration (`app/core/calibration.py`) with reference N-glycan and peptide calibrants, minimum calibrant count safeguards (`>= 2` linear, `>= 3` quadratic), span/residual warnings, and explicit calibration provenance metadata (`never silently calibrate`).
  - Implemented MassyTools-inspired single-spectrum and batch quality control engine (`app/analysis/quality_control.py`) with explicit threshold provenance (`USER_DEFINED`, `INSTRUMENT_SPECIFIC`, `LITERATURE_DERIVED`, `PROVISIONAL_DEFAULT`), profile S/N availability guards, batch missingness, replicate CVs, and outlier flagging.
  - Verified 4 comprehensive calibration and QC test suites in `tests/test_calibration_qc.py` (55 total tests passing).
- [x] **Stage 9 — Statistics, Visualization, and Optional Machine Learning**
  - Implemented descriptive statistics, PCA with explicit imputation/scaling provenance, sample/feature hierarchical clustering, Pearson/Spearman correlations with Benjamini-Hochberg FDR, two-group (Welch t-test, Mann-Whitney U, $\log_2$ fold change) and multi-group (One-way ANOVA, Kruskal-Wallis) comparisons, univariate ROC curves, leakage-free stratified cross-validation (`scikit-learn` `Pipeline` fitted strictly inside each training fold), hypothesis-labeled derived glycan traits, static `matplotlib` PNG plot rendering, and R `mixOmics` / `DIABLO` multi-block export tables (`app/analysis/statistics.py`).
  - Verified 4 comprehensive statistical and ML test suites in `tests/test_statistics.py` (59 total tests passing).
- [x] **Stage 10 — Complete Export and Reporting Center**
  - Implemented all 12 standardized export categories (`annotated_peaks`, `all_glycan_candidates`, `quantified_glycans`, `replicate_abundance_matrix`, `derived_glycan_traits`, `glycopeptide_candidates`, `calibration_diagnostics`, `quality_control_summary`, `statistical_analysis`, `glycan_reference_library`, `analysis_parameters_and_reproducibility`, `processing_and_error_log`) across CSV, TSV, multi-sheet `.xlsx` (ECMA-376 OpenXML / `openpyxl`), structured `.json`, multi-page vector `.pdf` summary reports (`matplotlib.backends.backend_pdf.PdfPages`), and complete `.zip` bundles (`app/exports/`).
  - Verified 4 comprehensive export test suites in `tests/test_exports.py` (63 total tests passing).
- [x] **Stage 11 — User Interface and Workflow Integration**
  - Implemented responsive, accessible 12-page English web application interface (`app/templates/index.html`, `app/static/css/styles.css`, `app/static/js/app.js`) covering `Dashboard`, `Import Data`, `Glycan Library`, `MS1 Glycan Search`, `Glycopeptide Search`, `Batch Processing`, `Quantification`, `Quality Control`, `Statistics and Visualization`, `Export Results`, `Settings`, and `Help and Documentation`.
  - Verified 3 automated UI & workflow integration test suites in `tests/test_ui.py` (66 total tests passing).
- [x] **Stage 12 — API, Security, and Deployment**
  - Implemented full REST API backend (`app/main.py`, `app/api/compat.py`, `app/api/state.py`, `app/api/health.py`, `app/api/libraries.py`, `app/api/spectra.py`, `app/api/search.py`, `app/api/glycopeptides.py`, `app/api/batch.py`, `app/api/qc.py`, `app/api/exports.py`) and local CLI runner (`scripts/run_local.py`).
  - Implemented security middleware headers (`X-Content-Type-Options: nosniff`, `X-Frame-Options: DENY`, `Referrer-Policy`, `Content-Security-Policy`), upload size enforcement, and path-traversal protection.
  - Created deployment and documentation files (`requirements.txt`, `.env.example`, `.gitignore`, `.dockerignore`, `Dockerfile`, `docker-compose.yml`, `.github/workflows/tests.yml`, `LICENSE`, `README.md`, `docs/INSTALLATION.md`, `docs/USER_GUIDE.md`, `docs/SCIENTIFIC_METHODS.md`, `docs/LIBRARY_FORMAT.md`, `docs/VALIDATION.md`, `docs/DEPLOYMENT.md`).
  - Verified 4 comprehensive API, security, and deployment test suites in `tests/test_api.py` (70 total tests passing).
- [x] **Stage 13 — End-to-End Validation and Release**
  - Implemented end-to-end validation test suite (`tests/test_e2e_validation.py`) verifying the full synthetic acceptance workflow, multipart HTTP wire parsing, AST completeness verification (zero stub functions across all Python files), and 12-category export bundle round-trip.
  - Created `scripts/package_release.py` and packaged the complete integrated repository into `MALDI-GlycoID-v1.0.0.zip`.
  - Verified all 73 automated unit, scientific, parser, UI, API, and E2E validation tests pass (`Ran 73 tests in 4.568s — OK`).

---

## 2. Planned Stages & Roadmap

- All 14 stages (Stage 0 through Stage 13) are **COMPLETED & VALIDATED**.

---

## 3. Files Created

| File Path | Stage Created | Purpose | Status |
|---|---|---|---|
| `PROJECT_MANIFEST.md` | Stage 0 | Tracks project state, files, tests, bugs, and stage progression | Active |
| `IMPLEMENTATION_STATUS.md` | Stage 0 | Detailed feature-by-feature and module-by-module implementation status | Active |
| `CHANGELOG.md` | Stage 0 | Chronological log of stage completions and architectural changes | Active |
| `app/__init__.py` | Stage 1 | Root application package metadata | Active |
| `app/config.py` | Stage 1 | Global configuration and scientific/security defaults | Active |
| `app/core/__init__.py` | Stage 1 | Core scientific engine package init | Active |
| `app/core/scientific_constants.py` | Stage 1 | Authoritative atomic, residue, adduct, and labeling definitions | Validated |
| `app/core/formula_calculator.py` | Stage 1 | Elemental formula parser, stoichiometry, and Hill formatter | Validated |
| `app/core/mass_calculator.py` | Stage 1 | Released, labeled, permethylated, attached glycan & peptide mass engine | Validated |
| `app/core/ion_calculator.py` | Stage 1 | Positive/negative ion adduct and multi-charge m/z calculator | Validated |
| `app/core/mass_error.py` | Stage 1 | Signed/absolute Da & ppm error and tolerance filtering | Validated |
| `app/core/isotope_scoring.py` | Stage 1 | Theoretical isotope envelope convolution and `idotp` similarity | Validated |
| `app/core/validation.py` | Stage 1 | Search parameter and human N-glycan biosynthetic validators | Validated |
| `tests/test_scientific_constants.py` | Stage 1 | Unit tests for scientific constants (9 tests) | Passing |
| `tests/test_mass_calculator.py` | Stage 1 | Unit tests for mass, ion, error, and isotope calculators (11 tests) | Passing |
| `app/parsers/__init__.py` | Stage 2 | Unified spectrum parser dispatcher (`parse_spectrum_file`) | Validated |
| `app/parsers/base_models.py` | Stage 2 | `SpectrumData`, `SpectrumPeak`, `SpectrumDataType`, and deduplication models | Validated |
| `app/parsers/text_parser.py` | Stage 2 | Delimited TXT/TSV parser with encoding & column detection | Validated |
| `app/parsers/xy_parser.py` | Stage 2 | Two-column XY spectrum parser | Validated |
| `app/parsers/csv_parser.py` | Stage 2 | Configurable CSV peak/spectrum parser | Validated |
| `app/parsers/mzml_parser.py` | Stage 2 | Verified PSI-MS mzML 32/64-bit & zlib parser | Validated |
| `app/parsers/mzxml_parser.py` | Stage 2 | Verified mzXML big-endian 32/64-bit & zlib parser | Validated |
| `app/parsers/bruker_adapter.py` | Stage 2 | Bruker flexAnalysis export, 1D `.fid`+`acqus` parser, & conversion workflow | Validated |
| `data/examples/synthetic_native_nglycans_centroid.txt` | Stage 2 | Synthetic validation TXT peak list | Active |
| `data/examples/synthetic_native_nglycans.xy` | Stage 2 | Synthetic validation XY peak list | Active |
| `data/examples/synthetic_native_nglycans.csv` | Stage 2 | Synthetic validation CSV peak list with area & S/N | Active |
| `tests/test_parsers.py` | Stage 2 | Unit tests for all spectrum parsers and edge cases (9 tests) | Passing |
| `app/libraries/__init__.py` | Stage 3 | Public N-glycan library package interface | Validated |
| `app/libraries/composition_generator.py` | Stage 3 | Biosynthetically constrained human N-glycan composition generator | Validated |
| `app/libraries/library_search.py` | Stage 3 | Offline composition, mass-range, exact mass, and class search engine | Validated |
| `app/libraries/library_validation.py` | Stage 3 | Library schema, Stage 1 mass consistency, and isomer/duplicate detector | Validated |
| `app/libraries/library_io.py` | Stage 3 | JSON, CSV, and TSV library import/export with provenance retention | Validated |
| `scripts/validate_library.py` | Stage 3 | CLI script to generate default library and run validation | Validated |
| `data/libraries/human_nglycan_native_v1.json` | Stage 3 | Default 333-entry human N-glycan library (JSON) | Validated |
| `data/libraries/human_nglycan_native_v1.csv` | Stage 3 | Default 333-entry human N-glycan library (CSV) | Validated |
| `tests/test_library_generation.py` | Stage 3 | Unit tests for library generation, search, validation, and I/O (5 tests) | Passing |
| `app/analysis/__init__.py` | Stage 4 | Public analysis package interface | Validated |
| `app/analysis/peak_processing.py` | Stage 4 | Centroid filtering & profile baseline/peak-picking/area integration | Validated |
| `app/analysis/glycan_search.py` | Stage 4 | Mode A (targeted) & Mode B (library) MS1 search, ranking, and ambiguity engine | Validated |
| `tests/test_glycan_search.py` | Stage 4 | Unit tests for MS1 search modes, tolerances, adducts, and ambiguity (6 tests) | Passing |
| `app/analysis/quantification.py` | Stage 5 | Base-peak normalization, peak area guard, adduct aggregation, & replicate matrix engine | Validated |
| `tests/test_quantification.py` | Stage 5 | Unit tests for quantification, aggregation policies, and missingness (4 tests) | Passing |
| `app/analysis/glycopeptide_search.py` | Stage 6 | In-silico protease digestion, `N-X-S/T` sequons, attached glycan moiety & MS1/MS2 matching | Validated |
| `data/examples/synthetic_igg1_fc.fasta` | Stage 6 | Synthetic human IgG1 Fc region FASTA example containing `EEQYNSTYR` | Active |
| `tests/test_glycopeptide_search.py` | Stage 6 | Unit tests for digestion, sequons, attached moiety mass, and MS2 confidence categories (4 tests) | Passing |
| `app/analysis/batch_processing.py` | Stage 7 | Multi-sample batch execution, safe ZIP extraction, fault isolation, retry, & reproducibility manifest | Validated |
| `tests/test_batch_processing.py` | Stage 7 | Unit tests for mixed valid/invalid batch isolation, retry, ZIP security, and cancellation (3 tests) | Passing |
| `app/core/calibration.py` | Stage 8 | Linear & quadratic mass calibration, calibrant matching, & residual diagnostics | Validated |
| `app/analysis/quality_control.py` | Stage 8 | MassyTools-inspired spectrum & batch QC evaluation with explicit threshold provenance | Validated |
| `tests/test_calibration_qc.py` | Stage 8 | Unit tests for linear/quadratic calibration, safeguards, S/N guards, and batch QC (4 tests) | Passing |
| `app/analysis/statistics.py` | Stage 9 | Descriptive stats, PCA, clustering, group tests, ROC, leakage-free ML, derived traits, & PNG plots | Validated |
| `tests/test_statistics.py` | Stage 9 | Unit tests for PCA, clustering, FDR group comparisons, ROC, leakage-free CV, and derived traits (4 tests) | Passing |
| `app/exports/__init__.py` | Stage 10 | Public export and reporting package interface | Validated |
| `app/exports/csv_export.py` | Stage 10 | Standardized CSV and TSV table exporters | Validated |
| `app/exports/excel_export.py` | Stage 10 | Multi-sheet `.xlsx` ECMA-376 OpenXML / `openpyxl` exporter & reader | Validated |
| `app/exports/json_export.py` | Stage 10 | Structured JSON analytical manifest exporter | Validated |
| `app/exports/pdf_export.py` | Stage 10 | Multi-page vector PDF analytical summary report generator | Validated |
| `app/exports/report_builder.py` | Stage 10 | 12-category export table builder and complete ZIP bundle generator | Validated |
| `tests/test_exports.py` | Stage 10 | Unit tests for 12 export categories, CSV/TSV, XLSX, JSON, PDF, and ZIP bundle (4 tests) | Passing |
| `app/templates/index.html` | Stage 11 | 12-page responsive English scientific web application interface | Validated |
| `app/static/css/styles.css` | Stage 11 | Responsive, accessible CSS3 stylesheet | Validated |
| `app/static/js/app.js` | Stage 11 | Frontend workflow controller, SVG spectrum renderer, & API bindings | Validated |
| `tests/test_ui.py` | Stage 11 | Automated verification of 12 navigation pages, 12 export cards, and JS bindings (3 tests) | Passing |
| `app/api/__init__.py` | Stage 12 | Public REST API package interface | Validated |
| `app/api/compat.py` | Stage 12 | HTTP router, security headers middleware, `TestClient`, & local HTTP server | Validated |
| `app/api/state.py` | Stage 12 | In-memory analytical session state manager | Validated |
| `app/api/health.py` | Stage 12 | `/api/health` endpoint | Validated |
| `app/api/libraries.py` | Stage 12 | `/api/libraries/*` generation, query, import, & export endpoints | Validated |
| `app/api/spectra.py` | Stage 12 | `/api/spectra/upload` & `/api/spectra/example` endpoints with path-traversal guards | Validated |
| `app/api/search.py` | Stage 12 | `/api/search/ms1` & `/api/search/quantify` endpoints | Validated |
| `app/api/glycopeptides.py` | Stage 12 | `/api/glycopeptides/screen` endpoint | Validated |
| `app/api/batch.py` | Stage 12 | `/api/batch/demo`, `/api/batch/run`, & `/api/batch/retry` endpoints | Validated |
| `app/api/qc.py` | Stage 12 | `/api/qc/calibrate`, `/api/qc/evaluate`, & `/api/qc/statistics-demo` endpoints | Validated |
| `app/api/exports.py` | Stage 12 | `/api/exports/download` endpoint across ZIP, XLSX, PDF, JSON, CSV, & TSV | Validated |
| `app/main.py` | Stage 12 | Main web application entrypoint and router registration | Validated |
| `scripts/run_local.py` | Stage 12 | CLI launcher script (`--host`, `--port`, `--check-only`) | Validated |
| `requirements.txt`, `.env.example`, `.gitignore`, `.dockerignore`, `Dockerfile`, `docker-compose.yml`, `.github/workflows/tests.yml`, `LICENSE`, `README.md` | Stage 12 | Containerization, CI workflow, and repository configuration files | Validated |
| `docs/INSTALLATION.md`, `docs/USER_GUIDE.md`, `docs/SCIENTIFIC_METHODS.md`, `docs/LIBRARY_FORMAT.md`, `docs/VALIDATION.md`, `docs/DEPLOYMENT.md` | Stage 12 | Comprehensive scientific, user, and deployment documentation | Validated |
| `tests/test_api.py` | Stage 12 | Unit and integration tests for REST API, security guards, and deployment configs (4 tests) | Passing |
| `tests/test_e2e_validation.py` | Stage 13 | End-to-end pipeline acceptance, wire protocol, and AST completeness tests (3 tests) | Passing |
| `scripts/package_release.py` | Stage 13 | Release archive builder generating `MALDI-GlycoID-v1.0.0.zip` | Validated |
| `MALDI-GlycoID-v1.0.0.zip` | Stage 13 | Complete integrated release archive | Validated |

---

## 4. Files Modified

| File Path | Stage Modified | Description of Changes |
|---|---|---|
| `PROJECT_MANIFEST.md` | Stage 0–11 | Updated with Stage 11 files and test results |
| `IMPLEMENTATION_STATUS.md` | Stage 0–11 | Marked Stage 11 modules as `IMPLEMENTED & VALIDATED` |
| `CHANGELOG.md` | Stage 0–11 | Logged Stage 11 implementation and verification |
| `app/analysis/__init__.py` | Stage 4–9 | Exported Stage 9 statistical, derived trait, ML, and plotting functions |

---

## 5. Test Execution Summary

| Test Suite / File | Stage | Passed | Failed | Skipped / Not Run | Notes |
|---|---|---|---|---|---|
| `tests/test_scientific_constants.py` | Stage 1 | 9 | 0 | 0 | Verified residue masses, free vs residue, proton vs H atom, cations, and registries |
| `tests/test_mass_calculator.py` | Stage 1 | 11 | 0 | 0 | Verified native, 2-AB, 2-AA, Procainamide, permethylation, EEC sialic, ions, ppm, and `idotp` |
| `tests/test_parsers.py` | Stage 2 | 9 | 0 | 0 | Verified TXT, XY, CSV, TSV, decimal commas, ambiguous columns, profile vs centroid, mzML, mzXML, and Bruker adapter |
| `tests/test_library_generation.py` | Stage 3 | 5 | 0 | 0 | Verified biosynthetic rules, Neu5Gc & EEC expansion, isomer vs duplicate detection, offline search, and JSON/CSV/TSV round-trip |
| `tests/test_glycan_search.py` | Stage 4 | 6 | 0 | 0 | Verified Mode A & Mode B search, tolerance boundaries, `[M+Na]+`/`[M+K]+` & Neu5Gc ambiguity, 2-AB `[M+2H]2+`, and profile/centroid guard |
| `tests/test_quantification.py` | Stage 5 | 4 | 0 | 0 | Verified base-peak normalization, area vs intensity, adduct aggregation provenance, and replicate CV/missingness without silent imputation |
| `tests/test_glycopeptide_search.py` | Stage 6 | 4 | 0 | 0 | Verified `N-X-S/T` (`X != P`) sequons, tryptic digestion, attached moiety vs released H2O guard, and MS2 oxonium/Y-ion confidence categories |
| `tests/test_batch_processing.py` | Stage 7 | 3 | 0 | 0 | Verified mixed valid/corrupt sample isolation, retry & merge, safe ZIP path-traversal/nested rejection, cancellation, and reproducibility metadata |
| `tests/test_calibration_qc.py` | Stage 8 | 4 | 0 | 0 | Verified linear/quadratic drift correction, minimum calibrant safeguards, explicit threshold sources, S/N availability guard, and batch outlier QC |
| `tests/test_statistics.py` | Stage 9 | 4 | 0 | 0 | Verified PCA, clustering, Pearson/Spearman BH-FDR, Welch/Mann-Whitney/ANOVA/Kruskal, ROC, leakage-free CV pipeline, derived traits, and PNG plots |
| `tests/test_exports.py` | Stage 10 | 4 | 0 | 0 | Verified all 12 export categories, CSV/TSV round-trip, multi-sheet `.xlsx` round-trip, JSON schema, 2-page PDF via `PyPDF2`, and ZIP bundle |
| `tests/test_ui.py` | Stage 11 | 3 | 0 | 0 | Verified all 12 navigation views, all 12 export category cards, scientific disclosures, and JS API endpoint bindings |
| `tests/test_api.py` | Stage 12 | 4 | 0 | 0 | Verified all REST API endpoints, security headers, upload/path-traversal guards, deployment files, and CLI `--check-only` |
| `tests/test_e2e_validation.py` | Stage 13 | 3 | 0 | 0 | Verified complete E2E synthetic acceptance workflow, raw multipart wire round-trip, and AST repository completeness (73 total tests passing) |

---

## 6. Known Bugs & Environment Observations

- **Known Application Bugs:** None (Stage 0).
- **Environment Observation:** Verified Python 3.11.11 (`C:\FragPipe\FragPipe-24.0\python\python.exe`) is available on the host Windows environment for running automated unit tests and application validation across Stages 1–13.

---

## 7. Outstanding Requirements (Stages 1–13)

1. **Core Mass & Ion Engine (Stage 1):** Authoritative monoisotopic masses (IUPAC/NIST/Unimod), proton/cation adduct conventions, native vs. labeled vs. permethylated vs. linkage-specific sialic acid derivatization, attached glycan moiety convention, ppm/Da error metrics, and theoretical isotopic distribution calculation.
2. **File Parsers (Stage 2):** TXT, XY, CSV, TSV, mzML, mzXML, and explicit Bruker conversion/XY-export adapter; strict distinction between centroid peak lists and profile spectra.
3. **N-Glycan Library Engine (Stage 3):** Biosynthetically constrained human N-glycan composition generator (Hex, HexNAc, dHex/Fuc, Neu5Ac, optional Neu5Gc), distinction between composition and topology/isomers, offline search and JSON/CSV/TSV I/O.
4. **MS1 Search Engine (Stage 4):** Targeted (Mode A) and Library-based (Mode B) MS1 search, one-to-many ambiguity retention, transparent heuristic ranking vs. optional `idotp` isotope scoring.
5. **Quantification Engine (Stage 5):** Base-peak normalization, profile-only peak integration, configurable adduct/charge aggregation with full provenance, replicate CV and missingness tracking without silent imputation.
6. **Glycopeptide Screening Engine (Stage 6):** FASTA/peptide input, in-silico trypsin digestion, `N-X-S/T` (`X != P`) sequon detection, peptide + attached glycan moiety mass calculation, precursor matching, extensible MS/MS fragment hook, explicit confidence categories (`MS1_MASS_MATCH`, `MS1_COMPOSITION_CANDIDATE`, etc.).
7. **Batch Processing (Stage 7):** Multi-file and safe ZIP batch execution, per-sample isolation, reproducibility metadata logging, combined matrices.
8. **Mass Calibration & QC (Stage 8):** Linear/quadratic internal/external mass calibration with minimum calibrant safeguards, residual analysis, MassyTools-inspired QC metrics with explicit threshold provenance.
9. **Statistics, Visualization & Optional ML (Stage 9):** PCA, hierarchical clustering, group comparisons, ROC, leakage-free cross-validated classification, hypothesis-labeled derived glycan traits.
10. **Export & Reporting Center (Stage 10):** 12 standardized export artifacts across CSV, TSV, multi-sheet XLSX, JSON, and PDF summary report, plus ZIP bundle.
11. **Web User Interface (Stage 11):** 12-page responsive English HTML5/CSS3/JS interface wired to real backend endpoints.
12. **API, Security & Deployment (Stage 12):** FastAPI router integration, upload validation, path-traversal protection, `Dockerfile`, `docker-compose.yml`, GitHub Actions CI, and local/cloud documentation.
13. **End-to-End Validation & Packaging (Stage 13):** Full scientific and software validation suite, synthetic acceptance workflow, and `MALDI-GlycoID-v1.0.0.zip` deliverable archive.

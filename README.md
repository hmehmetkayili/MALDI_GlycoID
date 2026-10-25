# MALDI-GlycoID

**MALDI-GlycoID** is a modular, scientifically transparent software platform for **MALDI-TOF mass spectrometry-based human N-glycan compositional profiling, relative quantification, quality control, and glycopeptide screening**.

---

## Core Scientific Principles & Safeguards

1. **Authoritative Monoisotopic Mass Engine:** All elemental, isotopic, monosaccharide residue, adduct, and reducing-end derivatization calculations use IUPAC/NIST/Unimod reference masses (`app/core/scientific_constants.py`). Proton (`[M+H]+`, `1.007276466621 Da`) and neutral hydrogen atom (`1.00782503223 Da`) masses are strictly distinguished.
2. **Released vs. Peptide-Attached Glycan Moiety Separation:** Released glycans include the terminal reducing-end $\text{H}_2\text{O}$ (`18.0105646865 Da`), whereas peptide-attached N-glycan moieties (`calculate_attached_glycan_moiety_mass`) omit the reducing-end $\text{H}_2\text{O}$ lost during glycosidic bond formation to asparagine.
3. **Centroid vs. Profile Spectrum Integrity:** Centroid peak lists and profile spectra are explicitly classified. Peak area integration, baseline subtraction, and local S/N estimation are never fabricated on centroid-only peak lists (`Insufficient profile data for peak-area calculation`).
4. **Composition Candidates vs. Structural Isomers:** MS1 mass matches represent candidate monosaccharide compositions (`COMPOSITION_CANDIDATE`), not unique glycosidic linkages or branching topologies, unless linkage-specific sialic acid derivatization (`EEC_SIALIC`, `DMA_SIALIC`) or validated reference standards are used.
5. **No Silent Missing-Value Imputation in Replicate Summaries:** Missing glycan signals across biological/technical replicates are preserved as `None`/`NaN` in replicate summaries and only imputed when explicitly configured for multivariate PCA/clustering/classification (`HALF_MIN_OBSERVED`, `MEDIAN`, `MEAN`, `ZERO`, `KNN`).
6. **Leakage-Free Cross-Validated Machine Learning:** In cohort classification (`app/analysis/statistics.py`), missing-value imputation, feature scaling, and feature selection are fitted strictly inside each cross-validation training fold (`scikit-learn` `Pipeline`).

---

## Quick Start

### 1. Run Local Web Application
```bash
python scripts/run_local.py --host 127.0.0.1 --port 8000
```
Then open `http://127.0.0.1:8000` in your browser.

### 2. Verify Health & Default N-Glycan Library
```bash
python scripts/run_local.py --check-only
python scripts/validate_library.py
```

### 3. Run Full Automated Test Suite
```bash
python -m unittest discover -s tests -p "test_*.py" -v
```

### 4. Run via Docker Compose
```bash
docker compose up --build
```

---

## Documentation Index

- [Installation Guide](docs/INSTALLATION.md)
- [User Guide (12-Page Web Interface & Workflows)](docs/USER_GUIDE.md)
- [Scientific Methods & Mathematical Formulas](docs/SCIENTIFIC_METHODS.md)
- [N-Glycan Library Format & Schema](docs/LIBRARY_FORMAT.md)
- [Validation & Test Suite Reference](docs/VALIDATION.md)
- [Deployment & Security Guide](docs/DEPLOYMENT.md)
- [Project Manifest](PROJECT_MANIFEST.md)
- [Implementation Status](IMPLEMENTATION_STATUS.md)
- [Changelog](CHANGELOG.md)

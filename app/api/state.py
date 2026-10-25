"""In-memory analytical session state manager for MALDI-GlycoID API endpoints (Stage 12)."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from app.analysis.batch_processing import BatchAnalysisResult
from app.analysis.glycan_search import MS1SearchResult, MS1SearchSettings, search_ms1_library
from app.analysis.glycopeptide_search import GlycopeptideScreeningResult
from app.analysis.quality_control import BatchQCReport, SpectrumQCReport, evaluate_spectrum_qc
from app.analysis.quantification import (
    MultiSampleQuantificationMatrix,
    SampleQuantificationResult,
    quantify_ms1_search_result,
)
from app.analysis.statistics import (
    CrossValidatedClassificationResult,
    GroupComparisonRecord,
    PCAResult,
    SampleDerivedTraitsReport,
    compute_derived_glycan_traits,
)
from app.analysis.spectrum_corrections import CorrectionResult
from app.config import EXAMPLES_DIR, LIBRARIES_DIR
from app.core.calibration import MassCalibrationResult
from app.libraries.composition_generator import GlycanLibrary, generate_human_n_glycan_library
from app.libraries.library_io import import_library_from_json
from app.parsers import parse_spectrum_file
from app.parsers.base_models import SpectrumData


@dataclass
class AnalyticalSessionState:
    """Holds the active library, spectrum, and analytical workflow results."""

    active_library: Optional[GlycanLibrary] = None
    active_spectrum: Optional[SpectrumData] = None
    original_spectrum: Optional[SpectrumData] = None  # as imported
    base_spectrum: Optional[SpectrumData] = None  # after optional corrections, before calibration
    correction_result: Optional[CorrectionResult] = None
    ms1_search_result: Optional[MS1SearchResult] = None
    quantification_result: Optional[SampleQuantificationResult] = None
    derived_traits_report: Optional[SampleDerivedTraitsReport] = None
    glycopeptide_result: Optional[GlycopeptideScreeningResult] = None
    calibration_result: Optional[MassCalibrationResult] = None
    spectrum_qc_report: Optional[SpectrumQCReport] = None
    batch_result: Optional[BatchAnalysisResult] = None
    batch_qc_report: Optional[BatchQCReport] = None
    abundance_matrix: Optional[MultiSampleQuantificationMatrix] = None
    pca_result: Optional[PCAResult] = None
    cv_classification_result: Optional[CrossValidatedClassificationResult] = None
    group_comparisons: List[GroupComparisonRecord] = field(default_factory=list)

    def get_or_load_library(self) -> GlycanLibrary:
        """Return the active N-glycan library, loading the default 333-entry library if needed."""
        if self.active_library is not None:
            return self.active_library

        default_json = LIBRARIES_DIR / "human_nglycan_native_v1.json"
        if default_json.is_file():
            self.active_library = import_library_from_json(default_json)
        else:
            self.active_library = generate_human_n_glycan_library()
        return self.active_library

    def reset_to_default_library(self) -> GlycanLibrary:
        """Force reload of the default 333-entry human N-glycan library."""
        self.active_library = None
        return self.get_or_load_library()

    def get_or_load_spectrum(self) -> SpectrumData:
        """Return the active spectrum, or load the synthetic CSV example spectrum if none is active."""
        if self.active_spectrum is not None:
            return self.active_spectrum

        default_csv = EXAMPLES_DIR / "synthetic_native_nglycans.csv"
        self.set_imported_spectrum(parse_spectrum_file(default_csv, sample_id="SYNTH_NATIVE_01"))
        return self.active_spectrum

    # ------------------------------------------------------------------ #
    # Spectrum lifecycle: import -> (optional corrections) -> (optional calibration) -> search
    # ------------------------------------------------------------------ #
    def invalidate_results(self) -> None:
        """Drop every downstream result that depends on the active spectrum."""
        self.ms1_search_result = None
        self.quantification_result = None
        self.derived_traits_report = None
        self.spectrum_qc_report = None

    def set_imported_spectrum(self, spectrum: SpectrumData) -> None:
        """Register a freshly imported spectrum and reset the whole processing chain."""
        self.original_spectrum = spectrum
        self.base_spectrum = spectrum
        self.active_spectrum = spectrum
        self.correction_result = None
        self.calibration_result = None
        self.invalidate_results()

    def set_processed_spectrum(self, spectrum: SpectrumData, correction: Optional[CorrectionResult]) -> None:
        """Store the corrected (pre-calibration) spectrum; calibration must be re-run afterwards."""
        self.base_spectrum = spectrum
        self.active_spectrum = spectrum
        self.correction_result = correction
        self.calibration_result = None
        self.invalidate_results()

    def set_calibration(self, result: MassCalibrationResult) -> None:
        """Store a calibration result; the calibrated spectrum becomes active when it was applied."""
        self.calibration_result = result
        self.active_spectrum = result.calibrated_spectrum if result.applied else self.base_spectrum
        self.invalidate_results()

    def reset_calibration(self) -> None:
        """Undo calibration (back to the corrected/imported spectrum)."""
        self.calibration_result = None
        if self.base_spectrum is not None:
            self.active_spectrum = self.base_spectrum
        self.invalidate_results()

    def get_base_spectrum(self) -> SpectrumData:
        """Spectrum to calibrate: always the uncalibrated base, so repeated calibration never compounds."""
        self.get_or_load_spectrum()
        return self.base_spectrum if self.base_spectrum is not None else self.active_spectrum

    def ensure_baseline_analysis(self) -> None:
        """Ensure MS1 search, quantification, derived traits, and QC have been computed on active spectrum."""
        spectrum = self.get_or_load_spectrum()
        library = self.get_or_load_library()
        if self.ms1_search_result is None:
            self.ms1_search_result = search_ms1_library(
                spectrum=spectrum,
                library=library,
                settings=MS1SearchSettings(
                    tolerance=20.0,
                    adduct_codes=("[M+Na]+", "[M+H]+", "[M+K]+"),
                    charge_states=(1,),
                ),
            )
        if self.quantification_result is None:
            self.quantification_result = quantify_ms1_search_result(
                search_result=self.ms1_search_result,
            )
        if self.derived_traits_report is None:
            self.derived_traits_report = compute_derived_glycan_traits(self.quantification_result)
        if self.spectrum_qc_report is None:
            self.spectrum_qc_report = evaluate_spectrum_qc(
                spectrum=spectrum,
                search_result=self.ms1_search_result,
                quantification_result=self.quantification_result,
                calibration_result=self.calibration_result,
            )


session_state = AnalyticalSessionState()

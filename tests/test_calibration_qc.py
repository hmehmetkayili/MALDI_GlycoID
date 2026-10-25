"""Unit tests for Stage 8 Mass Calibration and Quality Control (`app/core/calibration.py`, `app/analysis/quality_control.py`)."""

from __future__ import annotations

from pathlib import Path
import unittest

from app.analysis.batch_processing import BatchSampleInput, execute_batch_analysis
from app.analysis.glycan_search import MS1SearchSettings, search_ms1_library
from app.analysis.quality_control import (
    QCStatus,
    QCThresholdSource,
    QCToleranceThresholds,
    evaluate_batch_qc,
    evaluate_spectrum_qc,
)
from app.analysis.quantification import quantify_ms1_search_result
from app.core.calibration import (
    CalibrationModelType,
    CalibrationSettings,
    calibrate_spectrum,
    get_default_nglycan_calibrants,
    get_default_peptide_calibrants,
)
from app.libraries.library_io import import_library_from_json
from app.parsers import parse_spectrum_file
from app.parsers.base_models import SpectrumData, SpectrumDataType, SpectrumPeak


class TestMassCalibrationAndQC(unittest.TestCase):
    """Test suite for Stage 8 mass calibration and MassyTools-inspired quality control."""

    @classmethod
    def setUpClass(cls) -> None:
        lib_path = Path("data/libraries/human_nglycan_native_v1.json")
        cls.library = import_library_from_json(lib_path)

    def test_linear_and_quadratic_mass_calibration_recovery(self) -> None:
        """Verify linear and quadratic calibration accurately correct systematic mass drift."""
        calibrants = get_default_nglycan_calibrants(adduct="[M+Na]+", labeling="NONE")
        self.assertGreaterEqual(len(calibrants), 6)

        # Also verify peptide calibrants compute valid m/z values
        pep_cals = get_default_peptide_calibrants()
        self.assertEqual(len(pep_cals), 6)

        # Simulate +65 ppm + 0.015 Da systematic drift across 6 calibrants
        shifted_peaks = []
        for idx, cal in enumerate(calibrants[:6]):
            obs_mz = cal.theoretical_mz * (1.0 + 65e-6) + 0.015
            shifted_peaks.append(SpectrumPeak(mz=obs_mz, intensity=5000.0 + idx * 1000.0))

        uncalibrated_spec = SpectrumData(
            sample_id="Shifted_Sample",
            filename="shifted.txt",
            source_format="TXT",
            peaks=tuple(shifted_peaks),
            data_type=SpectrumDataType.CENTROID,
            user_specified_data_type=None,
            valid_points_count=len(shifted_peaks),
            min_mz=min(p.mz for p in shifted_peaks),
            max_mz=max(p.mz for p in shifted_peaks),
            max_intensity=max(p.intensity for p in shifted_peaks),
        )

        # 1. Linear calibration
        lin_settings = CalibrationSettings(
            model_type=CalibrationModelType.LINEAR,
            match_tolerance=120.0,
        )
        lin_res = calibrate_spectrum(uncalibrated_spec, calibrants, lin_settings)
        self.assertTrue(lin_res.applied)
        self.assertEqual(len(lin_res.matched_calibrants), 6)
        self.assertGreater(lin_res.rms_ppm_before, 60.0)
        self.assertLess(lin_res.rms_ppm_after, 0.5)
        self.assertTrue(lin_res.calibrated_spectrum.metadata.get("is_calibrated"))

        # 2. Quadratic calibration with slight quadratic distortion
        quad_peaks = []
        for idx, cal in enumerate(calibrants[:6]):
            obs_mz = cal.theoretical_mz * (1.0 + 50e-6) + 2e-8 * (cal.theoretical_mz ** 2)
            quad_peaks.append(SpectrumPeak(mz=obs_mz, intensity=8000.0))
        quad_spec = SpectrumData(
            sample_id="Quad_Sample",
            filename="quad.txt",
            source_format="TXT",
            peaks=tuple(quad_peaks),
            data_type=SpectrumDataType.CENTROID,
            user_specified_data_type=None,
            valid_points_count=len(quad_peaks),
            min_mz=min(p.mz for p in quad_peaks),
            max_mz=max(p.mz for p in quad_peaks),
            max_intensity=max(p.intensity for p in quad_peaks),
        )
        quad_settings = CalibrationSettings(
            model_type=CalibrationModelType.QUADRATIC,
            match_tolerance=150.0,
        )
        quad_res = calibrate_spectrum(quad_spec, calibrants, quad_settings)
        self.assertTrue(quad_res.applied)
        self.assertEqual(len(quad_res.coefficients), 3)
        self.assertLess(quad_res.rms_ppm_after, 0.1)

    def test_insufficient_calibrants_safeguard_and_never_silently_calibrate(self) -> None:
        """Verify minimum calibrant count safeguard leaves spectrum uncalibrated with explicit warning."""
        calibrants = get_default_nglycan_calibrants(adduct="[M+Na]+", labeling="NONE")
        # Spectrum with only 1 matching calibrant peak
        single_peak_spec = SpectrumData(
            sample_id="Single_Cal_Sample",
            filename="single.txt",
            source_format="TXT",
            peaks=(SpectrumPeak(mz=calibrants[0].theoretical_mz + 0.01, intensity=5000.0),),
            data_type=SpectrumDataType.CENTROID,
            user_specified_data_type=None,
            valid_points_count=1,
            min_mz=calibrants[0].theoretical_mz + 0.01,
            max_mz=calibrants[0].theoretical_mz + 0.01,
            max_intensity=5000.0,
        )

        res = calibrate_spectrum(
            single_peak_spec,
            calibrants,
            CalibrationSettings(model_type=CalibrationModelType.LINEAR, require_calibration_success=False),
        )
        self.assertFalse(res.applied)
        self.assertEqual(res.calibrated_spectrum.peaks[0].mz, single_peak_spec.peaks[0].mz)
        self.assertTrue(
            any("Insufficient matched calibrants" in w for w in res.warnings),
            f"Expected insufficient calibrants warning in {res.warnings}",
        )

        # When require_calibration_success=True, must raise ValueError
        with self.assertRaises(ValueError):
            calibrate_spectrum(
                single_peak_spec,
                calibrants,
                CalibrationSettings(model_type=CalibrationModelType.LINEAR, require_calibration_success=True),
            )

    def test_spectrum_qc_threshold_sources_and_centroid_vs_profile_sn(self) -> None:
        """Verify SpectrumQCReport checks, explicit threshold sources, and S/N availability guard."""
        # 1. CSV example has S/N and area columns -> profile_metrics_available = True
        spec_csv = parse_spectrum_file(
            Path("data/examples/synthetic_native_nglycans.csv"),
            sample_id="CSV_QC",
        )
        ms1_csv = search_ms1_library(
            spec_csv,
            self.library,
            MS1SearchSettings(tolerance=25.0, adduct_codes=("[M+Na]+", "[M+K]+")),
        )
        quant_csv = quantify_ms1_search_result(ms1_csv)
        qc_csv = evaluate_spectrum_qc(spec_csv, ms1_csv, quant_csv)

        self.assertEqual(qc_csv.overall_status, QCStatus.PASS)
        self.assertTrue(qc_csv.profile_metrics_available)
        self.assertIsNotNone(qc_csv.median_matched_sn)
        self.assertGreater(qc_csv.median_matched_sn, 20.0)

        # Verify every check explicitly declares a valid QCThresholdSource
        sources = {c.threshold_source for c in qc_csv.criteria_checks}
        self.assertIn(QCThresholdSource.LITERATURE_DERIVED, sources)
        self.assertIn(QCThresholdSource.INSTRUMENT_SPECIFIC, sources)
        self.assertIn(QCThresholdSource.PROVISIONAL_DEFAULT, sources)

        # 2. XY centroid example lacks S/N -> profile_metrics_available = False, median_matched_sn = None
        spec_xy = parse_spectrum_file(
            Path("data/examples/synthetic_native_nglycans.xy"),
            sample_id="XY_QC",
        )
        ms1_xy = search_ms1_library(
            spec_xy,
            self.library,
            MS1SearchSettings(tolerance=25.0, adduct_codes=("[M+Na]+", "[M+K]+")),
        )
        qc_xy = evaluate_spectrum_qc(spec_xy, ms1_xy)
        self.assertFalse(qc_xy.profile_metrics_available)
        self.assertIsNone(qc_xy.median_matched_sn)

    def test_batch_qc_missingness_cv_and_outlier_flagging(self) -> None:
        """Verify BatchQCReport evaluates batch missingness, replicate CVs, and flags outlier samples."""
        txt_bytes = Path("data/examples/synthetic_native_nglycans_centroid.txt").read_bytes()
        xy_bytes = Path("data/examples/synthetic_native_nglycans.xy").read_bytes()
        # Noise-only spectrum with 0 glycan matches
        noise_text = "m/z\tintensity\n999.111\t5000.0\n1050.222\t4000.0\n"

        batch_samples = [
            BatchSampleInput(sample_id="Good_1", filename="g1.txt", content=txt_bytes, group="Control"),
            BatchSampleInput(sample_id="Good_2", filename="g2.xy", content=xy_bytes, group="Control"),
            BatchSampleInput(sample_id="Noise_Outlier", filename="noise.txt", content=noise_text, group="Control"),
        ]

        batch_res = execute_batch_analysis(
            samples=batch_samples,
            library=self.library,
            search_settings=MS1SearchSettings(tolerance=25.0, adduct_codes=("[M+Na]+", "[M+K]+")),
        )

        batch_qc = evaluate_batch_qc(
            batch_res,
            thresholds=QCToleranceThresholds(
                max_batch_missingness_source=QCThresholdSource.USER_DEFINED,
            ),
        )
        self.assertEqual(batch_qc.total_samples, 3)
        self.assertEqual(batch_qc.passed_samples_count, 2)
        self.assertEqual(batch_qc.failed_samples_count, 1)
        self.assertIn("Noise_Outlier", batch_qc.flagged_outlier_sample_ids)
        self.assertIsNotNone(batch_qc.median_replicate_cv_pct)
        self.assertLess(batch_qc.median_replicate_cv_pct, 5.0)


if __name__ == "__main__":
    unittest.main()

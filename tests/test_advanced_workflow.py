"""Comprehensive tests for speed optimizations, negative-mode library m/z, S/N & quality filtering,
calibration-before-search workflow, interactive idotp isotope comparison, optional MALDI corrections,
and complex realistic MALDI-TOF test spectra (positive & negative modes).
"""

from __future__ import annotations

import time
import unittest

from app.analysis.glycan_search import MS1SearchSettings, search_ms1_library
from app.analysis.spectrum_corrections import (
    CorrectionConfig,
    apply_corrections,
    list_available_corrections,
)
from app.analysis.synthetic_spectra import generate_complex_nglycan_spectrum
from app.api.compat import TestClient
from app.api.state import session_state
from app.core.calibration import (
    CalibrationModelType,
    CalibrationSettings,
    calibrate_spectrum,
    get_calibrants_for_panel,
    select_best_calibrant_panel,
)
from app.core.scientific_constants import IonPolarity
from app.main import app


class TestAdvancedMALDIWorkflow(unittest.TestCase):
    """End-to-end verification of all newly added capabilities and complex spectrum benchmarks."""

    def setUp(self) -> None:
        self.client = TestClient(app)
        session_state.reset_to_default_library()

    def test_library_table_and_exports_include_negative_mode_mz(self) -> None:
        """Verify default, generated, queried, and CSV-exported libraries include negative-mode m/z values."""
        resp_def = self.client.get("/api/libraries/default")
        self.assertEqual(resp_def.status_code, 200)
        data = resp_def.json()
        self.assertGreater(len(data["entries"]), 300)
        entry = next(e for e in data["entries"] if e["composition"] == "Hex5HexNAc2")
        self.assertIn("negative_mode_mz", entry)
        neg = entry["negative_mode_mz"]
        self.assertAlmostEqual(neg["[M-H]-"], 1233.42615, places=3)
        self.assertAlmostEqual(neg["[M+Cl]-"], 1269.40283, places=3)
        self.assertAlmostEqual(neg["[M-2H]2-"], 616.20944, places=3)

        # Verify custom generated library also includes negative mode m/z
        resp_gen = self.client.post(
            "/api/libraries/generate",
            json={
                "min_hex": 3,
                "max_hex": 6,
                "min_hexnac": 2,
                "max_hexnac": 4,
                "max_dhex": 1,
                "max_neu5ac": 1,
                "labeling": "2-AB",
            },
        )
        self.assertEqual(resp_gen.status_code, 200)
        gen_entry = resp_gen.json()["entries"][0]
        self.assertIn("[M-H]-", gen_entry["negative_mode_mz"])
        self.assertGreater(gen_entry["negative_mode_mz"]["[M-H]-"], 800.0)

        # Verify CSV export has negative-mode columns populated
        resp_csv = self.client.get("/api/libraries/export?format=csv")
        self.assertEqual(resp_csv.status_code, 200)
        header = resp_csv.text.splitlines()[0]
        self.assertIn("mz_M_minus_H_1neg", header)
        self.assertIn("mz_M_plus_Cl_1neg", header)
        self.assertIn("mz_M_minus_2H_2neg", header)

    def test_upload_snr_filter_and_optional_maldi_corrections(self) -> None:
        """Verify import-time S/N filtering and optional MALDI spectrum corrections catalogue & execution."""
        resp_cat = self.client.get("/api/spectra/corrections")
        self.assertEqual(resp_cat.status_code, 200)
        cat_keys = {c["key"] for c in resp_cat.json()["corrections"]}
        for expected in (
            "smoothing",
            "baseline_correction",
            "pick_peaks",
            "min_snr",
            "remove_matrix_peaks",
            "merge_close_peaks",
            "normalize",
            "calibration",
            "exclude_isotope_satellites",
        ):
            self.assertIn(expected, cat_keys)

        # Upload a centroid list with noise peaks and filter by min_snr at import time
        csv_lines = ["mz,intensity,snr"]
        csv_lines.append("1257.4226,10000.0,45.0")  # Hex5HexNAc2 [M+Na]+
        csv_lines.append("1258.4260,5100.0,22.0")   # M+1 isotope
        csv_lines.append("1485.5337,8000.0,35.0")   # Hex3HexNAc4dHex1 [M+Na]+
        csv_lines.append("1310.1234,90.0,1.8")      # low S/N noise peak
        csv_lines.append("1410.5678,70.0,2.1")      # low S/N noise peak
        csv_bytes = ("\n".join(csv_lines)).encode("utf-8")

        resp_up = self.client.post(
            "/api/spectra/upload",
            files={"file": ("snr_test.csv", csv_bytes, "text/csv")},
            data={"min_snr": "5.0"},
        )
        self.assertEqual(resp_up.status_code, 200)
        up_data = resp_up.json()
        self.assertEqual(up_data["valid_points_count"], 3)
        self.assertIsNotNone(up_data["corrections"])

    def test_complex_spectrum_calibration_first_and_search_quality_filters(self) -> None:
        """Verify calibration-before-search, quality filters, and speed on complex + and - mode spectra."""
        library = session_state.reset_to_default_library()

        for polarity in ("positive", "negative"):
            syn = generate_complex_nglycan_spectrum(polarity=polarity, seed=11)
            self.assertGreater(len(syn.spectrum.peaks), 100000)

            # 1. Apply optional MALDI profile corrections: smoothing + SNIP baseline + peak picking + matrix exclusion
            corr_res = apply_corrections(
                syn.spectrum,
                CorrectionConfig(
                    smoothing=True,
                    baseline_correction=True,
                    pick_peaks=True,
                    pick_min_snr=5.0,
                    remove_matrix_peaks=True,
                    matrix_polarity=polarity,
                    min_mz=850.0,
                ),
            )
            centroid_spec = corr_res.spectrum
            self.assertLess(len(centroid_spec.peaks), 4000)
            self.assertGreater(len(centroid_spec.peaks), 100)

            adducts = ("[M-H]-",) if polarity == "negative" else ("[M+Na]+", "[M+K]+")
            pol_enum = IonPolarity.NEGATIVE if polarity == "negative" else IonPolarity.POSITIVE
            search_cfg = MS1SearchSettings(
                tolerance=20.0,
                polarity=pol_enum,
                adduct_codes=adducts,
                min_snr=6.0,
                min_score=55.0,
                min_idotp=0.75,
                exclude_isotope_satellites=True,
            )

            # 2. Search BEFORE calibration (systematic +55 ppm drift causes poor recovery at 20 ppm tolerance)
            uncal_res = search_ms1_library(centroid_spec, library, search_cfg)

            # 3. Auto-detect calibrant panel and calibrate BEFORE MS1 search
            best_panel, counts = select_best_calibrant_panel(centroid_spec)
            expected_panel = "nglycan_neg" if polarity == "negative" else "nglycan_na"
            self.assertEqual(best_panel, expected_panel, f"Panel counts: {counts}")

            cal_res = calibrate_spectrum(
                centroid_spec,
                get_calibrants_for_panel(best_panel),
                CalibrationSettings(model_type=CalibrationModelType.LINEAR, match_tolerance=150.0),
            )
            self.assertTrue(cal_res.applied)
            self.assertGreaterEqual(len(cal_res.matched_calibrants), 6)
            self.assertGreater(cal_res.rms_ppm_before, 40.0)
            self.assertLess(cal_res.rms_ppm_after, 3.0)

            # 4. Fast MS1 search on the calibrated spectrum
            t0 = time.perf_counter()
            cal_search = search_ms1_library(cal_res.calibrated_spectrum, library, search_cfg)
            elapsed = time.perf_counter() - t0
            self.assertLess(elapsed, 0.5, f"Expected sub-0.5s MS1 search, took {elapsed:.3f}s")

            truth_comps = {t.composition for t in syn.truth}
            uncal_recovered = len(truth_comps & {m.glycan_composition for m in uncal_res.top_matches_per_peak})
            cal_recovered = len(truth_comps & {m.glycan_composition for m in cal_search.top_matches_per_peak})
            self.assertGreater(
                cal_recovered,
                uncal_recovered + 20,
                f"Calibration should dramatically improve recovery ({uncal_recovered} -> {cal_recovered})",
            )
            self.assertGreaterEqual(cal_recovered, 35)

    def test_pipeline_and_interactive_idotp_isotope_comparison_api(self) -> None:
        """Verify /api/spectra/example?format=complex_pos, /api/search/pipeline, and /api/search/isotope."""
        # Load complex positive profile spectrum via API
        resp_ex = self.client.get("/api/spectra/example?format=complex_pos")
        self.assertEqual(resp_ex.status_code, 200)
        self.assertFalse(resp_ex.json()["is_calibrated"])

        # Run full pipeline: corrections -> calibration -> MS1 search
        resp_pipe = self.client.post(
            "/api/search/pipeline",
            json={
                "corrections": {
                    "smoothing": True,
                    "baseline_correction": True,
                    "pick_peaks": True,
                    "pick_min_snr": 5.0,
                    "remove_matrix_peaks": True,
                    "min_mz": 850.0,
                },
                "calibration": {
                    "model_type": "linear",
                    "reference_set": "auto",
                    "match_tolerance_ppm": 150.0,
                },
                "search": {
                    "search_mode": "library",
                    "tolerance": 20.0,
                    "tolerance_unit": "ppm",
                    "polarity": "positive",
                    "adduct_codes": ["[M+Na]+", "[M+K]+"],
                    "min_snr": 5.0,
                    "min_score": 50.0,
                    "exclude_isotope_satellites": True,
                },
            },
        )
        self.assertEqual(resp_pipe.status_code, 200)
        p_data = resp_pipe.json()
        self.assertEqual(p_data["stages"], ["corrections", "calibration", "search"])
        self.assertTrue(p_data["calibration"]["applied"])
        self.assertTrue(p_data["search"]["spectrum_is_calibrated"])
        self.assertGreater(p_data["search"]["matched_peaks_count"], 35)

        # Pick a matched glycan with an available idotp and call /api/search/isotope
        top_matches = p_data["search"]["top_matches_per_peak"]
        with_idotp = [m for m in top_matches if m["isotope_score_idotp"] is not None]
        self.assertGreater(len(with_idotp), 20)
        sample_match = with_idotp[0]

        resp_iso = self.client.post(
            "/api/search/isotope",
            json={
                "library_id": sample_match["library_id"],
                "composition": sample_match["glycan_composition"],
                "adduct_code": sample_match["assigned_adduct"],
                "num_isotopes": 5,
            },
        )
        self.assertEqual(resp_iso.status_code, 200)
        iso_data = resp_iso.json()
        self.assertTrue(iso_data["available"])
        self.assertGreater(iso_data["idotp"], 0.85)
        self.assertEqual(len(iso_data["isotopologues"]), 5)
        self.assertEqual(iso_data["isotopologues"][0]["label"], "M+0")
        self.assertIsNotNone(iso_data["isotopologues"][0]["observed_mz"])
        self.assertGreater(len(iso_data["observed_trace"]), 0)

        # Verify calibration reset endpoint restores uncalibrated state
        resp_reset = self.client.post("/api/qc/calibrate/reset")
        self.assertEqual(resp_reset.status_code, 200)
        self.assertFalse(resp_reset.json()["is_calibrated"])


if __name__ == "__main__":
    unittest.main()

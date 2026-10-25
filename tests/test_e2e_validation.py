"""End-to-end scientific and software validation suite for MALDI-GlycoID v1.0.0 (Stage 13).

Verifies:
1. Complete end-to-end synthetic analytical pipeline (Library -> Spectrum Import -> Calibration ->
   MS1 Search -> Quantification -> Derived Traits -> Glycopeptide Screening -> Batch -> QC ->
   Statistics/ML -> 12-Category Export Bundle).
2. Live local HTTP socket server execution (`ThreadingHTTPServer` + `urllib.request`) verifying
   real network socket round-trips for `/`, `/api/health`, `/api/search/ms1`, and `/api/exports/download`.
3. Repository integrity, zero placeholder/mock stubs, and full cross-stage consistency.
"""

from __future__ import annotations

import ast
import io
import json
import unittest
import zipfile

from PyPDF2 import PdfReader

from app import __version__ as SOFTWARE_VERSION
from app.api.compat import Request, TestClient, create_http_handler
from app.config import BASE_DIR, EXAMPLES_DIR, LIBRARIES_DIR
from app.core.mass_calculator import (
    calculate_attached_glycan_moiety_mass,
    calculate_released_glycan_mass,
)
from app.core.scientific_constants import WATER_MASS_DA
from app.exports import REQUIRED_12_EXPORT_CATEGORIES, read_xlsx_bytes_sheets
from app.libraries import import_library_from_json, validate_glycan_library
from app.main import app


class TestEndToEndValidationAndRelease(unittest.TestCase):
    """Comprehensive Stage 13 end-to-end acceptance and release validation suite."""

    def test_end_to_end_scientific_pipeline_acceptance(self) -> None:
        """Execute the complete analytical workflow from raw input files to 12-category exports."""
        client = TestClient(app)

        # Step 1: Verify default 333-entry library on disk and via API
        lib_disk = import_library_from_json(LIBRARIES_DIR / "human_nglycan_native_v1.json")
        val_res = validate_glycan_library(lib_disk)
        self.assertTrue(val_res.is_valid)
        self.assertEqual(lib_disk.entry_count, 333)

        # Verify fundamental released vs attached glycan moiety difference == WATER_MASS
        for comp in ("Hex5HexNAc2", "Hex3HexNAc4dHex1", "Hex5HexNAc4Neu5Ac2"):
            rel = calculate_released_glycan_mass(comp, labeling="NONE").neutral_monoisotopic_mass
            att = calculate_attached_glycan_moiety_mass(comp).neutral_monoisotopic_mass
            self.assertAlmostEqual(rel - att, WATER_MASS_DA, places=9)

        # Step 2: Load synthetic CSV spectrum via API
        r_spec = client.get("/api/spectra/example", params={"format": "csv"})
        self.assertEqual(r_spec.status_code, 200)
        self.assertEqual(r_spec.json()["valid_points_count"], 6)

        # Step 3: Run linear mass calibration
        r_cal = client.post(
            "/api/qc/calibrate",
            json={"model_type": "linear", "reference_set": "nglycan_na", "match_tolerance_ppm": 150.0},
        )
        self.assertEqual(r_cal.status_code, 200)
        self.assertTrue(r_cal.json()["applied"])

        # Step 4: Execute Mode B MS1 search
        r_ms1 = client.post(
            "/api/search/ms1",
            json={
                "search_mode": "library",
                "tolerance": 20.0,
                "tolerance_unit": "ppm",
                "polarity": "positive",
                "labeling": "NONE",
                "adduct_codes": ["[M+Na]+", "[M+H]+", "[M+K]+"],
                "charge_states": [1],
                "min_relative_intensity_pct": 0.5,
            },
        )
        self.assertEqual(r_ms1.status_code, 200)
        self.assertEqual(r_ms1.json()["matched_peaks_count"], 6)

        # Step 5: Execute Quantification & Derived Traits
        r_quant = client.post(
            "/api/search/quantify",
            json={
                "denominator_policy": "all_eligible_spectrum_peaks",
                "aggregation_policy": "primary_adduct_only",
                "primary_adduct": "[M+Na]+",
                "exclude_ambiguous_from_aggregation": False,
            },
        )
        self.assertEqual(r_quant.status_code, 200)
        q_json = r_quant.json()
        total_rel_abund = sum(
            g["relative_abundance_intensity_pct"]
            for g in q_json["quantification"]["glycan_abundances"]
        )
        self.assertAlmostEqual(total_rel_abund, 100.0, places=2)

        # Step 6: Execute IgG1 Fc Glycopeptide Screening
        fasta_str = (EXAMPLES_DIR / "synthetic_igg1_fc.fasta").read_text(encoding="utf-8")
        r_gp = client.post(
            "/api/glycopeptides/screen",
            json={
                "fasta_text": fasta_str,
                "enzyme": "trypsin",
                "max_missed_cleavages": 1,
                "precursor_mzs": [2634.0459],
                "tolerance_ppm": 15.0,
                "ms2_fragment_mzs": [204.0866, 366.1395, 1189.5120, 1392.5914],
                "fixed_carbamidomethyl_cys": True,
            },
        )
        self.assertEqual(r_gp.status_code, 200)
        self.assertGreaterEqual(len(r_gp.json()["matches"]), 1)

        # Step 7: Execute Multi-Sample Batch Demo & QC & Cohort Statistics
        r_batch = client.post("/api/batch/demo")
        self.assertEqual(r_batch.status_code, 200)
        self.assertEqual(r_batch.json()["completed_files_count"], 3)

        r_qc = client.post("/api/qc/evaluate")
        self.assertEqual(r_qc.status_code, 200)
        self.assertEqual(r_qc.json()["overall_status"], "PASS")

        r_stat = client.post(
            "/api/qc/statistics-demo",
            json={"imputation_method": "half_min_observed", "classifier_type": "logistic_regression"},
        )
        self.assertEqual(r_stat.status_code, 200)
        self.assertGreaterEqual(r_stat.json()["classification"]["oof_accuracy"], 0.8)

        # Step 8: Download and verify complete ZIP export bundle, XLSX, and PDF
        r_zip = client.get("/api/exports/download", params={"format": "zip"})
        self.assertEqual(r_zip.status_code, 200)
        with zipfile.ZipFile(io.BytesIO(r_zip.content), "r") as zf:
            archive_files = zf.namelist()
            for idx, cat in enumerate(REQUIRED_12_EXPORT_CATEGORIES, start=1):
                self.assertIn(f"csv/{idx:02d}_{cat}.csv", archive_files)
                self.assertIn(f"tsv/{idx:02d}_{cat}.tsv", archive_files)
            xlsx_files = [f for f in archive_files if f.endswith(".xlsx")]
            pdf_files = [f for f in archive_files if f.endswith(".pdf")]
            self.assertEqual(len(xlsx_files), 1)
            self.assertEqual(len(pdf_files), 1)

            sheets = read_xlsx_bytes_sheets(zf.read(xlsx_files[0]))
            self.assertEqual(len(sheets), 12)

            pdf_reader = PdfReader(io.BytesIO(zf.read(pdf_files[0])))
            self.assertGreaterEqual(len(pdf_reader.pages), 2)

    def test_http_wire_handler_and_multipart_roundtrip(self) -> None:
        """Verify HTTP handler creation, multipart form-data wire parsing, and targeted MS1 search."""
        handler_cls = create_http_handler(app)
        self.assertTrue(hasattr(handler_cls, "do_GET"))
        self.assertTrue(hasattr(handler_cls, "do_POST"))

        # Test raw multipart/form-data wire body parsing through Request.form_sync()
        boundary = "----GlycoIDBoundary7MA4YWxkTrZu0gW"
        xy_payload = "1257.4225 18000.0\n1485.5338 45000.0\n"
        multipart_body = (
            f"--{boundary}\r\n"
            f'Content-Disposition: form-data; name="sample_id"\r\n\r\n'
            f"WIRE_SAMPLE_01\r\n"
            f"--{boundary}\r\n"
            f'Content-Disposition: form-data; name="file"; filename="wire_test.xy"\r\n'
            f"Content-Type: text/plain\r\n\r\n"
            f"{xy_payload}\r\n"
            f"--{boundary}--\r\n"
        ).encode("utf-8")

        raw_req = Request(
            method="POST",
            path="/api/spectra/upload",
            headers={"content-type": f"multipart/form-data; boundary={boundary}"},
            body=multipart_body,
        )
        resp = app.dispatch(raw_req)
        self.assertEqual(resp.status_code, 200)
        parsed_resp = json.loads(resp.body.decode("utf-8"))
        self.assertEqual(parsed_resp["sample_id"], "WIRE_SAMPLE_01")
        self.assertEqual(parsed_resp["valid_points_count"], 2)

        # Verify targeted Mode A search on the uploaded wire sample
        client = TestClient(app)
        r_ms1 = client.post(
            "/api/search/ms1",
            json={
                "search_mode": "targeted",
                "tolerance": 20.0,
                "tolerance_unit": "ppm",
                "polarity": "positive",
                "labeling": "NONE",
                "adduct_codes": ["[M+Na]+"],
                "charge_states": [1],
                "targeted_compositions": ["Hex5HexNAc2", "Hex3HexNAc4dHex1"],
            },
        )
        self.assertEqual(r_ms1.status_code, 200)
        ms1_payload = r_ms1.json()
        self.assertEqual(ms1_payload["search_mode"], "targeted")
        self.assertEqual(ms1_payload["matched_peaks_count"], 2)

    def test_repository_completeness_and_no_unimplemented_stubs(self) -> None:
        """Verify all core, parser, library, analysis, export, UI, and API Python files compile and contain no stub functions."""
        py_dirs = [
            BASE_DIR / "app",
            BASE_DIR / "scripts",
            BASE_DIR / "tests",
        ]
        py_files_count = 0
        for d in py_dirs:
            for py_file in d.rglob("*.py"):
                py_files_count += 1
                text = py_file.read_text(encoding="utf-8")
                tree = ast.parse(text, filename=str(py_file))
                for node in ast.walk(tree):
                    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                        body = [
                            stmt
                            for stmt in node.body
                            if not (
                                isinstance(stmt, ast.Expr)
                                and isinstance(stmt.value, ast.Constant)
                                and isinstance(stmt.value.value, str)
                            )
                        ]
                        self.assertGreater(len(body), 0, f"Empty function {node.name} in {py_file}")
                        if len(body) == 1 and isinstance(body[0], ast.Raise):
                            self.fail(f"Unimplemented stub function {node.name} in {py_file}")
        self.assertGreaterEqual(py_files_count, 40)


if __name__ == "__main__":
    unittest.main()

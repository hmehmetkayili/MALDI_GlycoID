"""Unit and integration tests for Stage 12: REST API, security middleware, and deployment configs."""

from __future__ import annotations

import io
import unittest
import zipfile
from pathlib import Path

from PyPDF2 import PdfReader

from app.api.compat import TestClient
from app.config import BASE_DIR
from app.exports import REQUIRED_12_EXPORT_CATEGORIES, read_xlsx_bytes_sheets
from app.main import app
from scripts.run_local import main as run_local_main


class TestStage12APIAndDeployment(unittest.TestCase):
    """Comprehensive verification of all Stage 12 REST API endpoints, security guards, and deployment files."""

    def setUp(self) -> None:
        self.client = TestClient(app)

    def test_health_index_and_static_serving_with_security_headers(self) -> None:
        """Verify /api/health, / HTML template, /static/* assets, security headers, and static path traversal guard."""
        resp_health = self.client.get("/api/health")
        self.assertEqual(resp_health.status_code, 200)
        h_data = resp_health.json()
        self.assertEqual(h_data["status"], "ok")
        self.assertEqual(h_data["app_name"], "MALDI-GlycoID")
        self.assertEqual(h_data["active_library_entries_count"], 333)

        # Verify security headers on API responses
        self.assertEqual(resp_health.headers.get("x-content-type-options"), "nosniff")
        self.assertEqual(resp_health.headers.get("x-frame-options"), "DENY")
        self.assertIn("default-src 'self'", resp_health.headers.get("content-security-policy", ""))

        # Verify root HTML page
        resp_root = self.client.get("/")
        self.assertEqual(resp_root.status_code, 200)
        self.assertIn("MALDI-GlycoID", resp_root.text)
        self.assertIn("view-dashboard", resp_root.text)

        # Verify static CSS & JS assets
        resp_css = self.client.get("/static/css/styles.css")
        self.assertEqual(resp_css.status_code, 200)
        self.assertIn("--bg-body", resp_css.text)

        resp_js = self.client.get("/static/js/app.js")
        self.assertEqual(resp_js.status_code, 200)
        self.assertIn("/api/search/ms1", resp_js.text)

        # Verify static path traversal attempt is rejected
        resp_traversal = self.client.get("/static/../app/config.py")
        self.assertEqual(resp_traversal.status_code, 400)

    def test_libraries_and_spectra_api_endpoints_and_upload_security(self) -> None:
        """Verify library generation/query/export/import and spectrum upload/security validation."""
        # 1. Default library
        resp_def = self.client.get("/api/libraries/default")
        self.assertEqual(resp_def.status_code, 200)
        self.assertEqual(resp_def.json()["entry_count"], 333)

        # 2. Custom library generation
        resp_gen = self.client.post(
            "/api/libraries/generate",
            json={
                "min_hex": 3,
                "max_hex": 6,
                "min_hexnac": 2,
                "max_hexnac": 4,
                "max_dhex": 1,
                "max_neu5ac": 1,
                "labeling": "NONE",
                "include_neu5gc": False,
                "expand_sialic_linkages": False,
            },
        )
        self.assertEqual(resp_gen.status_code, 200)
        gen_data = resp_gen.json()
        self.assertTrue(gen_data["validation_valid"])
        self.assertGreater(gen_data["entry_count"], 10)

        # 3. Query library
        resp_q = self.client.post(
            "/api/libraries/query",
            json={"composition_query": "Hex5HexNAc2", "glycan_class": "high_mannose"},
        )
        self.assertEqual(resp_q.status_code, 200)
        self.assertEqual(resp_q.json()["matched_count"], 1)
        self.assertEqual(resp_q.json()["entries"][0]["composition"], "Hex5HexNAc2")

        # 4. Export library to CSV and re-import
        resp_exp = self.client.get("/api/libraries/export", params={"format": "csv"})
        self.assertEqual(resp_exp.status_code, 200)
        csv_bytes = resp_exp.content
        self.assertIn(b"Hex5HexNAc2", csv_bytes)

        resp_imp = self.client.post(
            "/api/libraries/import",
            files={"file": ("custom_lib.csv", csv_bytes, "text/csv")},
        )
        self.assertEqual(resp_imp.status_code, 200)
        self.assertEqual(resp_imp.json()["entry_count"], gen_data["entry_count"])

        # Reset back to default 333-entry library
        self.client.get("/api/libraries/default")

        # 5. Load synthetic example spectrum
        resp_ex = self.client.get("/api/spectra/example", params={"format": "csv"})
        self.assertEqual(resp_ex.status_code, 200)
        self.assertEqual(resp_ex.json()["valid_points_count"], 6)

        # 6. Upload valid XY spectrum
        xy_text = "1257.4225 15000.0\n1485.5338 42000.0\n1647.5866 38000.0\n"
        resp_up = self.client.post(
            "/api/spectra/upload",
            data={"sample_id": "UPLOAD_TEST_01"},
            files={"file": ("test_spec.xy", xy_text.encode("utf-8"), "text/plain")},
        )
        self.assertEqual(resp_up.status_code, 200)
        self.assertEqual(resp_up.json()["sample_id"], "UPLOAD_TEST_01")
        self.assertEqual(resp_up.json()["valid_points_count"], 3)

        # 7. Security: reject path traversal in filename
        resp_bad_path = self.client.post(
            "/api/spectra/upload",
            files={"file": ("../../etc/passwd.xy", xy_text.encode("utf-8"), "text/plain")},
        )
        self.assertEqual(resp_bad_path.status_code, 400)
        self.assertIn("path traversal", resp_bad_path.json()["detail"].lower())

        # 8. Security: reject unsupported extension & proprietary Bruker .baf container
        resp_bad_ext = self.client.post(
            "/api/spectra/upload",
            files={"file": ("malware.exe", b"MZ1234", "application/octet-stream")},
        )
        self.assertEqual(resp_bad_ext.status_code, 400)
        self.assertIn("Unsupported file format", resp_bad_ext.json()["detail"])

        resp_baf = self.client.post(
            "/api/spectra/upload",
            files={"file": ("analysis.baf", b"\x00\x01\x02\x03", "application/octet-stream")},
        )
        self.assertEqual(resp_baf.status_code, 400)
        self.assertIn("msconvert", resp_baf.json()["detail"])

    def test_search_quantify_glycopeptides_batch_qc_statistics_and_exports_api(self) -> None:
        """Verify full workflow across MS1 search, quantification, glycopeptides, batch, QC, stats, and exports."""
        # Load synthetic CSV spectrum
        self.client.get("/api/spectra/example", params={"format": "csv"})

        # 1. Mode B (Library) MS1 Search
        resp_ms1 = self.client.post(
            "/api/search/ms1",
            json={
                "search_mode": "library",
                "tolerance": 25.0,
                "tolerance_unit": "ppm",
                "polarity": "positive",
                "labeling": "NONE",
                "adduct_codes": ["[M+Na]+", "[M+H]+", "[M+K]+"],
                "charge_states": [1],
                "min_relative_intensity_pct": 0.5,
            },
        )
        self.assertEqual(resp_ms1.status_code, 200)
        ms1_data = resp_ms1.json()
        self.assertEqual(ms1_data["matched_peaks_count"], 6)

        # 2. Quantification & derived traits
        resp_quant = self.client.post(
            "/api/search/quantify",
            json={
                "denominator_policy": "all_eligible_spectrum_peaks",
                "aggregation_policy": "primary_adduct_only",
                "primary_adduct": "[M+Na]+",
                "exclude_ambiguous_from_aggregation": False,
            },
        )
        self.assertEqual(resp_quant.status_code, 200)
        q_data = resp_quant.json()
        self.assertGreaterEqual(len(q_data["quantification"]["glycan_abundances"]), 5)
        self.assertGreaterEqual(len(q_data["derived_traits"]), 5)

        # 3. Glycopeptide screening with IgG1 EEQYNSTYR + Hex3HexNAc4dHex1 ([M+H]+ = 2634.0459) and MS2 ions
        resp_gp = self.client.post(
            "/api/glycopeptides/screen",
            json={
                "fasta_text": ">sp|P01857|IGHG1_HUMAN\nTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYK",
                "enzyme": "trypsin",
                "max_missed_cleavages": 1,
                "precursor_mzs": [2634.0459],
                "tolerance_ppm": 15.0,
                "ms2_fragment_mzs": [204.0866, 366.1395, 1189.5120, 1392.5914],
                "fixed_carbamidomethyl_cys": True,
            },
        )
        self.assertEqual(resp_gp.status_code, 200)
        gp_data = resp_gp.json()
        self.assertGreaterEqual(len(gp_data["matches"]), 1)
        top_gp = gp_data["matches"][0]
        self.assertEqual(top_gp["peptide_sequence"], "EEQYNSTYR")
        self.assertEqual(top_gp["glycan_composition"], "Hex3HexNAc4dHex1")
        self.assertIn(top_gp["confidence_category"], ("MSMS_SUPPORTED", "SITE_LOCALIZED"))

        # 4. Demo Batch Processing
        resp_batch = self.client.post("/api/batch/demo")
        self.assertEqual(resp_batch.status_code, 200)
        b_data = resp_batch.json()
        self.assertEqual(b_data["total_files"], 3)
        self.assertEqual(b_data["completed_files_count"], 3)

        # 5. Calibration & QC
        resp_cal = self.client.post(
            "/api/qc/calibrate",
            json={"model_type": "linear", "reference_set": "nglycan_na", "match_tolerance_ppm": 150.0},
        )
        self.assertEqual(resp_cal.status_code, 200)
        self.assertTrue(resp_cal.json()["applied"])

        resp_qc = self.client.post("/api/qc/evaluate")
        self.assertEqual(resp_qc.status_code, 200)
        self.assertIn(resp_qc.json()["overall_status"], ("PASS", "WARN", "FAIL"))
        self.assertGreaterEqual(len(resp_qc.json()["criteria_checks"]), 4)

        # 6. Cohort Statistics & Leakage-Free ML Demo
        resp_stat = self.client.post(
            "/api/qc/statistics-demo",
            json={"imputation_method": "half_min_observed", "classifier_type": "logistic_regression"},
        )
        self.assertEqual(resp_stat.status_code, 200)
        stat_data = resp_stat.json()
        self.assertGreaterEqual(stat_data["classification"]["oof_accuracy"], 0.8)
        self.assertTrue(len(stat_data["pca_plot_base64"]) > 100)

        # 7. Export Downloads (ZIP, XLSX, PDF, JSON, CSV)
        resp_zip = self.client.get("/api/exports/download", params={"format": "zip"})
        self.assertEqual(resp_zip.status_code, 200)
        with zipfile.ZipFile(io.BytesIO(resp_zip.content), "r") as zf:
            names = zf.namelist()
            for idx, cat in enumerate(REQUIRED_12_EXPORT_CATEGORIES, start=1):
                self.assertIn(f"csv/{idx:02d}_{cat}.csv", names)

        resp_xlsx = self.client.get("/api/exports/download", params={"format": "xlsx"})
        self.assertEqual(resp_xlsx.status_code, 200)
        sheets = read_xlsx_bytes_sheets(resp_xlsx.content)
        self.assertEqual(len(sheets), 12)

        resp_pdf = self.client.get("/api/exports/download", params={"format": "pdf"})
        self.assertEqual(resp_pdf.status_code, 200)
        pdf_reader = PdfReader(io.BytesIO(resp_pdf.content))
        self.assertGreaterEqual(len(pdf_reader.pages), 2)

        resp_csv = self.client.get(
            "/api/exports/download",
            params={"format": "csv", "category": "annotated_peaks"},
        )
        self.assertEqual(resp_csv.status_code, 200)
        self.assertIn("glycan_composition", resp_csv.text)

    def test_deployment_configs_and_run_local_cli(self) -> None:
        """Verify Dockerfile, docker-compose.yml, CI workflow, docs, and scripts/run_local.py --check-only."""
        required_files = [
            "requirements.txt",
            ".env.example",
            ".gitignore",
            ".dockerignore",
            "Dockerfile",
            "docker-compose.yml",
            "render.yaml",
            "Procfile",
            "runtime.txt",
            ".github/workflows/tests.yml",
            "LICENSE",
            "README.md",
            "docs/INSTALLATION.md",
            "docs/USER_GUIDE.md",
            "docs/SCIENTIFIC_METHODS.md",
            "docs/LIBRARY_FORMAT.md",
            "docs/VALIDATION.md",
            "docs/DEPLOYMENT.md",
        ]
        for rel_path in required_files:
            full_path = BASE_DIR / rel_path
            self.assertTrue(full_path.is_file(), f"Missing deployment/doc file: {rel_path}")
            self.assertGreater(full_path.stat().st_size, 10, f"Empty file: {rel_path}")

        # Verify TXT visual summary enrichment (plot_peaks, top_peaks, base_peak_mz, raw_text_lines)
        resp_txt = self.client.get("/api/spectra/example", params={"format": "txt"})
        self.assertEqual(resp_txt.status_code, 200)
        txt_summary = resp_txt.json()
        self.assertIn("plot_peaks", txt_summary)
        self.assertIn("top_peaks", txt_summary)
        self.assertIn("base_peak_mz", txt_summary)
        self.assertIn("raw_text_lines", txt_summary)
        self.assertGreater(len(txt_summary["raw_text_lines"]), 3)
        self.assertIn("relative_intensity_pct", txt_summary["preview"][0])

        # Verify CLI check-only mode
        import sys
        old_argv = sys.argv
        try:
            sys.argv = ["run_local.py", "--check-only"]
            exit_code = run_local_main()
            self.assertEqual(exit_code, 0)
        finally:
            sys.argv = old_argv


if __name__ == "__main__":
    unittest.main()

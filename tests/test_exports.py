"""Unit tests for Stage 10 Complete Export and Reporting Center (`app/exports/`)."""

from __future__ import annotations

import csv
import io
import json
from pathlib import Path
import unittest
import zipfile

from PyPDF2 import PdfReader

from app.analysis.glycan_search import MS1SearchSettings, search_ms1_library
from app.analysis.glycopeptide_search import PrecursorObservation, screen_glycopeptides
from app.analysis.quality_control import evaluate_spectrum_qc
from app.analysis.quantification import quantify_ms1_search_result, summarize_replicates_and_build_matrix
from app.analysis.statistics import compute_descriptive_statistics, run_pca_analysis
from app.core.calibration import calibrate_spectrum, get_default_nglycan_calibrants
from app.exports import (
    REQUIRED_12_EXPORT_CATEGORIES,
    AnalysisExportBundleInput,
    build_standardized_export_tables,
    export_analysis_to_json_bytes,
    export_rows_to_csv_bytes,
    export_rows_to_tsv_bytes,
    export_summary_report_to_pdf_bytes,
    export_tables_to_xlsx_bytes,
    generate_complete_export_zip_bundle,
    read_xlsx_bytes_sheets,
)
from app.libraries.library_io import import_library_from_json
from app.parsers import parse_spectrum_file


class TestExportsAndReportingCenter(unittest.TestCase):
    """Test suite for Stage 10 CSV, TSV, multi-sheet XLSX, JSON, PDF, and ZIP bundle exports."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.library = import_library_from_json(Path("data/libraries/human_nglycan_native_v1.json"))
        cls.spec1 = parse_spectrum_file(
            Path("data/examples/synthetic_native_nglycans.csv"),
            sample_id="Sample_A",
        )
        cls.spec2 = parse_spectrum_file(
            Path("data/examples/synthetic_native_nglycans.xy"),
            sample_id="Sample_B",
        )
        cls.cal_res = calibrate_spectrum(
            cls.spec1,
            get_default_nglycan_calibrants(adduct="[M+Na]+"),
        )
        cls.ms1_res1 = search_ms1_library(
            cls.spec1,
            cls.library,
            MS1SearchSettings(tolerance=25.0, adduct_codes=("[M+Na]+", "[M+K]+")),
        )
        cls.ms1_res2 = search_ms1_library(
            cls.spec2,
            cls.library,
            MS1SearchSettings(tolerance=25.0, adduct_codes=("[M+Na]+", "[M+K]+")),
        )
        cls.quant1 = quantify_ms1_search_result(cls.ms1_res1)
        cls.quant2 = quantify_ms1_search_result(cls.ms1_res2)
        cls.matrix = summarize_replicates_and_build_matrix([cls.quant1, cls.quant2])
        cls.qc1 = evaluate_spectrum_qc(cls.spec1, cls.ms1_res1, cls.quant1, cls.cal_res)

        # Glycopeptide screening example (EEQYNSTYR + Hex3HexNAc4dHex1 [M+2H]2+ = 1317.5274)
        cls.gp_res = screen_glycopeptides(
            precursors=[PrecursorObservation(sample_id="Sample_A", observed_mz=1317.5274, intensity=9000.0)],
            fasta_or_proteins=Path("data/examples/synthetic_igg1_fc.fasta"),
            glycan_compositions=["Hex3HexNAc4dHex1"],
        )

        cls.desc_stats = compute_descriptive_statistics(cls.matrix)
        cls.pca_res = run_pca_analysis(cls.matrix)

        cls.bundle_input = AnalysisExportBundleInput(
            sample_id="Sample_A",
            spectrum=cls.spec1,
            ms1_search_result=cls.ms1_res1,
            quantification_result=cls.quant1,
            abundance_matrix=cls.matrix,
            glycopeptide_result=cls.gp_res,
            calibration_result=cls.cal_res,
            spectrum_qc_report=cls.qc1,
            descriptive_stats=cls.desc_stats,
            pca_result=cls.pca_res,
            reference_library=cls.library,
        )
        cls.tables = build_standardized_export_tables(cls.bundle_input)

    def test_all_12_export_categories_built_and_non_empty(self) -> None:
        """Verify all 12 standardized export categories are constructed and populated."""
        self.assertEqual(len(REQUIRED_12_EXPORT_CATEGORIES), 12)
        for cat in REQUIRED_12_EXPORT_CATEGORIES:
            self.assertIn(cat, self.tables)
            self.assertGreater(
                len(self.tables[cat]),
                0,
                f"Expected category '{cat}' to be non-empty in complete export bundle",
            )

    def test_csv_and_tsv_roundtrip_consistency(self) -> None:
        """Verify CSV and TSV exporters preserve headers and rows cleanly."""
        quant_rows = self.tables["quantified_glycans"]
        csv_bytes = export_rows_to_csv_bytes(quant_rows)
        tsv_bytes = export_rows_to_tsv_bytes(quant_rows)

        csv_parsed = list(csv.DictReader(io.StringIO(csv_bytes.decode("utf-8"))))
        tsv_parsed = list(csv.DictReader(io.StringIO(tsv_bytes.decode("utf-8")), delimiter="\t"))

        self.assertEqual(len(csv_parsed), len(quant_rows))
        self.assertEqual(len(tsv_parsed), len(quant_rows))
        self.assertEqual(csv_parsed[0]["glycan_composition"], quant_rows[0]["glycan_composition"])
        self.assertEqual(tsv_parsed[0]["glycan_composition"], quant_rows[0]["glycan_composition"])

    def test_multisheet_xlsx_export_and_sheet_roundtrip(self) -> None:
        """Verify multi-sheet .xlsx workbook generation and round-trip sheet parsing."""
        xlsx_bytes = export_tables_to_xlsx_bytes(self.tables)
        self.assertTrue(xlsx_bytes.startswith(b"PK\x03\x04"))

        parsed_sheets = read_xlsx_bytes_sheets(xlsx_bytes)
        self.assertEqual(len(parsed_sheets), 12)
        self.assertIn("annotated_peaks", parsed_sheets)
        self.assertIn("quantified_glycans", parsed_sheets)
        self.assertEqual(
            len(parsed_sheets["quantified_glycans"]),
            len(self.tables["quantified_glycans"]),
        )

    def test_json_pdf_and_complete_zip_bundle_verification(self) -> None:
        """Verify JSON payload, multi-page PDF report via PyPDF2, and complete ZIP export archive."""
        # 1. JSON verification
        json_bytes = export_analysis_to_json_bytes(self.tables)
        payload = json.loads(json_bytes.decode("utf-8"))
        self.assertEqual(payload["schema_version"], "1.0.0")
        self.assertEqual(len(payload["scientific_disclaimers"]), 3)
        self.assertEqual(len(payload["tables"]), 12)

        # 2. PDF verification using PyPDF2
        pdf_bytes = export_summary_report_to_pdf_bytes(
            standardized_tables=self.tables,
            spectrum=self.spec1,
            search_result=self.ms1_res1,
            quantification_result=self.quant1,
            qc_report=self.qc1,
        )
        self.assertTrue(pdf_bytes.startswith(b"%PDF-"))
        pdf_reader = PdfReader(io.BytesIO(pdf_bytes))
        self.assertEqual(len(pdf_reader.pages), 2)

        # 3. Complete ZIP bundle verification
        zip_bytes = generate_complete_export_zip_bundle(self.bundle_input)
        with zipfile.ZipFile(io.BytesIO(zip_bytes), "r") as zf:
            names = set(zf.namelist())
            self.assertIn("csv/01_annotated_peaks.csv", names)
            self.assertIn("tsv/03_quantified_glycans.tsv", names)
            self.assertIn("MALDI_GlycoID_Sample_A_Workbook.xlsx", names)
            self.assertIn("MALDI_GlycoID_Sample_A_Results.json", names)
            self.assertIn("MALDI_GlycoID_Sample_A_Summary_Report.pdf", names)
            self.assertIn("plots/annotated_ms1_spectrum.png", names)
            self.assertIn("plots/pca_scores_plot.png", names)


if __name__ == "__main__":
    unittest.main()

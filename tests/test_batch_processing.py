"""Unit tests for Stage 7: Multi-sample batch processing, fault isolation, ZIP safety, cancellation, and retry."""

from __future__ import annotations

import io
import unittest
import zipfile

from app.analysis import (
    BatchCancellationToken,
    BatchItemStatus,
    BatchJobStatus,
    BatchSampleInput,
    MS1SearchSettings,
    execute_batch_analysis,
    extract_batch_zip_safely,
    retry_failed_batch_items,
)
from app.config import EXAMPLES_DIR
from app.libraries import CompositionConstraints, generate_human_n_glycan_library


class TestBatchProcessing(unittest.TestCase):
    """Comprehensive verification of Stage 7 batch processing, fault isolation, and reproducibility."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.library = generate_human_n_glycan_library(
            constraints=CompositionConstraints(
                min_hex=3,
                max_hex=9,
                min_hexnac=2,
                max_hexnac=6,
                min_dhex=0,
                max_dhex=2,
                min_neu5ac=0,
                max_neu5ac=3,
                include_neu5gc=False,
                labeling="NONE",
            ),
            library_name="Stage 7 Batch Test Library",
            library_version="1.0.0",
        )

    def test_mixed_valid_and_invalid_batch_fault_isolation_and_retry(self) -> None:
        """Verify failed samples do not corrupt valid samples and can be retried cleanly."""
        samples = [
            BatchSampleInput(
                sample_id="Sample_A",
                filename="sample_a.txt",
                content="m/z intensity\n1257.4226 5000\n1485.5336 8000\n",
                group="Control",
                replicate_id="Rep1",
                metadata={"cohort": "Healthy"},
            ),
            BatchSampleInput(
                sample_id="Sample_B_Corrupt",
                filename="sample_b_corrupt.txt",
                content="# Empty or corrupted file with no valid peaks\n",
                group="Control",
                replicate_id="Rep2",
                metadata={"cohort": "Healthy"},
            ),
            BatchSampleInput(
                sample_id="Sample_C",
                filename="sample_c.xy",
                content="1257.4226 6000\n1809.6393 12000\n",
                group="Case",
                replicate_id="Rep1",
                metadata={"cohort": "Disease"},
            ),
        ]

        res = execute_batch_analysis(
            samples=samples,
            library=self.library,
            search_settings=MS1SearchSettings(tolerance=15.0, adduct_codes=("[M+Na]+",)),
        )

        self.assertEqual(res.status, BatchJobStatus.COMPLETED_WITH_FAILURES)
        self.assertEqual(res.total_files, 3)
        self.assertEqual(res.completed_files_count, 2)
        self.assertEqual(res.failed_files_count, 1)
        self.assertIn("Batch item failed; other samples completed successfully", res.warnings)
        self.assertEqual(len(res.failed_files_report), 1)
        self.assertEqual(res.failed_files_report[0]["sample_id"], "Sample_B_Corrupt")

        # Verify Sample_A and Sample_C results and sample identities are intact
        self.assertIsNotNone(res.abundance_matrix)
        assert res.abundance_matrix is not None
        self.assertEqual(res.abundance_matrix.sample_ids, ("Sample_A", "Sample_C"))
        self.assertEqual(len(res.combined_annotated_peaks), 4)  # 2 peaks in Sample_A + 2 peaks in Sample_C

        # Now retry the failed item with corrected spectrum data
        retried_res = retry_failed_batch_items(
            previous_result=res,
            corrected_inputs={
                "Sample_B_Corrupt": BatchSampleInput(
                    sample_id="Sample_B_Corrupt",
                    filename="sample_b_fixed.txt",
                    content="m/z intensity\n1257.4226 5500\n1485.5336 7500\n",
                    group="Control",
                    replicate_id="Rep2",
                )
            },
            library=self.library,
            search_settings=MS1SearchSettings(tolerance=15.0, adduct_codes=("[M+Na]+",)),
        )
        self.assertEqual(retried_res.status, BatchJobStatus.COMPLETED)
        self.assertEqual(retried_res.completed_files_count, 3)
        self.assertEqual(retried_res.failed_files_count, 0)
        assert retried_res.abundance_matrix is not None
        self.assertEqual(len(retried_res.abundance_matrix.sample_ids), 3)

    def test_safe_zip_extraction_and_path_traversal_rejection(self) -> None:
        """Verify safe in-memory ZIP extraction and strict rejection of path traversal and nested ZIPs."""
        # 1. Valid ZIP containing 2 spectra
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as zf:
            zf.writestr("sub/spec1.txt", "m/z intensity\n1257.4226 5000\n")
            zf.writestr("sub/spec2.xy", "1485.5336 8000\n")
            zf.writestr("sub/readme.md", "Ignored non-spectrum file")
        extracted = extract_batch_zip_safely(buf.getvalue())
        self.assertEqual(len(extracted), 2)
        self.assertEqual({s.sample_id for s in extracted}, {"spec1", "spec2"})

        # 2. Malicious path-traversal ZIP ('../../etc/passwd.txt') -> MUST raise ValueError
        bad_buf = io.BytesIO()
        with zipfile.ZipFile(bad_buf, "w") as zf:
            zf.writestr("../../evil.txt", "1257.4226 5000\n")
        with self.assertRaises(ValueError) as ctx1:
            extract_batch_zip_safely(bad_buf.getvalue())
        self.assertIn("Unsafe path traversal", str(ctx1.exception))

        # 3. Nested archive inside ZIP -> MUST raise ValueError
        nested_buf = io.BytesIO()
        with zipfile.ZipFile(nested_buf, "w") as zf:
            zf.writestr("inner.zip", b"PK\x03\x04")
        with self.assertRaises(ValueError) as ctx2:
            extract_batch_zip_safely(nested_buf.getvalue())
        self.assertIn("Nested archives are not permitted", str(ctx2.exception))

    def test_batch_cancellation_and_reproducibility_metadata(self) -> None:
        """Verify cancellation token halts remaining items and reproducibility manifest is populated."""
        token = BatchCancellationToken()

        def _cancel_after_first(done_count: int, total_count: int, current_sid: str) -> None:
            if done_count >= 1:
                token.cancel()

        samples = [
            BatchSampleInput(sample_id="S1", filename="s1.txt", content="1257.4226 5000\n"),
            BatchSampleInput(sample_id="S2", filename="s2.txt", content="1485.5336 8000\n"),
            BatchSampleInput(sample_id="S3", filename="s3.txt", content="1809.6393 9000\n"),
        ]
        res = execute_batch_analysis(
            samples=samples,
            library=self.library,
            cancellation_token=token,
            progress_callback=_cancel_after_first,
        )
        self.assertEqual(res.status, BatchJobStatus.CANCELLED)
        self.assertEqual(res.completed_files_count, 1)
        self.assertEqual(res.cancelled_files_count, 2)
        self.assertEqual(res.sample_outcomes[0].status, BatchItemStatus.COMPLETED)
        self.assertEqual(res.sample_outcomes[1].status, BatchItemStatus.CANCELLED)

        # Check reproducibility record
        self.assertEqual(res.reproducibility.library_version, "1.0.0")
        self.assertEqual(res.reproducibility.mass_tolerance, 20.0)


if __name__ == "__main__":
    unittest.main()

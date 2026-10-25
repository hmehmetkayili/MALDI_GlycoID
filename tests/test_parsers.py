"""Unit tests for Stage 2: Input files and peak-list parsing (TXT, XY, CSV, TSV, mzML, mzXML, Bruker adapter)."""

from __future__ import annotations

import base64
import struct
import unittest
import zlib
from pathlib import Path

import numpy as np

from app.config import EXAMPLES_DIR
from app.parsers import (
    AmbiguousColumnError,
    BrukerFormatConversionRequiredError,
    DuplicateMzPolicy,
    SpectrumDataType,
    SpectrumParseError,
    parse_bruker_fid_with_acqus,
    parse_bruker_input,
    parse_csv_spectrum,
    parse_mzml_spectrum,
    parse_mzxml_spectrum,
    parse_spectrum_file,
    parse_txt_spectrum,
    parse_xy_spectrum,
)


def _build_synthetic_mzml(
    mzs: list[float],
    intensities: list[float],
    is_64bit: bool = True,
    use_zlib: bool = True,
    is_profile: bool = False,
) -> str:
    """Construct a valid PSI-MS mzML XML document with encoded binaryDataArrays."""
    fmt = f"<{len(mzs)}{'d' if is_64bit else 'f'}"
    mz_bytes = struct.pack(fmt, *mzs)
    int_bytes = struct.pack(fmt, *intensities)
    if use_zlib:
        mz_bytes = zlib.compress(mz_bytes)
        int_bytes = zlib.compress(int_bytes)

    mz_b64 = base64.b64encode(mz_bytes).decode("ascii")
    int_b64 = base64.b64encode(int_bytes).decode("ascii")

    bit_cv = (
        '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>'
        if is_64bit
        else '<cvParam cvRef="MS" accession="MS:1000521" name="32-bit float"/>'
    )
    comp_cv = (
        '<cvParam cvRef="MS" accession="MS:1000574" name="zlib compression"/>'
        if use_zlib
        else '<cvParam cvRef="MS" accession="MS:1000576" name="no compression"/>'
    )
    spec_type_cv = (
        '<cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum"/>'
        if is_profile
        else '<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum"/>'
    )

    return f"""<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">
  <run id="synthetic_run">
    <spectrumList count="1">
      <spectrum index="0" id="scan=1" defaultArrayLength="{len(mzs)}">
        <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>
        {spec_type_cv}
        <binaryDataArrayList count="2">
          <binaryDataArray>
            {bit_cv}
            {comp_cv}
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array"/>
            <binary>{mz_b64}</binary>
          </binaryDataArray>
          <binaryDataArray>
            {bit_cv}
            {comp_cv}
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array"/>
            <binary>{int_b64}</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
    </spectrumList>
  </run>
</mzML>"""


def _build_synthetic_mzxml(
    mzs: list[float],
    intensities: list[float],
    precision: int = 64,
    use_zlib: bool = True,
) -> str:
    """Construct a valid mzXML document with big-endian (network order) interleaved <peaks>."""
    interleaved: list[float] = []
    for m, i in zip(mzs, intensities):
        interleaved.extend([m, i])
    fmt = f">{len(interleaved)}{'d' if precision == 64 else 'f'}"
    raw = struct.pack(fmt, *interleaved)
    if use_zlib:
        raw = zlib.compress(raw)
    b64 = base64.b64encode(raw).decode("ascii")
    comp_attr = 'compressionType="zlib"' if use_zlib else 'compressionType="none"'
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2">
  <msRun scanCount="1">
    <scan num="1" msLevel="1" peaksCount="{len(mzs)}" centroided="1" retentionTime="PT45.0S">
      <peaks precision="{precision}" byteOrder="network" {comp_attr}>{b64}</peaks>
    </scan>
  </msRun>
</mzXML>"""


class TestSpectrumParsers(unittest.TestCase):
    """Comprehensive tests for TXT, XY, CSV, TSV, mzML, mzXML, and Bruker parsers."""

    def test_example_files_parsing(self) -> None:
        """Verify synthetic example files in data/examples/ parse cleanly into SpectrumData."""
        txt_spec = parse_spectrum_file(EXAMPLES_DIR / "synthetic_native_nglycans_centroid.txt")
        self.assertEqual(txt_spec.source_format, "TXT")
        self.assertEqual(txt_spec.data_type, SpectrumDataType.CENTROID)
        self.assertEqual(txt_spec.valid_points_count, 10)
        self.assertAlmostEqual(txt_spec.peaks[0].mz, 1257.4226, places=4)

        xy_spec = parse_spectrum_file(EXAMPLES_DIR / "synthetic_native_nglycans.xy")
        self.assertEqual(xy_spec.source_format, "XY")
        self.assertEqual(xy_spec.valid_points_count, 6)

        csv_spec = parse_spectrum_file(EXAMPLES_DIR / "synthetic_native_nglycans.csv")
        self.assertEqual(csv_spec.source_format, "CSV")
        self.assertEqual(csv_spec.valid_points_count, 6)
        self.assertAlmostEqual(csv_spec.peaks[0].area or 0.0, 1240.5, places=2)
        self.assertAlmostEqual(csv_spec.peaks[0].snr or 0.0, 42.3, places=2)

    def test_delimited_formats_and_decimal_comma(self) -> None:
        """Test semicolon-delimited files with European decimal commas and UTF-16 encoding."""
        eu_text = "# Comment line\nm/z;intensity\n1257,4226;8450,5\n1485,5336;15200,25\n"
        spec = parse_csv_spectrum(eu_text.encode("utf-16"), filename="eu_sample.csv")
        self.assertEqual(spec.valid_points_count, 2)
        self.assertAlmostEqual(spec.peaks[0].mz, 1257.4226, places=4)
        self.assertAlmostEqual(spec.peaks[0].intensity, 8450.5, places=2)
        self.assertEqual(spec.metadata["encoding"], "utf-16")

    def test_malformed_rows_duplicates_and_negative_intensities(self) -> None:
        """Test handling of malformed lines, negative intensities, zero intensities, and duplicate m/z."""
        raw = """m/z\tintensity
1000.1234\t1500
1000.1234\t4200
1001.1268\t-300
1002.1301\t0.0
bad_mz\t1200
1005.5000\tnot_a_number
1010.2500\t2800
"""
        spec_max = parse_txt_spectrum(
            raw,
            filename="dirty.tsv",
            duplicate_policy=DuplicateMzPolicy.KEEP_MAX_INTENSITY,
        )
        self.assertEqual(spec_max.source_format, "TSV")
        self.assertEqual(spec_max.valid_points_count, 2)  # 1000.1234 (deduped) and 1010.2500
        self.assertEqual(spec_max.duplicate_mz_count, 1)
        self.assertEqual(spec_max.negative_intensity_count, 1)
        self.assertEqual(spec_max.skipped_invalid_rows, 2)
        self.assertAlmostEqual(spec_max.peaks[0].intensity, 4200.0, places=2)

        spec_sum = parse_txt_spectrum(
            raw,
            filename="dirty.tsv",
            duplicate_policy=DuplicateMzPolicy.SUM_INTENSITIES,
        )
        self.assertAlmostEqual(spec_sum.peaks[0].intensity, 5700.0, places=2)

    def test_ambiguous_columns_never_silently_reinterpreted(self) -> None:
        """Verify that a multi-column file with ambiguous headers raises AmbiguousColumnError until user specifies columns."""
        ambiguous_csv = "ColA,ColB,ColC,ColD\n1,1257.4226,8450.0,0.12\n2,1485.5336,15200.0,0.14\n"
        with self.assertRaises(AmbiguousColumnError) as ctx:
            parse_csv_spectrum(ambiguous_csv, filename="ambiguous.csv")
        self.assertEqual(ctx.exception.detected_columns, ("ColA", "ColB", "ColC", "ColD"))

        # Now provide explicit user-selected column names
        resolved = parse_csv_spectrum(
            ambiguous_csv,
            filename="ambiguous.csv",
            mz_column="ColB",
            intensity_column="ColC",
        )
        self.assertEqual(resolved.valid_points_count, 2)
        self.assertEqual(resolved.selected_mz_column, "ColB")
        self.assertEqual(resolved.selected_intensity_column, "ColC")
        self.assertAlmostEqual(resolved.peaks[0].mz, 1257.4226, places=4)

    def test_profile_vs_centroid_classification_and_override(self) -> None:
        """Verify auto-detection of dense profile spectra and explicit user data-type override."""
        # Create a 300-point dense profile spectrum with step 0.02 Da
        mzs = [1257.0 + i * 0.02 for i in range(300)]
        ints = [0.0 if (i < 20 or i > 280) else 1000.0 for i in range(300)]
        xy_lines = "\n".join(f"{m:.4f} {y:.2f}" for m, y in zip(mzs, ints))

        auto_prof = parse_xy_spectrum(xy_lines, filename="dense_profile.xy")
        self.assertEqual(auto_prof.data_type, SpectrumDataType.PROFILE)
        # Profile spectra retain 0-intensity baseline points
        self.assertEqual(auto_prof.valid_points_count, 300)

        # User explicitly overrides as centroid
        user_cent = parse_xy_spectrum(
            xy_lines,
            filename="dense_profile.xy",
            user_specified_data_type=SpectrumDataType.CENTROID,
        )
        self.assertEqual(user_cent.data_type, SpectrumDataType.CENTROID)
        self.assertEqual(user_cent.user_specified_data_type, SpectrumDataType.CENTROID)

    def test_mzml_parser_32_64_zlib(self) -> None:
        """Test verified mzML parser with 64-bit zlib and 32-bit uncompressed arrays."""
        mzs = [1257.4226, 1485.5336, 1647.5865]
        ints = [8450.0, 15200.0, 24800.0]

        mzml_64_zlib = _build_synthetic_mzml(mzs, ints, is_64bit=True, use_zlib=True, is_profile=False)
        spec1 = parse_mzml_spectrum(mzml_64_zlib, filename="test64.mzML")
        self.assertEqual(spec1.source_format, "MZML")
        self.assertEqual(spec1.data_type, SpectrumDataType.CENTROID)
        self.assertEqual(spec1.valid_points_count, 3)
        self.assertAlmostEqual(spec1.peaks[0].mz, 1257.4226, places=6)

        mzml_32_raw = _build_synthetic_mzml(mzs, ints, is_64bit=False, use_zlib=False, is_profile=True)
        spec2 = parse_mzml_spectrum(mzml_32_raw, filename="test32.mzML")
        self.assertEqual(spec2.data_type, SpectrumDataType.PROFILE)
        self.assertAlmostEqual(spec2.peaks[1].mz, 1485.5336, places=2)

    def test_mzxml_parser(self) -> None:
        """Test verified mzXML parser with network-order 64-bit zlib-compressed peaks."""
        mzs = [1257.4226, 1485.5336, 1809.6393]
        ints = [8450.0, 15200.0, 19600.0]
        mzxml_str = _build_synthetic_mzxml(mzs, ints, precision=64, use_zlib=True)
        spec = parse_mzxml_spectrum(mzxml_str, filename="test.mzXML")
        self.assertEqual(spec.source_format, "MZXML")
        self.assertEqual(spec.data_type, SpectrumDataType.CENTROID)
        self.assertEqual(spec.valid_points_count, 3)
        self.assertAlmostEqual(spec.retention_time_sec or 0.0, 45.0, places=2)
        self.assertAlmostEqual(spec.peaks[2].mz, 1809.6393, places=5)

    def test_bruker_adapter_and_conversion_workflow(self) -> None:
        """Verify Bruker .fid + acqus parsing and explicit conversion error on proprietary .baf/.tdf."""
        # 1. Proprietary .baf or .tdf must raise BrukerFormatConversionRequiredError with instructions
        with self.assertRaises(BrukerFormatConversionRequiredError) as ctx:
            parse_bruker_input(b"BINARY_BAF_HEADER", filename="analysis.baf")
        self.assertIn("flexAnalysis", str(ctx.exception))
        self.assertIn("msconvert", str(ctx.exception))

        # 2. Valid 1D Bruker .fid + acqus pair
        acqus_text = """##TITLE= Synthetic Bruker Flex TOF Parameter File
##$TD= 10
##$DELAY= 20000
##$DW= 100
##$ML1= 250000.0
##$ML2= 100.0
##$ML3= 0.0
##$BYTORDP= 0
"""
        fid_ints = np.array([10, 20, 500, 1200, 300, 50, 10, 5, 0, 0], dtype="<i4").tobytes()
        fid_spec = parse_bruker_fid_with_acqus(fid_ints, acqus_text, filename="fid")
        self.assertEqual(fid_spec.source_format, "BRUKER_FID")
        self.assertEqual(fid_spec.data_type, SpectrumDataType.PROFILE)
        self.assertEqual(fid_spec.valid_points_count, 10)
        self.assertGreater(fid_spec.max_mz, fid_spec.min_mz)

    def test_invalid_and_unsupported_files(self) -> None:
        """Verify clear error messages on unsupported file formats and empty files."""
        with self.assertRaises(SpectrumParseError) as ctx1:
            parse_spectrum_file("1000 2000", filename="document.pdf")
        self.assertIn("Unsupported file format", str(ctx1.exception))

        with self.assertRaises(SpectrumParseError) as ctx2:
            parse_txt_spectrum("# Only comments\n\n", filename="empty.txt")
        self.assertIn("No valid m/z-intensity pairs found", str(ctx2.exception))


if __name__ == "__main__":
    unittest.main()

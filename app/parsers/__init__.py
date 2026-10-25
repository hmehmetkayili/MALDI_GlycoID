"""Unified spectrum parsing dispatcher and public interface for Stage 2."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

from app.parsers.base_models import (
    AmbiguousColumnError,
    DuplicateMzPolicy,
    SpectrumData,
    SpectrumDataType,
    SpectrumParseError,
    SpectrumPeak,
)
from app.parsers.bruker_adapter import (
    BRUKER_CONVERSION_GUIDE,
    BrukerFormatConversionRequiredError,
    parse_bruker_fid_with_acqus,
    parse_bruker_input,
)
from app.parsers.csv_parser import parse_csv_spectrum
from app.parsers.mzml_parser import parse_mzml_spectrum
from app.parsers.mzxml_parser import parse_mzxml_spectrum
from app.parsers.text_parser import parse_delimited_PeakList, parse_txt_spectrum
from app.parsers.xy_parser import parse_xy_spectrum


def parse_spectrum_file(
    content: Union[str, bytes, Path],
    filename: Optional[str] = None,
    sample_id: Optional[str] = None,
    delimiter: Optional[str] = None,
    decimal_separator: Optional[str] = None,
    mz_column: Optional[Union[str, int]] = None,
    intensity_column: Optional[Union[str, int]] = None,
    area_column: Optional[Union[str, int]] = None,
    snr_column: Optional[Union[str, int]] = None,
    user_specified_data_type: Optional[SpectrumDataType] = None,
    duplicate_policy: DuplicateMzPolicy = DuplicateMzPolicy.KEEP_MAX_INTENSITY,
    spectrum_index: int = 0,
    acqus_content: Optional[Union[str, bytes, Path]] = None,
) -> SpectrumData:
    """Auto-dispatch and parse any supported spectrum or peak list file into `SpectrumData`.

    Raises:
      - `BrukerFormatConversionRequiredError` for proprietary binary Bruker containers (`.baf`, `.tdf`, uncalibrated `.fid`).
      - `AmbiguousColumnError` when tabular files have uncertain columns requiring user assignment.
      - `SpectrumParseError` with `'Unsupported file format'` or `'No valid m/z-intensity pairs found'`.
    """
    if filename is None:
        if isinstance(content, Path):
            filename = content.name
        else:
            raise SpectrumParseError("Filename must be provided when parsing from raw string or bytes")

    ext = Path(filename).suffix.lower()
    lower_name = Path(filename).name.lower()

    if ext in (".baf", ".tdf", ".tsf", ".mcf", ".yep", ".fid") or lower_name in ("fid", "analysis.baf", "analysis.tdf"):
        return parse_bruker_input(
            content=content,
            filename=filename,
            sample_id=sample_id,
            acqus_content=acqus_content,
            user_specified_data_type=user_specified_data_type,
            duplicate_policy=duplicate_policy,
        )

    if ext in (".txt", ".tsv"):
        return parse_txt_spectrum(
            content=content,
            filename=filename,
            sample_id=sample_id,
            delimiter=delimiter,
            decimal_separator=decimal_separator,
            mz_column=mz_column,
            intensity_column=intensity_column,
            area_column=area_column,
            snr_column=snr_column,
            user_specified_data_type=user_specified_data_type,
            duplicate_policy=duplicate_policy,
        )

    if ext == ".xy":
        return parse_xy_spectrum(
            content=content,
            filename=filename,
            sample_id=sample_id,
            delimiter=delimiter,
            decimal_separator=decimal_separator,
            mz_column=mz_column,
            intensity_column=intensity_column,
            user_specified_data_type=user_specified_data_type,
            duplicate_policy=duplicate_policy,
        )

    if ext == ".csv":
        return parse_csv_spectrum(
            content=content,
            filename=filename,
            sample_id=sample_id,
            delimiter=delimiter,
            decimal_separator=decimal_separator,
            mz_column=mz_column,
            intensity_column=intensity_column,
            area_column=area_column,
            snr_column=snr_column,
            user_specified_data_type=user_specified_data_type,
            duplicate_policy=duplicate_policy,
        )

    if ext == ".mzml":
        return parse_mzml_spectrum(
            content=content,
            filename=filename,
            sample_id=sample_id,
            spectrum_index=spectrum_index,
            user_specified_data_type=user_specified_data_type,
            duplicate_policy=duplicate_policy,
        )

    if ext == ".mzxml":
        return parse_mzxml_spectrum(
            content=content,
            filename=filename,
            sample_id=sample_id,
            scan_index=spectrum_index,
            user_specified_data_type=user_specified_data_type,
            duplicate_policy=duplicate_policy,
        )

    raise SpectrumParseError(f"Unsupported file format: '{ext or filename}'")


__all__ = [
    "AmbiguousColumnError",
    "BRUKER_CONVERSION_GUIDE",
    "BrukerFormatConversionRequiredError",
    "DuplicateMzPolicy",
    "SpectrumData",
    "SpectrumDataType",
    "SpectrumParseError",
    "SpectrumPeak",
    "parse_bruker_fid_with_acqus",
    "parse_bruker_input",
    "parse_csv_spectrum",
    "parse_delimited_PeakList",
    "parse_mzml_spectrum",
    "parse_mzxml_spectrum",
    "parse_spectrum_file",
    "parse_txt_spectrum",
    "parse_xy_spectrum",
]

"""Parser for two-column .xy mass spectrum files.

Supports Bruker flexAnalysis and generic instrument XY exports where the first
column is m/z (x-axis) and the second column is intensity (y-axis).
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

from app.parsers.base_models import DuplicateMzPolicy, SpectrumData, SpectrumDataType
from app.parsers.text_parser import parse_delimited_PeakList


def parse_xy_spectrum(
    content: Union[str, bytes, Path],
    filename: str = "spectrum.xy",
    sample_id: Optional[str] = None,
    delimiter: Optional[str] = None,
    decimal_separator: Optional[str] = None,
    mz_column: Optional[Union[str, int]] = None,
    intensity_column: Optional[Union[str, int]] = None,
    user_specified_data_type: Optional[SpectrumDataType] = None,
    duplicate_policy: DuplicateMzPolicy = DuplicateMzPolicy.KEEP_MAX_INTENSITY,
) -> SpectrumData:
    """Parse a two-column .xy spectrum file into `SpectrumData`."""
    return parse_delimited_PeakList(
        content=content,
        filename=filename,
        sample_id=sample_id,
        source_format="XY",
        delimiter=delimiter,
        decimal_separator=decimal_separator,
        mz_column=mz_column,
        intensity_column=intensity_column,
        user_specified_data_type=user_specified_data_type,
        duplicate_policy=duplicate_policy,
    )

"""Parser for delimited text (.txt, .tsv) and generic tabular peak lists.

Supports whitespace, tab, semicolon, and comma delimiters, decimal separator
options ('.' or ','), header detection, comment lines ('#', ';', '//'),
and user-selectable column mapping without silently guessing ambiguous multi-column layouts.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

from app.parsers.base_models import (
    AmbiguousColumnError,
    DuplicateMzPolicy,
    SpectrumData,
    SpectrumDataType,
    SpectrumParseError,
    SpectrumPeak,
    decode_text_bytes,
    finalize_peaks_and_build_spectrum,
)

# Recognized column header synonyms (normalized lowercase without spaces/punctuation)
MZ_HEADER_ALIASES = {
    "mz",
    "m/z",
    "mass/charge",
    "masstocharge",
    "mass",
    "centroidmz",
    "observedmz",
    "moverz",
    "x",
    "xaxis",
}

INTENSITY_HEADER_ALIASES = {
    "intensity",
    "int",
    "intens",
    "abund",
    "abundance",
    "height",
    "peakheight",
    "signal",
    "counts",
    "absintensity",
    "relativeintensity",
    "relintensity",
    "y",
    "yaxis",
}

AREA_HEADER_ALIASES = {
    "area",
    "peakarea",
    "integratedarea",
    "auc",
    "volume",
}

SNR_HEADER_ALIASES = {
    "snr",
    "s/n",
    "signaltonoise",
    "sn",
}

FWHM_HEADER_ALIASES = {
    "fwhm",
    "width",
    "peakwidth",
}

CHARGE_HEADER_ALIASES = {
    "z",
    "charge",
    "chargestate",
}


def _normalize_header_token(token: str) -> str:
    return re.sub(r"[^a-z0-9/]", "", token.strip().lower())


def _try_parse_float(val_str: str, decimal_separator: str = ".") -> Optional[float]:
    s = val_str.strip().strip('"').strip("'")
    if not s or s.lower() in ("na", "n/a", "nan", "null", "none", "-"):
        return None
    if decimal_separator == ",":
        s = s.replace(".", "").replace(",", ".")
    try:
        v = float(s)
        if math.isnan(v) or math.isinf(v):
            return None
        return v
    except ValueError:
        return None


def _split_line_tokens(line: str, delimiter: Optional[str]) -> List[str]:
    if delimiter == "whitespace" or delimiter is None:
        return [tok for tok in re.split(r"\s+", line.strip()) if tok]
    return [tok.strip().strip('"').strip("'") for tok in line.split(delimiter)]


def detect_delimiter_and_decimal(
    non_empty_lines: Sequence[str],
    explicit_delimiter: Optional[str] = None,
    explicit_decimal: Optional[str] = None,
) -> Tuple[str, str]:
    """Detect column delimiter (tab, semicolon, comma, or whitespace) and decimal separator."""
    sample_lines = list(non_empty_lines[:20])
    if explicit_delimiter is not None:
        delim = explicit_delimiter
    else:
        # Check tab, semicolon, comma consistency across first lines
        if any("\t" in ln for ln in sample_lines):
            delim = "\t"
        elif any(";" in ln for ln in sample_lines):
            delim = ";"
        elif any("," in ln for ln in sample_lines):
            # Distinguish comma-delimited vs European decimal comma with whitespace delimiter
            comma_counts = [ln.count(",") for ln in sample_lines[:5]]
            space_tokens = [len(re.split(r"\s+", ln.strip())) for ln in sample_lines[:5]]
            if all(sp >= 2 for sp in space_tokens) and all(
                all(re.match(r"^[+-]?\d+,\d+([eE][+-]?\d+)?$", tok) for tok in re.split(r"\s+", ln.strip()))
                for ln in sample_lines[1:4]
                if ln.strip()
            ):
                delim = "whitespace"
            elif comma_counts and comma_counts[0] >= 1:
                delim = ","
            else:
                delim = "whitespace"
        else:
            delim = "whitespace"

    if explicit_decimal is not None:
        dec = explicit_decimal
    else:
        if delim in ("\t", ";", "whitespace"):
            # Check if numbers use comma as decimal separator (e.g., 1234,5678)
            tokens_to_test: List[str] = []
            for ln in sample_lines[1:6] or sample_lines[:1]:
                tokens_to_test.extend(_split_line_tokens(ln, delim))
            if tokens_to_test and any(re.match(r"^[+-]?\d+,\d+$", t) for t in tokens_to_test):
                dec = ","
            else:
                dec = "."
        else:
            dec = "."

    return delim, dec


def resolve_column_indices(
    header_tokens: Optional[Sequence[str]],
    num_columns: int,
    mz_column: Optional[Union[str, int]] = None,
    intensity_column: Optional[Union[str, int]] = None,
    area_column: Optional[Union[str, int]] = None,
    snr_column: Optional[Union[str, int]] = None,
) -> Tuple[int, int, Optional[int], Optional[int], Optional[int], Optional[int], List[str]]:
    """Resolve column indices for m/z, intensity, area, S/N, FWHM, and charge.

    CRITICAL RULE: Never silently reinterpret ambiguous columns!
    If a file has > 2 columns and no header (or an unrecognizable header) and the user
    has not explicitly specified `mz_column` and `intensity_column`, raise `AmbiguousColumnError`.
    """
    if num_columns < 2:
        raise SpectrumParseError(
            f"Spectrum file must contain at least 2 columns (m/z and intensity); found {num_columns}."
        )

    if header_tokens is not None:
        col_names = [h if h else f"Column_{i + 1}" for i, h in enumerate(header_tokens)]
    else:
        col_names = [f"Column_{i + 1}" for i in range(num_columns)]

    def _find_user_col(spec: Optional[Union[str, int]], role_name: str) -> Optional[int]:
        if spec is None:
            return None
        if isinstance(spec, int):
            if 0 <= spec < num_columns:
                return spec
            if 1 <= spec <= num_columns:
                return spec - 1
            raise SpectrumParseError(f"Selected {role_name} column index {spec} is out of range (1..{num_columns}).")
        spec_str = str(spec).strip()
        if not spec_str:
            return None
        # Exact match
        for i, name in enumerate(col_names):
            if name == spec_str or name.lower() == spec_str.lower():
                return i
        # Try 1-based numeric string
        if spec_str.isdigit():
            idx1 = int(spec_str)
            if 1 <= idx1 <= num_columns:
                return idx1 - 1
        raise SpectrumParseError(
            f"Selected {role_name} column '{spec}' not found among detected columns: {col_names}"
        )

    mz_idx = _find_user_col(mz_column, "m/z")
    int_idx = _find_user_col(intensity_column, "intensity")
    area_idx = _find_user_col(area_column, "area")
    snr_idx = _find_user_col(snr_column, "S/N")
    fwhm_idx: Optional[int] = None
    charge_idx: Optional[int] = None

    # Auto-detect from header names when available
    if header_tokens is not None:
        norm_headers = [_normalize_header_token(h) for h in col_names]
        for i, nh in enumerate(norm_headers):
            if mz_idx is None and nh in MZ_HEADER_ALIASES:
                mz_idx = i
            elif int_idx is None and nh in INTENSITY_HEADER_ALIASES:
                int_idx = i
            elif area_idx is None and nh in AREA_HEADER_ALIASES:
                area_idx = i
            elif snr_idx is None and nh in SNR_HEADER_ALIASES:
                snr_idx = i
            elif fwhm_idx is None and nh in FWHM_HEADER_ALIASES:
                fwhm_idx = i
            elif charge_idx is None and nh in CHARGE_HEADER_ALIASES:
                charge_idx = i

    # Standard 2-column layout without headers (or with generic two columns)
    if mz_idx is None and int_idx is None and num_columns == 2 and header_tokens is None:
        mz_idx = 0
        int_idx = 1

    # If still unresolved, refuse to silently guess on multi-column or ambiguous header files
    if mz_idx is None or int_idx is None:
        raise AmbiguousColumnError(
            f"Ambiguous column layout in file (detected columns: {col_names}). "
            f"Please explicitly specify the m/z and intensity columns.",
            detected_columns=col_names,
        )

    if mz_idx == int_idx:
        raise SpectrumParseError("m/z column and intensity column cannot be the same column.")

    return mz_idx, int_idx, area_idx, snr_idx, fwhm_idx, charge_idx, col_names


def parse_delimited_PeakList(
    content: Union[str, bytes, Path],
    filename: str = "spectrum.txt",
    sample_id: Optional[str] = None,
    source_format: str = "TXT",
    delimiter: Optional[str] = None,
    decimal_separator: Optional[str] = None,
    mz_column: Optional[Union[str, int]] = None,
    intensity_column: Optional[Union[str, int]] = None,
    area_column: Optional[Union[str, int]] = None,
    snr_column: Optional[Union[str, int]] = None,
    user_specified_data_type: Optional[SpectrumDataType] = None,
    duplicate_policy: DuplicateMzPolicy = DuplicateMzPolicy.KEEP_MAX_INTENSITY,
) -> SpectrumData:
    """Parse a delimited text/TSV/CSV spectrum or peak list into `SpectrumData`."""
    warnings: List[str] = []
    encoding_used = "utf-8"

    if isinstance(content, Path):
        filename = content.name
        raw_bytes = content.read_bytes()
        text, encoding_used = decode_text_bytes(raw_bytes)
    elif isinstance(content, bytes):
        text, encoding_used = decode_text_bytes(content)
    elif isinstance(content, str):
        text = content
    else:
        raise SpectrumParseError(f"Unsupported input content type: {type(content).__name__}")

    if sample_id is None:
        sample_id = Path(filename).stem

    # Filter out empty lines and comment lines (#, //, ;)
    raw_lines = text.splitlines()
    data_lines: List[Tuple[int, str]] = []
    comment_lines_count = 0

    for line_num, raw_ln in enumerate(raw_lines, start=1):
        stripped = raw_ln.strip()
        if not stripped:
            continue
        if stripped.startswith("#") or stripped.startswith("//"):
            comment_lines_count += 1
            continue
        data_lines.append((line_num, stripped))

    if not data_lines:
        raise SpectrumParseError("No valid m/z-intensity pairs found (file is empty or contains only comments)")

    delim, dec = detect_delimiter_and_decimal(
        [ln for _, ln in data_lines],
        explicit_delimiter=delimiter,
        explicit_decimal=decimal_separator,
    )

    # Inspect the first non-comment line to determine if it is a header row
    first_line_num, first_line_str = data_lines[0]
    first_tokens = _split_line_tokens(first_line_str, delim)
    if len(first_tokens) < 2 and len(data_lines) > 1:
        # Check if the first line is a title/metadata banner and second line is the actual header/data
        second_tokens = _split_line_tokens(data_lines[1][1], delim)
        if len(second_tokens) >= 2:
            warnings.append(f"Skipped non-tabular title banner on line {first_line_num}: '{first_line_str[:60]}'")
            data_lines = data_lines[1:]
            first_line_num, first_line_str = data_lines[0]
            first_tokens = second_tokens

    if len(first_tokens) < 2:
        raise SpectrumParseError("No valid m/z-intensity pairs found (fewer than 2 columns detected)")

    is_first_row_numeric = all(_try_parse_float(tok, dec) is not None for tok in first_tokens[:2])
    header_tokens: Optional[List[str]] = None
    start_idx = 0
    if not is_first_row_numeric:
        header_tokens = first_tokens
        start_idx = 1

    mz_idx, int_idx, area_idx, snr_idx, fwhm_idx, charge_idx, col_names = resolve_column_indices(
        header_tokens=header_tokens,
        num_columns=len(first_tokens),
        mz_column=mz_column,
        intensity_column=intensity_column,
        area_column=area_column,
        snr_column=snr_column,
    )

    peaks: List[SpectrumPeak] = []
    total_rows_read = 0
    skipped_invalid_rows = 0
    negative_intensity_count = 0
    zero_intensity_count = 0

    for line_num, ln_str in data_lines[start_idx:]:
        total_rows_read += 1
        tokens = _split_line_tokens(ln_str, delim)
        if len(tokens) <= max(mz_idx, int_idx):
            skipped_invalid_rows += 1
            continue

        mz_val = _try_parse_float(tokens[mz_idx], dec)
        int_val = _try_parse_float(tokens[int_idx], dec)

        if mz_val is None or int_val is None or mz_val <= 0.0:
            skipped_invalid_rows += 1
            continue

        if int_val < 0.0:
            negative_intensity_count += 1
            continue

        if int_val == 0.0:
            zero_intensity_count += 1

        area_val = (
            _try_parse_float(tokens[area_idx], dec)
            if (area_idx is not None and area_idx < len(tokens))
            else None
        )
        snr_val = (
            _try_parse_float(tokens[snr_idx], dec)
            if (snr_idx is not None and snr_idx < len(tokens))
            else None
        )
        fwhm_val = (
            _try_parse_float(tokens[fwhm_idx], dec)
            if (fwhm_idx is not None and fwhm_idx < len(tokens))
            else None
        )
        charge_val: Optional[int] = None
        if charge_idx is not None and charge_idx < len(tokens):
            c_float = _try_parse_float(tokens[charge_idx], dec)
            if c_float is not None and c_float.is_integer() and int(c_float) != 0:
                charge_val = int(c_float)

        peaks.append(
            SpectrumPeak(
                mz=mz_val,
                intensity=int_val,
                area=area_val,
                snr=snr_val,
                fwhm=fwhm_val,
                charge=charge_val,
                row_index=line_num,
            )
        )

    return finalize_peaks_and_build_spectrum(
        raw_peaks=peaks,
        sample_id=sample_id,
        filename=filename,
        source_format=source_format,
        user_specified_data_type=user_specified_data_type,
        detected_columns=col_names,
        selected_mz_column=col_names[mz_idx],
        selected_intensity_column=col_names[int_idx],
        selected_area_column=col_names[area_idx] if area_idx is not None else None,
        selected_snr_column=col_names[snr_idx] if snr_idx is not None else None,
        total_rows_read=total_rows_read,
        skipped_invalid_rows=skipped_invalid_rows,
        negative_intensity_count=negative_intensity_count,
        zero_intensity_count=zero_intensity_count,
        duplicate_policy=duplicate_policy,
        warnings=warnings,
        metadata={
            "encoding": encoding_used,
            "delimiter": "whitespace" if delim == "whitespace" else repr(delim),
            "decimal_separator": dec,
            "comment_lines_skipped": comment_lines_count,
        },
    )


def parse_txt_spectrum(
    content: Union[str, bytes, Path],
    filename: str = "spectrum.txt",
    sample_id: Optional[str] = None,
    delimiter: Optional[str] = None,
    decimal_separator: Optional[str] = None,
    mz_column: Optional[Union[str, int]] = None,
    intensity_column: Optional[Union[str, int]] = None,
    area_column: Optional[Union[str, int]] = None,
    snr_column: Optional[Union[str, int]] = None,
    user_specified_data_type: Optional[SpectrumDataType] = None,
    duplicate_policy: DuplicateMzPolicy = DuplicateMzPolicy.KEEP_MAX_INTENSITY,
) -> SpectrumData:
    """Parse a .txt or .tsv mass spectrum file."""
    ext = Path(filename).suffix.lower()
    fmt = "TSV" if ext == ".tsv" else "TXT"
    default_delim = "\t" if (ext == ".tsv" and delimiter is None) else delimiter
    return parse_delimited_PeakList(
        content=content,
        filename=filename,
        sample_id=sample_id,
        source_format=fmt,
        delimiter=default_delim,
        decimal_separator=decimal_separator,
        mz_column=mz_column,
        intensity_column=intensity_column,
        area_column=area_column,
        snr_column=snr_column,
        user_specified_data_type=user_specified_data_type,
        duplicate_policy=duplicate_policy,
    )

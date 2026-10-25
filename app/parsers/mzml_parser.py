"""Verified PSI-MS mzML parser supporting 32-bit/64-bit IEEE-754 and zlib-compressed binary data arrays.

Extracts MS level, centroid/profile CV terms (MS:1000127 / MS:1000128), precursor m/z
and charge (for MS2 spectra), and binary m/z (MS:1000514) and intensity (MS:1000515) arrays.
"""

from __future__ import annotations

import base64
import math
import xml.etree.ElementTree as ET
import zlib
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from app.parsers.base_models import (
    DuplicateMzPolicy,
    SpectrumData,
    SpectrumDataType,
    SpectrumParseError,
    SpectrumPeak,
    finalize_peaks_and_build_spectrum,
)


def _strip_ns(tag: str) -> str:
    if "}" in tag:
        return tag.split("}", 1)[1]
    return tag


def _decode_mzml_binary_array(
    encoded_b64: str,
    is_64bit: bool,
    is_zlib: bool,
) -> np.ndarray:
    """Decode base64 (and optional zlib) mzML binary array into a 1D float64 numpy array."""
    cleaned = encoded_b64.strip()
    if not cleaned:
        return np.array([], dtype=np.float64)
    try:
        raw_bytes = base64.b64decode(cleaned)
        if is_zlib:
            raw_bytes = zlib.decompress(raw_bytes)
        dtype = "<f8" if is_64bit else "<f4"
        arr = np.frombuffer(raw_bytes, dtype=dtype).astype(np.float64)
        return arr
    except Exception as exc:
        raise SpectrumParseError(f"Failed to decode mzML binary data array: {exc}") from exc


def parse_mzml_spectrum(
    content: Union[str, bytes, Path],
    filename: str = "spectrum.mzML",
    sample_id: Optional[str] = None,
    spectrum_index: int = 0,
    target_ms_level: Optional[int] = None,
    user_specified_data_type: Optional[SpectrumDataType] = None,
    duplicate_policy: DuplicateMzPolicy = DuplicateMzPolicy.KEEP_MAX_INTENSITY,
) -> SpectrumData:
    """Parse an mzML file and extract the requested spectrum into `SpectrumData`."""
    if isinstance(content, Path):
        filename = content.name
        xml_bytes = content.read_bytes()
    elif isinstance(content, bytes):
        xml_bytes = content
    elif isinstance(content, str):
        xml_bytes = content.encode("utf-8")
    else:
        raise SpectrumParseError(f"Unsupported mzML input type: {type(content).__name__}")

    if sample_id is None:
        sample_id = Path(filename).stem

    try:
        root = ET.fromstring(xml_bytes)
    except ET.ParseError as exc:
        raise SpectrumParseError(f"Invalid XML syntax in mzML file '{filename}': {exc}") from exc

    # Collect all <spectrum> elements
    spectrum_elements = [el for el in root.iter() if _strip_ns(el.tag) == "spectrum"]
    if not spectrum_elements:
        raise SpectrumParseError(f"No <spectrum> elements found in mzML file '{filename}'")

    selected_el: Optional[ET.Element] = None
    matched_idx = 0

    for spec_el in spectrum_elements:
        # Check MS level if target_ms_level is specified
        ms_level = 1
        for cv in spec_el.iter():
            if _strip_ns(cv.tag) == "cvParam":
                acc = cv.attrib.get("accession", "")
                name = cv.attrib.get("name", "").lower()
                if acc == "MS:1000511" or name == "ms level":
                    try:
                        ms_level = int(cv.attrib.get("value", "1"))
                    except ValueError:
                        ms_level = 1
        if target_ms_level is not None and ms_level != target_ms_level:
            continue
        if matched_idx == spectrum_index:
            selected_el = spec_el
            break
        matched_idx += 1

    if selected_el is None:
        raise SpectrumParseError(
            f"Requested spectrum index {spectrum_index} (target_ms_level={target_ms_level}) not found in '{filename}'."
        )

    spectrum_id = selected_el.attrib.get("id", f"index={spectrum_index}")
    ms_level = 1
    cv_data_type: Optional[SpectrumDataType] = None
    precursor_mz: Optional[float] = None
    precursor_charge: Optional[int] = None
    retention_time_sec: Optional[float] = None

    # Inspect spectrum-level cvParams
    for el in selected_el:
        tag = _strip_ns(el.tag)
        if tag == "cvParam":
            acc = el.attrib.get("accession", "")
            name = el.attrib.get("name", "").lower()
            val = el.attrib.get("value", "")
            if acc == "MS:1000511" or name == "ms level":
                try:
                    ms_level = int(val)
                except ValueError:
                    pass
            elif acc == "MS:1000127" or "centroid spectrum" in name:
                cv_data_type = SpectrumDataType.CENTROID
            elif acc == "MS:1000128" or "profile spectrum" in name:
                cv_data_type = SpectrumDataType.PROFILE
        elif tag == "scanList":
            for cv in el.iter():
                if _strip_ns(cv.tag) == "cvParam":
                    acc = cv.attrib.get("accession", "")
                    name = cv.attrib.get("name", "").lower()
                    if acc == "MS:1000016" or "scan start time" in name:
                        try:
                            rt_val = float(cv.attrib.get("value", "0"))
                            unit_name = cv.attrib.get("unitName", "").lower()
                            retention_time_sec = rt_val * 60.0 if "minute" in unit_name else rt_val
                        except ValueError:
                            pass
        elif tag == "precursorList":
            for cv in el.iter():
                if _strip_ns(cv.tag) == "cvParam":
                    acc = cv.attrib.get("accession", "")
                    name = cv.attrib.get("name", "").lower()
                    if acc in ("MS:1000744", "MS:1000040") or "selected ion m/z" in name:
                        try:
                            precursor_mz = float(cv.attrib.get("value", "0"))
                        except ValueError:
                            pass
                    elif acc == "MS:1000041" or "charge state" in name:
                        try:
                            precursor_charge = int(cv.attrib.get("value", "0"))
                        except ValueError:
                            pass

    # Parse binaryDataArrayList
    mz_arr: Optional[np.ndarray] = None
    int_arr: Optional[np.ndarray] = None

    for bda in selected_el.iter():
        if _strip_ns(bda.tag) != "binaryDataArray":
            continue
        is_64bit = True
        is_zlib = False
        array_type: Optional[str] = None
        b64_text = ""

        for child in bda:
            ctag = _strip_ns(child.tag)
            if ctag == "cvParam":
                acc = child.attrib.get("accession", "")
                name = child.attrib.get("name", "").lower()
                if acc == "MS:1000523" or "64-bit float" in name:
                    is_64bit = True
                elif acc == "MS:1000521" or "32-bit float" in name:
                    is_64bit = False
                elif acc == "MS:1000574" or "zlib" in name:
                    is_zlib = True
                elif acc == "MS:1000576" or "no compression" in name:
                    is_zlib = False
                elif acc == "MS:1000514" or "m/z array" in name:
                    array_type = "mz"
                elif acc == "MS:1000515" or "intensity array" in name:
                    array_type = "intensity"
            elif ctag == "binary":
                b64_text = child.text or ""

        if array_type == "mz":
            mz_arr = _decode_mzml_binary_array(b64_text, is_64bit=is_64bit, is_zlib=is_zlib)
        elif array_type == "intensity":
            int_arr = _decode_mzml_binary_array(b64_text, is_64bit=is_64bit, is_zlib=is_zlib)

    if mz_arr is None or int_arr is None:
        raise SpectrumParseError(f"mzML spectrum '{spectrum_id}' is missing m/z or intensity binaryDataArray")

    if len(mz_arr) != len(int_arr):
        raise SpectrumParseError(
            f"mzML binary array length mismatch in '{spectrum_id}': m/z ({len(mz_arr)}) != intensity ({len(int_arr)})"
        )

    raw_peaks: List[SpectrumPeak] = []
    skipped = 0
    neg_count = 0
    zero_count = 0

    for idx, (mz_v, int_v) in enumerate(zip(mz_arr.tolist(), int_arr.tolist()), start=1):
        if math.isnan(mz_v) or math.isinf(mz_v) or math.isnan(int_v) or math.isinf(int_v) or mz_v <= 0.0:
            skipped += 1
            continue
        if int_v < 0.0:
            neg_count += 1
            continue
        if int_v == 0.0:
            zero_count += 1
        raw_peaks.append(SpectrumPeak(mz=float(mz_v), intensity=float(int_v), row_index=idx))

    return finalize_peaks_and_build_spectrum(
        raw_peaks=raw_peaks,
        sample_id=sample_id,
        filename=filename,
        source_format="MZML",
        user_specified_data_type=user_specified_data_type,
        file_metadata_data_type=cv_data_type,
        detected_columns=("m/z", "intensity"),
        selected_mz_column="m/z",
        selected_intensity_column="intensity",
        total_rows_read=len(mz_arr),
        skipped_invalid_rows=skipped,
        negative_intensity_count=neg_count,
        zero_intensity_count=zero_count,
        duplicate_policy=duplicate_policy,
        ms_level=ms_level,
        precursor_mz=precursor_mz,
        precursor_charge=precursor_charge,
        retention_time_sec=retention_time_sec,
        metadata={
            "spectrum_id": spectrum_id,
            "total_spectra_in_file": len(spectrum_elements),
        },
    )

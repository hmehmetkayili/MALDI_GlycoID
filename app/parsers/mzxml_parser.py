"""Verified mzXML spectrum parser supporting 32-bit/64-bit big-endian (network order) and zlib-compressed <peaks> elements."""

from __future__ import annotations

import base64
import math
import re
import xml.etree.ElementTree as ET
import zlib
from pathlib import Path
from typing import List, Optional, Tuple, Union

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


def _parse_iso8601_duration_seconds(rt_str: str) -> Optional[float]:
    """Parse mzXML retentionTime attribute (e.g. 'PT125.4S' or 'PT2M5.4S') into seconds."""
    m = re.match(r"^PT(?:(\d+(?:\.\d+)?)M)?(?:(\d+(?:\.\d+)?)S)?$", rt_str.strip())
    if not m:
        return None
    mins = float(m.group(1)) if m.group(1) else 0.0
    secs = float(m.group(2)) if m.group(2) else 0.0
    return mins * 60.0 + secs


def parse_mzxml_spectrum(
    content: Union[str, bytes, Path],
    filename: str = "spectrum.mzXML",
    sample_id: Optional[str] = None,
    scan_index: int = 0,
    target_ms_level: Optional[int] = None,
    user_specified_data_type: Optional[SpectrumDataType] = None,
    duplicate_policy: DuplicateMzPolicy = DuplicateMzPolicy.KEEP_MAX_INTENSITY,
) -> SpectrumData:
    """Parse an mzXML file and extract the requested <scan> into `SpectrumData`."""
    if isinstance(content, Path):
        filename = content.name
        xml_bytes = content.read_bytes()
    elif isinstance(content, bytes):
        xml_bytes = content
    elif isinstance(content, str):
        xml_bytes = content.encode("utf-8")
    else:
        raise SpectrumParseError(f"Unsupported mzXML input type: {type(content).__name__}")

    if sample_id is None:
        sample_id = Path(filename).stem

    try:
        root = ET.fromstring(xml_bytes)
    except ET.ParseError as exc:
        raise SpectrumParseError(f"Invalid XML syntax in mzXML file '{filename}': {exc}") from exc

    scan_elements = [el for el in root.iter() if _strip_ns(el.tag) == "scan"]
    if not scan_elements:
        raise SpectrumParseError(f"No <scan> elements found in mzXML file '{filename}'")

    selected_scan: Optional[ET.Element] = None
    matched_idx = 0
    for scan_el in scan_elements:
        try:
            ms_level = int(scan_el.attrib.get("msLevel", "1"))
        except ValueError:
            ms_level = 1
        if target_ms_level is not None and ms_level != target_ms_level:
            continue
        if matched_idx == scan_index:
            selected_scan = scan_el
            break
        matched_idx += 1

    if selected_scan is None:
        raise SpectrumParseError(
            f"Requested scan index {scan_index} (target_ms_level={target_ms_level}) not found in '{filename}'."
        )

    scan_num = selected_scan.attrib.get("num", str(scan_index + 1))
    ms_level = int(selected_scan.attrib.get("msLevel", "1"))
    centroided_attr = selected_scan.attrib.get("centroided", "").strip()
    file_data_type: Optional[SpectrumDataType] = None
    if centroided_attr == "1":
        file_data_type = SpectrumDataType.CENTROID
    elif centroided_attr == "0":
        file_data_type = SpectrumDataType.PROFILE

    rt_sec: Optional[float] = None
    if "retentionTime" in selected_scan.attrib:
        rt_sec = _parse_iso8601_duration_seconds(selected_scan.attrib["retentionTime"])

    precursor_mz: Optional[float] = None
    precursor_charge: Optional[int] = None
    peaks_el: Optional[ET.Element] = None

    for child in selected_scan:
        ctag = _strip_ns(child.tag)
        if ctag == "precursorMz":
            try:
                if child.text:
                    precursor_mz = float(child.text.strip())
            except ValueError:
                pass
            if "precursorCharge" in child.attrib:
                try:
                    precursor_charge = int(child.attrib["precursorCharge"])
                except ValueError:
                    pass
        elif ctag == "peaks":
            peaks_el = child

    if peaks_el is None:
        raise SpectrumParseError(f"No <peaks> element found inside <scan num='{scan_num}'> of '{filename}'")

    precision = peaks_el.attrib.get("precision", "32").strip()
    compression = peaks_el.attrib.get("compressionType", "none").strip().lower()
    byte_order = peaks_el.attrib.get("byteOrder", "network").strip().lower()

    b64_payload = (peaks_el.text or "").strip()
    if not b64_payload:
        raise SpectrumParseError("No valid m/z-intensity pairs found (empty <peaks> array in mzXML)")

    try:
        raw_bytes = base64.b64decode(b64_payload)
        if compression == "zlib":
            raw_bytes = zlib.decompress(raw_bytes)
        endian = ">" if byte_order in ("network", "big") else "<"
        float_code = "f8" if precision == "64" else "f4"
        interleaved = np.frombuffer(raw_bytes, dtype=f"{endian}{float_code}").astype(np.float64)
    except Exception as exc:
        raise SpectrumParseError(f"Failed to decode mzXML <peaks> binary data: {exc}") from exc

    if len(interleaved) % 2 != 0:
        raise SpectrumParseError(f"Corrupted mzXML <peaks> array: odd number of float values ({len(interleaved)})")

    mz_arr = interleaved[0::2]
    int_arr = interleaved[1::2]

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
        source_format="MZXML",
        user_specified_data_type=user_specified_data_type,
        file_metadata_data_type=file_data_type,
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
        retention_time_sec=rt_sec,
        metadata={
            "scan_num": scan_num,
            "precision": precision,
            "compression": compression,
            "total_scans_in_file": len(scan_elements),
        },
    )

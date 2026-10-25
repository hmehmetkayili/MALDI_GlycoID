"""Bruker MALDI-TOF format adapter and explicit conversion workflow.

Scientific and technical honesty policy:
1. Bruker flexAnalysis exported text (`.txt`), `.xy`, and `.csv` peak lists or profile spectra
   are natively parsed and validated.
2. Unencrypted 1D Bruker XMASS/Flex `.fid` + `acqus` pairs (where `acqus` provides JCAMP-DX
   calibration constants `$ML1`, `$ML2`, `$ML3`, `$DELAY`, `$DW`, `$TD`, `$BYTORDP`) can be read
   if both files are provided and calibration constants are non-zero.
3. Proprietary encrypted or vendor-SDK-locked Bruker containers (`.baf`, `.tdf`, `.tsf`, `.mcf`,
   or `.fid` without `acqus` calibration metadata) CANNOT be reliably parsed without proprietary
   vendor DLLs. When encountered, `bruker_adapter` raises a descriptive `BrukerFormatConversionRequiredError`
   with step-by-step instructions to export `.xy`/`.txt` from flexAnalysis or convert to `.mzML`
   via ProteoWizard `msconvert` / CompassXport.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from app.parsers.base_models import (
    DuplicateMzPolicy,
    SpectrumData,
    SpectrumDataType,
    SpectrumParseError,
    SpectrumPeak,
    decode_text_bytes,
    finalize_peaks_and_build_spectrum,
)
from app.parsers.text_parser import parse_delimited_PeakList


class BrukerFormatConversionRequiredError(SpectrumParseError):
    """Raised when a proprietary binary Bruker file requires vendor export or msconvert conversion."""

    def __init__(self, filename: str, format_type: str, conversion_instructions: Tuple[str, ...]) -> None:
        msg = (
            f"Unsupported proprietary Bruker binary format ({format_type}) in '{filename}'. "
            f"{' '.join(conversion_instructions)}"
        )
        super().__init__(msg)
        self.filename = filename
        self.format_type = format_type
        self.conversion_instructions = conversion_instructions


BRUKER_CONVERSION_GUIDE: Tuple[str, ...] = (
    "Option 1 (Bruker flexAnalysis): Open the spectrum in flexAnalysis -> File -> Export -> Mass List (.txt / .csv) or Spectrum (.xy).",
    "Option 2 (ProteoWizard msconvert): Run `msconvert <input.d> --mzML --64 --zlib` to produce a standard PSI-MS .mzML file.",
    "Option 3 (Bruker CompassXport): Convert the .d / .baf / .fid directory to .mzXML or .mzML.",
)


def parse_jcamp_acqus_parameters(acqus_content: Union[str, bytes, Path]) -> Dict[str, float]:
    """Parse Bruker JCAMP-DX `acqus` or `acqu` file for TOF calibration parameters.

    Extracts:
      $TD (number of time-domain points)
      $DELAY (initial TOF delay in ns)
      $DW (dwell time per point in ns)
      $ML1, $ML2, $ML3 (Bruker TOF mass calibration constants)
      $BYTORDP (0 = little endian, 1 = big endian)
    """
    if isinstance(acqus_content, Path):
        text, _ = decode_text_bytes(acqus_content.read_bytes())
    elif isinstance(acqus_content, bytes):
        text, _ = decode_text_bytes(acqus_content)
    else:
        text = acqus_content

    params: Dict[str, float] = {}
    pattern = re.compile(r"^##\$([A-Za-z0-9_]+)=\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", re.MULTILINE)
    for m in pattern.finditer(text):
        key = m.group(1).upper()
        try:
            params[key] = float(m.group(2))
        except ValueError:
            continue
    return params


def convert_bruker_tof_time_to_mz(
    index_array: np.ndarray,
    delay_ns: float,
    dw_ns: float,
    ml1: float,
    ml2: float,
    ml3: float,
) -> np.ndarray:
    """Convert Bruker TOF time-of-flight bin indices to m/z using Bruker XMASS calibration equation.

    Bruker Flex TOF equation:
      t_i = DELAY + i * DW  (in nanoseconds; Bruker calibration uses time in microseconds or 100-ns units:
      specifically, A*t^2 + B*t + C = 0 where:
      sqrt(m/z) satisfies: t_tof = DELAY + i * DW, and for linear/quadratic ML3=0:
      m/z = ((t_tof - ML2)^2 * ML1) / 1e12).
    """
    if ml1 <= 0.0:
        raise SpectrumParseError(f"Invalid Bruker calibration constant ML1 ({ml1}) in acqus file.")

    t_ns = delay_ns + index_array * dw_ns
    dt = np.maximum(0.0, t_ns - ml2)

    # Standard Bruker linear sqrt(m/z) calibration (ML3 == 0 or negligible):
    # sqrt(m/z) = sqrt(ML1 / 1e12) * (t_ns - ML2) => m/z = (ML1 * (t_ns - ML2)^2) / 1e12
    mz = (ml1 * (dt ** 2)) / 1e12
    if abs(ml3) > 1e-15:
        # Apply small quadratic perturbation when ML3 is non-zero
        mz = mz + ml3 * ((dt / 1e6) ** 4)
    return mz


def parse_bruker_fid_with_acqus(
    fid_content: Union[bytes, Path],
    acqus_content: Union[str, bytes, Path],
    filename: str = "fid",
    sample_id: Optional[str] = None,
    user_specified_data_type: Optional[SpectrumDataType] = None,
) -> SpectrumData:
    """Parse a 1D Bruker `.fid` 32-bit integer array paired with its JCAMP-DX `acqus` parameter file."""
    if isinstance(fid_content, Path):
        filename = fid_content.name
        if sample_id is None:
            sample_id = fid_content.parent.name or "bruker_fid"
        raw_fid = fid_content.read_bytes()
    else:
        raw_fid = fid_content
        if sample_id is None:
            sample_id = Path(filename).stem

    params = parse_jcamp_acqus_parameters(acqus_content)
    required = ("ML1", "ML2", "DELAY", "DW")
    missing = [k for k in required if k not in params]
    if missing or params.get("ML1", 0.0) <= 0.0 or params.get("DW", 0.0) <= 0.0:
        raise BrukerFormatConversionRequiredError(
            filename=filename,
            format_type="Bruker .fid with incomplete/missing acqus calibration parameters",
            conversion_instructions=BRUKER_CONVERSION_GUIDE,
        )

    bytordp = int(params.get("BYTORDP", 0))
    dtype = ">i4" if bytordp == 1 else "<i4"
    if len(raw_fid) < 8 or len(raw_fid) % 4 != 0:
        raise SpectrumParseError(f"Invalid Bruker .fid byte length ({len(raw_fid)} bytes)")

    intensities = np.frombuffer(raw_fid, dtype=dtype).astype(np.float64)
    td = int(params.get("TD", len(intensities)))
    if 0 < td < len(intensities):
        intensities = intensities[:td]

    indices = np.arange(len(intensities), dtype=np.float64)
    mzs = convert_bruker_tof_time_to_mz(
        index_array=indices,
        delay_ns=params["DELAY"],
        dw_ns=params["DW"],
        ml1=params["ML1"],
        ml2=params["ML2"],
        ml3=params.get("ML3", 0.0),
    )

    raw_peaks: List[SpectrumPeak] = []
    neg_count = 0
    zero_count = 0
    for idx, (mz_v, int_v) in enumerate(zip(mzs.tolist(), intensities.tolist()), start=1):
        if mz_v <= 0.0 or math.isnan(mz_v):
            continue
        if int_v < 0.0:
            neg_count += 1
            int_v = 0.0
        if int_v == 0.0:
            zero_count += 1
        raw_peaks.append(SpectrumPeak(mz=mz_v, intensity=int_v, row_index=idx))

    return finalize_peaks_and_build_spectrum(
        raw_peaks=raw_peaks,
        sample_id=sample_id,
        filename=filename,
        source_format="BRUKER_FID",
        user_specified_data_type=user_specified_data_type,
        file_metadata_data_type=SpectrumDataType.PROFILE,
        detected_columns=("m/z (calibrated from TOF)", "intensity"),
        selected_mz_column="m/z (calibrated from TOF)",
        selected_intensity_column="intensity",
        total_rows_read=len(intensities),
        negative_intensity_count=neg_count,
        zero_intensity_count=zero_count,
        keep_zero_intensities=True,
        metadata={
            "ML1": params["ML1"],
            "ML2": params["ML2"],
            "ML3": params.get("ML3", 0.0),
            "DELAY": params["DELAY"],
            "DW": params["DW"],
        },
    )


def parse_bruker_input(
    content: Union[str, bytes, Path],
    filename: str = "bruker_export.txt",
    sample_id: Optional[str] = None,
    acqus_content: Optional[Union[str, bytes, Path]] = None,
    user_specified_data_type: Optional[SpectrumDataType] = None,
    duplicate_policy: DuplicateMzPolicy = DuplicateMzPolicy.KEEP_MAX_INTENSITY,
) -> SpectrumData:
    """Parse Bruker flexAnalysis text/XY/CSV exports or `.fid`+`acqus`, or raise conversion guidance for `.baf`/`.tdf`."""
    ext = Path(filename).suffix.lower()
    lower_name = Path(filename).name.lower()

    if ext in (".baf", ".tdf", ".tsf", ".mcf", ".yep") or lower_name in ("analysis.baf", "analysis.tdf"):
        raise BrukerFormatConversionRequiredError(
            filename=filename,
            format_type=f"Bruker {ext or lower_name} binary container",
            conversion_instructions=BRUKER_CONVERSION_GUIDE,
        )

    if lower_name == "fid" or ext == ".fid":
        if acqus_content is None and isinstance(content, Path):
            for candidate_name in ("acqus", "acqu"):
                candidate_path = content.parent / candidate_name
                if candidate_path.is_file():
                    acqus_content = candidate_path
                    break
        if acqus_content is None:
            raise BrukerFormatConversionRequiredError(
                filename=filename,
                format_type="Raw Bruker .fid without accompanying 'acqus' calibration parameter file",
                conversion_instructions=BRUKER_CONVERSION_GUIDE,
            )
        if isinstance(content, str):
            raise SpectrumParseError("Bruker .fid content must be binary bytes or a file Path")
        return parse_bruker_fid_with_acqus(
            fid_content=content,
            acqus_content=acqus_content,
            filename=filename,
            sample_id=sample_id,
            user_specified_data_type=user_specified_data_type,
        )

    # Otherwise, parse as a Bruker flexAnalysis text/XY/CSV peak list or profile export
    return parse_delimited_PeakList(
        content=content,
        filename=filename,
        sample_id=sample_id,
        source_format="BRUKER_EXPORT",
        user_specified_data_type=user_specified_data_type,
        duplicate_policy=duplicate_policy,
    )

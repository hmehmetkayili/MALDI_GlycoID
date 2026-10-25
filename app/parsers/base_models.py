"""Standardized data models and parser utilities for mass spectrometry inputs.

Every successful spectrum import produces a `SpectrumData` object with documented
fields, explicit units (`Da/e` for m/z, arbitrary instrument counts for intensity),
and an explicit `SpectrumDataType` (`CENTROID`, `PROFILE`, or `UNKNOWN`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Sequence, Tuple


class SpectrumDataType(str, Enum):
    """Classification of imported mass spectrum signal representation."""

    CENTROID = "centroid"
    PROFILE = "profile"
    UNKNOWN = "unknown"


class DuplicateMzPolicy(str, Enum):
    """Policy for handling duplicate m/z values in a peak list."""

    KEEP_ALL = "keep_all"
    KEEP_MAX_INTENSITY = "keep_max_intensity"
    SUM_INTENSITIES = "sum_intensities"
    AVERAGE_INTENSITIES = "average_intensities"


class SpectrumParseError(ValueError):
    """Raised when an input file cannot be parsed into a valid spectrum."""


class AmbiguousColumnError(SpectrumParseError):
    """Raised when a tabular file has multiple numeric columns and no clear m/z or intensity header."""

    def __init__(self, message: str, detected_columns: Sequence[str]) -> None:
        super().__init__(message)
        self.detected_columns = tuple(detected_columns)


@dataclass(frozen=True)
class SpectrumPeak:
    """Single data point or centroided peak in a mass spectrum."""

    mz: float  # Mass-to-charge ratio (units: m/z, Thomson / Da per elementary charge)
    intensity: float  # Raw observed signal intensity (instrument units / counts)
    area: Optional[float] = None  # Integrated peak area if present in file or integrated from profile
    snr: Optional[float] = None  # Signal-to-noise ratio if available
    fwhm: Optional[float] = None  # Full width at half maximum (Da) if available
    charge: Optional[int] = None  # Assigned charge state if present in peak list
    row_index: Optional[int] = None  # 1-based source row or point index


@dataclass(frozen=True)
class SpectrumData:
    """Standardized internal representation of an imported mass spectrum or peak list."""

    sample_id: str
    filename: str
    source_format: str  # "TXT", "XY", "CSV", "TSV", "MZML", "MZXML", "BRUKER_EXPORT", "BRUKER_FID"
    data_type: SpectrumDataType
    user_specified_data_type: Optional[SpectrumDataType]
    peaks: Tuple[SpectrumPeak, ...]
    mz_unit: str = "m/z (Da/e)"
    intensity_unit: str = "intensity (a.u.)"
    detected_columns: Tuple[str, ...] = field(default_factory=tuple)
    selected_mz_column: Optional[str] = None
    selected_intensity_column: Optional[str] = None
    selected_area_column: Optional[str] = None
    selected_snr_column: Optional[str] = None
    total_rows_read: int = 0
    valid_points_count: int = 0
    skipped_invalid_rows: int = 0
    duplicate_mz_count: int = 0
    negative_intensity_count: int = 0
    zero_intensity_count: int = 0
    min_mz: float = 0.0
    max_mz: float = 0.0
    max_intensity: float = 0.0
    ms_level: int = 1
    precursor_mz: Optional[float] = None
    precursor_charge: Optional[int] = None
    retention_time_sec: Optional[float] = None
    warnings: Tuple[str, ...] = field(default_factory=tuple)
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def mz_array(self) -> List[float]:
        """Return list of m/z values sorted in ascending order."""
        return [p.mz for p in self.peaks]

    @property
    def intensity_array(self) -> List[float]:
        """Return list of intensity values corresponding to `peaks`."""
        return [p.intensity for p in self.peaks]

    def _serialize_peak(self, idx_1based: int, p: SpectrumPeak, max_int: float) -> Dict[str, Any]:
        rel_pct = (p.intensity / max_int * 100.0) if max_int > 0 else 0.0
        return {
            "index": idx_1based,
            "mz": round(p.mz, 6),
            "intensity": round(p.intensity, 4),
            "relative_intensity_pct": round(rel_pct, 2),
            "area": round(p.area, 4) if p.area is not None else None,
            "snr": round(p.snr, 2) if p.snr is not None else None,
            "fwhm": round(p.fwhm, 5) if p.fwhm is not None else None,
            "charge": p.charge,
        }

    def get_preview_rows(self, limit: int = 20) -> List[Dict[str, Any]]:
        """Return the first `limit` peaks as dictionaries for UI preview."""
        max_int = self.max_intensity if self.max_intensity > 0 else max((pk.intensity for pk in self.peaks), default=0.0)
        return [
            self._serialize_peak(idx + 1, p, max_int)
            for idx, p in enumerate(self.peaks[:limit])
        ]

    def get_plot_rows(self, max_points: int = 600) -> List[Dict[str, Any]]:
        """Return representative peaks across the full m/z span (always preserving top-intensity peaks)."""
        n = len(self.peaks)
        max_int = self.max_intensity if self.max_intensity > 0 else max((pk.intensity for pk in self.peaks), default=0.0)
        if n <= max_points:
            return [self._serialize_peak(idx + 1, p, max_int) for idx, p in enumerate(self.peaks)]

        # Always include the top 120 highest-intensity peaks plus evenly spaced points across m/z
        indexed = list(enumerate(self.peaks))
        top_indices = {
            idx for idx, _ in sorted(indexed, key=lambda item: -item[1].intensity)[:120]
        }
        step = max(1, n // max(1, max_points - len(top_indices)))
        sampled_indices = sorted(top_indices.union(range(0, n, step)))
        return [self._serialize_peak(idx + 1, self.peaks[idx], max_int) for idx in sampled_indices]

    def get_top_intensity_rows(self, limit: int = 50) -> List[Dict[str, Any]]:
        """Return the highest-intensity peaks sorted by intensity descending."""
        max_int = self.max_intensity if self.max_intensity > 0 else max((pk.intensity for pk in self.peaks), default=0.0)
        indexed = list(enumerate(self.peaks))
        indexed.sort(key=lambda item: (-item[1].intensity, item[1].mz))
        return [self._serialize_peak(idx + 1, p, max_int) for idx, p in indexed[:limit]]

    def to_summary_dict(self, preview_limit: int = 20) -> Dict[str, Any]:
        """Serialize spectrum metadata, warnings, and preview points for API/UI inspection."""
        base_peak_mz = 0.0
        total_ion_current = 0.0
        if self.peaks:
            best_pk = max(self.peaks, key=lambda pk: pk.intensity)
            base_peak_mz = round(best_pk.mz, 6)
            total_ion_current = round(sum(pk.intensity for pk in self.peaks), 4)

        return {
            "sample_id": self.sample_id,
            "filename": self.filename,
            "source_format": self.source_format,
            "data_type": self.data_type.value,
            "user_specified_data_type": (
                self.user_specified_data_type.value if self.user_specified_data_type else None
            ),
            "mz_unit": self.mz_unit,
            "intensity_unit": self.intensity_unit,
            "detected_columns": list(self.detected_columns),
            "selected_mz_column": self.selected_mz_column,
            "selected_intensity_column": self.selected_intensity_column,
            "selected_area_column": self.selected_area_column,
            "selected_snr_column": self.selected_snr_column,
            "total_rows_read": self.total_rows_read,
            "valid_points_count": self.valid_points_count,
            "skipped_invalid_rows": self.skipped_invalid_rows,
            "duplicate_mz_count": self.duplicate_mz_count,
            "negative_intensity_count": self.negative_intensity_count,
            "zero_intensity_count": self.zero_intensity_count,
            "min_mz": round(self.min_mz, 6),
            "max_mz": round(self.max_mz, 6),
            "max_intensity": round(self.max_intensity, 4),
            "base_peak_mz": base_peak_mz,
            "total_ion_current": total_ion_current,
            "ms_level": self.ms_level,
            "precursor_mz": self.precursor_mz,
            "precursor_charge": self.precursor_charge,
            "retention_time_sec": self.retention_time_sec,
            "warnings": list(self.warnings),
            "metadata": self.metadata,
            "preview": self.get_preview_rows(limit=preview_limit),
            "top_peaks": self.get_top_intensity_rows(limit=min(preview_limit, 100)),
            "plot_peaks": self.get_plot_rows(max_points=600),
        }


def decode_text_bytes(raw_bytes: bytes) -> Tuple[str, str]:
    """Decode raw bytes trying UTF-8 (with BOM), UTF-16, and Latin-1/CP1252 safely."""
    for enc in ("utf-8-sig", "utf-8", "utf-16", "cp1252", "latin-1"):
        try:
            return raw_bytes.decode(enc), enc
        except UnicodeDecodeError:
            continue
    return raw_bytes.decode("utf-8", errors="replace"), "utf-8-replace"


def infer_spectrum_data_type(
    peaks: Sequence[SpectrumPeak],
    user_specified: Optional[SpectrumDataType] = None,
    file_metadata_hint: Optional[SpectrumDataType] = None,
) -> Tuple[SpectrumDataType, List[str]]:
    """Determine or infer whether a spectrum is CENTROID, PROFILE, or UNKNOWN.

    Respects explicit user specification first, then file metadata (e.g. mzML CV terms),
    and finally inspects point count, median delta m/z spacing, and zero-intensity baseline points.
    """
    notes: List[str] = []
    if user_specified is not None and user_specified != SpectrumDataType.UNKNOWN:
        if file_metadata_hint is not None and file_metadata_hint != SpectrumDataType.UNKNOWN:
            if user_specified != file_metadata_hint:
                notes.append(
                    f"User-specified data type '{user_specified.value}' overrides file metadata hint '{file_metadata_hint.value}'."
                )
        return user_specified, notes

    if file_metadata_hint is not None and file_metadata_hint != SpectrumDataType.UNKNOWN:
        return file_metadata_hint, notes

    n = len(peaks)
    if n < 50:
        return SpectrumDataType.CENTROID, notes

    # Inspect m/z spacing across consecutive sorted peaks
    mzs = [p.mz for p in peaks]
    diffs = [mzs[i + 1] - mzs[i] for i in range(n - 1) if mzs[i + 1] > mzs[i]]
    if not diffs:
        return SpectrumDataType.UNKNOWN, notes

    diffs_sorted = sorted(diffs)
    median_step = diffs_sorted[len(diffs_sorted) // 2]

    # Dense, regular spacing (< 0.08 Da median step with >= 250 points) is characteristic of profile spectra
    if n >= 250 and median_step <= 0.08:
        notes.append(
            f"Auto-classified as 'profile' based on high sampling density ({n} points, median delta m/z = {median_step:.5f} Da)."
        )
        return SpectrumDataType.PROFILE, notes

    if median_step > 0.15:
        return SpectrumDataType.CENTROID, notes

    notes.append(
        f"Data type could not be unambiguously inferred ({n} points, median delta m/z = {median_step:.4f} Da); classified as 'unknown' (treated conservatively as centroid list unless user specifies profile)."
    )
    return SpectrumDataType.UNKNOWN, notes


def finalize_peaks_and_build_spectrum(
    raw_peaks: List[SpectrumPeak],
    sample_id: str,
    filename: str,
    source_format: str,
    user_specified_data_type: Optional[SpectrumDataType] = None,
    file_metadata_data_type: Optional[SpectrumDataType] = None,
    detected_columns: Sequence[str] = (),
    selected_mz_column: Optional[str] = None,
    selected_intensity_column: Optional[str] = None,
    selected_area_column: Optional[str] = None,
    selected_snr_column: Optional[str] = None,
    total_rows_read: int = 0,
    skipped_invalid_rows: int = 0,
    negative_intensity_count: int = 0,
    zero_intensity_count: int = 0,
    duplicate_policy: DuplicateMzPolicy = DuplicateMzPolicy.KEEP_MAX_INTENSITY,
    keep_zero_intensities: Optional[bool] = None,
    ms_level: int = 1,
    precursor_mz: Optional[float] = None,
    precursor_charge: Optional[int] = None,
    retention_time_sec: Optional[float] = None,
    warnings: Optional[List[str]] = None,
    metadata: Optional[Dict[str, Any]] = None,
) -> SpectrumData:
    """Sort peaks, resolve duplicate m/z entries, infer data type, and build `SpectrumData`."""
    warn_list: List[str] = list(warnings) if warnings else []

    if not raw_peaks:
        raise SpectrumParseError("No valid m/z-intensity pairs found")

    # Sort by m/z ascending
    sorted_peaks = sorted(raw_peaks, key=lambda p: p.mz)

    # Determine spectrum data type before deciding whether to retain 0-intensity baseline points
    inferred_type, type_notes = infer_spectrum_data_type(
        sorted_peaks,
        user_specified=user_specified_data_type,
        file_metadata_hint=file_metadata_data_type,
    )
    warn_list.extend(type_notes)

    if keep_zero_intensities is None:
        # Profile spectra need zero/near-zero baseline points for baseline & area integration;
        # centroid peak lists filter out zero-intensity peaks unless all peaks are zero.
        keep_zeros = inferred_type == SpectrumDataType.PROFILE
    else:
        keep_zeros = keep_zero_intensities

    if not keep_zeros:
        non_zero = [p for p in sorted_peaks if p.intensity > 0.0]
        if non_zero:
            if len(non_zero) < len(sorted_peaks):
                dropped_zeros = len(sorted_peaks) - len(non_zero)
                warn_list.append(f"Excluded {dropped_zeros} zero-intensity entry/entries from centroid peak list.")
            sorted_peaks = non_zero
        else:
            raise SpectrumParseError("No valid m/z-intensity pairs found (all intensities are zero or negative)")

    # Handle exact duplicate m/z values (within 1e-9 Da)
    deduped: List[SpectrumPeak] = []
    duplicate_count = 0

    if duplicate_policy == DuplicateMzPolicy.KEEP_ALL:
        # Count duplicates for warning reporting
        for i in range(1, len(sorted_peaks)):
            if abs(sorted_peaks[i].mz - sorted_peaks[i - 1].mz) <= 1e-9:
                duplicate_count += 1
        deduped = sorted_peaks
    else:
        idx = 0
        n = len(sorted_peaks)
        while idx < n:
            group = [sorted_peaks[idx]]
            j = idx + 1
            while j < n and abs(sorted_peaks[j].mz - sorted_peaks[idx].mz) <= 1e-9:
                group.append(sorted_peaks[j])
                j += 1
            if len(group) > 1:
                duplicate_count += len(group) - 1
                if duplicate_policy == DuplicateMzPolicy.KEEP_MAX_INTENSITY:
                    best = max(group, key=lambda p: p.intensity)
                    deduped.append(best)
                elif duplicate_policy == DuplicateMzPolicy.SUM_INTENSITIES:
                    base = group[0]
                    sum_int = sum(p.intensity for p in group)
                    sum_area = (
                        sum(p.area for p in group if p.area is not None)
                        if any(p.area is not None for p in group)
                        else None
                    )
                    deduped.append(
                        SpectrumPeak(
                            mz=base.mz,
                            intensity=sum_int,
                            area=sum_area,
                            snr=base.snr,
                            fwhm=base.fwhm,
                            charge=base.charge,
                            row_index=base.row_index,
                        )
                    )
                elif duplicate_policy == DuplicateMzPolicy.AVERAGE_INTENSITIES:
                    base = group[0]
                    avg_int = sum(p.intensity for p in group) / len(group)
                    deduped.append(
                        SpectrumPeak(
                            mz=base.mz,
                            intensity=avg_int,
                            area=base.area,
                            snr=base.snr,
                            fwhm=base.fwhm,
                            charge=base.charge,
                            row_index=base.row_index,
                        )
                    )
            else:
                deduped.append(group[0])
            idx = j

    if duplicate_count > 0:
        warn_list.append(
            f"Detected {duplicate_count} duplicate m/z value(s); resolved using policy '{duplicate_policy.value}'."
        )

    if negative_intensity_count > 0:
        warn_list.append(f"Skipped {negative_intensity_count} row(s) with negative intensity values.")

    if skipped_invalid_rows > 0:
        warn_list.append(f"Skipped {skipped_invalid_rows} malformed or non-numeric data row(s).")

    min_mz = min(p.mz for p in deduped)
    max_mz = max(p.mz for p in deduped)
    max_int = max(p.intensity for p in deduped)

    return SpectrumData(
        sample_id=sample_id,
        filename=filename,
        source_format=source_format,
        data_type=inferred_type,
        user_specified_data_type=user_specified_data_type,
        peaks=tuple(deduped),
        detected_columns=tuple(detected_columns),
        selected_mz_column=selected_mz_column,
        selected_intensity_column=selected_intensity_column,
        selected_area_column=selected_area_column,
        selected_snr_column=selected_snr_column,
        total_rows_read=total_rows_read,
        valid_points_count=len(deduped),
        skipped_invalid_rows=skipped_invalid_rows,
        duplicate_mz_count=duplicate_count,
        negative_intensity_count=negative_intensity_count,
        zero_intensity_count=zero_intensity_count,
        min_mz=min_mz,
        max_mz=max_mz,
        max_intensity=max_int,
        ms_level=ms_level,
        precursor_mz=precursor_mz,
        precursor_charge=precursor_charge,
        retention_time_sec=retention_time_sec,
        warnings=tuple(warn_list),
        metadata=dict(metadata) if metadata else {},
    )

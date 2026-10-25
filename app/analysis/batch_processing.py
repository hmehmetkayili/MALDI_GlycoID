"""Multi-sample batch processing, safe ZIP extraction, fault isolation, and reproducibility engine.

Implements:
- Safe flat ZIP archive extraction (`extract_batch_zip_safely`) with path-traversal
  protection, extension whitelisting, file count limits, and decompression-bomb guards.
- Fault-isolated multi-sample execution (`execute_batch_analysis`): a failed sample
  never corrupts successfully processed samples (`Batch item failed; other samples completed successfully`).
- Real progress tracking, cancellation support (`BatchCancellationToken`), and retry
  of failed items (`retry_failed_batch_items`).
- Complete combined outputs:
  1. Combined annotated peak table (rank-1 per peak across samples).
  2. Combined glycan candidate table (all ranked candidates across samples).
  3. Per-sample glycan quantification.
  4. Sample-by-glycan abundance matrix & replicate summary.
  5. Sample metadata table.
  6. Batch QC summary table.
  7. Failed-file report.
  8. Timestamped batch processing log.
  9. Combined glycopeptide candidate table (when glycopeptide screening is enabled).
- Full reproducibility manifest (`BatchReproducibilityRecord`).
"""

from __future__ import annotations

import io
import os
import re
import threading
import zipfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple, Union

from app import __version__ as SOFTWARE_VERSION
from app.analysis.glycan_search import (
    MS1CandidateMatch,
    MS1SearchResult,
    MS1SearchSettings,
    SearchMode,
    TargetedQueryItem,
    search_ms1_library,
    search_ms1_targeted,
)
from app.analysis.glycopeptide_search import (
    GlycopeptideCandidateMatch,
    GlycopeptideScreeningResult,
    GlycopeptideSearchSettings,
    PrecursorObservation,
    screen_glycopeptides,
)
from app.analysis.quantification import (
    DenominatorSetPolicy,
    GlycanAggregationPolicy,
    MultiSampleQuantificationMatrix,
    SampleQuantificationResult,
    quantify_ms1_search_result,
    summarize_replicates_and_build_matrix,
)
from app.config import settings as app_settings
from app.libraries.composition_generator import GlycanLibrary
from app.parsers import (
    DuplicateMzPolicy,
    SpectrumData,
    SpectrumDataType,
    parse_spectrum_file,
)


class BatchItemStatus(str, Enum):
    """Execution status of a single sample within a batch job."""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class BatchJobStatus(str, Enum):
    """Overall status of a batch job."""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    COMPLETED_WITH_FAILURES = "completed_with_failures"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass(frozen=True)
class BatchSampleInput:
    """Single sample input specification for a batch analysis run."""

    sample_id: str
    filename: str
    content: Union[str, bytes, Path]
    group: str = "Ungrouped"
    replicate_id: str = ""
    metadata: Mapping[str, Any] = field(default_factory=dict)
    delimiter: Optional[str] = None
    decimal_separator: Optional[str] = None
    mz_column: Optional[Union[str, int]] = None
    intensity_column: Optional[Union[str, int]] = None
    user_specified_data_type: Optional[SpectrumDataType] = None
    duplicate_policy: DuplicateMzPolicy = DuplicateMzPolicy.KEEP_MAX_INTENSITY


@dataclass(frozen=True)
class BatchReproducibilityRecord:
    """Complete provenance and parameter snapshot for a batch analysis run."""

    software_version: str
    library_name: str
    library_version: str
    search_mode: str
    mass_tolerance: float
    mass_tolerance_unit: str
    ionization_polarity: str
    adduct_codes: Tuple[str, ...]
    charge_states: Tuple[int, ...]
    labeling_method: str
    denominator_policy: str
    aggregation_policy: str
    primary_adduct: str
    glycopeptide_screening_enabled: bool
    processing_started_utc: str
    processing_completed_utc: str
    algorithm_settings: Mapping[str, Any]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "software_version": self.software_version,
            "library_name": self.library_name,
            "library_version": self.library_version,
            "search_mode": self.search_mode,
            "mass_tolerance": self.mass_tolerance,
            "mass_tolerance_unit": self.mass_tolerance_unit,
            "ionization_polarity": self.ionization_polarity,
            "adduct_codes": list(self.adduct_codes),
            "charge_states": list(self.charge_states),
            "labeling_method": self.labeling_method,
            "denominator_policy": self.denominator_policy,
            "aggregation_policy": self.aggregation_policy,
            "primary_adduct": self.primary_adduct,
            "glycopeptide_screening_enabled": self.glycopeptide_screening_enabled,
            "processing_started_utc": self.processing_started_utc,
            "processing_completed_utc": self.processing_completed_utc,
            "algorithm_settings": dict(self.algorithm_settings),
        }


@dataclass(frozen=True)
class BatchLogEntry:
    """Timestamped log entry produced during batch execution."""

    timestamp_utc: str
    level: str  # "INFO", "WARNING", "ERROR"
    sample_id: Optional[str]
    filename: Optional[str]
    message: str

    def to_dict(self) -> Dict[str, Any]:
        return {
            "timestamp_utc": self.timestamp_utc,
            "level": self.level,
            "sample_id": self.sample_id or "",
            "filename": self.filename or "",
            "message": self.message,
        }


@dataclass(frozen=True)
class BatchSampleOutcome:
    """Execution outcome and analytical outputs for a single batch item."""

    sample_id: str
    filename: str
    group: str
    replicate_id: str
    status: BatchItemStatus
    error_message: Optional[str] = None
    spectrum: Optional[SpectrumData] = None
    ms1_result: Optional[MS1SearchResult] = None
    quantification: Optional[SampleQuantificationResult] = None
    glycopeptide_result: Optional[GlycopeptideScreeningResult] = None
    warnings: Tuple[str, ...] = field(default_factory=tuple)

    def to_summary_dict(self) -> Dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "filename": self.filename,
            "group": self.group,
            "replicate_id": self.replicate_id,
            "status": self.status.value,
            "error_message": self.error_message,
            "imported_peaks": self.spectrum.valid_points_count if self.spectrum else 0,
            "matched_peaks": self.ms1_result.matched_peaks_count if self.ms1_result else 0,
            "quantified_glycans": len(self.quantification.glycan_abundances) if self.quantification else 0,
            "glycopeptide_matches": (
                self.glycopeptide_result.total_candidate_matches if self.glycopeptide_result else 0
            ),
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class BatchAnalysisResult:
    """Integrated multi-sample batch analysis result."""

    batch_id: str
    status: BatchJobStatus
    total_files: int
    completed_files_count: int
    failed_files_count: int
    cancelled_files_count: int
    progress_pct: float
    sample_outcomes: Tuple[BatchSampleOutcome, ...]
    combined_annotated_peaks: Tuple[MS1CandidateMatch, ...]  # Top rank-1 matches across all samples
    combined_glycan_candidates: Tuple[MS1CandidateMatch, ...]  # All ranked candidates across all samples
    per_sample_quantifications: Tuple[SampleQuantificationResult, ...]
    abundance_matrix: Optional[MultiSampleQuantificationMatrix]
    sample_metadata_table: Tuple[Mapping[str, Any], ...]
    qc_summary_table: Tuple[Mapping[str, Any], ...]
    failed_files_report: Tuple[Mapping[str, Any], ...]
    processing_log: Tuple[BatchLogEntry, ...]
    combined_glycopeptide_candidates: Tuple[GlycopeptideCandidateMatch, ...]
    reproducibility: BatchReproducibilityRecord
    warnings: Tuple[str, ...]

    def to_summary_dict(self) -> Dict[str, Any]:
        return {
            "batch_id": self.batch_id,
            "status": self.status.value,
            "total_files": self.total_files,
            "completed_files_count": self.completed_files_count,
            "failed_files_count": self.failed_files_count,
            "cancelled_files_count": self.cancelled_files_count,
            "progress_pct": round(self.progress_pct, 2),
            "sample_outcomes": [o.to_summary_dict() for o in self.sample_outcomes],
            "combined_annotated_peaks_count": len(self.combined_annotated_peaks),
            "combined_glycan_candidates_count": len(self.combined_glycan_candidates),
            "combined_glycopeptide_candidates_count": len(self.combined_glycopeptide_candidates),
            "abundance_matrix": self.abundance_matrix.to_summary_dict() if self.abundance_matrix else None,
            "sample_metadata_table": [dict(r) for r in self.sample_metadata_table],
            "qc_summary_table": [dict(r) for r in self.qc_summary_table],
            "failed_files_report": [dict(r) for r in self.failed_files_report],
            "processing_log": [l.to_dict() for l in self.processing_log],
            "reproducibility": self.reproducibility.to_dict(),
            "warnings": list(self.warnings),
        }


class BatchCancellationToken:
    """Thread-safe cancellation token for stopping an active batch run."""

    def __init__(self) -> None:
        self._event = threading.Event()

    def cancel(self) -> None:
        self._event.set()

    @property
    def is_cancelled(self) -> bool:
        return self._event.is_set()


_SAFE_FILENAME_RE = re.compile(r"[^A-Za-z0-9._\-]")


def sanitize_filename(raw_name: str) -> str:
    """Sanitize a filename by stripping directory components and unsafe characters."""
    base = Path(raw_name.replace("\\", "/")).name
    sanitized = _SAFE_FILENAME_RE.sub("_", base).strip("._")
    return sanitized or "sample.txt"


def extract_batch_zip_safely(
    zip_content: Union[bytes, Path],
    max_files: int = app_settings.max_batch_files,
    max_total_uncompressed_bytes: int = app_settings.max_upload_size_bytes * 2,
    allowed_extensions: Tuple[str, ...] = app_settings.allowed_spectrum_extensions,
) -> List[BatchSampleInput]:
    """Safely inspect and extract spectrum files from a ZIP archive in memory without writing arbitrary paths.

    Security protections:
    1. Rejects path traversal (`..`, leading `/` or `\\`, Windows drive letters `C:`).
    2. Rejects nested `.zip` archives (no recursive ZIP bomb extraction).
    3. Enforces `max_files` and `max_total_uncompressed_bytes` limits.
    4. Validates file extensions against `allowed_extensions` and sanitizes internal filenames.
    """
    if isinstance(zip_content, Path):
        raw_bytes = zip_content.read_bytes()
    else:
        raw_bytes = zip_content

    try:
        zf = zipfile.ZipFile(io.BytesIO(raw_bytes), mode="r")
    except zipfile.BadZipFile as exc:
        raise ValueError(f"Invalid or corrupted ZIP archive: {exc}") from exc

    extracted_samples: List[BatchSampleInput] = []
    total_bytes = 0
    seen_sample_ids: Dict[str, int] = {}

    with zf:
        for info in zf.infolist():
            if info.is_dir():
                continue

            raw_entry_name = info.filename
            norm_parts = raw_entry_name.replace("\\", "/").split("/")

            # Security check 1: Path traversal or drive prefix
            if (
                ".." in norm_parts
                or raw_entry_name.startswith("/")
                or raw_entry_name.startswith("\\")
                or ":" in raw_entry_name
            ):
                raise ValueError(f"Unsafe path traversal detected in ZIP archive entry: '{raw_entry_name}'")

            base_name = norm_parts[-1]
            # Skip macOS metadata resource forks
            if base_name.startswith("._") or "__MACOSX" in norm_parts:
                continue

            ext = Path(base_name).suffix.lower()
            # Security check 2: Reject nested archives
            if ext in (".zip", ".tar", ".gz", ".7z", ".rar"):
                raise ValueError(f"Nested archives are not permitted inside batch ZIP upload: '{base_name}'")

            if ext not in allowed_extensions:
                continue

            if len(extracted_samples) >= max_files:
                raise ValueError(f"ZIP archive exceeds maximum allowed batch file count ({max_files}).")

            total_bytes += info.file_size
            if total_bytes > max_total_uncompressed_bytes:
                raise ValueError(
                    f"ZIP archive exceeds maximum uncompressed size limit ({max_total_uncompressed_bytes} bytes)."
                )

            file_bytes = zf.read(info.filename)
            if len(file_bytes) > max_total_uncompressed_bytes:
                raise ValueError(f"Extracted size of '{base_name}' exceeds safety limit.")

            safe_name = sanitize_filename(base_name)
            stem = Path(safe_name).stem
            if stem in seen_sample_ids:
                seen_sample_ids[stem] += 1
                unique_id = f"{stem}_{seen_sample_ids[stem]}"
            else:
                seen_sample_ids[stem] = 1
                unique_id = stem

            extracted_samples.append(
                BatchSampleInput(
                    sample_id=unique_id,
                    filename=safe_name,
                    content=file_bytes,
                )
            )

    if not extracted_samples:
        raise ValueError("No supported spectrum files found inside the uploaded ZIP archive.")

    return extracted_samples


def _now_utc_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def execute_batch_analysis(
    samples: Sequence[BatchSampleInput],
    library: GlycanLibrary,
    search_settings: Optional[MS1SearchSettings] = None,
    targeted_items: Optional[Sequence[TargetedQueryItem]] = None,
    denominator_policy: DenominatorSetPolicy = DenominatorSetPolicy.ALL_ELIGIBLE_SPECTRUM_PEAKS,
    aggregation_policy: GlycanAggregationPolicy = GlycanAggregationPolicy.PRIMARY_ADDUCT_ONLY,
    primary_adduct: str = "[M+Na]+",
    exclude_ambiguous_from_aggregation: bool = False,
    enable_glycopeptide_screening: bool = False,
    fasta_or_proteins: Optional[Union[str, bytes, Path]] = None,
    glycopeptide_settings: Optional[GlycopeptideSearchSettings] = None,
    batch_id: Optional[str] = None,
    cancellation_token: Optional[BatchCancellationToken] = None,
    progress_callback: Optional[Callable[[int, int, str], None]] = None,
) -> BatchAnalysisResult:
    """Execute fault-isolated multi-sample batch processing using the validated scientific engine."""
    if not samples:
        raise ValueError("Batch analysis requires at least one sample input.")

    if len(samples) > app_settings.max_batch_files:
        raise ValueError(f"Batch size ({len(samples)}) exceeds maximum limit ({app_settings.max_batch_files}).")

    if search_settings is None:
        search_settings = MS1SearchSettings(search_mode=SearchMode.LIBRARY)

    if batch_id is None:
        batch_id = f"BATCH_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}"

    started_utc = _now_utc_iso()
    logs: List[BatchLogEntry] = [
        BatchLogEntry(
            timestamp_utc=started_utc,
            level="INFO",
            sample_id=None,
            filename=None,
            message=f"Started batch job '{batch_id}' with {len(samples)} sample(s) using library '{library.library_name}' (v{library.library_version}).",
        )
    ]

    outcomes: List[BatchSampleOutcome] = []
    combined_top_peaks: List[MS1CandidateMatch] = []
    combined_all_candidates: List[MS1CandidateMatch] = []
    per_sample_quants: List[SampleQuantificationResult] = []
    combined_glycopeptides: List[GlycopeptideCandidateMatch] = []
    sample_metadata_rows: List[Dict[str, Any]] = []
    qc_summary_rows: List[Dict[str, Any]] = []
    failed_rows: List[Dict[str, Any]] = []
    batch_warnings: List[str] = []

    seen_ids: Dict[str, int] = {}
    total_count = len(samples)

    for idx, item in enumerate(samples, start=1):
        # Ensure unique sample_id across the batch
        base_sid = item.sample_id.strip() or Path(item.filename).stem or f"Sample_{idx}"
        if base_sid in seen_ids:
            seen_ids[base_sid] += 1
            sid = f"{base_sid}_{seen_ids[base_sid]}"
        else:
            seen_ids[base_sid] = 1
            sid = base_sid

        safe_fname = sanitize_filename(item.filename)

        if cancellation_token is not None and cancellation_token.is_cancelled:
            logs.append(
                BatchLogEntry(
                    timestamp_utc=_now_utc_iso(),
                    level="WARNING",
                    sample_id=sid,
                    filename=safe_fname,
                    message="Sample skipped because batch job was cancelled by user.",
                )
            )
            outcomes.append(
                BatchSampleOutcome(
                    sample_id=sid,
                    filename=safe_fname,
                    group=item.group,
                    replicate_id=item.replicate_id,
                    status=BatchItemStatus.CANCELLED,
                    error_message="Cancelled before execution",
                )
            )
            continue

        try:
            # Step 1: Parse spectrum
            spec = parse_spectrum_file(
                content=item.content,
                filename=safe_fname,
                sample_id=sid,
                delimiter=item.delimiter,
                decimal_separator=item.decimal_separator,
                mz_column=item.mz_column,
                intensity_column=item.intensity_column,
                user_specified_data_type=item.user_specified_data_type,
                duplicate_policy=item.duplicate_policy,
            )

            # Step 2: Run MS1 Search (Mode A or Mode B)
            if search_settings.search_mode == SearchMode.TARGETED and targeted_items:
                ms1_res = search_ms1_targeted(
                    spectrum=spec,
                    targets=targeted_items,
                    settings=search_settings,
                    reference_library=library,
                )
            else:
                ms1_res = search_ms1_library(
                    spectrum=spec,
                    library=library,
                    settings=search_settings,
                )

            # Step 3: Run Quantification
            quant_res = quantify_ms1_search_result(
                search_result=ms1_res,
                denominator_policy=denominator_policy,
                aggregation_policy=aggregation_policy,
                primary_adduct=primary_adduct,
                exclude_ambiguous_from_aggregation=exclude_ambiguous_from_aggregation,
            )

            # Step 4: Optional Glycopeptide Screening on unmatched or all peaks
            gp_res: Optional[GlycopeptideScreeningResult] = None
            if enable_glycopeptide_screening and fasta_or_proteins is not None:
                gp_precursors = [
                    PrecursorObservation(
                        sample_id=sid,
                        observed_mz=p.mz,
                        intensity=p.intensity,
                        known_charge=p.charge,
                    )
                    for p in spec.peaks
                ]
                gp_res = screen_glycopeptides(
                    precursors=gp_precursors,
                    fasta_or_proteins=fasta_or_proteins,
                    glycan_compositions=library,
                    settings=glycopeptide_settings,
                )
                combined_glycopeptides.extend(gp_res.matches)

            combined_top_peaks.extend(ms1_res.top_matches_per_peak)
            combined_all_candidates.extend(ms1_res.matches)
            per_sample_quants.append(quant_res)

            # Compute per-sample QC metrics for the batch QC summary table
            ppm_errors = [m.abs_error_ppm for m in ms1_res.top_matches_per_peak]
            signed_ppms = [m.signed_error_ppm for m in ms1_res.top_matches_per_peak]
            import numpy as np

            median_abs_ppm = float(np.median(ppm_errors)) if ppm_errors else None
            mean_signed_ppm = float(np.mean(signed_ppms)) if signed_ppms else None
            match_rate_pct = (
                (ms1_res.matched_peaks_count / ms1_res.eligible_peaks_searched) * 100.0
                if ms1_res.eligible_peaks_searched > 0
                else 0.0
            )

            qc_summary_rows.append(
                {
                    "sample_id": sid,
                    "filename": safe_fname,
                    "group": item.group,
                    "data_type": spec.data_type.value,
                    "imported_peaks": spec.valid_points_count,
                    "eligible_peaks": ms1_res.eligible_peaks_searched,
                    "matched_peaks": ms1_res.matched_peaks_count,
                    "match_rate_pct": round(match_rate_pct, 2),
                    "quantified_glycans": len(quant_res.glycan_abundances),
                    "median_abs_ppm_error": round(median_abs_ppm, 4) if median_abs_ppm is not None else None,
                    "mean_signed_ppm_error": round(mean_signed_ppm, 4) if mean_signed_ppm is not None else None,
                    "status": BatchItemStatus.COMPLETED.value,
                }
            )

            sample_metadata_rows.append(
                {
                    "sample_id": sid,
                    "filename": safe_fname,
                    "group": item.group,
                    "replicate_id": item.replicate_id,
                    "source_format": spec.source_format,
                    "data_type": spec.data_type.value,
                    "status": BatchItemStatus.COMPLETED.value,
                    **dict(item.metadata),
                }
            )

            outcomes.append(
                BatchSampleOutcome(
                    sample_id=sid,
                    filename=safe_fname,
                    group=item.group,
                    replicate_id=item.replicate_id,
                    status=BatchItemStatus.COMPLETED,
                    spectrum=spec,
                    ms1_result=ms1_res,
                    quantification=quant_res,
                    glycopeptide_result=gp_res,
                    warnings=spec.warnings + ms1_res.warnings,
                )
            )
            logs.append(
                BatchLogEntry(
                    timestamp_utc=_now_utc_iso(),
                    level="INFO",
                    sample_id=sid,
                    filename=safe_fname,
                    message=(
                        f"Processed sample '{sid}': {spec.valid_points_count} peaks imported, "
                        f"{ms1_res.matched_peaks_count} matched, {len(quant_res.glycan_abundances)} glycans quantified."
                    ),
                )
            )

        except Exception as exc:
            err_msg = str(exc)
            failed_rows.append(
                {
                    "sample_id": sid,
                    "filename": safe_fname,
                    "group": item.group,
                    "error_message": err_msg,
                    "timestamp_utc": _now_utc_iso(),
                }
            )
            sample_metadata_rows.append(
                {
                    "sample_id": sid,
                    "filename": safe_fname,
                    "group": item.group,
                    "replicate_id": item.replicate_id,
                    "source_format": "FAILED",
                    "data_type": "unknown",
                    "status": BatchItemStatus.FAILED.value,
                    **dict(item.metadata),
                }
            )
            outcomes.append(
                BatchSampleOutcome(
                    sample_id=sid,
                    filename=safe_fname,
                    group=item.group,
                    replicate_id=item.replicate_id,
                    status=BatchItemStatus.FAILED,
                    error_message=err_msg,
                )
            )
            logs.append(
                BatchLogEntry(
                    timestamp_utc=_now_utc_iso(),
                    level="ERROR",
                    sample_id=sid,
                    filename=safe_fname,
                    message=f"Sample failed: {err_msg}",
                )
            )

        if progress_callback is not None:
            progress_callback(idx, total_count, sid)

    # Build combined sample-by-glycan abundance matrix across all successfully quantified samples
    abundance_matrix: Optional[MultiSampleQuantificationMatrix] = None
    if per_sample_quants:
        abundance_matrix = summarize_replicates_and_build_matrix(per_sample_quants)

    completed_cnt = sum(1 for o in outcomes if o.status == BatchItemStatus.COMPLETED)
    failed_cnt = sum(1 for o in outcomes if o.status == BatchItemStatus.FAILED)
    cancelled_cnt = sum(1 for o in outcomes if o.status == BatchItemStatus.CANCELLED)

    if cancelled_cnt > 0:
        overall_status = BatchJobStatus.CANCELLED
        batch_warnings.append(f"Batch execution was cancelled ({completed_cnt} completed, {cancelled_cnt} cancelled).")
    elif failed_cnt > 0 and completed_cnt > 0:
        overall_status = BatchJobStatus.COMPLETED_WITH_FAILURES
        batch_warnings.append("Batch item failed; other samples completed successfully")
    elif failed_cnt > 0 and completed_cnt == 0:
        overall_status = BatchJobStatus.FAILED
        batch_warnings.append("All samples in the batch failed to process.")
    else:
        overall_status = BatchJobStatus.COMPLETED

    completed_utc = _now_utc_iso()
    logs.append(
        BatchLogEntry(
            timestamp_utc=completed_utc,
            level="INFO",
            sample_id=None,
            filename=None,
            message=(
                f"Batch '{batch_id}' finished with status '{overall_status.value}' "
                f"(completed={completed_cnt}, failed={failed_cnt}, cancelled={cancelled_cnt})."
            ),
        )
    )

    repro = BatchReproducibilityRecord(
        software_version=SOFTWARE_VERSION,
        library_name=library.library_name,
        library_version=library.library_version,
        search_mode=search_settings.search_mode.value,
        mass_tolerance=search_settings.tolerance,
        mass_tolerance_unit=search_settings.tolerance_unit.value,
        ionization_polarity=search_settings.polarity.value,
        adduct_codes=search_settings.adduct_codes,
        charge_states=search_settings.charge_states,
        labeling_method=search_settings.labeling,
        denominator_policy=denominator_policy.value,
        aggregation_policy=aggregation_policy.value,
        primary_adduct=primary_adduct,
        glycopeptide_screening_enabled=enable_glycopeptide_screening,
        processing_started_utc=started_utc,
        processing_completed_utc=completed_utc,
        algorithm_settings={
            "min_intensity": search_settings.min_intensity,
            "min_relative_intensity_pct": search_settings.min_relative_intensity_pct,
            "candidate_limit_per_peak": search_settings.candidate_limit_per_peak,
            "enable_isotope_scoring": search_settings.enable_isotope_scoring,
            "exclude_ambiguous_from_aggregation": exclude_ambiguous_from_aggregation,
        },
    )

    prog_pct = ((completed_cnt + failed_cnt) / total_count) * 100.0 if total_count > 0 else 100.0

    return BatchAnalysisResult(
        batch_id=batch_id,
        status=overall_status,
        total_files=total_count,
        completed_files_count=completed_cnt,
        failed_files_count=failed_cnt,
        cancelled_files_count=cancelled_cnt,
        progress_pct=prog_pct,
        sample_outcomes=tuple(outcomes),
        combined_annotated_peaks=tuple(combined_top_peaks),
        combined_glycan_candidates=tuple(combined_all_candidates),
        per_sample_quantifications=tuple(per_sample_quants),
        abundance_matrix=abundance_matrix,
        sample_metadata_table=tuple(sample_metadata_rows),
        qc_summary_table=tuple(qc_summary_rows),
        failed_files_report=tuple(failed_rows),
        processing_log=tuple(logs),
        combined_glycopeptide_candidates=tuple(combined_glycopeptides),
        reproducibility=repro,
        warnings=tuple(batch_warnings),
    )


def retry_failed_batch_items(
    previous_result: BatchAnalysisResult,
    corrected_inputs: Mapping[str, BatchSampleInput],
    library: GlycanLibrary,
    search_settings: Optional[MS1SearchSettings] = None,
    targeted_items: Optional[Sequence[TargetedQueryItem]] = None,
    denominator_policy: DenominatorSetPolicy = DenominatorSetPolicy.ALL_ELIGIBLE_SPECTRUM_PEAKS,
    aggregation_policy: GlycanAggregationPolicy = GlycanAggregationPolicy.PRIMARY_ADDUCT_ONLY,
    primary_adduct: str = "[M+Na]+",
    exclude_ambiguous_from_aggregation: bool = False,
) -> BatchAnalysisResult:
    """Retry failed or cancelled samples from `previous_result` using `corrected_inputs` and merge with existing successes."""
    retry_list: List[BatchSampleInput] = []
    for outcome in previous_result.sample_outcomes:
        if outcome.status in (BatchItemStatus.FAILED, BatchItemStatus.CANCELLED):
            if outcome.sample_id in corrected_inputs:
                retry_list.append(corrected_inputs[outcome.sample_id])

    if not retry_list:
        return previous_result

    retry_res = execute_batch_analysis(
        samples=retry_list,
        library=library,
        search_settings=search_settings,
        targeted_items=targeted_items,
        denominator_policy=denominator_policy,
        aggregation_policy=aggregation_policy,
        primary_adduct=primary_adduct,
        exclude_ambiguous_from_aggregation=exclude_ambiguous_from_aggregation,
        batch_id=f"{previous_result.batch_id}_RETRY",
    )

    retry_by_id = {o.sample_id: o for o in retry_res.sample_outcomes}
    merged_outcomes: List[BatchSampleOutcome] = []
    for old_o in previous_result.sample_outcomes:
        if old_o.sample_id in retry_by_id:
            merged_outcomes.append(retry_by_id[old_o.sample_id])
        else:
            merged_outcomes.append(old_o)

    all_top_peaks: List[MS1CandidateMatch] = []
    all_candidates: List[MS1CandidateMatch] = []
    all_quants: List[SampleQuantificationResult] = []
    all_gps: List[GlycopeptideCandidateMatch] = []
    failed_report: List[Mapping[str, Any]] = []
    meta_rows: List[Mapping[str, Any]] = []
    qc_rows: List[Mapping[str, Any]] = []

    for o in merged_outcomes:
        meta_rows.append(
            {
                "sample_id": o.sample_id,
                "filename": o.filename,
                "group": o.group,
                "replicate_id": o.replicate_id,
                "source_format": o.spectrum.source_format if o.spectrum else "FAILED",
                "data_type": o.spectrum.data_type.value if o.spectrum else "unknown",
                "status": o.status.value,
            }
        )
        if o.status == BatchItemStatus.COMPLETED and o.ms1_result and o.quantification and o.spectrum:
            all_top_peaks.extend(o.ms1_result.top_matches_per_peak)
            all_candidates.extend(o.ms1_result.matches)
            all_quants.append(o.quantification)
            if o.glycopeptide_result:
                all_gps.extend(o.glycopeptide_result.matches)
            qc_rows.append(
                {
                    "sample_id": o.sample_id,
                    "filename": o.filename,
                    "group": o.group,
                    "data_type": o.spectrum.data_type.value,
                    "imported_peaks": o.spectrum.valid_points_count,
                    "eligible_peaks": o.ms1_result.eligible_peaks_searched,
                    "matched_peaks": o.ms1_result.matched_peaks_count,
                    "quantified_glycans": len(o.quantification.glycan_abundances),
                    "status": o.status.value,
                }
            )
        elif o.status == BatchItemStatus.FAILED:
            failed_report.append(
                {
                    "sample_id": o.sample_id,
                    "filename": o.filename,
                    "group": o.group,
                    "error_message": o.error_message or "Unknown error",
                }
            )

    matrix = summarize_replicates_and_build_matrix(all_quants) if all_quants else None
    comp_cnt = sum(1 for o in merged_outcomes if o.status == BatchItemStatus.COMPLETED)
    fail_cnt = sum(1 for o in merged_outcomes if o.status == BatchItemStatus.FAILED)
    canc_cnt = sum(1 for o in merged_outcomes if o.status == BatchItemStatus.CANCELLED)

    status = (
        BatchJobStatus.COMPLETED
        if fail_cnt == 0 and canc_cnt == 0
        else BatchJobStatus.COMPLETED_WITH_FAILURES
    )

    return BatchAnalysisResult(
        batch_id=previous_result.batch_id,
        status=status,
        total_files=len(merged_outcomes),
        completed_files_count=comp_cnt,
        failed_files_count=fail_cnt,
        cancelled_files_count=canc_cnt,
        progress_pct=100.0,
        sample_outcomes=tuple(merged_outcomes),
        combined_annotated_peaks=tuple(all_top_peaks),
        combined_glycan_candidates=tuple(all_candidates),
        per_sample_quantifications=tuple(all_quants),
        abundance_matrix=matrix,
        sample_metadata_table=tuple(meta_rows),
        qc_summary_table=tuple(qc_rows),
        failed_files_report=tuple(failed_report),
        processing_log=previous_result.processing_log + retry_res.processing_log,
        combined_glycopeptide_candidates=tuple(all_gps),
        reproducibility=retry_res.reproducibility,
        warnings=tuple(retry_res.warnings),
    )

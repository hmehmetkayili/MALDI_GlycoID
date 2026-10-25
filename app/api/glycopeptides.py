"""Glycopeptide in-silico digestion, sequon screening, and MS1/MS2 matching API endpoint (Stage 12)."""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from app.analysis.glycopeptide_search import (
    GlycopeptideSearchSettings,
    PrecursorObservation,
    ProteaseEnzyme,
    screen_glycopeptides,
)
from app.api.compat import APIRouter, HTTPException, Request
from app.api.state import session_state
from app.config import EXAMPLES_DIR
from app.core.mass_error import MassToleranceUnit
from app.parsers.base_models import SpectrumData, SpectrumDataType, SpectrumPeak

router = APIRouter(prefix="/api/glycopeptides", tags=["glycopeptides"])


@router.post("/screen")
def screen_glycopeptides_endpoint(request: Request) -> Dict[str, Any]:
    """Execute in-silico protease digestion, N-X-S/T (X!=P) sequon detection, and precursor/MS2 matching."""
    payload = request.json_sync()
    library = session_state.get_or_load_library()

    fasta_text = str(payload.get("fasta_text", "") or "").strip()
    if not fasta_text:
        default_fasta = EXAMPLES_DIR / "synthetic_igg1_fc.fasta"
        fasta_text = default_fasta.read_text(encoding="utf-8")

    try:
        enzyme_str = str(payload.get("enzyme", "trypsin")).strip().lower()
        enzyme = ProteaseEnzyme(enzyme_str)
        max_mc = int(payload.get("max_missed_cleavages", 2))
        tol_ppm = float(payload.get("tolerance_ppm", 20.0))
        fixed_cam = bool(payload.get("fixed_carbamidomethyl_cys", True))

        raw_precursors = payload.get("precursor_mzs", []) or []
        raw_ms2_mzs = payload.get("ms2_fragment_mzs", []) or []

        ms2_spec: Optional[SpectrumData] = None
        if raw_ms2_mzs:
            ms2_peaks = tuple(
                SpectrumPeak(mz=float(mz_val), intensity=1000.0)
                for mz_val in sorted(float(x) for x in raw_ms2_mzs)
            )
            ms2_spec = SpectrumData(
                sample_id="MS2_FRAGMENTS",
                filename="user_ms2_fragments.txt",
                source_format="TXT",
                data_type=SpectrumDataType.CENTROID,
                user_specified_data_type=SpectrumDataType.CENTROID,
                peaks=ms2_peaks,
                valid_points_count=len(ms2_peaks),
                ms_level=2,
            )

        precursors: List[PrecursorObservation] = []
        for idx, mz_val in enumerate(raw_precursors, start=1):
            precursors.append(
                PrecursorObservation(
                    sample_id=f"GP_PREC_{idx}",
                    observed_mz=float(mz_val),
                    intensity=10000.0,
                    ms2_spectrum=ms2_spec,
                )
            )

        gp_settings = GlycopeptideSearchSettings(
            enzyme=enzyme,
            max_missed_cleavages=max_mc,
            fixed_carbamidomethyl_cys=fixed_cam,
            tolerance=tol_ppm,
            tolerance_unit=MassToleranceUnit.PPM,
        )

        result = screen_glycopeptides(
            precursors=precursors,
            fasta_or_proteins=fasta_text,
            glycan_compositions=library,
            settings=gp_settings,
        )
        session_state.glycopeptide_result = result
        out = result.to_summary_dict()
        out["digested_sequon_peptides_count"] = result.sequon_peptides_count
        return out
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Glycopeptide screening error: {exc}") from exc

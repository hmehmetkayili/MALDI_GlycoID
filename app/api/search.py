"""MS1 glycan search and relative abundance quantification API endpoints (Stage 12)."""

from __future__ import annotations

from typing import Any, Dict, List

from app.analysis.glycan_search import (
    MS1SearchSettings,
    SearchMode,
    TargetedQueryItem,
    search_ms1_library,
    search_ms1_targeted,
)
from app.analysis.quantification import (
    DenominatorSetPolicy,
    GlycanAggregationPolicy,
    quantify_ms1_search_result,
)
from app.analysis.statistics import compute_derived_glycan_traits
from app.api.compat import APIRouter, HTTPException, Request
from app.api.state import session_state
from app.core.mass_error import MassToleranceUnit
from app.core.scientific_constants import IonPolarity

router = APIRouter(prefix="/api/search", tags=["search"])


@router.post("/ms1")
def run_ms1_search_endpoint(request: Request) -> Dict[str, Any]:
    """Execute Mode A (Targeted) or Mode B (Library) MS1 glycan search on the active spectrum."""
    payload = request.json_sync()
    spectrum = session_state.get_or_load_spectrum()
    library = session_state.get_or_load_library()

    try:
        mode_str = str(payload.get("search_mode", "library")).strip().lower()
        search_mode = SearchMode.TARGETED if mode_str == "targeted" else SearchMode.LIBRARY

        tol = float(payload.get("tolerance", 20.0))
        unit_str = str(payload.get("tolerance_unit", "ppm")).strip().lower()
        tol_unit = MassToleranceUnit.PPM if unit_str == "ppm" else MassToleranceUnit.DA

        pol_str = str(payload.get("polarity", "positive")).strip().lower()
        polarity = IonPolarity.NEGATIVE if pol_str == "negative" else IonPolarity.POSITIVE

        labeling = str(payload.get("labeling", "NONE")).strip().upper()
        adducts = tuple(str(a).strip() for a in payload.get("adduct_codes", ["[M+Na]+"]) if str(a).strip())
        if not adducts:
            adducts = ("[M+Na]+",)

        charges = tuple(int(z) for z in payload.get("charge_states", [1]))
        if not charges:
            charges = (1,)

        min_rel_int = float(payload.get("min_relative_intensity_pct", 0.0))

        search_settings = MS1SearchSettings(
            search_mode=search_mode,
            tolerance=tol,
            tolerance_unit=tol_unit,
            polarity=polarity,
            adduct_codes=adducts,
            charge_states=charges,
            labeling=labeling,
            min_relative_intensity_pct=min_rel_int,
        )

        if search_mode == SearchMode.TARGETED:
            raw_targets = payload.get("targeted_compositions", []) or []
            if not raw_targets:
                raw_targets = ["Hex5HexNAc2", "Hex3HexNAc4dHex1", "Hex5HexNAc4dHex1", "Hex5HexNAc4Neu5Ac2"]
            targets: List[TargetedQueryItem] = []
            for item in raw_targets:
                s = str(item).strip()
                if not s:
                    continue
                if s.upper().startswith("NG_"):
                    targets.append(TargetedQueryItem(query_label=s, target_type="library_id", library_id=s))
                else:
                    targets.append(TargetedQueryItem(query_label=s, target_type="composition", composition=s))
            result = search_ms1_targeted(
                spectrum=spectrum,
                targets=targets,
                settings=search_settings,
                reference_library=library,
            )
        else:
            result = search_ms1_library(
                spectrum=spectrum,
                library=library,
                settings=search_settings,
            )

        session_state.ms1_search_result = result
        session_state.quantification_result = None
        session_state.derived_traits_report = None
        return result.to_summary_dict()
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"MS1 search error: {exc}") from exc


@router.post("/quantify")
def run_quantification_endpoint(request: Request) -> Dict[str, Any]:
    """Compute base-peak normalization, relative glycan abundance, and derived glycosylation traits."""
    payload = request.json_sync()
    session_state.ensure_baseline_analysis()
    if session_state.ms1_search_result is None:
        raise HTTPException(status_code=400, detail="Run MS1 search before quantification.")

    try:
        denom_str = str(payload.get("denominator_policy", "all_eligible_spectrum_peaks")).strip().lower()
        denom_policy = DenominatorSetPolicy(denom_str)

        agg_str = str(payload.get("aggregation_policy", "primary_adduct_only")).strip().lower()
        agg_policy = GlycanAggregationPolicy(agg_str)

        primary_adduct = str(payload.get("primary_adduct", "[M+Na]+")).strip() or "[M+Na]+"
        exclude_ambig = bool(payload.get("exclude_ambiguous_from_aggregation", False))

        quant_res = quantify_ms1_search_result(
            search_result=session_state.ms1_search_result,
            denominator_policy=denom_policy,
            aggregation_policy=agg_policy,
            primary_adduct=primary_adduct,
            exclude_ambiguous_from_aggregation=exclude_ambig,
        )
        traits_report = compute_derived_glycan_traits(quant_res)

        session_state.quantification_result = quant_res
        session_state.derived_traits_report = traits_report

        return {
            "quantification": quant_res.to_summary_dict(),
            "derived_traits": [
                {
                    "trait_code": t.trait_code,
                    "trait_name": t.trait_name,
                    "value": t.value,
                    "unit": t.unit,
                    "is_composition_inferred_hypothesis": t.is_composition_inferred_hypothesis,
                    "structural_disclaimer": t.structural_disclaimer,
                }
                for t in traits_report.traits
            ],
        }
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Quantification error: {exc}") from exc

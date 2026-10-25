"""MS1 glycan search, isotope-envelope comparison and relative abundance quantification API (Stage 12)."""

from __future__ import annotations

import bisect
import time
from typing import Any, Dict, List, Optional

import numpy as np

from app.analysis.glycan_search import (
    MS1SearchSettings,
    SearchMode,
    TargetedQueryItem,
    build_preprocessing_config,
    search_ms1_library,
    search_ms1_targeted,
)
from app.analysis.peak_processing import preprocess_spectrum
from app.analysis.quantification import (
    DenominatorSetPolicy,
    GlycanAggregationPolicy,
    quantify_ms1_search_result,
)
from app.analysis.statistics import compute_derived_glycan_traits
from app.api.compat import APIRouter, HTTPException, Request
from app.api.qc import run_calibration
from app.api.spectra import run_spectrum_processing
from app.api.state import session_state
from app.core.ion_calculator import calculate_mz_from_neutral_mass
from app.core.isotope_scoring import (
    generate_theoretical_isotope_envelope,
    score_envelope_against_peaks,
)
from app.core.mass_calculator import (
    calculate_released_glycan_mass,
    parse_glycan_composition,
)
from app.core.mass_error import MassToleranceUnit
from app.core.scientific_constants import IonPolarity

router = APIRouter(prefix="/api/search", tags=["search"])

DEFAULT_ADDUCTS = {
    IonPolarity.POSITIVE: ["[M+Na]+"],
    IonPolarity.NEGATIVE: ["[M-H]-"],
}


def _truthy(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "on")
    return bool(value)


def _float_or_zero(value: Any) -> float:
    if value is None or value == "":
        return 0.0
    return float(value)


def build_search_settings(payload: Dict[str, Any], search_mode: SearchMode) -> MS1SearchSettings:
    """Translate an API payload into `MS1SearchSettings` (including polarity/adduct consistency)."""
    tol = float(payload.get("tolerance", 20.0))
    unit_str = str(payload.get("tolerance_unit", "ppm")).strip().lower()
    tol_unit = MassToleranceUnit.PPM if unit_str == "ppm" else MassToleranceUnit.DA

    pol_str = str(payload.get("polarity", "positive")).strip().lower()
    polarity = IonPolarity.NEGATIVE if pol_str == "negative" else IonPolarity.POSITIVE

    labeling = str(payload.get("labeling", "NONE")).strip().upper()
    adducts = tuple(str(a).strip() for a in (payload.get("adduct_codes") or []) if str(a).strip())
    if not adducts:
        adducts = tuple(DEFAULT_ADDUCTS[polarity])

    charges = tuple(int(z) for z in (payload.get("charge_states") or [1]))
    if not charges:
        charges = (1,)

    return MS1SearchSettings(
        search_mode=search_mode,
        tolerance=tol,
        tolerance_unit=tol_unit,
        polarity=polarity,
        adduct_codes=adducts,
        charge_states=charges,
        labeling=labeling,
        min_relative_intensity_pct=float(payload.get("min_relative_intensity_pct", 0.0) or 0.0),
        min_snr=_float_or_zero(payload.get("min_snr")),
        min_score=_float_or_zero(payload.get("min_score")),
        min_idotp=_float_or_zero(payload.get("min_idotp")),
        require_isotope_support=_truthy(payload.get("require_isotope_support", False)),
        exclude_isotope_satellites=_truthy(payload.get("exclude_isotope_satellites", False)),
    )


def _adduct_polarity_warnings(settings: MS1SearchSettings) -> List[str]:
    """Warn when none of the requested adducts belongs to the selected polarity (search would be empty)."""
    from app.core.ion_calculator import resolve_adduct_definition

    notes: List[str] = []
    usable = [a for a in settings.adduct_codes if resolve_adduct_definition(a).polarity == settings.polarity]
    ignored = [a for a in settings.adduct_codes if a not in usable]
    if ignored:
        notes.append(
            f"Adducts {ignored} do not match the selected {settings.polarity.value} polarity and were ignored."
        )
    if not usable:
        notes.append(
            f"No requested adduct matches {settings.polarity.value} polarity; the search cannot return candidates. "
            f"Use e.g. {DEFAULT_ADDUCTS[settings.polarity]}."
        )
    return notes


def run_ms1_search(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Run Mode A (targeted) or Mode B (library) MS1 search on the active spectrum."""
    spectrum = session_state.get_or_load_spectrum()
    library = session_state.get_or_load_library()

    mode_str = str(payload.get("search_mode", "library")).strip().lower()
    search_mode = SearchMode.TARGETED if mode_str == "targeted" else SearchMode.LIBRARY
    search_settings = build_search_settings(payload, search_mode)
    extra_notes = _adduct_polarity_warnings(search_settings)

    t0 = time.perf_counter()
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
    elapsed = time.perf_counter() - t0

    session_state.ms1_search_result = result
    session_state.quantification_result = None
    session_state.derived_traits_report = None

    summary = result.to_summary_dict()
    summary["spectrum_is_calibrated"] = bool(spectrum.metadata.get("is_calibrated", False))
    summary["spectrum_corrections_applied"] = list(spectrum.metadata.get("corrections_applied", []))
    summary["polarity_notes"] = extra_notes
    summary["elapsed_seconds"] = round(elapsed, 4)
    summary["warnings"] = list(summary["warnings"]) + extra_notes
    return summary


@router.post("/ms1")
def run_ms1_search_endpoint(request: Request) -> Dict[str, Any]:
    """Execute Mode A (Targeted) or Mode B (Library) MS1 glycan search on the active spectrum."""
    payload = request.json_sync()
    try:
        return run_ms1_search(payload)
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"MS1 search error: {exc}") from exc


@router.post("/pipeline")
def run_full_pipeline_endpoint(request: Request) -> Dict[str, Any]:
    """Run the recommended workflow in order: optional corrections -> optional calibration -> MS1 search.

    Payload keys: `corrections` (dict, optional), `calibration` (dict, optional), `search` (dict).
    Each stage is skipped when its key is absent / null.
    """
    payload = request.json_sync()
    out: Dict[str, Any] = {"stages": []}
    try:
        t0 = time.perf_counter()
        if payload.get("corrections"):
            proc = run_spectrum_processing(payload["corrections"])
            out["stages"].append("corrections")
            out["corrections"] = proc.get("corrections")
        if payload.get("calibration"):
            out["calibration"] = run_calibration(payload["calibration"])
            out["stages"].append("calibration")
        out["search"] = run_ms1_search(payload.get("search") or {})
        out["stages"].append("search")
        out["total_elapsed_seconds"] = round(time.perf_counter() - t0, 4)
        return out
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Pipeline error: {exc}") from exc


def _formula_string(formula: Dict[str, int]) -> str:
    order = ["C", "H", "N", "O", "Na", "K", "Cl", "P", "S"]
    keys = [k for k in order if k in formula] + sorted(k for k in formula if k not in order)
    return "".join(f"{k}{formula[k] if formula[k] != 1 else ''}" for k in keys if formula[k] != 0)


@router.post("/isotope")
def isotope_comparison_endpoint(request: Request) -> Dict[str, Any]:
    """Theoretical vs experimental isotope envelope for one candidate (opened from the idotp cell).

    Payload: `library_id` or `composition`, `adduct_code`, optional `labeling`, `num_isotopes`,
    `tolerance`, `tolerance_unit`. Search settings of the last MS1 search are used as defaults, so the
    experimental peaks are exactly those that the search saw (same S/N / intensity filters).
    """
    payload = request.json_sync()
    spectrum = session_state.get_or_load_spectrum()
    library = session_state.get_or_load_library()
    last = session_state.ms1_search_result
    last_settings: Optional[MS1SearchSettings] = last.settings if last is not None else None

    try:
        adduct_code = str(payload.get("adduct_code") or (last_settings.adduct_codes[0] if last_settings else "[M+Na]+"))
        labeling = str(payload.get("labeling") or (last_settings.labeling if last_settings else "NONE")).upper()
        tol = float(payload.get("tolerance", last_settings.tolerance if last_settings else 20.0))
        unit_str = str(
            payload.get("tolerance_unit", last_settings.tolerance_unit.value if last_settings else "ppm")
        ).lower()
        tol_unit = MassToleranceUnit.PPM if unit_str == "ppm" else MassToleranceUnit.DA
        n_iso = max(3, min(8, int(payload.get("num_isotopes", 5))))

        comp_str = payload.get("composition")
        lib_id = str(payload.get("library_id") or "").strip().upper()
        if lib_id and not str(lib_id).startswith("TARGET_"):
            entry = next((e for e in library.entries if e.library_id.upper() == lib_id), None)
            if entry is None:
                raise ValueError(f"library_id '{lib_id}' not found in the active library.")
            counts = dict(entry.residue_counts)
            comp_label = entry.composition
        elif comp_str:
            counts = parse_glycan_composition(str(comp_str))
            comp_label = str(comp_str)
        else:
            raise ValueError("Provide 'library_id' or 'composition'.")

        mass_res = calculate_released_glycan_mass(counts, labeling=labeling)
        ion_res = calculate_mz_from_neutral_mass(
            neutral_mass=mass_res.neutral_monoisotopic_mass,
            adduct_code=adduct_code,
            neutral_formula=mass_res.elemental_formula,
        )
        if ion_res.ion_formula is None:
            raise ValueError("Ion elemental formula unavailable; isotope envelope cannot be computed.")

        env = generate_theoretical_isotope_envelope(
            ion_formula=ion_res.ion_formula,
            monoisotopic_mz=ion_res.theoretical_mz,
            charge=ion_res.charge,
            num_peaks=n_iso,
        )

        # Experimental peaks exactly as seen by the search (same preprocessing / S/N filters)
        prep_cfg = build_preprocessing_config(last_settings) if last_settings else None
        proc = preprocess_spectrum(spectrum, config=prep_cfg)
        pool = sorted(proc.isotope_pool_peaks or proc.eligible_peaks, key=lambda p: p.mz)
        pool_mz = np.array([p.mz for p in pool], dtype=float)
        pool_int = np.array([p.intensity for p in pool], dtype=float)

        score = score_envelope_against_peaks(
            theo_envelope=env,
            obs_mz_sorted=pool_mz,
            obs_int_sorted=pool_int,
            tolerance=tol,
            tolerance_unit=tol_unit,
        )
        comparison = score.to_comparison_dict(ion_formula_str=_formula_string(dict(ion_res.ion_formula)))

        # Raw observed trace around the envelope (profile curve or centroid sticks)
        lo_mz = env[0].mz - 1.5 / max(1, abs(ion_res.charge))
        hi_mz = env[-1].mz + 1.5 / max(1, abs(ion_res.charge))
        raw_mz = [p.mz for p in spectrum.peaks]
        i0 = bisect.bisect_left(raw_mz, lo_mz)
        i1 = bisect.bisect_right(raw_mz, hi_mz)
        trace = [
            {"mz": round(p.mz, 5), "intensity": round(p.intensity, 3)} for p in spectrum.peaks[i0:i1]
        ]
        if len(trace) > 1500:
            step = len(trace) // 1500 + 1
            trace = trace[::step]

        comparison.update(
            {
                "composition": comp_label,
                "library_id": lib_id or None,
                "adduct_code": ion_res.adduct_code,
                "charge": ion_res.charge,
                "labeling": labeling,
                "neutral_mass": round(mass_res.neutral_monoisotopic_mass, 6),
                "theoretical_mz": round(ion_res.theoretical_mz, 6),
                "tolerance": tol,
                "tolerance_unit": tol_unit.value,
                "spectrum_data_type": spectrum.data_type.value,
                "spectrum_is_calibrated": bool(spectrum.metadata.get("is_calibrated", False)),
                "observed_trace": trace,
                "observed_trace_type": "profile" if spectrum.data_type.value == "profile" else "centroid",
            }
        )
        return comparison
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Isotope comparison error: {exc}") from exc


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

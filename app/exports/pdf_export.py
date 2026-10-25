"""Multi-page vector PDF analytical summary report generator for MALDI-GlycoID (Stage 10).

Uses `matplotlib.backends.backend_pdf.PdfPages` to produce a self-contained,
publication-grade PDF report containing:
- Page 1: Executive Summary, Reproducibility Parameters, Quality Control Metrics,
  and Mandatory Scientific Disclaimers.
- Page 2: Annotated MS1 Mass Spectrum and Top Quantified N-Glycan Relative
  Abundance Bar Chart.
"""

from __future__ import annotations

from datetime import datetime, timezone
import io
from typing import Any, Mapping, Optional, Sequence

import matplotlib
matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

from app.analysis.glycan_search import MS1SearchResult
from app.analysis.quality_control import SpectrumQCReport
from app.analysis.quantification import SampleQuantificationResult
from app.parsers.base_models import SpectrumData


def export_summary_report_to_pdf_bytes(
    standardized_tables: Mapping[str, Sequence[Mapping[str, Any]]],
    spectrum: Optional[SpectrumData] = None,
    search_result: Optional[MS1SearchResult] = None,
    quantification_result: Optional[SampleQuantificationResult] = None,
    qc_report: Optional[SpectrumQCReport] = None,
    report_title: str = "MALDI-GlycoID — Analytical Summary Report",
) -> bytes:
    """Render a multi-page PDF analytical summary report and return raw PDF bytes."""
    buf = io.BytesIO()
    now_str = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    with PdfPages(buf) as pdf:
        # =====================================================================
        # PAGE 1: Executive Summary, Parameters, QC & Scientific Disclaimers
        # =====================================================================
        fig1 = plt.figure(figsize=(8.5, 11.0), dpi=150)
        fig1.patch.set_facecolor("white")
        ax1 = fig1.add_axes([0.06, 0.05, 0.88, 0.90])
        ax1.axis("off")

        ax1.text(0.0, 0.97, report_title, fontsize=15, fontweight="bold", color="#0f172a")
        ax1.text(0.0, 0.94, f"Generated: {now_str} | Platform: MALDI-GlycoID v1.0.0", fontsize=9, color="#475569")
        ax1.axhline(0.925, color="#cbd5e1", linewidth=1.2)

        # Section 1: Reproducibility Parameters
        ax1.text(0.0, 0.895, "1. Analysis Parameters & Reproducibility", fontsize=11, fontweight="bold", color="#0284c7")
        param_rows = standardized_tables.get("analysis_parameters_and_reproducibility", [])
        y_cursor = 0.865
        if param_rows:
            for k, v in list(param_rows[0].items())[:12]:
                ax1.text(0.02, y_cursor, f"{k}:", fontsize=8.5, fontweight="bold", color="#1e293b")
                ax1.text(0.38, y_cursor, f"{v}", fontsize=8.5, color="#334155")
                y_cursor -= 0.022
        else:
            ax1.text(0.02, y_cursor, "Default MS1 Library Search Parameters applied.", fontsize=8.5, color="#334155")
            y_cursor -= 0.025

        # Section 2: Quality Control Summary
        y_cursor -= 0.015
        ax1.text(0.0, y_cursor, "2. Quality Control Summary", fontsize=11, fontweight="bold", color="#0284c7")
        y_cursor -= 0.028
        qc_rows = standardized_tables.get("quality_control_summary", [])
        if qc_report is not None:
            qc_items = [
                ("Sample ID", qc_report.sample_id),
                ("Overall QC Status", qc_report.overall_status.value),
                ("Data Type", qc_report.data_type),
                ("Imported / Eligible / Matched Peaks", f"{qc_report.imported_peaks_count} / {qc_report.eligible_peaks_count} / {qc_report.matched_peaks_count} ({qc_report.match_rate_pct:.1f}%)"),
                ("Quantified Glycan Compositions", str(qc_report.quantified_glycans_count)),
                ("Median Abs Error (ppm)", f"{qc_report.median_abs_ppm_error:.3f}" if qc_report.median_abs_ppm_error is not None else "N/A"),
                ("Mean Signed Error (ppm)", f"{qc_report.mean_signed_ppm_error:+.3f}" if qc_report.mean_signed_ppm_error is not None else "N/A"),
                ("Annotated Intensity Fraction (%)", f"{qc_report.annotated_intensity_fraction_pct:.2f}%"),
                ("Calibrated", "Yes" if qc_report.is_calibrated else "No (Uncalibrated)"),
            ]
            for label, val in qc_items:
                ax1.text(0.02, y_cursor, f"{label}:", fontsize=8.5, fontweight="bold", color="#1e293b")
                ax1.text(0.38, y_cursor, str(val), fontsize=8.5, color="#334155")
                y_cursor -= 0.022
        elif qc_rows:
            for k, v in list(qc_rows[0].items())[:9]:
                ax1.text(0.02, y_cursor, f"{k}:", fontsize=8.5, fontweight="bold", color="#1e293b")
                ax1.text(0.38, y_cursor, f"{v}", fontsize=8.5, color="#334155")
                y_cursor -= 0.022
        else:
            ax1.text(0.02, y_cursor, "No QC summary rows recorded.", fontsize=8.5, color="#475569")
            y_cursor -= 0.022

        # Section 3: Top Quantified N-Glycan Candidates Table
        y_cursor -= 0.018
        ax1.text(0.0, y_cursor, "3. Top Quantified N-Glycan Compositions", fontsize=11, fontweight="bold", color="#0284c7")
        y_cursor -= 0.026
        quant_rows = standardized_tables.get("quantified_glycans", [])
        if quant_rows:
            header_str = f"{'Composition':<28} {'Class':<15} {'Adducts':<16} {'BasePeak %':>12} {'Rel.Abund %':>13}"
            ax1.text(0.02, y_cursor, header_str, fontsize=8, fontweight="bold", family="monospace", color="#0f172a")
            y_cursor -= 0.018
            for r in quant_rows[:12]:
                comp_s = str(r.get("glycan_composition", ""))[:26]
                cls_s = str(r.get("candidate_class", ""))[:14]
                add_s = str(r.get("contributing_adducts", ""))[:15]
                bp_s = f"{float(r.get('base_peak_relative_intensity_pct', 0.0)):.2f}%"
                ra_s = f"{float(r.get('relative_abundance_intensity_pct', 0.0)):.2f}%"
                line_s = f"{comp_s:<28} {cls_s:<15} {add_s:<16} {bp_s:>12} {ra_s:>13}"
                ax1.text(0.02, y_cursor, line_s, fontsize=7.8, family="monospace", color="#334155")
                y_cursor -= 0.017
        else:
            ax1.text(0.02, y_cursor, "No quantified glycans available.", fontsize=8.5, color="#475569")
            y_cursor -= 0.022

        # Section 4: Mandatory Scientific Interpretation Disclaimers
        ax1.axhline(0.155, color="#cbd5e1", linewidth=0.9)
        ax1.text(0.0, 0.13, "4. Scientific Scope & Interpretation Disclaimers", fontsize=9.5, fontweight="bold", color="#b91c1c")
        disclaimers = [
            "• MS1 mass matching identifies candidate monosaccharide compositions and adducts; structural isomers (branching,",
            "  core vs. antennary fucosylation, and glycosidic linkages without linkage-specific derivatization) require MS/MS.",
            "• Relative abundances (% of total annotated glycan signal) are compositional metrics, not absolute molar concentrations.",
            "• Candidate ranking scores are heuristic combinations of mass accuracy, intensity, and idotp (not FDR probabilities).",
        ]
        d_y = 0.105
        for d_line in disclaimers:
            ax1.text(0.01, d_y, d_line, fontsize=7.8, color="#334155")
            d_y -= 0.020

        pdf.savefig(fig1)
        plt.close(fig1)

        # =====================================================================
        # PAGE 2: Annotated Spectrum & Relative Abundance Bar Chart
        # =====================================================================
        fig2, (ax_spec, ax_bar) = plt.subplots(2, 1, figsize=(8.5, 11.0), dpi=150)
        fig2.patch.set_facecolor("white")

        # Top Plot: Annotated MS1 Spectrum
        if spectrum is not None and spectrum.peaks:
            mzs = [p.mz for p in spectrum.peaks]
            ints = [p.intensity for p in spectrum.peaks]
            ax_spec.vlines(mzs, 0, ints, color="#94a3b8", linewidth=0.9, label="All Observed Peaks")
            if search_result is not None and search_result.top_matches_per_peak:
                m_mzs = [m.observed_mz for m in search_result.top_matches_per_peak]
                m_ints = [m.observed_intensity for m in search_result.top_matches_per_peak]
                ax_spec.vlines(m_mzs, 0, m_ints, color="#0284c7", linewidth=1.5, label="Annotated N-Glycans")
                top_peaks = sorted(search_result.top_matches_per_peak, key=lambda m: -m.observed_intensity)[:6]
                for m in top_peaks:
                    ax_spec.annotate(
                        f"{m.glycan_composition}\n({m.observed_mz:.2f})",
                        xy=(m.observed_mz, m.observed_intensity),
                        xytext=(0, 5),
                        textcoords="offset points",
                        ha="center",
                        fontsize=7,
                    )
            ax_spec.set_ylim(bottom=0, top=max(ints) * 1.25)
            ax_spec.legend(loc="upper right", fontsize=8)
        else:
            ax_spec.text(0.5, 0.5, "Spectrum plot not provided in export bundle", ha="center", va="center")
        ax_spec.set_title("Annotated MALDI-TOF MS1 Spectrum", fontsize=11, fontweight="bold")
        ax_spec.set_xlabel("m/z (Da/e)")
        ax_spec.set_ylabel("Intensity (a.u.)")
        ax_spec.grid(True, linestyle="--", alpha=0.3)

        # Bottom Plot: Relative Abundance Bar Chart
        quant_rows_sorted = sorted(
            quant_rows,
            key=lambda r: float(r.get("relative_abundance_intensity_pct", 0.0)),
            reverse=True,
        )[:12]
        if quant_rows_sorted:
            labels = [str(r.get("glycan_composition", "")) for r in quant_rows_sorted]
            vals = [float(r.get("relative_abundance_intensity_pct", 0.0)) for r in quant_rows_sorted]
            y_pos = list(range(len(labels)))
            ax_bar.barh(y_pos, vals, color="#0284c7", alpha=0.85)
            ax_bar.set_yticks(y_pos)
            ax_bar.set_yticklabels(labels, fontsize=8)
            ax_bar.invert_yaxis()
            ax_bar.set_xlabel("Relative Abundance (% of Total Annotated Glycan Signal)")
            ax_bar.set_title("Top Quantified N-Glycan Relative Abundances", fontsize=11, fontweight="bold")
            ax_bar.grid(True, axis="x", linestyle="--", alpha=0.3)
        else:
            ax_bar.text(0.5, 0.5, "No quantified glycans available for bar chart", ha="center", va="center")

        fig2.tight_layout(pad=2.5)
        pdf.savefig(fig2)
        plt.close(fig2)

    return buf.getvalue()

/**
 * MALDI-GlycoID — Frontend Workflow Controller (Stage 11)
 *
 * Connects all 12 navigation views to the FastAPI scientific backend endpoints:
 * - /api/health
 * - /api/spectra/upload, /api/spectra/example
 * - /api/libraries/default, /api/libraries/generate, /api/libraries/query
 * - /api/search/ms1, /api/search/quantify
 * - /api/glycopeptides/screen
 * - /api/batch/run, /api/batch/demo
 * - /api/qc/calibrate, /api/qc/evaluate, /api/qc/statistics-demo
 * - /api/exports/download
 */

(function () {
  "use strict";

  const appState = {
    activeSampleId: null,
    spectrumSummary: null,
    libraryEntries: [],
    ms1SearchResult: null,
    quantificationResult: null,
    derivedTraits: [],
    glycopeptideResult: null,
    calibrationResult: null,
    qcReport: null,
    batchResult: null,
  };

  function showAlert(message, isError = false) {
    const banner = document.getElementById("global-alert-banner");
    const textEl = document.getElementById("global-alert-text");
    if (!banner || !textEl) return;
    textEl.textContent = message;
    banner.classList.remove("hidden", "error");
    if (isError) {
      banner.classList.add("error");
    }
  }

  function hideAlert() {
    const banner = document.getElementById("global-alert-banner");
    if (banner) banner.classList.add("hidden");
  }

  function switchView(viewName) {
    document.querySelectorAll(".nav-item").forEach((btn) => {
      btn.classList.toggle("active", btn.getAttribute("data-view") === viewName);
    });
    document.querySelectorAll(".view-section").forEach((sec) => {
      sec.classList.toggle("active", sec.id === `view-${viewName}`);
    });
  }

  async function apiRequest(url, options = {}) {
    const resp = await fetch(url, options);
    const contentType = resp.headers.get("content-type") || "";
    if (!resp.ok) {
      let errMsg = `HTTP ${resp.status}`;
      if (contentType.includes("application/json")) {
        const errData = await resp.json();
        errMsg = errData.detail || JSON.stringify(errData);
      } else {
        errMsg = await resp.text();
      }
      throw new Error(errMsg);
    }
    if (contentType.includes("application/json")) {
      return await resp.json();
    }
    return resp;
  }

  // Detect isotopic [M+1] satellite peaks (approx +1.003 Da from a stronger preceding peak)
  function classifyPeaks(peaks, matchedPeaks = []) {
    const sorted = [...peaks].sort((a, b) => Number(a.mz) - Number(b.mz));
    const maxInt = Math.max(...sorted.map((p) => Number(p.intensity)), 1);
    const matchedMap = new Map();
    (matchedPeaks || []).forEach((m) => {
      matchedMap.set(Number(m.observed_mz).toFixed(3), m);
    });

    const classified = [];
    for (let i = 0; i < sorted.length; i++) {
      const pk = sorted[i];
      const mz = Number(pk.mz);
      const intensity = Number(pk.intensity);
      const relPct =
        pk.relative_intensity_pct !== undefined && pk.relative_intensity_pct !== null
          ? Number(pk.relative_intensity_pct)
          : (intensity / maxInt) * 100.0;

      let isIsotopeSatellite = false;
      let parentMz = null;
      for (let j = i - 1; j >= 0; j--) {
        const prevMz = Number(sorted[j].mz);
        const diff = mz - prevMz;
        if (diff > 1.15) break;
        if (diff >= 0.92 && diff <= 1.09 && Number(sorted[j].intensity) > intensity * 0.35) {
          isIsotopeSatellite = true;
          parentMz = prevMz;
          break;
        }
      }

      const matchObj = matchedMap.get(mz.toFixed(3)) || null;
      const isBasePeak = Math.abs(intensity - maxInt) < 1e-6 && maxInt > 0;

      classified.push({
        ...pk,
        mz,
        intensity,
        relPct,
        isBasePeak,
        isIsotopeSatellite,
        parentMz,
        matchObj,
      });
    }
    return classified;
  }

  // Render interactive publication-grade SVG spectrum viewer
  function renderSpectrumSvg(rawPlotRows, matchedPeaks = []) {
    const svg = document.getElementById("spectrum-preview-svg");
    const tooltip = document.getElementById("spectrum-hover-tooltip");
    if (!svg) return;
    svg.innerHTML = "";

    const sourceRows =
      rawPlotRows && rawPlotRows.length > 0
        ? rawPlotRows
        : appState.spectrumSummary
        ? appState.spectrumSummary.plot_peaks || appState.spectrumSummary.preview || []
        : [];
    if (!sourceRows || sourceRows.length === 0) return;

    const plotStyleSel = document.getElementById("spec-plot-style")?.value || "auto";
    const labelMode = document.getElementById("spec-label-mode")?.value || "smart";
    const zoomMinVal = parseFloat(document.getElementById("spec-zoom-min")?.value || "");
    const zoomMaxVal = parseFloat(document.getElementById("spec-zoom-max")?.value || "");

    const allClassified = classifyPeaks(sourceRows, matchedPeaks);
    const filteredRows = allClassified.filter((p) => {
      if (!Number.isNaN(zoomMinVal) && p.mz < zoomMinVal) return false;
      if (!Number.isNaN(zoomMaxVal) && p.mz > zoomMaxVal) return false;
      return true;
    });

    const plotRows = filteredRows.length > 0 ? filteredRows : allClassified;

    const width = 880;
    const height = 340;
    const padLeft = 58;
    const padRight = 26;
    const padTop = 36;
    const padBottom = 42;
    const plotW = width - padLeft - padRight;
    const plotH = height - padTop - padBottom;

    const mzs = plotRows.map((p) => p.mz);
    const ints = plotRows.map((p) => p.intensity);
    const rawMinMz = Math.min(...mzs);
    const rawMaxMz = Math.max(...mzs);
    const spanMz = Math.max(rawMaxMz - rawMinMz, 20);
    const minMz = !Number.isNaN(zoomMinVal) ? zoomMinVal : Math.max(0, rawMinMz - spanMz * 0.05);
    const maxMz = !Number.isNaN(zoomMaxVal) ? zoomMaxVal : rawMaxMz + spanMz * 0.05;
    const maxInt = Math.max(...ints, 1);

    const ns = "http://www.w3.org/2000/svg";

    // Background plot area
    const bgRect = document.createElementNS(ns, "rect");
    bgRect.setAttribute("x", String(padLeft));
    bgRect.setAttribute("y", String(padTop));
    bgRect.setAttribute("width", String(plotW));
    bgRect.setAttribute("height", String(plotH));
    bgRect.setAttribute("fill", "#f8fafc");
    bgRect.setAttribute("rx", "4");
    svg.appendChild(bgRect);

    // Horizontal Y-axis gridlines & percentage ticks (0%, 25%, 50%, 75%, 100%)
    [0, 25, 50, 75, 100].forEach((pct) => {
      const y = padTop + plotH - (pct / 100) * plotH;
      const gLine = document.createElementNS(ns, "line");
      gLine.setAttribute("x1", String(padLeft));
      gLine.setAttribute("y1", String(y));
      gLine.setAttribute("x2", String(width - padRight));
      gLine.setAttribute("y2", String(y));
      gLine.setAttribute("stroke", pct === 0 ? "#64748b" : "#e2e8f0");
      gLine.setAttribute("stroke-width", pct === 0 ? "1.4" : "1");
      if (pct > 0 && pct < 100) {
        gLine.setAttribute("stroke-dasharray", "3,3");
      }
      svg.appendChild(gLine);

      const yLbl = document.createElementNS(ns, "text");
      yLbl.setAttribute("x", String(padLeft - 8));
      yLbl.setAttribute("y", String(y + 3.5));
      yLbl.setAttribute("text-anchor", "end");
      yLbl.setAttribute("font-size", "10");
      yLbl.setAttribute("fill", "#64748b");
      yLbl.textContent = `${pct}%`;
      svg.appendChild(yLbl);
    });

    // Vertical X-axis m/z ticks
    const numXTicks = 7;
    for (let i = 0; i <= numXTicks; i++) {
      const tickMz = minMz + (i / numXTicks) * (maxMz - minMz);
      const x = padLeft + (i / numXTicks) * plotW;

      const tickLine = document.createElementNS(ns, "line");
      tickLine.setAttribute("x1", String(x));
      tickLine.setAttribute("y1", String(padTop + plotH));
      tickLine.setAttribute("x2", String(x));
      tickLine.setAttribute("y2", String(padTop + plotH + 5));
      tickLine.setAttribute("stroke", "#64748b");
      tickLine.setAttribute("stroke-width", "1");
      svg.appendChild(tickLine);

      const xLbl = document.createElementNS(ns, "text");
      xLbl.setAttribute("x", String(x));
      xLbl.setAttribute("y", String(padTop + plotH + 17));
      xLbl.setAttribute("text-anchor", "middle");
      xLbl.setAttribute("font-size", "10");
      xLbl.setAttribute("fill", "#475569");
      xLbl.textContent = tickMz.toFixed(0);
      svg.appendChild(xLbl);
    }

    // Axis titles
    const xAxisTitle = document.createElementNS(ns, "text");
    xAxisTitle.setAttribute("x", String(padLeft + plotW / 2));
    xAxisTitle.setAttribute("y", String(height - 6));
    xAxisTitle.setAttribute("text-anchor", "middle");
    xAxisTitle.setAttribute("font-size", "11");
    xAxisTitle.setAttribute("font-weight", "600");
    xAxisTitle.setAttribute("fill", "#334155");
    xAxisTitle.textContent = "m/z (Da/e)";
    svg.appendChild(xAxisTitle);

    const yAxisTitle = document.createElementNS(ns, "text");
    yAxisTitle.setAttribute("transform", `translate(14, ${padTop + plotH / 2}) rotate(-90)`);
    yAxisTitle.setAttribute("text-anchor", "middle");
    yAxisTitle.setAttribute("font-size", "11");
    yAxisTitle.setAttribute("font-weight", "600");
    yAxisTitle.setAttribute("fill", "#334155");
    yAxisTitle.textContent = "Relative Intensity (%)";
    svg.appendChild(yAxisTitle);

    // Determine effective plot style
    const isProfileMode =
      plotStyleSel === "profile" ||
      (plotStyleSel === "auto" &&
        appState.spectrumSummary &&
        appState.spectrumSummary.data_type === "profile" &&
        plotRows.length > 60);

    const mapX = (mz) => padLeft + ((mz - minMz) / Math.max(maxMz - minMz, 1e-6)) * plotW;
    const mapY = (intVal) => padTop + plotH - (intVal / maxInt) * plotH;

    if (isProfileMode && plotRows.length >= 2) {
      // Shaded area + continuous polyline for profile trace
      const polyPoints = plotRows.map((p) => `${mapX(p.mz).toFixed(1)},${mapY(p.intensity).toFixed(1)}`).join(" ");
      const areaPoints = `${mapX(plotRows[0].mz).toFixed(1)},${padTop + plotH} ${polyPoints} ${mapX(
        plotRows[plotRows.length - 1].mz
      ).toFixed(1)},${padTop + plotH}`;

      const areaPoly = document.createElementNS(ns, "polygon");
      areaPoly.setAttribute("points", areaPoints);
      areaPoly.setAttribute("fill", "rgba(37, 99, 235, 0.12)");
      svg.appendChild(areaPoly);

      const linePoly = document.createElementNS(ns, "polyline");
      linePoly.setAttribute("points", polyPoints);
      linePoly.setAttribute("fill", "none");
      linePoly.setAttribute("stroke", "#2563eb");
      linePoly.setAttribute("stroke-width", "1.6");
      svg.appendChild(linePoly);
    }

    // Draw centroid sticks & apex lollipop circles
    plotRows.forEach((pk) => {
      const x = mapX(pk.mz);
      const y = mapY(pk.intensity);

      let strokeColor = "#1e293b";
      let strokeW = "1.7";
      let dash = "";
      if (pk.isBasePeak) {
        strokeColor = "#059669";
        strokeW = "2.5";
      } else if (pk.matchObj) {
        strokeColor = "#2563eb";
        strokeW = "2.3";
      } else if (pk.isIsotopeSatellite) {
        strokeColor = "#94a3b8";
        strokeW = "1.3";
        dash = "2,2";
      }

      if (!isProfileMode || pk.matchObj || pk.isBasePeak) {
        const stick = document.createElementNS(ns, "line");
        stick.setAttribute("x1", String(x));
        stick.setAttribute("y1", String(padTop + plotH));
        stick.setAttribute("x2", String(x));
        stick.setAttribute("y2", String(y));
        stick.setAttribute("stroke", strokeColor);
        stick.setAttribute("stroke-width", strokeW);
        if (dash) stick.setAttribute("stroke-dasharray", dash);
        svg.appendChild(stick);
      }

      // Apex circle for monoisotopic/matched/base peaks (or all sticks when count <= 100)
      if (!pk.isIsotopeSatellite || plotRows.length <= 40) {
        const dot = document.createElementNS(ns, "circle");
        dot.setAttribute("cx", String(x));
        dot.setAttribute("cy", String(y));
        dot.setAttribute("r", pk.isBasePeak || pk.matchObj ? "3.8" : pk.isIsotopeSatellite ? "2.2" : "3.0");
        dot.setAttribute("fill", strokeColor);
        dot.setAttribute("stroke", "#ffffff");
        dot.setAttribute("stroke-width", "1");
        dot.style.cursor = "pointer";

        // Hover tooltip interaction
        const showTip = (evt) => {
          if (!tooltip) return;
          const clsLabel = pk.isBasePeak
            ? "Base Peak (100%)"
            : pk.matchObj
            ? `Annotated: ${pk.matchObj.glycan_composition} (${pk.matchObj.assigned_adduct})`
            : pk.isIsotopeSatellite
            ? `Isotopic [M+1] Satellite (of ${pk.parentMz.toFixed(4)})`
            : "Monoisotopic Peak";

          tooltip.innerHTML = `
            <div><strong>Peak #${pk.index}</strong> &mdash; <code>${pk.mz.toFixed(4)} Da/e</code></div>
            <div>Intensity: <strong>${pk.intensity.toLocaleString(undefined, { maximumFractionDigits: 1 })}</strong> (${pk.relPct.toFixed(2)}%)</div>
            <div style="color:#93c5fd;margin-top:2px;">${clsLabel}</div>
          `;
          tooltip.classList.remove("hidden");
          const containerRect = svg.parentElement.getBoundingClientRect();
          const leftPx = Math.min(
            Math.max(10, evt.clientX - containerRect.left + 12),
            Math.max(20, containerRect.width - 230)
          );
          const topPx = Math.max(10, evt.clientY - containerRect.top - 58);
          tooltip.style.left = `${leftPx}px`;
          tooltip.style.top = `${topPx}px`;
        };

        const hideTip = () => {
          if (tooltip) tooltip.classList.add("hidden");
        };

        dot.addEventListener("mouseenter", showTip);
        dot.addEventListener("mousemove", showTip);
        dot.addEventListener("mouseleave", hideTip);
        svg.appendChild(dot);
      }
    });

    // Smart non-overlapping peak labels (sorted by intensity descending, skipping [M+1] satellites)
    if (labelMode !== "none") {
      const candidatesForLabel = [...plotRows]
        .filter((p) => {
          if (labelMode === "matched") return Boolean(p.matchObj || p.isBasePeak);
          if (p.isIsotopeSatellite) return false;
          return p.relPct >= 6.0 || Boolean(p.matchObj);
        })
        .sort((a, b) => b.intensity - a.intensity);

      const placedLabelPositions = [];
      const minHorizontalGapPx = 54;

      candidatesForLabel.forEach((pk) => {
        const x = mapX(pk.mz);
        const y = mapY(pk.intensity);

        const collides = placedLabelPositions.some((px) => Math.abs(px - x) < minHorizontalGapPx);
        if (collides && !pk.isBasePeak) return;

        placedLabelPositions.push(x);

        const labelY = Math.max(padTop + 10, y - 8);
        const mzText = document.createElementNS(ns, "text");
        mzText.setAttribute("x", String(x));
        mzText.setAttribute("y", String(labelY));
        mzText.setAttribute("text-anchor", "middle");
        mzText.setAttribute("font-size", pk.isBasePeak ? "10.5" : "9.5");
        mzText.setAttribute("font-weight", pk.isBasePeak || pk.matchObj ? "700" : "600");
        mzText.setAttribute("fill", pk.isBasePeak ? "#047857" : pk.matchObj ? "#1d4ed8" : "#0f172a");
        mzText.textContent = pk.mz.toFixed(2);
        svg.appendChild(mzText);

        if (pk.matchObj) {
          const compText = document.createElementNS(ns, "text");
          compText.setAttribute("x", String(x));
          compText.setAttribute("y", String(Math.max(padTop - 2, labelY - 11)));
          compText.setAttribute("text-anchor", "middle");
          compText.setAttribute("font-size", "8.5");
          compText.setAttribute("font-weight", "700");
          compText.setAttribute("fill", "#2563eb");
          compText.textContent = pk.matchObj.glycan_composition;
          svg.appendChild(compText);
        }
      });
    }
  }

  // Render Raw .TXT / .XY / .CSV File Structure Inspector with syntax highlighting
  function renderRawTxtInspector(summary) {
    const panel = document.getElementById("raw-txt-preview-panel");
    const codeEl = document.getElementById("raw-txt-code");
    const badgeEl = document.getElementById("raw-txt-format-badge");
    if (!panel || !codeEl) return;

    const rawLines = summary.raw_text_lines || [];
    if (rawLines.length === 0) {
      panel.style.display = "none";
      return;
    }

    panel.style.display = "block";
    if (badgeEl) {
      badgeEl.textContent = `${summary.filename} (${summary.source_format.toUpperCase()} — ${summary.data_type.toUpperCase()})`;
    }

    const escapeHtml = (str) =>
      String(str)
        .replace(/&/g, "&amp;")
        .replace(/</g, "&lt;")
        .replace(/>/g, "&gt;");

    let headerSeen = false;
    const formattedLines = rawLines.map((line, idx) => {
      const lineNo = `<span class="txt-line-num">${String(idx + 1).padStart(2, " ")}</span>`;
      const trimmed = line.trim();
      if (!trimmed) {
        return `${lineNo}`;
      }
      if (trimmed.startsWith("#") || trimmed.startsWith("//") || trimmed.startsWith(";")) {
        return `${lineNo}<span class="txt-line-comment">${escapeHtml(line)}</span>`;
      }
      const hasLetters = /[a-zA-Z]/.test(trimmed) && !/^[0-9.,eE+-\s\t;]+$/.test(trimmed);
      if (hasLetters && !headerSeen) {
        headerSeen = true;
        return `${lineNo}<span class="txt-line-header">${escapeHtml(line)}</span>`;
      }
      const parts = trimmed.split(/[\s\t,;]+/);
      if (parts.length >= 2) {
        const mzPart = escapeHtml(parts[0]);
        const intPart = escapeHtml(parts[1]);
        const restPart = parts.length > 2 ? "   " + escapeHtml(parts.slice(2).join("  ")) : "";
        return `${lineNo}<span class="txt-line-mz">${mzPart.padEnd(14, " ")}</span><span class="txt-line-int">${intPart}</span>${restPart}`;
      }
      return `${lineNo}${escapeHtml(line)}`;
    });

    codeEl.innerHTML = formattedLines.join("\n");
  }

  // Render enriched Peak List Table with inline relative intensity bar & classification badges
  function renderImportPeakTable() {
    const summary = appState.spectrumSummary;
    const tbody = document.querySelector("#table-import-preview tbody");
    const countBadge = document.getElementById("import-table-count-badge");
    if (!tbody || !summary) return;

    const topMatches = appState.ms1SearchResult ? appState.ms1SearchResult.top_matches_per_peak || [] : [];
    const classified = classifyPeaks(summary.preview || [], topMatches);

    const sortMode = document.getElementById("import-table-sort")?.value || "intensity_desc";
    const filterQuery = (document.getElementById("import-table-filter")?.value || "").trim().toLowerCase();

    let rows = classified.filter((r) => {
      if (!filterQuery) return true;
      const mzStr = r.mz.toFixed(4);
      const glyStr = r.matchObj ? r.matchObj.glycan_composition.toLowerCase() : "";
      const clsStr = r.isIsotopeSatellite ? "isotopic m+1" : "monoisotopic";
      return mzStr.includes(filterQuery) || glyStr.includes(filterQuery) || clsStr.includes(filterQuery);
    });

    if (sortMode === "intensity_desc") {
      rows.sort((a, b) => b.intensity - a.intensity);
    } else {
      rows.sort((a, b) => a.mz - b.mz);
    }

    if (countBadge) {
      countBadge.textContent = `${rows.length} of ${summary.valid_points_count} peaks`;
    }

    if (rows.length === 0) {
      tbody.innerHTML = `<tr><td colspan="8" class="muted">No peaks match the current filter.</td></tr>`;
      return;
    }

    tbody.innerHTML = rows
      .map((r) => {
        const barWidth = Math.min(100, Math.max(1, r.relPct)).toFixed(1);
        const barClass = r.isBasePeak ? "base-peak" : r.isIsotopeSatellite ? "isotope-peak" : "";

        let clsBadge = `<span class="badge badge-mono">Monoisotopic</span>`;
        if (r.isBasePeak) {
          clsBadge = `<span class="badge badge-pass">Base Peak (100%)</span>`;
        }
        if (r.isIsotopeSatellite) {
          clsBadge = `<span class="badge badge-iso">Isotopic [M+1] (${r.parentMz.toFixed(2)})</span>`;
        }
        if (r.matchObj) {
          clsBadge += ` <span class="badge badge-pass" style="margin-left:4px;">${r.matchObj.glycan_composition} (${r.matchObj.assigned_adduct})</span>`;
        }

        return `
      <tr>
        <td><code>#${r.index}</code></td>
        <td><strong>${r.mz.toFixed(4)}</strong></td>
        <td>${r.intensity.toLocaleString(undefined, { minimumFractionDigits: 1, maximumFractionDigits: 1 })}</td>
        <td class="rel-int-cell">
          <div class="rel-int-wrapper">
            <span class="rel-int-pct">${r.relPct.toFixed(2)}%</span>
            <div class="rel-int-bar-bg">
              <div class="rel-int-bar-fill ${barClass}" style="width: ${barWidth}%;"></div>
            </div>
          </div>
        </td>
        <td>${clsBadge}</td>
        <td>${r.area !== null && r.area !== undefined ? Number(r.area).toFixed(1) : '<span class="muted">2-col TXT</span>'}</td>
        <td>${r.snr !== null && r.snr !== undefined ? Number(r.snr).toFixed(1) : '<span class="muted">—</span>'}</td>
        <td>${r.fwhm !== null && r.fwhm !== undefined ? Number(r.fwhm).toFixed(4) : '<span class="muted">—</span>'}</td>
      </tr>
    `;
      })
      .join("");
  }

  function updateSpectrumUI(summary) {
    appState.spectrumSummary = summary;
    appState.activeSampleId = summary.sample_id;

    document.getElementById("active-sample-badge").textContent = `${summary.sample_id} (${summary.filename})`;
    document.getElementById("active-datatype-badge").textContent = summary.data_type.toUpperCase();
    document.getElementById("kpi-peak-count").textContent = String(summary.valid_points_count);
    document.getElementById("kpi-spectrum-type").textContent = `${summary.source_format} | ${summary.data_type}`;

    const box = document.getElementById("import-summary-box");
    const warningsHtml =
      summary.warnings && summary.warnings.length > 0
        ? `<ul class="safeguard-list" style="margin-top:8px;">${summary.warnings.map((w) => `<li>${w}</li>`).join("")}</ul>`
        : `<p class="muted" style="margin-top:8px;margin-bottom:0;">No parser warnings. Ready for MS1 N-glycan search.</p>`;

    const basePeakStr =
      summary.base_peak_mz !== undefined && summary.base_peak_mz !== null
        ? `${Number(summary.base_peak_mz).toFixed(4)} Da/e`
        : "—";

    box.innerHTML = `
      <div class="summary-metrics-grid">
        <div class="summary-metric-item">
          <span class="sm-label">Sample &amp; Format</span>
          <span class="sm-val">${summary.sample_id} <small class="muted">(${summary.source_format.toUpperCase()})</small></span>
        </div>
        <div class="summary-metric-item">
          <span class="sm-label">Classified Type</span>
          <span class="sm-val"><code>${summary.data_type.toUpperCase()}</code></span>
        </div>
        <div class="summary-metric-item">
          <span class="sm-label">Valid Peaks</span>
          <span class="sm-val">${summary.valid_points_count}</span>
        </div>
        <div class="summary-metric-item">
          <span class="sm-label">m/z Span</span>
          <span class="sm-val">${Number(summary.min_mz).toFixed(1)} &ndash; ${Number(summary.max_mz).toFixed(1)}</span>
        </div>
        <div class="summary-metric-item">
          <span class="sm-label">Base Peak (100%)</span>
          <span class="sm-val" style="color:#059669;">${basePeakStr}</span>
        </div>
        <div class="summary-metric-item">
          <span class="sm-label">Max Intensity</span>
          <span class="sm-val">${Number(summary.max_intensity).toLocaleString()}</span>
        </div>
      </div>
      ${warningsHtml}
    `;

    renderRawTxtInspector(summary);
    renderImportPeakTable();
    renderSpectrumSvg(summary.plot_peaks || summary.preview || []);
    updateCalibrationStatusBanner(summary.is_calibrated ? appState.calibrationResult : null);
  }

  function updateCalibrationStatusBanner(calData) {
    const statusText = document.getElementById("search-cal-status-text");
    if (!statusText) return;
    if (calData && calData.applied) {
      const pre = calData.rms_ppm_before !== null ? Number(calData.rms_ppm_before).toFixed(2) : "—";
      const post = calData.rms_ppm_after !== null ? Number(calData.rms_ppm_after).toFixed(2) : "—";
      statusText.innerHTML = `<span class="badge badge-pass">CALIBRATED</span> Model: <code>${calData.model_type}</code> (${calData.matched_calibrants_count} calibrants, RMS error ${pre} &rarr; <strong>${post} ppm</strong>). Ready for MS1 search.`;
    } else {
      statusText.innerHTML = `<span class="badge badge-warn">UNCALIBRATED</span> Calibrate on imported data first (Import Data or button below) before running MS1 search.`;
    }
  }

  function syncPolarityAdducts(polarity) {
    const polEl = document.getElementById("search-polarity");
    const addEl = document.getElementById("search-adducts");
    if (polEl) polEl.value = polarity;
    if (addEl) {
      const cur = addEl.value.trim();
      if (polarity === "negative" && (!cur || cur.includes("+"))) {
        addEl.value = "[M-H]-, [M+Cl]-";
      } else if (polarity === "positive" && (!cur || cur.includes("-"))) {
        addEl.value = "[M+Na]+, [M+K]+, [M+H]+";
      }
    }
  }

  async function loadSyntheticExample(format = "csv", runPipeline = false) {
    try {
      hideAlert();
      const data = await apiRequest(`/api/spectra/example?format=${encodeURIComponent(format)}`);
      appState.calibrationResult = null;
      appState.ms1SearchResult = null;
      if (format === "complex_neg") {
        syncPolarityAdducts("negative");
      } else if (format === "complex_pos") {
        syncPolarityAdducts("positive");
      }
      updateSpectrumUI(data);
      showAlert(`Loaded synthetic validation spectrum '${data.filename}' (${data.valid_points_count} peaks).`);
      if (runPipeline) {
        await executeMS1Search();
        await runQuantification();
        await evaluateSpectrumQC();
      }
    } catch (err) {
      showAlert(`Failed to load example spectrum: ${err.message}`, true);
    }
  }

  async function uploadSpectrumFile() {
    const fileInput = document.getElementById("input-spectrum-file");
    if (!fileInput || !fileInput.files || fileInput.files.length === 0) {
      showAlert("Please select a spectrum file to upload first.", true);
      return;
    }
    const formData = new FormData();
    formData.append("file", fileInput.files[0]);
    const sid = document.getElementById("import-sample-id").value.trim();
    const dtype = document.getElementById("import-data-type").value;
    const delim = document.getElementById("import-delimiter").value;
    const dec = document.getElementById("import-decimal").value;
    const mzCol = document.getElementById("import-mz-col").value.trim();
    const intCol = document.getElementById("import-int-col").value.trim();
    const dupPol = document.getElementById("import-dup-policy").value;
    const minSnrEl = document.getElementById("import-min-snr");
    const minSnr = minSnrEl ? minSnrEl.value.trim() : "";

    if (sid) formData.append("sample_id", sid);
    if (dtype) formData.append("user_specified_data_type", dtype);
    if (delim) formData.append("delimiter", delim);
    if (dec) formData.append("decimal_separator", dec);
    if (mzCol) formData.append("mz_column", mzCol);
    if (intCol) formData.append("intensity_column", intCol);
    if (dupPol) formData.append("duplicate_policy", dupPol);
    if (minSnr) formData.append("min_snr", minSnr);

    try {
      hideAlert();
      const data = await apiRequest("/api/spectra/upload", {
        method: "POST",
        body: formData,
      });
      appState.calibrationResult = null;
      appState.ms1SearchResult = null;
      updateSpectrumUI(data);
      showAlert(`Successfully parsed '${data.filename}' (${data.valid_points_count} valid peaks).`);
    } catch (err) {
      showAlert(`Spectrum import error: ${err.message}`, true);
    }
  }

  async function applySpectrumCorrections(resetOnly = false) {
    try {
      hideAlert();
      const pol = document.getElementById("search-polarity")?.value || "positive";
      const payload = resetOnly
        ? {}
        : {
            min_snr: document.getElementById("corr-min-snr")?.value || null,
            min_relative_intensity_pct: Number(document.getElementById("corr-min-rel-int")?.value || 0),
            min_mz: document.getElementById("corr-min-mz")?.value || null,
            max_mz: document.getElementById("corr-max-mz")?.value || null,
            normalize: document.getElementById("corr-normalize")?.value || null,
            smoothing: Boolean(document.getElementById("corr-smoothing")?.checked),
            baseline_correction: Boolean(document.getElementById("corr-baseline")?.checked),
            pick_peaks: Boolean(document.getElementById("corr-pick-peaks")?.checked),
            pick_min_snr: Number(document.getElementById("corr-min-snr")?.value || 5.0) || 5.0,
            remove_matrix_peaks: Boolean(document.getElementById("corr-remove-matrix")?.checked),
            matrix_polarity: pol,
            merge_close_peaks: Boolean(document.getElementById("corr-merge-close")?.checked),
          };
      const data = await apiRequest("/api/spectra/process", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload),
      });
      appState.calibrationResult = null;
      appState.ms1SearchResult = null;
      updateSpectrumUI(data);
      const box = document.getElementById("corrections-result-box");
      if (box && data.corrections) {
        const steps = (data.corrections.steps || []).map((s) => `<code>${s.step}</code> (&rarr; ${s.peaks_after})`).join(", ");
        box.innerHTML = `<p><strong>Peaks (Before &rarr; After):</strong> ${data.corrections.peaks_before} &rarr; <strong>${data.corrections.peaks_after}</strong> | <strong>Steps:</strong> ${steps || "none (raw spectrum restored)"}</p>`;
      }
      showAlert(
        resetOnly
          ? `Restored raw imported spectrum (${data.valid_points_count} peaks).`
          : `Applied MALDI corrections: ${data.valid_points_count} peaks active.`
      );
    } catch (err) {
      showAlert(`MALDI spectrum processing error: ${err.message}`, true);
    }
  }

  function renderLibraryTable(entries, totalCount) {
    appState.libraryEntries = entries;
    document.getElementById("lib-shown-count").textContent = `${entries.length} of ${totalCount}`;
    document.getElementById("kpi-lib-count").textContent = String(totalCount);

    const tbody = document.querySelector("#table-library-entries tbody");
    if (!entries || entries.length === 0) {
      tbody.innerHTML = `<tr><td colspan="12" class="muted">No matching library entries found.</td></tr>`;
      return;
    }
    tbody.innerHTML = entries
      .slice(0, 100)
      .map((e) => {
        const adducts = e.adduct_mzs || {};
        const neg = e.negative_mode_mz || {};
        const mhNeg = neg["[M-H]-"] ?? adducts["[M-H]-"];
        const mclNeg = neg["[M+Cl]-"] ?? adducts["[M+Cl]-"];
        const m2hNeg = neg["[M-2H]2-"] ?? adducts["[M-2H]2-"];
        return `
        <tr>
          <td><code>${e.library_id}</code></td>
          <td><strong>${e.composition}</strong></td>
          <td>${e.glycan_class}</td>
          <td>${e.antennary_hypothesis}</td>
          <td>${Number(e.theoretical_neutral_mass).toFixed(4)}</td>
          <td>${adducts["[M+Na]+"] ? Number(adducts["[M+Na]+"]).toFixed(4) : "—"}</td>
          <td>${adducts["[M+H]+"] ? Number(adducts["[M+H]+"]).toFixed(4) : "—"}</td>
          <td>${adducts["[M+K]+"] ? Number(adducts["[M+K]+"]).toFixed(4) : "—"}</td>
          <td>${mhNeg ? Number(mhNeg).toFixed(4) : "—"}</td>
          <td>${mclNeg ? Number(mclNeg).toFixed(4) : "—"}</td>
          <td>${m2hNeg ? Number(m2hNeg).toFixed(4) : "—"}</td>
          <td><span class="badge badge-pass">${e.structural_annotation_status}</span></td>
        </tr>
      `;
      })
      .join("");
  }

  async function loadDefaultLibrary() {
    try {
      const data = await apiRequest("/api/libraries/default");
      document.getElementById("active-library-badge").textContent = `${data.library_name} v${data.library_version} (${data.entry_count} entries)`;
      renderLibraryTable(data.entries, data.entry_count);
    } catch (err) {
      showAlert(`Could not load default library: ${err.message}`, true);
    }
  }

  async function generateCustomLibrary() {
    try {
      hideAlert();
      const payload = {
        min_hex: Number(document.getElementById("lib-min-hex").value),
        max_hex: Number(document.getElementById("lib-max-hex").value),
        min_hexnac: Number(document.getElementById("lib-min-hexnac").value),
        max_hexnac: Number(document.getElementById("lib-max-hexnac").value),
        max_dhex: Number(document.getElementById("lib-max-dhex").value),
        max_neu5ac: Number(document.getElementById("lib-max-neu5ac").value),
        labeling: document.getElementById("lib-labeling").value,
        include_neu5gc: document.getElementById("lib-include-neu5gc").checked,
        expand_sialic_linkages: document.getElementById("lib-expand-linkages").checked,
      };
      const data = await apiRequest("/api/libraries/generate", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload),
      });
      document.getElementById("active-library-badge").textContent = `${data.library_name} (${data.entry_count} entries)`;
      renderLibraryTable(data.entries, data.entry_count);
      showAlert(`Generated & validated library with ${data.entry_count} N-glycan compositions.`);
    } catch (err) {
      showAlert(`Library generation failed: ${err.message}`, true);
    }
  }

  async function filterLibrary() {
    try {
      const payload = {
        composition_query: document.getElementById("lib-filter-comp").value.trim() || null,
        glycan_class: document.getElementById("lib-filter-class").value || null,
        min_neutral_mass: document.getElementById("lib-filter-min-mass").value
          ? Number(document.getElementById("lib-filter-min-mass").value)
          : null,
        max_neutral_mass: document.getElementById("lib-filter-max-mass").value
          ? Number(document.getElementById("lib-filter-max-mass").value)
          : null,
      };
      const data = await apiRequest("/api/libraries/query", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload),
      });
      renderLibraryTable(data.entries, data.total_in_library);
    } catch (err) {
      showAlert(`Library filter error: ${err.message}`, true);
    }
  }

  async function executeMS1Search() {
    try {
      hideAlert();
      const adducts = document
        .getElementById("search-adducts")
        .value.split(",")
        .map((s) => s.trim())
        .filter(Boolean);
      const charges = document
        .getElementById("search-charges")
        .value.split(",")
        .map((s) => parseInt(s.trim(), 10))
        .filter((n) => !Number.isNaN(n));
      const targets = document
        .getElementById("search-targeted-query")
        .value.split(",")
        .map((s) => s.trim())
        .filter(Boolean);

      const payload = {
        search_mode: document.getElementById("search-mode").value,
        tolerance: Number(document.getElementById("search-tolerance").value),
        tolerance_unit: document.getElementById("search-tolerance-unit").value,
        polarity: document.getElementById("search-polarity").value,
        labeling: document.getElementById("search-labeling").value,
        adduct_codes: adducts,
        charge_states: charges,
        min_relative_intensity_pct: Number(document.getElementById("search-min-rel-int").value || 0),
        min_snr: Number(document.getElementById("search-min-snr")?.value || 0),
        min_score: Number(document.getElementById("search-min-score")?.value || 0),
        min_idotp: Number(document.getElementById("search-min-idotp")?.value || 0),
        require_isotope_support: Boolean(document.getElementById("search-require-isotope")?.checked),
        exclude_isotope_satellites: Boolean(document.getElementById("search-exclude-satellites")?.checked),
        targeted_compositions: targets,
      };

      const data = await apiRequest("/api/search/ms1", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload),
      });

      appState.ms1SearchResult = data;
      document.getElementById("search-matches-count").textContent = String(data.matched_peaks_count);
      document.getElementById("kpi-matched-count").textContent = String(data.matched_peaks_count);
      const rate =
        data.eligible_peaks_searched > 0
          ? ((data.matched_peaks_count / data.eligible_peaks_searched) * 100).toFixed(1)
          : "0.0";
      document.getElementById("kpi-match-rate").textContent = `Match rate: ${rate}%`;

      const tbody = document.querySelector("#table-ms1-matches tbody");
      const topMatches = data.top_matches_per_peak || [];
      if (topMatches.length === 0) {
        tbody.innerHTML = `<tr><td colspan="12" class="muted">No candidates found within tolerance.</td></tr>`;
      } else {
        tbody.innerHTML = topMatches
          .map((m) => {
            const ambigBadge = m.is_ambiguous
              ? `<span class="badge badge-warn">AMBIGUOUS (${m.total_candidates_for_peak})</span> <small>${(m.alternative_candidates || []).slice(0, 2).join("; ")}</small>`
              : `<span class="badge badge-pass">Unique</span>`;
            const idotpLabel =
              m.isotope_score_idotp !== null && m.isotope_score_idotp !== undefined
                ? Number(m.isotope_score_idotp).toFixed(3)
                : "Plot M0..M+4";
            const idotpBtn = `<button type="button" class="btn btn-outline btn-xs btn-idotp-plot" data-lib-id="${m.library_id || ""}" data-comp="${m.glycan_composition}" data-adduct="${m.assigned_adduct}" title="Click to compare Theoretical vs. Experimental Isotope Envelope">${idotpLabel}</button>`;
            return `
            <tr>
              <td>${m.peak_index}</td>
              <td>${Number(m.observed_mz).toFixed(4)}</td>
              <td>${Number(m.observed_intensity).toFixed(1)}</td>
              <td>${Number(m.relative_intensity_pct).toFixed(2)}%</td>
              <td><strong>${m.glycan_composition}</strong></td>
              <td>${m.candidate_class}</td>
              <td><code>${m.assigned_adduct}</code></td>
              <td>${Number(m.theoretical_mz).toFixed(4)}</td>
              <td>${Number(m.signed_error_ppm).toFixed(2)}</td>
              <td><strong>${Number(m.composite_heuristic_score).toFixed(1)}</strong></td>
              <td>${idotpBtn}</td>
              <td>${ambigBadge}</td>
            </tr>
          `;
          })
          .join("");

        tbody.querySelectorAll(".btn-idotp-plot").forEach((btn) => {
          btn.addEventListener("click", () => {
            openIsotopeComparison({
              library_id: btn.getAttribute("data-lib-id"),
              composition: btn.getAttribute("data-comp"),
              adduct_code: btn.getAttribute("data-adduct"),
            });
          });
        });
      }

      if (appState.spectrumSummary) {
        renderSpectrumSvg(appState.spectrumSummary.plot_peaks || appState.spectrumSummary.preview || [], topMatches);
        renderImportPeakTable();
      }
      const elapsedNote = data.elapsed_seconds !== undefined ? ` in ${(Number(data.elapsed_seconds) * 1000).toFixed(0)} ms` : "";
      showAlert(`MS1 search completed${elapsedNote}: ${data.matched_peaks_count} of ${data.eligible_peaks_searched} peaks annotated.`);
    } catch (err) {
      showAlert(`MS1 search failed: ${err.message}`, true);
    }
  }

  async function runQuantification() {
    try {
      const payload = {
        denominator_policy: document.getElementById("quant-denom-policy").value,
        aggregation_policy: document.getElementById("quant-agg-policy").value,
        primary_adduct: document.getElementById("quant-primary-adduct").value.trim(),
        exclude_ambiguous_from_aggregation: document.getElementById("quant-exclude-ambiguous").checked,
      };
      const data = await apiRequest("/api/search/quantify", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload),
      });
      appState.quantificationResult = data.quantification;
      appState.derivedTraits = data.derived_traits || [];

      const qBody = document.querySelector("#table-quantified-glycans tbody");
      const abunds = data.quantification.glycan_abundances || [];
      qBody.innerHTML =
        abunds.length === 0
          ? `<tr><td colspan="7" class="muted">No glycans quantified.</td></tr>`
          : abunds
              .map(
                (g) => `
          <tr>
            <td><strong>${g.glycan_composition}</strong></td>
            <td>${g.candidate_class}</td>
            <td>${Number(g.aggregated_intensity).toFixed(1)}</td>
            <td>${g.aggregated_area !== null && g.aggregated_area !== undefined ? Number(g.aggregated_area).toFixed(1) : "—"}</td>
            <td>${Number(g.base_peak_relative_intensity_pct).toFixed(2)}%</td>
            <td><strong>${Number(g.relative_abundance_intensity_pct).toFixed(2)}%</strong></td>
            <td><small>${(g.contributing_adducts || []).join(", ")}</small></td>
          </tr>
        `
              )
              .join("");

      const tBody = document.querySelector("#table-derived-traits tbody");
      tBody.innerHTML = appState.derivedTraits
        .map(
          (t) => `
        <tr>
          <td><code>${t.trait_code}</code></td>
          <td>${t.trait_name}</td>
          <td><strong>${t.value !== null && t.value !== undefined ? Number(t.value).toFixed(2) + "%" : "—"}</strong></td>
          <td>${t.is_composition_inferred_hypothesis ? '<span class="badge badge-warn">Composition Hypothesis</span>' : '<span class="badge badge-pass">Direct Stoichiometry</span>'}</td>
          <td><small>${t.structural_disclaimer}</small></td>
        </tr>
      `
        )
        .join("");
    } catch (err) {
      showAlert(`Quantification error: ${err.message}`, true);
    }
  }

  async function runGlycopeptideScreening() {
    try {
      hideAlert();
      const precursors = document
        .getElementById("gp-precursor-mzs")
        .value.split(",")
        .map((s) => parseFloat(s.trim()))
        .filter((n) => !Number.isNaN(n));
      const oxoniums = document
        .getElementById("gp-oxonium-mzs")
        .value.split(",")
        .map((s) => parseFloat(s.trim()))
        .filter((n) => !Number.isNaN(n));

      const payload = {
        fasta_text: document.getElementById("gp-fasta-input").value,
        enzyme: document.getElementById("gp-protease").value,
        max_missed_cleavages: Number(document.getElementById("gp-missed-cleavages").value),
        precursor_mzs: precursors,
        tolerance_ppm: Number(document.getElementById("gp-tolerance-ppm").value),
        ms2_fragment_mzs: oxoniums,
        fixed_carbamidomethyl_cys: document.getElementById("gp-cys-cam").checked,
      };

      const data = await apiRequest("/api/glycopeptides/screen", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload),
      });

      appState.glycopeptideResult = data;
      document.getElementById("gp-matches-count").textContent = String((data.matches || []).length);

      const seqBox = document.getElementById("gp-sequon-summary");
      seqBox.innerHTML = `<p><strong>Digested Sequon-Containing Peptides:</strong> ${data.digested_sequon_peptides_count} | <strong>Matches Found:</strong> ${(data.matches || []).length}</p>`;

      const tbody = document.querySelector("#table-glycopeptide-matches tbody");
      const matches = data.matches || [];
      tbody.innerHTML =
        matches.length === 0
          ? `<tr><td colspan="11" class="muted">No glycopeptide candidates matched within tolerance.</td></tr>`
          : matches
              .map(
                (m) => `
          <tr>
            <td><code>${m.protein_id}</code></td>
            <td><strong>${m.peptide_sequence}</strong></td>
            <td>${(m.candidate_glycosylation_sites || []).join(", ")}</td>
            <td><strong>${m.glycan_composition}</strong></td>
            <td>${Number(m.peptide_neutral_mass).toFixed(4)}</td>
            <td>${Number(m.attached_glycan_moiety_mass).toFixed(4)}</td>
            <td>${Number(m.observed_precursor_mz).toFixed(4)}</td>
            <td>${Number(m.theoretical_precursor_mz).toFixed(4)}</td>
            <td><code>${m.adduct} (z=${m.charge_state})</code></td>
            <td>${Number(m.signed_error_ppm).toFixed(2)}</td>
            <td><span class="badge badge-pass">${m.confidence_category}</span></td>
          </tr>
        `
              )
              .join("");
      showAlert(`Glycopeptide screening finished: ${matches.length} candidate(s) identified.`);
    } catch (err) {
      showAlert(`Glycopeptide screening failed: ${err.message}`, true);
    }
  }

  async function runDemoBatch() {
    try {
      hideAlert();
      const data = await apiRequest("/api/batch/demo", { method: "POST" });
      renderBatchResult(data);
      showAlert(`Batch job '${data.batch_id}' completed (${data.completed_files_count}/${data.total_files} samples).`);
    } catch (err) {
      showAlert(`Batch execution failed: ${err.message}`, true);
    }
  }

  async function runUploadBatch() {
    const input = document.getElementById("input-batch-files");
    if (!input || !input.files || input.files.length === 0) {
      showAlert("Please select one or more spectrum files or a .zip archive for batch processing.", true);
      return;
    }
    const formData = new FormData();
    Array.from(input.files).forEach((f) => formData.append("files", f));
    formData.append("default_group", document.getElementById("batch-default-group").value.trim() || "Cohort_A");
    formData.append("tolerance_ppm", document.getElementById("batch-tolerance-ppm").value);

    try {
      hideAlert();
      const data = await apiRequest("/api/batch/run", {
        method: "POST",
        body: formData,
      });
      renderBatchResult(data);
      showAlert(`Batch job '${data.batch_id}' finished with status '${data.status}'.`);
    } catch (err) {
      showAlert(`Batch upload error: ${err.message}`, true);
    }
  }

  function renderBatchResult(data) {
    appState.batchResult = data;
    const box = document.getElementById("batch-summary-box");
    box.innerHTML = `
      <p><strong>Batch ID:</strong> <code>${data.batch_id}</code> | <strong>Status:</strong> <strong>${data.status}</strong></p>
      <p><strong>Completed:</strong> ${data.completed_files_count} | <strong>Failed:</strong> ${data.failed_files_count} | <strong>Total Files:</strong> ${data.total_files}</p>
      <p><strong>Reproducibility Library:</strong> ${data.reproducibility.library_name} (v${data.reproducibility.library_version}) | <strong>Tolerance:</strong> ${data.reproducibility.mass_tolerance} ${data.reproducibility.mass_tolerance_unit}</p>
    `;

    const tbody = document.querySelector("#table-batch-outcomes tbody");
    const rows = data.qc_summary_table || [];
    tbody.innerHTML = rows
      .map(
        (r) => `
      <tr>
        <td><strong>${r.sample_id}</strong></td>
        <td>${r.filename}</td>
        <td>${r.group || "—"}</td>
        <td>${r.data_type}</td>
        <td>${r.imported_peaks}</td>
        <td>${r.matched_peaks} (${r.match_rate_pct}%)</td>
        <td>${r.quantified_glycans}</td>
        <td>${r.median_abs_ppm_error !== null ? r.median_abs_ppm_error : "—"}</td>
        <td><span class="badge badge-pass">${r.status}</span></td>
      </tr>
    `
      )
      .join("");
  }

  async function openIsotopeComparison({ library_id, composition, adduct_code }) {
    try {
      hideAlert();
      const data = await apiRequest("/api/search/isotope", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          library_id: library_id || null,
          composition: composition || null,
          adduct_code: adduct_code || "[M+Na]+",
          num_isotopes: 5,
        }),
      });
      const panel = document.getElementById("isotope-comparison-panel");
      if (!panel) return;
      panel.classList.remove("hidden");
      const titleEl = document.getElementById("isotope-panel-title");
      if (titleEl) {
        titleEl.textContent = `Isotope Envelope Comparison — ${data.composition} ${data.adduct_code} (${data.ion_formula || ""})`;
      }
      const metaEl = document.getElementById("isotope-panel-meta");
      if (metaEl) {
        const idotpStr =
          data.idotp !== null && data.idotp !== undefined
            ? `<span class="badge badge-pass">idotp = ${Number(data.idotp).toFixed(4)}</span>`
            : `<span class="badge badge-warn">idotp = N/A (single peak / insufficient satellites)</span>`;
        metaEl.innerHTML = `
          <p><strong>Candidate:</strong> <code>${data.composition}</code> (${data.adduct_code}, z = ${data.charge}) | <strong>Monoisotopic m/z:</strong> ${Number(data.theoretical_mz).toFixed(4)} | <strong>Matched Isotopes:</strong> ${data.matched_isotope_peaks} / ${data.expected_isotope_peaks} | ${idotpStr}</p>
          <p class="muted" style="margin:0;">${data.status_message} | Tolerance: &plusmn;${data.tolerance} ${data.tolerance_unit}</p>
        `;
      }

      const tbody = document.querySelector("#table-isotope-envelope tbody");
      if (tbody) {
        tbody.innerHTML = (data.isotopologues || [])
          .map((iso) => {
            const obsMz = iso.observed_mz !== null && iso.observed_mz !== undefined ? Number(iso.observed_mz).toFixed(4) : '<span class="muted">Not detected</span>';
            const errPpm = iso.error_ppm !== null && iso.error_ppm !== undefined ? `${Number(iso.error_ppm).toFixed(2)} ppm` : "—";
            const dRatio =
              iso.delta_base_ratio_pct !== null && iso.delta_base_ratio_pct !== undefined
                ? `${Number(iso.delta_base_ratio_pct) >= 0 ? "+" : ""}${Number(iso.delta_base_ratio_pct).toFixed(1)}%`
                : "—";
            return `
              <tr>
                <td><strong>${iso.label}</strong></td>
                <td>${Number(iso.theoretical_mz).toFixed(4)}</td>
                <td>${obsMz}</td>
                <td>${errPpm}</td>
                <td>${Number(iso.theoretical_base_ratio_pct).toFixed(1)}%</td>
                <td>${Number(iso.observed_base_ratio_pct).toFixed(1)}%</td>
                <td>${dRatio}</td>
              </tr>
            `;
          })
          .join("");
      }

      renderIsotopeEnvelopeSvg(data);
      panel.scrollIntoView({ behavior: "smooth", block: "nearest" });
    } catch (err) {
      showAlert(`Could not load isotope envelope comparison: ${err.message}`, true);
    }
  }

  function renderIsotopeEnvelopeSvg(data) {
    const svg = document.getElementById("isotope-comparison-svg");
    if (!svg) return;
    const isos = data.isotopologues || [];
    if (isos.length === 0) {
      svg.innerHTML = "";
      return;
    }
    const W = 640;
    const H = 310;
    const padL = 56;
    const padR = 24;
    const padT = 26;
    const padB = 38;
    const plotW = W - padL - padR;
    const plotH = H - padT - padB;
    const midY = padT + plotH * 0.52;
    const topH = plotH * 0.46;
    const botH = plotH * 0.42;

    const minMz = Number(isos[0].theoretical_mz) - 0.65;
    const maxMz = Number(isos[isos.length - 1].theoretical_mz) + 0.65;
    const xOf = (mz) => padL + ((Number(mz) - minMz) / Math.max(0.1, maxMz - minMz)) * plotW;

    const parts = [];
    parts.push(`<rect x="${padL}" y="${padT}" width="${plotW}" height="${plotH}" fill="#f8fafc" stroke="#cbd5e1" rx="4"/>`);
    parts.push(`<line x1="${padL}" y1="${midY}" x2="${W - padR}" y2="${midY}" stroke="#475569" stroke-width="1.2"/>`);
    parts.push(`<text x="${padL + 6}" y="${padT + 14}" font-size="10" fill="#1e40af" font-weight="600">Upper: Side-by-Side Base Ratio (%) — Theoretical (Blue) vs. Experimental (Green)</text>`);
    parts.push(`<text x="${padL + 6}" y="${midY + botH + 14}" font-size="10" fill="#047857" font-weight="600">Lower Mirror: Experimental Observed Envelope &amp; Local Spectrum Trace</text>`);

    // Local raw spectrum trace in the mirror region if present
    const trace = data.observed_trace || [];
    if (trace.length > 1) {
      const maxTraceInt = Math.max(...trace.map((p) => Number(p.intensity)), 1);
      const pts = trace
        .filter((p) => Number(p.mz) >= minMz && Number(p.mz) <= maxMz)
        .map((p) => {
          const x = xOf(p.mz).toFixed(1);
          const y = (midY + (Number(p.intensity) / maxTraceInt) * botH).toFixed(1);
          return `${x},${y}`;
        });
      if (pts.length > 1) {
        parts.push(`<polyline fill="none" stroke="#94a3b8" stroke-width="1.1" opacity="0.65" points="${pts.join(" ")}"/>`);
      }
    }

    const barW = Math.min(18, Math.max(8, plotW / (isos.length * 4.5)));
    isos.forEach((iso) => {
      const cx = xOf(iso.theoretical_mz);
      const theoRatio = Math.max(0, Math.min(100, Number(iso.theoretical_base_ratio_pct || 0)));
      const obsRatio = Math.max(0, Math.min(100, Number(iso.observed_base_ratio_pct || 0)));
      const hTheo = (theoRatio / 100.0) * topH;
      const hObs = (obsRatio / 100.0) * topH;
      const hMirror = (obsRatio / 100.0) * botH;

      // Upper side-by-side bars (Theoretical left, Observed right)
      parts.push(`<rect x="${(cx - barW - 1).toFixed(1)}" y="${(midY - hTheo).toFixed(1)}" width="${barW}" height="${Math.max(1, hTheo).toFixed(1)}" fill="#2563eb" opacity="0.88" rx="2"><title>${iso.label} Theoretical: ${theoRatio.toFixed(1)}% (${Number(iso.theoretical_mz).toFixed(4)})</title></rect>`);
      if (obsRatio > 0) {
        parts.push(`<rect x="${(cx + 1).toFixed(1)}" y="${(midY - hObs).toFixed(1)}" width="${barW}" height="${Math.max(1, hObs).toFixed(1)}" fill="#059669" opacity="0.9" rx="2"><title>${iso.label} Experimental: ${obsRatio.toFixed(1)}% (${iso.observed_mz})</title></rect>`);
        // Lower mirror stick for observed isotopologue
        const ox = iso.observed_mz ? xOf(iso.observed_mz) : cx;
        parts.push(`<line x1="${ox.toFixed(1)}" y1="${midY}" x2="${ox.toFixed(1)}" y2="${(midY + hMirror).toFixed(1)}" stroke="#059669" stroke-width="3"/>`);
        parts.push(`<circle cx="${ox.toFixed(1)}" cy="${(midY + hMirror).toFixed(1)}" r="3.5" fill="#059669"/>`);
      } else {
        parts.push(`<line x1="${(cx + 1).toFixed(1)}" y1="${midY - 4}" x2="${(cx + 1 + barW).toFixed(1)}" y2="${midY - 4}" stroke="#ef4444" stroke-width="1.5" stroke-dasharray="2,2"/>`);
      }

      // Label above bar pair
      const topY = midY - Math.max(hTheo, hObs) - 5;
      parts.push(`<text x="${cx.toFixed(1)}" y="${Math.max(padT + 24, topY).toFixed(1)}" text-anchor="middle" font-size="10" font-weight="600" fill="#1e293b">${iso.label}</text>`);
      // m/z tick at bottom
      parts.push(`<text x="${cx.toFixed(1)}" y="${H - 12}" text-anchor="middle" font-size="10" fill="#475569">${Number(iso.theoretical_mz).toFixed(2)}</text>`);
    });

    svg.innerHTML = parts.join("");
  }

  function renderCalibrationSummaryBoxes(data) {
    const html = `
      <p><strong>Applied:</strong> ${data.applied ? '<span class="badge badge-pass">YES</span>' : '<span class="badge badge-warn">NO (Uncalibrated)</span>'} | <strong>Model:</strong> ${data.model_type} | <strong>Panel:</strong> <code>${data.reference_set || "nglycan_na"}</code></p>
      <p><strong>Matched Calibrants:</strong> ${data.matched_calibrants_count} | <strong>RMS ppm (Before &rarr; After):</strong> ${data.rms_ppm_before !== null ? Number(data.rms_ppm_before).toFixed(3) : "—"} &rarr; <strong>${data.rms_ppm_after !== null ? Number(data.rms_ppm_after).toFixed(3) : "—"} ppm</strong></p>
      ${data.warnings && data.warnings.length ? `<p class="muted">${data.warnings.join("; ")}</p>` : ""}
    `;
    ["calibration-result-box", "import-calibration-box"].forEach((id) => {
      const box = document.getElementById(id);
      if (box) box.innerHTML = html;
    });
    updateCalibrationStatusBanner(data);
  }

  async function runCalibration(source = "qc") {
    try {
      hideAlert();
      const useImport = source === "import";
      const model_type = document.getElementById(useImport ? "import-cal-model" : "cal-model-type")?.value || "linear";
      const reference_set = document.getElementById(useImport ? "import-cal-panel" : "cal-reference-set")?.value || "auto";
      const match_tolerance_ppm = Number(
        document.getElementById(useImport ? "import-cal-tol" : "cal-tolerance-ppm")?.value || 150.0
      );
      const payload = { model_type, reference_set, match_tolerance_ppm };
      const data = await apiRequest("/api/qc/calibrate", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload),
      });
      appState.calibrationResult = data;
      renderCalibrationSummaryBoxes(data);
      if (data.reference_set === "nglycan_neg") {
        syncPolarityAdducts("negative");
      } else if (data.reference_set && data.reference_set !== "auto") {
        syncPolarityAdducts("positive");
      }
      showAlert(
        data.applied
          ? `Mass calibration applied (${data.matched_calibrants_count} calibrants, RMS ${Number(data.rms_ppm_before).toFixed(1)} → ${Number(data.rms_ppm_after).toFixed(2)} ppm).`
          : "Calibration safeguard triggered; spectrum left uncalibrated."
      );
      return data;
    } catch (err) {
      showAlert(`Mass calibration failed: ${err.message}`, true);
      return null;
    }
  }

  async function resetCalibration() {
    try {
      await apiRequest("/api/qc/calibrate/reset", { method: "POST" });
      appState.calibrationResult = null;
      ["calibration-result-box", "import-calibration-box"].forEach((id) => {
        const box = document.getElementById(id);
        if (box) box.innerHTML = `<p class="muted">Calibration removed; spectrum is uncalibrated.</p>`;
      });
      updateCalibrationStatusBanner(null);
      showAlert("Calibration reset; active spectrum restored to uncalibrated m/z axis.");
    } catch (err) {
      showAlert(`Failed to reset calibration: ${err.message}`, true);
    }
  }

  async function evaluateSpectrumQC() {
    try {
      const data = await apiRequest("/api/qc/evaluate", { method: "POST" });
      appState.qcReport = data;
      document.getElementById("kpi-qc-status").textContent = data.overall_status;
      document.getElementById("kpi-median-ppm").textContent =
        data.median_abs_ppm_error !== null ? `Median |ppm|: ${Number(data.median_abs_ppm_error).toFixed(2)}` : "Median |ppm|: —";

      const tbody = document.querySelector("#table-qc-checks tbody");
      tbody.innerHTML = (data.criteria_checks || [])
        .map((c) => {
          const cls = c.status === "PASS" ? "badge-pass" : c.status === "WARN" ? "badge-warn" : "badge-fail";
          return `
          <tr>
            <td><code>${c.metric_name}</code></td>
            <td><strong>${c.observed_value !== null ? c.observed_value : "—"}</strong></td>
            <td>${c.pass_threshold}</td>
            <td><span class="badge ${cls}">${c.status}</span></td>
            <td><code>${c.threshold_source}</code></td>
            <td><small>${c.rationale}</small></td>
          </tr>
        `;
        })
        .join("");
    } catch (err) {
      showAlert(`QC evaluation error: ${err.message}`, true);
    }
  }

  async function runCohortStatisticsDemo() {
    try {
      hideAlert();
      const payload = {
        imputation_method: document.getElementById("stat-imputation").value,
        classifier_type: document.getElementById("stat-classifier").value,
      };
      const data = await apiRequest("/api/qc/statistics-demo", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload),
      });

      const box = document.getElementById("stat-ml-summary-box");
      box.innerHTML = `
        <p><strong>PCA Variance Explained:</strong> PC1 = ${(data.pca.explained_variance_ratio[0] * 100).toFixed(1)}%, PC2 = ${((data.pca.explained_variance_ratio[1] || 0) * 100).toFixed(1)}% (Imputation: <code>${data.pca.imputation_method_used}</code>)</p>
        <p><strong>Leakage-Free Cross-Validated Classification (${data.classification.classifier_name}):</strong> Out-of-Fold Accuracy = ${(data.classification.oof_accuracy * 100).toFixed(1)}%, ROC AUC = ${data.classification.oof_roc_auc}</p>
        <p><small>${data.classification.leakage_prevention_statement}</small></p>
      `;

      if (data.pca_plot_base64) {
        const plotBox = document.getElementById("stat-plot-container");
        plotBox.innerHTML = `<img src="data:image/png;base64,${data.pca_plot_base64}" alt="PCA Score Plot" style="max-width:100%;height:auto;border-radius:6px;" />`;
      }

      const tbody = document.querySelector("#table-group-comparisons tbody");
      tbody.innerHTML = (data.group_comparisons || [])
        .map(
          (g) => `
        <tr>
          <td><strong>${g.glycan_composition}</strong></td>
          <td>${(g.groups_compared || []).join(" vs ")}</td>
          <td>${JSON.stringify(g.group_means)}</td>
          <td>${g.log2_fold_change !== null ? Number(g.log2_fold_change).toFixed(3) : "—"}</td>
          <td>${g.parametric_test_name}</td>
          <td>${g.parametric_p_value !== null ? Number(g.parametric_p_value).toExponential(2) : "—"}</td>
          <td><strong>${g.parametric_fdr_q_value !== null ? Number(g.parametric_fdr_q_value).toExponential(2) : "—"}</strong></td>
          <td>${g.nonparametric_p_value !== null ? Number(g.nonparametric_p_value).toExponential(2) : "—"}</td>
          <td>${g.significant_fdr_0_05 ? '<span class="badge badge-pass">YES</span>' : '<span class="badge badge-warn">NO</span>'}</td>
        </tr>
      `
        )
        .join("");

      showAlert("Completed multivariate PCA, differential BH-FDR tests, and leakage-free CV classification.");
    } catch (err) {
      showAlert(`Statistical suite error: ${err.message}`, true);
    }
  }

  function triggerDownload(format, category = "") {
    const url = `/api/exports/download?format=${encodeURIComponent(format)}${
      category ? `&category=${encodeURIComponent(category)}` : ""
    }`;
    window.location.href = url;
  }

  // Wire event listeners on DOMContentLoaded
  document.addEventListener("DOMContentLoaded", () => {
    // Sidebar navigation
    document.querySelectorAll(".nav-item").forEach((btn) => {
      btn.addEventListener("click", () => switchView(btn.getAttribute("data-view")));
    });
    document.querySelectorAll("[data-goto]").forEach((btn) => {
      btn.addEventListener("click", () => switchView(btn.getAttribute("data-goto")));
    });

    const closeAlertBtn = document.getElementById("btn-close-alert");
    if (closeAlertBtn) closeAlertBtn.addEventListener("click", hideAlert);

    // Dropzone & file input
    const dropzone = document.getElementById("dropzone-spectrum");
    const fileInput = document.getElementById("input-spectrum-file");
    if (dropzone && fileInput) {
      dropzone.addEventListener("click", () => fileInput.click());
      dropzone.addEventListener("dragover", (e) => {
        e.preventDefault();
        dropzone.classList.add("dragover");
      });
      dropzone.addEventListener("dragleave", () => dropzone.classList.remove("dragover"));
      dropzone.addEventListener("drop", (e) => {
        e.preventDefault();
        dropzone.classList.remove("dragover");
        if (e.dataTransfer && e.dataTransfer.files.length > 0) {
          fileInput.files = e.dataTransfer.files;
          uploadSpectrumFile();
        }
      });
    }

    // Buttons across views
    document.getElementById("btn-load-example-Top")?.addEventListener("click", () => loadSyntheticExample("csv", true));
    document.getElementById("btn-dashboard-load-example")?.addEventListener("click", () => loadSyntheticExample("csv", true));
    document.getElementById("btn-quick-run-search")?.addEventListener("click", () => {
      switchView("search");
      executeMS1Search();
    });

    document.getElementById("btn-upload-spectrum")?.addEventListener("click", uploadSpectrumFile);
    document.getElementById("btn-load-synthetic-csv")?.addEventListener("click", () => loadSyntheticExample("csv", false));
    document.getElementById("btn-load-synthetic-txt")?.addEventListener("click", () => loadSyntheticExample("txt", false));
    document.getElementById("btn-load-synthetic-xy")?.addEventListener("click", () => loadSyntheticExample("xy", false));
    document.getElementById("btn-load-complex-pos")?.addEventListener("click", () => loadSyntheticExample("complex_pos", false));
    document.getElementById("btn-load-complex-neg")?.addEventListener("click", () => loadSyntheticExample("complex_neg", false));

    document.getElementById("btn-apply-corrections")?.addEventListener("click", () => applySpectrumCorrections(false));
    document.getElementById("btn-reset-corrections")?.addEventListener("click", () => applySpectrumCorrections(true));
    document.getElementById("btn-import-calibrate")?.addEventListener("click", () => runCalibration("import"));
    document.getElementById("btn-import-reset-calibration")?.addEventListener("click", resetCalibration);
    document.getElementById("btn-reset-calibration")?.addEventListener("click", resetCalibration);
    document.getElementById("btn-calibrate-and-search")?.addEventListener("click", async () => {
      await runCalibration("import");
      switchView("search");
      await executeMS1Search();
    });
    document.getElementById("btn-search-calibrate-first")?.addEventListener("click", async () => {
      await runCalibration("import");
      await executeMS1Search();
    });
    document.getElementById("search-polarity")?.addEventListener("change", (e) => {
      syncPolarityAdducts(e.target.value);
    });
    document.getElementById("btn-close-isotope-panel")?.addEventListener("click", () => {
      document.getElementById("isotope-comparison-panel")?.classList.add("hidden");
    });

    // Interactive Spectrum Viewer Toolbar controls
    const redrawCurrentSpectrum = () => {
      if (!appState.spectrumSummary) return;
      const topMatches = appState.ms1SearchResult ? appState.ms1SearchResult.top_matches_per_peak || [] : [];
      renderSpectrumSvg(appState.spectrumSummary.plot_peaks || appState.spectrumSummary.preview || [], topMatches);
    };
    document.getElementById("spec-plot-style")?.addEventListener("change", redrawCurrentSpectrum);
    document.getElementById("spec-label-mode")?.addEventListener("change", redrawCurrentSpectrum);
    document.getElementById("btn-spec-zoom-apply")?.addEventListener("click", redrawCurrentSpectrum);
    document.getElementById("btn-spec-zoom-reset")?.addEventListener("click", () => {
      const minEl = document.getElementById("spec-zoom-min");
      const maxEl = document.getElementById("spec-zoom-max");
      if (minEl) minEl.value = "";
      if (maxEl) maxEl.value = "";
      redrawCurrentSpectrum();
    });

    // Peak list table sort & filter controls
    document.getElementById("import-table-sort")?.addEventListener("change", renderImportPeakTable);
    document.getElementById("import-table-filter")?.addEventListener("input", renderImportPeakTable);

    document.getElementById("btn-generate-library")?.addEventListener("click", generateCustomLibrary);
    document.getElementById("btn-reset-default-library")?.addEventListener("click", loadDefaultLibrary);
    document.getElementById("btn-query-library")?.addEventListener("click", filterLibrary);

    document.getElementById("btn-execute-ms1-search")?.addEventListener("click", executeMS1Search);
    document.getElementById("btn-run-glycopeptide-search")?.addEventListener("click", runGlycopeptideScreening);
    document.getElementById("btn-run-batch-upload")?.addEventListener("click", runUploadBatch);
    document.getElementById("btn-run-demo-batch")?.addEventListener("click", runDemoBatch);
    document.getElementById("btn-run-quantification")?.addEventListener("click", runQuantification);
    document.getElementById("btn-run-calibration")?.addEventListener("click", () => runCalibration("qc"));
    document.getElementById("btn-evaluate-qc")?.addEventListener("click", evaluateSpectrumQC);
    document.getElementById("btn-run-demo-statistics")?.addEventListener("click", runCohortStatisticsDemo);

    // Export buttons
    document.getElementById("btn-export-zip")?.addEventListener("click", () => triggerDownload("zip"));
    document.getElementById("btn-export-xlsx")?.addEventListener("click", () => triggerDownload("xlsx"));
    document.getElementById("btn-export-pdf")?.addEventListener("click", () => triggerDownload("pdf"));
    document.getElementById("btn-export-json")?.addEventListener("click", () => triggerDownload("json"));
    document.querySelectorAll(".btn-dl-cat").forEach((btn) => {
      btn.addEventListener("click", () => {
        triggerDownload(btn.getAttribute("data-fmt"), btn.getAttribute("data-cat"));
      });
    });

    // Settings buttons
    document.getElementById("btn-save-settings")?.addEventListener("click", () => {
      const ppm = document.getElementById("set-default-ppm").value;
      document.getElementById("search-tolerance").value = ppm;
      showAlert(`Applied default MS1 tolerance (${ppm} ppm) to active session.`);
    });
    document.getElementById("btn-clear-session")?.addEventListener("click", () => {
      window.location.reload();
    });

    // Load default library on startup
    loadDefaultLibrary();
  });
})();

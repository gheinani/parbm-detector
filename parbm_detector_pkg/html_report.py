"""
Interactive HTML report generator using Plotly.

Produces a single self-contained HTML file with:
  • Protein identity header cards
  • Tabbed layout: Overview | Protein Canvas | PARylation Sites |
                   Disorder | Pathways & Functions | Literature
  • Zoomable, hoverable Plotly figures
  • Filterable site table
  • Summary statistics cards
"""

from __future__ import annotations
import re
import json
from datetime import datetime
from typing import Dict, List, Optional

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from plotly.offline import plot as _plotly_plot
    import plotly.io as pio
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False


# ── colour palette (matches matplotlib figures) ───────────────────────────────
C_PARBM   = "#E63946"
C_DOMAIN  = "#457B9D"
C_HC      = "#F4A261"
C_BG      = "#F8F9FA"
C_DARK    = "#1D3557"
C_TEAL    = "#2A9D8F"
C_PURPLE  = "#9B5DE5"

MOTIF_COLORS = {
    "SxxE motif":  C_PARBM,
    "SxxD motif":  C_HC,
    "TxxE motif":  C_TEAL,
    "ExxE motif":  C_DOMAIN,
    "Acidic region": C_PURPLE,
}


def _esc(s: str) -> str:
    """HTML-escape a string."""
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


class HTMLReportGenerator:
    """Generate a full interactive HTML report for one protein analysis result."""

    def generate(self, result: Dict, output_file: str) -> Optional[str]:
        if not PLOTLY_AVAILABLE:
            print("  ✗ plotly not installed. Run: pip install plotly")
            return None
        if "error" in result:
            print(f"  ✗ {result['error']}")
            return None

        html = self._build_html(result)
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(html)
        print(f"  ✓ HTML report → {output_file}")
        return output_file

    # ── top-level builder ─────────────────────────────────────────────────────

    def _build_html(self, result: Dict) -> str:
        uid   = result.get("uniprot_id") or "N/A"
        gene  = result.get("gene_name") or "N/A"
        pname = result.get("protein_name", "Unknown")
        org   = result.get("organism", "Unknown")
        seqlen = result.get("sequence_length", 0)
        s      = result.get("summary", {})
        enrich = result.get("enrichment", {})

        canvas_div    = self._protein_canvas(result)
        conf_div      = self._confidence_figure(result)
        residue_div   = self._residue_breakdown(result)
        disorder_div  = self._disorder_figure(result)
        pathways_html = self._pathways_section(enrich)
        literature_html = self._literature_section(enrich)
        site_table    = self._site_table(result)
        domain_table  = self._domain_table(result)
        stats_cards   = self._stats_cards(result, enrich)
        quality_html  = self._quality_metrics_section(result)
        motif_div     = self._motif_context_section(enrich)

        hc    = s.get("high_confidence_sites_count", 0)
        total = s.get("total_paryation_sites", 0)
        top   = s.get("top_site", "N/A")
        parbm = s.get("parbm_domains_label", "—")

        now = datetime.now().strftime("%Y-%m-%d %H:%M")

        return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>PARBM Detector – {_esc(pname)}</title>
<link rel="stylesheet"
      href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
<style>
  body {{ background:{C_BG}; font-family:'Segoe UI',sans-serif; color:{C_DARK}; }}
  .header-bar {{ background:{C_DARK}; color:#fff; padding:22px 32px 14px; }}
  .header-bar h1 {{ font-size:1.6rem; font-weight:700; margin:0 0 4px; }}
  .header-bar .sub {{ color:#A8DADC; font-size:0.9rem; }}
  .id-badge {{ background:rgba(255,255,255,.12); border-radius:6px;
               padding:3px 10px; font-size:0.82rem; color:#A8DADC; display:inline-block; margin:2px; }}
  .stat-card {{ background:#fff; border-left:4px solid {C_PARBM};
                border-radius:8px; padding:14px 18px; box-shadow:0 1px 4px rgba(0,0,0,.08); }}
  .stat-card.blue {{ border-left-color:{C_DOMAIN}; }}
  .stat-card.teal {{ border-left-color:{C_TEAL}; }}
  .stat-card.purple {{ border-left-color:{C_PURPLE}; }}
  .stat-card .val {{ font-size:1.8rem; font-weight:700; color:{C_DARK}; }}
  .stat-card .lbl {{ font-size:0.8rem; color:#6c757d; text-transform:uppercase; letter-spacing:.04em; }}
  .section-title {{ font-size:1.05rem; font-weight:700; color:{C_DARK};
                    border-bottom:2px solid {C_PARBM}; padding-bottom:6px; margin:28px 0 14px; }}
  .nav-tabs .nav-link.active {{ border-bottom:3px solid {C_PARBM}; font-weight:600; color:{C_DARK}; }}
  .nav-tabs .nav-link {{ color:#6c757d; }}
  .pub-card {{ background:#fff; border-radius:8px; padding:12px 16px;
               margin-bottom:8px; box-shadow:0 1px 3px rgba(0,0,0,.07); }}
  .pub-card .title {{ font-weight:600; font-size:0.92rem; }}
  .pub-card .meta {{ font-size:0.8rem; color:#6c757d; }}
  .pw-badge {{ display:inline-block; background:#e3f2fd; color:#1565c0;
               border-radius:4px; padding:2px 8px; font-size:0.8rem; margin:2px; }}
  .pw-badge.kegg {{ background:#e8f5e9; color:#2e7d32; }}
  table.site-tbl {{ font-size:0.82rem; }}
  table.site-tbl th {{ background:{C_DARK}; color:#fff; position:sticky; top:0; }}
  .badge-high {{ background:{C_PARBM}; color:#fff; border-radius:4px; padding:1px 7px; font-size:0.78rem; }}
  .badge-med  {{ background:{C_HC};   color:#fff; border-radius:4px; padding:1px 7px; font-size:0.78rem; }}
  .badge-low  {{ background:#adb5bd; color:#fff; border-radius:4px; padding:1px 7px; font-size:0.78rem; }}
  .context-seq {{ font-family:monospace; background:#f1f3f5; padding:1px 6px; border-radius:3px; font-size:0.85rem; }}
  .footer {{ text-align:center; font-size:0.78rem; color:#adb5bd; margin:40px 0 20px; }}
</style>
</head>
<body>

<!-- ── Header ─────────────────────────────────────────────────────────── -->
<div class="header-bar">
  <h1>{_esc(pname)}</h1>
  <div class="sub">PARylation &amp; PAR-binding Motif Analysis &nbsp;·&nbsp; Generated {now}</div>
  <div class="mt-2">
    <span class="id-badge">UniProt: {_esc(uid)}</span>
    <span class="id-badge">Gene: {_esc(gene)}</span>
    <span class="id-badge">Organism: {_esc(org)}</span>
    <span class="id-badge">Length: {seqlen} aa</span>
    <span class="id-badge">PARBM Detector v3.0</span>
  </div>
</div>

<div class="container-fluid px-4 py-3">

<!-- ── Summary cards ───────────────────────────────────────────────────── -->
{stats_cards}

<!-- ── PAR-binding domains ─────────────────────────────────────────────── -->
<div class="section-title">PAR-binding Domains</div>
<p>{_esc(parbm)}</p>
{domain_table}

<!-- ── Tabs ────────────────────────────────────────────────────────────── -->
<ul class="nav nav-tabs mt-4" id="mainTabs">
  <li class="nav-item"><a class="nav-link active" data-bs-toggle="tab" href="#tab-canvas">Protein Canvas</a></li>
  <li class="nav-item"><a class="nav-link" data-bs-toggle="tab" href="#tab-conf">Confidence Landscape</a></li>
  <li class="nav-item"><a class="nav-link" data-bs-toggle="tab" href="#tab-residue">Site Breakdown</a></li>
  <li class="nav-item"><a class="nav-link" data-bs-toggle="tab" href="#tab-disorder">Disorder</a></li>
  <li class="nav-item"><a class="nav-link" data-bs-toggle="tab" href="#tab-motif">Motif Context</a></li>
  <li class="nav-item"><a class="nav-link" data-bs-toggle="tab" href="#tab-sites">All Sites</a></li>
  <li class="nav-item"><a class="nav-link" data-bs-toggle="tab" href="#tab-pathways">Pathways</a></li>
  <li class="nav-item"><a class="nav-link" data-bs-toggle="tab" href="#tab-lit">Literature</a></li>
  <li class="nav-item"><a class="nav-link" data-bs-toggle="tab" href="#tab-quality">Quality Metrics</a></li>
</ul>

<div class="tab-content pt-3" id="mainTabContent">

  <div class="tab-pane fade show active" id="tab-canvas">
    <p class="text-muted small">Linear protein map showing domains and all predicted PARylation sites.
       Hover for details. Drag to zoom, double-click to reset.</p>
    {canvas_div}
  </div>

  <div class="tab-pane fade" id="tab-conf">
    <p class="text-muted small">Gaussian-smoothed confidence landscape. Dashed line = high-confidence
       threshold (0.70). Domain regions shaded.</p>
    {conf_div}
  </div>

  <div class="tab-pane fade" id="tab-residue">
    {residue_div}
  </div>

  <div class="tab-pane fade" id="tab-disorder">
    <p class="text-muted small">Per-residue intrinsic disorder score. Scores &gt;0.5 indicate
       disordered regions (typically preferred for PARylation). PARylation sites overlaid as markers.</p>
    {disorder_div}
  </div>

  <div class="tab-pane fade" id="tab-motif">
    {motif_div}
  </div>

  <div class="tab-pane fade" id="tab-sites">
    <p class="text-muted small">All predicted PARylation sites. Use browser search (Ctrl+F) to filter.</p>
    {site_table}
  </div>

  <div class="tab-pane fade" id="tab-pathways">
    {pathways_html}
  </div>

  <div class="tab-pane fade" id="tab-lit">
    {literature_html}
  </div>

  <div class="tab-pane fade" id="tab-quality">
    {quality_html}
  </div>

</div><!-- tab-content -->

<div class="footer">
  PARBM Detector v3.0 &nbsp;·&nbsp; Data sources: UniProt, InterPro, PHOSPHO.ELM, IUPred2A, EuropePMC, Reactome, KEGG
</div>
</div><!-- container -->

<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>"""

    # ── Summary stat cards ────────────────────────────────────────────────────

    def _stats_cards(self, result: Dict, enrich: Dict) -> str:
        s      = result.get("summary", {})
        hc     = s.get("high_confidence_sites_count", 0)
        total  = s.get("total_paryation_sites", 0)
        parbm  = s.get("parbm_domains_count", 0)
        top    = s.get("top_site", "N/A")
        dis    = enrich.get("disorder", {})
        dis_f  = f"{dis.get('disordered_fraction', 0)*100:.0f}%" if dis else "—"
        lit_n  = len(enrich.get("literature", []))
        ptm_n  = len(enrich.get("known_ptm_sites", []))
        return f"""
<div class="row g-3 mt-2">
  <div class="col-6 col-md-3">
    <div class="stat-card">
      <div class="val">{total}</div>
      <div class="lbl">Predicted PARylation Sites</div>
    </div>
  </div>
  <div class="col-6 col-md-3">
    <div class="stat-card blue">
      <div class="val">{hc}</div>
      <div class="lbl">High-Confidence Sites (≥0.70)</div>
    </div>
  </div>
  <div class="col-6 col-md-3">
    <div class="stat-card teal">
      <div class="val">{parbm}</div>
      <div class="lbl">PAR-Binding Domains</div>
    </div>
  </div>
  <div class="col-6 col-md-3">
    <div class="stat-card purple">
      <div class="val">{dis_f}</div>
      <div class="lbl">Disordered Residues</div>
    </div>
  </div>
  <div class="col-6 col-md-3">
    <div class="stat-card">
      <div class="val" style="font-size:1rem;padding-top:4px">{_esc(top)}</div>
      <div class="lbl">Top PARylation Site</div>
    </div>
  </div>
  <div class="col-6 col-md-3">
    <div class="stat-card blue">
      <div class="val">{ptm_n}</div>
      <div class="lbl">Known PTM Sites (PHOSPHO.ELM)</div>
    </div>
  </div>
  <div class="col-6 col-md-3">
    <div class="stat-card teal">
      <div class="val">{lit_n}</div>
      <div class="lbl">Related Publications</div>
    </div>
  </div>
  <div class="col-6 col-md-3">
    <div class="stat-card purple">
      <div class="val">{result.get('sequence_length', 0)}</div>
      <div class="lbl">Sequence Length (aa)</div>
    </div>
  </div>
</div>"""

    # ── Protein canvas (interactive linear map) ───────────────────────────────

    def _protein_canvas(self, result: Dict) -> str:
        seqlen  = result.get("sequence_length", 1)
        domains = result.get("all_domains", [])
        parbm   = result.get("parbm_domains", [])
        sites   = result.get("paryation_sites", [])
        enrich  = result.get("enrichment", {})
        known   = enrich.get("known_ptm_sites", [])

        fig = go.Figure()

        # backbone
        fig.add_shape(type="line", x0=1, x1=seqlen, y0=0, y1=0,
                      line=dict(color="#868E96", width=5))

        # domain boxes
        y_track = {}
        track_ends: list = []
        for dom in domains:
            placed = False
            for ti, te in enumerate(track_ends):
                if dom.start > te + 2:
                    track_ends[ti] = dom.end
                    y_track[id(dom)] = ti
                    placed = True
                    break
            if not placed:
                track_ends.append(dom.end)
                y_track[id(dom)] = len(track_ends) - 1

        for dom in domains:
            is_p = dom in parbm
            col  = C_PARBM if is_p else C_DOMAIN
            ti   = y_track[id(dom)]
            dy   = (ti + 1) * 0.25 * (1 if ti % 2 == 0 else -1)
            label = dom.name if len(dom.name) <= 22 else dom.name[:20] + "…"
            fig.add_shape(type="rect",
                          x0=dom.start, x1=max(dom.end, dom.start + 3),
                          y0=dy - 0.1, y1=dy + 0.1,
                          fillcolor=col, line=dict(color="white", width=1),
                          opacity=0.88)
            fig.add_annotation(x=(dom.start + dom.end) / 2, y=dy,
                                text=label, showarrow=False,
                                font=dict(size=10, color="white", family="Arial"),
                                bgcolor=col, borderpad=2, opacity=0.0)

        # predicted sites
        for site in sites:
            col = C_PARBM if site.score >= 0.7 else C_HC
            fig.add_trace(go.Scatter(
                x=[site.position], y=[0],
                mode="markers",
                marker=dict(symbol="triangle-down", size=12 if site.score >= 0.7 else 7,
                            color=col, line=dict(color="white", width=0.5)),
                name=f"{site.residue}{site.position} ({site.motif_type}, {site.score:.2f})",
                hovertemplate=(
                    f"<b>{site.residue}{site.position}</b><br>"
                    f"Motif: {site.motif_type}<br>"
                    f"Score: {site.score:.3f}<br>"
                    f"Category: {'High-conf' if site.score >= 0.7 else 'Moderate'}"
                    "<extra></extra>"
                ),
                showlegend=False,
            ))

        # known PTM sites (if any)
        for ks in known:
            if ks.get("position"):
                fig.add_trace(go.Scatter(
                    x=[ks["position"]], y=[0],
                    mode="markers",
                    marker=dict(symbol="star", size=11, color="#FFD700",
                                line=dict(color="black", width=0.8)),
                    hovertemplate=(
                        f"<b>Known PTM: {ks.get('residue','?')}{ks['position']}</b><br>"
                        f"Modification: {ks.get('modification','')}<br>"
                        f"Source: PHOSPHO.ELM"
                        "<extra></extra>"
                    ),
                    showlegend=False,
                ))

        fig.update_layout(
            height=300,
            margin=dict(l=40, r=40, t=30, b=50),
            plot_bgcolor="white",
            paper_bgcolor="white",
            xaxis=dict(title="Residue position", showgrid=True,
                       gridcolor="#DEE2E6", zeroline=False),
            yaxis=dict(showticklabels=False, zeroline=False,
                       showgrid=False, range=[-0.5, 0.7]),
            font=dict(family="Segoe UI", color=C_DARK),
        )
        return pio.to_html(fig, full_html=False, include_plotlyjs=False)

    # ── Confidence landscape ──────────────────────────────────────────────────

    def _confidence_figure(self, result: Dict) -> str:
        import numpy as np
        seq    = result.get("sequence", "")
        seqlen = result.get("sequence_length", 1)
        sites  = result.get("paryation_sites", [])
        domains = result.get("all_domains", [])
        parbm   = result.get("parbm_domains", [])

        x = list(range(1, seqlen + 1))
        score_arr = np.zeros(seqlen)
        sigma = max(seqlen / 120, 3)
        for site in sites:
            idx = site.position - 1
            for j in range(max(0, idx - int(4 * sigma)),
                           min(seqlen, idx + int(4 * sigma) + 1)):
                score_arr[j] = max(score_arr[j],
                                   site.score * np.exp(-0.5 * ((j - idx) / sigma) ** 2))

        fig = go.Figure()

        # domain shading
        for dom in domains:
            col = "rgba(230,57,70,0.07)" if dom in parbm else "rgba(69,123,157,0.07)"
            fig.add_vrect(x0=dom.start, x1=dom.end, fillcolor=col, line_width=0)

        fig.add_trace(go.Scatter(
            x=x, y=list(score_arr),
            fill="tozeroy", fillcolor="rgba(230,57,70,0.15)",
            line=dict(color=C_PARBM, width=1.5),
            name="Confidence envelope",
            hovertemplate="Position %{x}<br>Score: %{y:.3f}<extra></extra>",
        ))

        # individual site stems
        for site in sites:
            col = C_PARBM if site.score >= 0.7 else C_HC
            fig.add_trace(go.Scatter(
                x=[site.position, site.position], y=[0, site.score],
                mode="lines", line=dict(color=col, width=1),
                showlegend=False,
                hoverinfo="skip",
            ))
            fig.add_trace(go.Scatter(
                x=[site.position], y=[site.score],
                mode="markers",
                marker=dict(color=col, size=5 if site.score >= 0.7 else 3,
                            line=dict(color="white", width=0.4)),
                showlegend=False,
                hovertemplate=(
                    f"<b>{site.residue}{site.position}</b><br>"
                    f"{site.motif_type}<br>Score: {site.score:.3f}"
                    "<extra></extra>"
                ),
            ))

        fig.add_hline(y=0.7, line_dash="dash", line_color=C_PARBM,
                      annotation_text="High-confidence (0.70)",
                      annotation_position="top right", line_width=1.2)
        fig.add_hline(y=0.5, line_dash="dot", line_color="#ADB5BD",
                      annotation_text="Moderate (0.50)",
                      annotation_position="bottom right", line_width=0.8)

        fig.update_layout(
            height=320, margin=dict(l=50, r=40, t=20, b=50),
            plot_bgcolor="white", paper_bgcolor="white",
            xaxis=dict(title="Residue position", showgrid=True, gridcolor="#DEE2E6"),
            yaxis=dict(title="Confidence score", range=[0, 1.05],
                       showgrid=True, gridcolor="#DEE2E6"),
            showlegend=False,
            font=dict(family="Segoe UI", color=C_DARK),
        )
        return pio.to_html(fig, full_html=False, include_plotlyjs=False)

    # ── Residue type breakdown ────────────────────────────────────────────────

    def _residue_breakdown(self, result: Dict) -> str:
        sites  = result.get("paryation_sites", [])
        if not sites:
            return "<p class='text-muted'>No predicted sites.</p>"

        from collections import Counter
        aa_count  = Counter(s.residue for s in sites)
        aa_hc     = Counter(s.residue for s in sites if s.score >= 0.7)
        aa_order  = ["S", "T", "E", "D", "Y", "K", "R"]
        colors    = [C_PARBM, C_HC, C_TEAL, C_DOMAIN, C_PURPLE, "#20B2AA", "#FF8C00"]

        fig = make_subplots(rows=1, cols=2,
                            subplot_titles=("All Sites by Residue Type",
                                            "High-Confidence Sites Only"))
        for col_idx, counter in enumerate([aa_count, aa_hc], start=1):
            labels = [aa for aa in aa_order if counter.get(aa, 0) > 0]
            values = [counter.get(aa, 0) for aa in labels]
            cols_  = [colors[aa_order.index(aa)] for aa in labels]
            fig.add_trace(
                go.Bar(x=labels, y=values, marker_color=cols_,
                       text=values, textposition="outside",
                       hovertemplate="%{x}: %{y} sites<extra></extra>"),
                row=1, col=col_idx,
            )

        fig.update_layout(
            height=320, margin=dict(l=40, r=40, t=50, b=40),
            plot_bgcolor="white", paper_bgcolor="white",
            showlegend=False,
            font=dict(family="Segoe UI", color=C_DARK),
        )
        fig.update_yaxes(title_text="Site count", showgrid=True, gridcolor="#DEE2E6")
        fig.update_xaxes(title_text="Residue")

        # confidence tier table
        tiers = {"High (≥0.70)": sum(1 for s in sites if s.score >= 0.7),
                 "Medium (0.50–0.69)": sum(1 for s in sites if 0.5 <= s.score < 0.7),
                 "Low (<0.50)": sum(1 for s in sites if s.score < 0.5)}
        tier_rows = "".join(
            f"<tr><td>{k}</td><td>{v}</td>"
            f"<td>{v/len(sites)*100:.1f}%</td></tr>"
            for k, v in tiers.items()
        )
        tier_table = f"""
<div class="mt-3">
  <h6 class="fw-bold">Confidence Tiers</h6>
  <table class="table table-sm table-bordered" style="max-width:400px">
    <thead class="table-dark"><tr><th>Tier</th><th>Count</th><th>%</th></tr></thead>
    <tbody>{tier_rows}</tbody>
  </table>
</div>"""

        return pio.to_html(fig, full_html=False, include_plotlyjs=False) + tier_table

    # ── Disorder figure ───────────────────────────────────────────────────────

    def _disorder_figure(self, result: Dict) -> str:
        enrich = result.get("enrichment", {})
        dis    = enrich.get("disorder", {})
        seq    = result.get("sequence", "")
        sites  = result.get("high_confidence_sites", [])
        seqlen = result.get("sequence_length", 1)

        scores = dis.get("scores", [])
        method = dis.get("method", "unknown")
        mean_d = dis.get("mean_disorder", 0)
        frac_d = dis.get("disordered_fraction", 0)

        if not scores:
            return ("<p class='text-muted'>Disorder data unavailable "
                    "(enrichment not run or IUPred2A unreachable).</p>")

        x = list(range(1, len(scores) + 1))
        fig = go.Figure()

        # disorder score
        fig.add_trace(go.Scatter(
            x=x, y=scores,
            fill="tozeroy", fillcolor="rgba(155,93,229,0.15)",
            line=dict(color=C_PURPLE, width=1.5),
            name=f"Disorder ({method})",
            hovertemplate="Position %{x}<br>Disorder: %{y:.3f}<extra></extra>",
        ))

        # threshold
        fig.add_hline(y=0.5, line_dash="dash", line_color="#6C757D",
                      annotation_text="Disorder threshold (0.5)", line_width=1)

        # site markers
        site_x = [s.position for s in sites if s.position <= len(scores)]
        site_y = [scores[s.position - 1] for s in sites if s.position <= len(scores)]
        fig.add_trace(go.Scatter(
            x=site_x, y=site_y,
            mode="markers",
            marker=dict(symbol="triangle-down", size=10, color=C_PARBM,
                        line=dict(color="white", width=0.5)),
            name="High-conf PARylation site",
            hovertemplate="Site at %{x}<br>Disorder: %{y:.3f}<extra></extra>",
        ))

        fig.update_layout(
            height=320, margin=dict(l=50, r=40, t=20, b=50),
            plot_bgcolor="white", paper_bgcolor="white",
            xaxis=dict(title="Residue position", showgrid=True, gridcolor="#DEE2E6"),
            yaxis=dict(title="Disorder score", range=[0, 1.05],
                       showgrid=True, gridcolor="#DEE2E6"),
            legend=dict(orientation="h", yanchor="bottom", y=1.02),
            font=dict(family="Segoe UI", color=C_DARK),
            annotations=[dict(
                x=0.01, y=0.98, xref="paper", yref="paper",
                text=f"Mean disorder: {mean_d:.2f} | Disordered residues: {frac_d*100:.0f}%",
                showarrow=False, font=dict(size=10, color="#6C757D"),
                align="left",
            )],
        )
        return pio.to_html(fig, full_html=False, include_plotlyjs=False)

    # ── Motif context section ─────────────────────────────────────────────────

    def _motif_context_section(self, enrich: Dict) -> str:
        contexts = enrich.get("site_contexts", [])
        if not contexts:
            return "<p class='text-muted'>Run enrichment to see motif context analysis.</p>"

        rows = ""
        for ctx in contexts[:20]:
            seq  = ctx["context"]
            mid  = len(seq) // 2
            coloured = (
                f'<span class="text-muted">{_esc(seq[:mid])}</span>'
                f'<span style="color:{C_PARBM};font-weight:700">{_esc(seq[mid:mid+1])}</span>'
                f'<span class="text-muted">{_esc(seq[mid+1:])}</span>'
            )
            tier = ("High" if ctx["score"] >= 0.7 else "Med")
            badge = (f'<span class="badge-high">{tier}</span>'
                     if tier == "High" else f'<span class="badge-med">{tier}</span>')
            rows += (f"<tr><td>{ctx['position']}</td>"
                     f"<td>{ctx['residue']}</td>"
                     f"<td><span class='context-seq'>{coloured}</span></td>"
                     f"<td>{ctx['score']:.3f}</td>"
                     f"<td>{badge}</td></tr>")

        return f"""
<h6 class="fw-bold mt-2">±4 aa Sequence Context Around High-Confidence Sites</h6>
<p class="text-muted small">Central residue highlighted in red.</p>
<div style="max-height:380px;overflow-y:auto">
<table class="table table-sm table-hover site-tbl">
  <thead><tr><th>Position</th><th>Residue</th><th>Context (±4 aa)</th>
  <th>Score</th><th>Tier</th></tr></thead>
  <tbody>{rows}</tbody>
</table></div>"""

    # ── All sites table ───────────────────────────────────────────────────────

    def _site_table(self, result: Dict) -> str:
        sites  = result.get("paryation_sites", [])
        enrich = result.get("enrichment", {})
        known  = {s["position"] for s in enrich.get("known_ptm_sites", [])
                  if s.get("position")}
        dis_sc = enrich.get("disorder", {}).get("scores", [])

        rows = ""
        for site in sites:
            tier = ("High" if site.score >= 0.7 else
                    "Medium" if site.score >= 0.5 else "Low")
            badge = (f'<span class="badge-high">High</span>' if tier == "High"
                     else f'<span class="badge-med">Med</span>' if tier == "Medium"
                     else f'<span class="badge-low">Low</span>')
            exp = "✓" if site.position in known else "—"
            dis_val = (f"{dis_sc[site.position-1]:.2f}"
                       if dis_sc and site.position <= len(dis_sc) else "—")
            rows += (f"<tr><td>{site.position}</td><td><b>{site.residue}</b></td>"
                     f"<td>{_esc(site.motif_type)}</td>"
                     f"<td>{site.score:.3f}</td><td>{badge}</td>"
                     f"<td>{dis_val}</td><td>{exp}</td></tr>")

        return f"""
<div style="max-height:500px;overflow-y:auto">
<table class="table table-sm table-hover site-tbl">
  <thead><tr><th>Position</th><th>Residue</th><th>Motif</th>
  <th>Score</th><th>Tier</th><th>Disorder</th><th>Known PTM</th></tr></thead>
  <tbody>{rows}</tbody>
</table></div>"""

    # ── Domain table ──────────────────────────────────────────────────────────

    def _domain_table(self, result: Dict) -> str:
        domains = result.get("all_domains", [])
        parbm   = result.get("parbm_domains", [])
        if not domains:
            return "<p class='text-muted small'>No domain annotations (sequence input or API unavailable).</p>"

        rows = ""
        for dom in domains:
            is_p = dom in parbm
            badge = (f'<span style="background:{C_PARBM};color:white;'
                     f'border-radius:4px;padding:1px 7px;font-size:0.78rem">'
                     f'PAR-binding ✓</span>' if is_p else "—")
            rows += (f"<tr><td>{_esc(dom.name)}</td>"
                     f"<td><code>{dom.interpro_id}</code></td>"
                     f"<td>{dom.start}–{dom.end}</td><td>{badge}</td></tr>")

        return f"""
<table class="table table-sm table-bordered mt-2" style="max-width:750px">
  <thead class="table-dark"><tr><th>Domain</th><th>InterPro ID</th>
  <th>Position</th><th>PAR-binding?</th></tr></thead>
  <tbody>{rows}</tbody>
</table>"""

    # ── Pathways ──────────────────────────────────────────────────────────────

    def _pathways_section(self, enrich: Dict) -> str:
        pws    = enrich.get("pathways", {})
        react  = pws.get("reactome", [])
        kegg   = pws.get("kegg", [])
        if not react and not kegg:
            return "<p class='text-muted'>No pathway data (enrichment not run or unavailable).</p>"

        html = '<div class="section-title">Biological Pathways</div>'
        if react:
            html += "<h6 class='fw-bold'>Reactome</h6>"
            html += "".join(
                f'<span class="pw-badge">{_esc(p["name"])}'
                f' <small>({_esc(p["id"])})</small></span>'
                for p in react
            )
        if kegg:
            html += "<h6 class='fw-bold mt-3'>KEGG</h6>"
            html += "".join(
                f'<span class="pw-badge kegg">{_esc(p["name"])}'
                f' <small>({_esc(p["id"])})</small></span>'
                for p in kegg
            )
        return html

    # ── Literature ────────────────────────────────────────────────────────────

    def _literature_section(self, enrich: Dict) -> str:
        pubs = enrich.get("literature", [])
        if not pubs:
            return "<p class='text-muted'>No literature data (enrichment not run or no results).</p>"

        html = '<div class="section-title">Related Publications</div>'
        for pub in pubs:
            pmid = pub.get("pmid", "")
            pm_link = (f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}" '
                       f'target="_blank">PubMed ↗</a>' if pmid else "")
            doi = pub.get("doi", "")
            doi_link = (f' &nbsp;·&nbsp; <a href="https://doi.org/{doi}" '
                        f'target="_blank">DOI ↗</a>' if doi else "")
            html += f"""
<div class="pub-card">
  <div class="title">{_esc(pub.get('title',''))}</div>
  <div class="meta mt-1">
    {_esc(pub.get('authors','')[:100])} &nbsp;·&nbsp;
    <i>{_esc(pub.get('journal',''))}</i> {_esc(pub.get('year',''))}
    &nbsp;·&nbsp; Cited by: {pub.get('cited_by',0)}
    &nbsp;·&nbsp; {pm_link}{doi_link}
  </div>
</div>"""
        return html

    # ── Quality metrics ───────────────────────────────────────────────────────

    def _quality_metrics_section(self, result: Dict) -> str:
        from .enrichment import EnrichmentEngine
        metrics = EnrichmentEngine.quality_metrics(result)

        if not metrics.get("available"):
            return f"""
<div class="section-title">Prediction Quality Metrics</div>
<p class="text-muted">{metrics.get('note','')}</p>
<p class="text-muted small">
  Run <code>detector.enrich_result(result)</code> and ensure PHOSPHO.ELM
  returned known sites for this protein to compute sensitivity / precision.
</p>"""

        s   = metrics["sensitivity"]
        pr  = metrics["precision"]
        f1  = metrics["f1_score"]
        tp, fp, fn = metrics["tp"], metrics["fp"], metrics["fn"]

        def bar(pct, col):
            return (f'<div style="background:#e9ecef;border-radius:4px;height:12px;'
                    f'width:200px;display:inline-block;vertical-align:middle">'
                    f'<div style="width:{min(pct,100):.0f}%;background:{col};'
                    f'height:100%;border-radius:4px"></div></div>')

        return f"""
<div class="section-title">Prediction Quality Metrics</div>
<p class="text-muted small">{metrics.get('note','')}</p>
<table class="table table-sm" style="max-width:500px">
<thead class="table-dark"><tr><th>Metric</th><th>Value</th><th></th></tr></thead>
<tbody>
  <tr><td>Sensitivity (Recall)</td><td><b>{s:.1f}%</b></td>
      <td>{bar(s, C_TEAL)}</td></tr>
  <tr><td>Precision (PPV)</td><td><b>{pr:.1f}%</b></td>
      <td>{bar(pr, C_DOMAIN)}</td></tr>
  <tr><td>F1 Score</td><td><b>{f1:.3f}</b></td>
      <td>{bar(f1*100, C_PARBM)}</td></tr>
  <tr><td>True Positives</td><td>{tp}</td><td></td></tr>
  <tr><td>False Positives</td><td>{fp}</td><td></td></tr>
  <tr><td>False Negatives</td><td>{fn}</td><td></td></tr>
  <tr><td>Known sites (PHOSPHO.ELM)</td><td>{metrics['known_sites']}</td><td></td></tr>
  <tr><td>Predicted sites</td><td>{metrics['pred_sites']}</td><td></td></tr>
</tbody>
</table>"""

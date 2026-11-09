"""Build the five portfolio diagrams (SVG with editable text, plus PNG).

Usage: python portfolio/src/make_diagrams.py
Writes portfolio/diagrams/D1..D5_*.{svg,png}. Reads values only through values.py
(frozen outputs, read-only). Fits no model and computes no statistic.
"""
import math
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

import values  # noqa: E402

OUT = values.ROOT / "portfolio" / "diagrams"

INK, MUTED, FAINT, GRID, PANEL = "#1f2933", "#52606d", "#8a94a3", "#e3e8ee", "#f5f7fa"
PRIMARY, S7, FIVE, APPARENT = "#2a78d6", "#eb6834", "#1baf7a", "#9aa5b1"
TAG = {"Observed": "#2a78d6", "Literature": "#7b5ea7", "Hypothesis": "#b7791f",
       "Not tested": "#7b8794", "Textbook": "#52606d", "Conceptual": "#7b8794"}

plt.rcParams.update({
    "font.family": ["Arial", "DejaVu Sans"], "font.size": 11, "svg.fonttype": "none",
    "svg.hashsalt": "gse39582-portfolio", "axes.edgecolor": FAINT, "axes.labelcolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
})
W, H = 13.333, 7.5


def canvas(title, subtitle=None):
    fig = plt.figure(figsize=(W, H), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    ax.text(0.45, H - 0.55, title, fontsize=22, fontweight="bold", va="top")
    if subtitle:
        ax.text(0.45, H - 1.05, subtitle, fontsize=12.5, color=MUTED, va="top")
    return fig, ax


def wrap(text, width_in, fs):
    chars = max(8, int(width_in * 72 / (fs * 0.53)))
    return "\n".join(textwrap.fill(p, chars) for p in text.split("\n"))


def box(ax, x, y, w, h, fc="white", ec=GRID, lw=1.2, ls="-", r=0.08):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                fc=fc, ec=ec, lw=lw, ls=ls))


def tag(ax, x, y, label, fs=8.5):
    c = TAG[label]
    t = ax.text(x, y, label.upper(), fontsize=fs, color="white", fontweight="bold", va="center", ha="left",
                bbox=dict(boxstyle="round,pad=0.25,rounding_size=0.15", fc=c, ec=c))
    return t


def footer(ax, text):
    ax.text(0.45, 0.28, wrap(text, W - 0.9, 9), fontsize=9, color=MUTED, va="bottom")


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    meta = {"Date": None} if True else {}
    fig.savefig(OUT / f"{name}.svg", metadata=meta)
    fig.savefig(OUT / f"{name}.png", dpi=200, metadata={"Software": None})
    plt.close(fig)


def arrow(ax, x0, y0, x1, y1, color=FAINT):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.6, shrinkA=0, shrinkB=0))


# ---------------------------------------------------------------- D1
def d1(V):
    fig, ax = canvas("Four biological axes of colorectal cancer, one question",
                     "Do four predefined pathway-expression scores add recurrence-discrimination information beyond clinical variables?")
    lanes = [
        ("Wnt/β-catenin", "Initiating driver pathway",
         "Loss of APC (or activating CTNNB1 mutation) stabilises β-catenin, which with TCF/LEF drives target genes such as MYC and AXIN2.",
         "Wnt-pathway genes are genomically altered in 93 % of tumours (R6). Marked WNT/MYC activation characterises CMS2 (R1). Alteration is not measured activation.",
         f"Mean standardised mRNA of {V['panel']['sizes']['Wnt/β-catenin']} member genes"),
        ("E2F targets", "Proliferation",
         "CDK–cyclin phosphorylation of RB releases E2F factors, which activate G1/S and DNA-replication genes.",
         "Evidence linking proliferation to outcome is mixed: lower in relapsed primaries (R4), but high Ki-67 with poorer survival (R5).",
         f"Mean standardised mRNA of {V['panel']['sizes']['E2F targets']} member genes"),
        ("p53 pathway", "DNA-damage checkpoint",
         "Cellular stress stabilises p53, which induces arrest (CDKN1A), feedback (MDM2) and apoptosis genes.",
         "TP53 is mutated in about 60 % of non-hypermutated and about 20 % of hypermutated tumours (R6).",
         f"Mean standardised mRNA of {V['panel']['sizes']['p53 pathway']} member genes; not a TP53 mutation or function assay"),
        ("EMT", "Mesenchymal / stromal state",
         "Epithelial cells lose E-cadherin junctions and gain mesenchymal genes (vimentin, fibronectin, collagens) under SNAIL, ZEB and TWIST.",
         "Mesenchymal CMS4 had worse relapse-free survival (R1; GSE39582 contributed, so not independent). Bulk mesenchymal signal is mostly stromal (R2, R3).",
         f"Mean standardised mRNA of {V['panel']['sizes']['EMT']} member genes; not a direct measure of EMT in tumour cells"),
    ]
    x0, gap, top = 0.45, 0.18, H - 1.5
    lw = (W - 0.9 - 3 * gap) / 4
    rows = [("Textbook", 1.25), ("Literature", 1.45), ("Observed", 1.0)]
    for i, (name, axis_, mech, lit, score) in enumerate(lanes):
        x = x0 + i * (lw + gap)
        box(ax, x, top - 0.75, lw, 0.75, fc=INK, ec=INK)
        ax.text(x + 0.18, top - 0.27, name, fontsize=15, color="white", fontweight="bold", va="center")
        ax.text(x + 0.18, top - 0.56, axis_, fontsize=10.5, color="#d9e2ec", va="center")
        y = top - 0.85
        for (lab, h), txt in zip(rows, (mech, lit, score)):
            box(ax, x, y - h, lw, h, fc="white", ec=GRID)
            tag(ax, x + 0.15, y - 0.22, lab if lab != "Observed" else "Observed", fs=7.5)
            if lab == "Observed":
                ax.text(x + 1.25, y - 0.22, "project score", fontsize=8, color=MUTED, va="center")
            ax.text(x + 0.15, y - 0.45, wrap(txt, lw - 0.3, 10), fontsize=10, va="top", linespacing=1.3)
            y -= h + 0.1
    box(ax, 0.45, 0.66, W - 0.9, 0.52, fc=PANEL, ec=PANEL)
    ax.text(0.65, 0.92, wrap("Every project score is the mean training-fold-standardised mRNA of a fixed gene set in bulk tumour tissue. "
                            "It summarises expression. It does not measure pathway activation, protein activity, mutation status or the cell type the mRNA came from.",
                            W - 1.3, 11), fontsize=11, va="center")
    footer(ax, "Conceptual schematic: lanes are parallel axes, not a causal sequence. Mechanisms are textbook background. "
               f"Context: verified references R1–R8 (docs/BIOLOGICAL_EVIDENCE_AND_CITATION_AUDIT.md). Gene counts: MSigDB Hallmark {V['panel']['version']}, config/pathways/m2_panel_4set.json.")
    save(fig, "D1_biological_rationale")


# ---------------------------------------------------------------- D2
def d2(V):
    P, pn = V["probes"], V["panel"]
    fig, ax = canvas("From probes to pathway scores",
                     "Expression processing and feature engineering, fixed before any outcome was examined")
    steps = [
        (f"{P['total']:,}", "Affymetrix probe sets", "HG-U133 Plus 2.0, deposited RMA + ComBat"),
        (f"{P['single_gene']:,}", "probes map to one gene", f"{P['unmapped']:,} unmapped and {P['ambiguous']:,} ambiguous excluded"),
        (f"{P['genes']:,}", "unique genes", f"median of probes per gene ({P['multi_probe_genes']:,} genes have several)"),
        (f"{pn['union4']}", "pathway-member genes", "EMT {EMT} · E2F {E2F} · Wnt {W} · p53 {P53}".format(
            EMT=pn["sizes"]["EMT"], E2F=pn["sizes"]["E2F targets"], W=pn["sizes"]["Wnt/β-catenin"], P53=pn["sizes"]["p53 pathway"])),
        ("4", "pathway scores per patient", "mean of member-gene z-scores, standardised inside each training fold"),
    ]
    n = len(steps)
    bw, gap, y, bh = 2.25, 0.32, 3.55, 2.2
    x0 = (W - (n * bw + (n - 1) * gap)) / 2
    for i, (num, lab, note) in enumerate(steps):
        x = x0 + i * (bw + gap)
        last = i == n - 1
        box(ax, x, y, bw, bh, fc=PRIMARY if last else "white", ec=PRIMARY, lw=1.6)
        c = "white" if last else INK
        ax.text(x + bw / 2, y + bh - 0.55, num, fontsize=28, fontweight="bold", ha="center", va="center", color=c)
        ax.text(x + bw / 2, y + bh - 1.08, wrap(lab, bw - 0.2, 12), fontsize=12, ha="center", va="center", color=c, fontweight="bold")
        ax.text(x + bw / 2, y + 0.55, wrap(note, bw - 0.25, 9.5), fontsize=9.5, ha="center", va="center",
                color="white" if last else MUTED, linespacing=1.25)
        if not last:
            arrow(ax, x + bw + 0.04, y + bh / 2, x + bw + gap - 0.04, y + bh / 2, color=PRIMARY)
    # panel selection strip
    box(ax, 0.45, 0.95, W - 0.9, 2.15, fc=PANEL, ec=PANEL)
    ax.text(0.7, 2.83, "How the four sets were chosen (outcome-blind)", fontsize=13, fontweight="bold", va="center")
    sel = [
        ("10 candidate Hallmark sets", "score-matrix condition number 77"),
        ("Correlation clusters", "six sets on one mesenchymal axis (r 0.60–0.91 with EMT)"),
        ("Primary: 4 sets", f"EMT, E2F, Wnt, p53 · condition number {pn['cond4']:.2f}, max |r| {pn['maxr4']:.2f}"),
        ("Sensitivity only: 5 sets", f"+ PI3K/AKT/mTOR · condition number {pn['cond5']:.2f}, max |r| {pn['maxr5']:.2f}"),
    ]
    sw, sg = 2.85, 0.27
    sx0 = 0.7
    for i, (a, b) in enumerate(sel):
        x = sx0 + i * (sw + sg)
        box(ax, x, 1.12, sw, 1.35, fc="white", ec=FIVE if i == 3 else (PRIMARY if i == 2 else GRID), lw=1.4)
        ax.text(x + 0.15, 2.2, a, fontsize=11.5, fontweight="bold", va="center")
        ax.text(x + 0.15, 1.95, wrap(b, sw - 0.3, 9.5), fontsize=9.5, color=MUTED, va="top", linespacing=1.25)
        if i < 2:
            arrow(ax, x + sw + 0.03, 1.8, x + sw + sg - 0.03, 1.8)
    footer(ax, "Condition numbers are descriptive full-cohort score-matrix values (Phase 6C), not the Phase 9 model-fit values. "
               "Selection used gene-set metadata, literature and correlation structure only; no outcome consulted; frozen 2 Nov 2026, 21:13 UTC. "
               "Sources: reports/figures/final/figure_values.json (fig2.probes); config/pathways/m2_panel_4set.json and m2_panel_5set.json; docs/PHASE6C_PATHWAY_PANEL_REVIEW.md §5.")
    save(fig, "D2_processing_feature_engineering")


# ---------------------------------------------------------------- D3
def fmt_q(q):
    if q < 1e-3:
        m, ex = f"{q:.1e}".split("e")
        return rf"${m}\times10^{{{int(ex)}}}$"
    return f"{q:.4f}"


def d3(V):
    fig, ax = canvas("Pathway by pathway: estimate, uncertainty and what follows",
                     "Primary four-set Cox model (median representation) with representation and batch checks")
    lo, hi = math.log(0.7), math.log(1.35)
    interp = {
        "EMT": [("Observed", "HR above 1 in every reported fit; interval includes 1"),
                ("Literature", "Compatible with stromal/mesenchymal poor-prognosis reports (R1–R3)"),
                ("Not tested", "Tumour-cell EMT versus stromal expression")],
        "E2F targets": [("Observed", "HR below 1 in every reported fit; no primary interval excludes 1"),
                        ("Literature", "Proliferation–outcome evidence is mixed (R4 vs R5)"),
                        ("Not tested", "Any mechanism; the direction is not established biology")],
        "Wnt/β-catenin": [("Observed", "Near zero; least representation-stable score (r 0.858)"),
                          ("Literature", "Alteration in 93 % of tumours (R6) is not measured activation"),
                          ("Not tested", "No explanation is offered for the null estimate")],
        "p53 pathway": [("Observed", "Near zero; strongest MMR association of the four scores"),
                        ("Literature", "TP53 mutation rates differ by hypermutation status (R6)"),
                        ("Not tested", "TP53 status; direction or mechanism of the MMR link")],
    }
    gap = 0.18
    cw = (W - 0.9 - 3 * gap) / 4
    T, bottom = H - 1.45, 0.95
    for i, p in enumerate(V["pathways"]):
        x = 0.45 + i * (cw + gap)
        box(ax, x, bottom, cw, T - bottom, fc="white", ec=GRID)
        ax.text(x + 0.15, T - 0.3, p["label"], fontsize=15, fontweight="bold", va="center")
        ax.text(x + 0.15, T - 0.6, f"{p['genes']} genes · coherence r̄ {p['coherence']:.3f}", fontsize=9.5, color=MUTED, va="center")
        (h, l, u), (h7, l7, u7) = p["hr"], p["hr_s7"]
        ax.text(x + 0.15, T - 0.9, f"● Primary HR {h:.3f} ({l:.3f}–{u:.3f})", fontsize=9, color=PRIMARY, va="center")
        ax.text(x + 0.15, T - 1.13, f"■ S-7 HR {h7:.3f} ({l7:.3f}–{u7:.3f})", fontsize=9, color=S7, va="center")
        fx0, fx1 = x + 0.3, x + cw - 0.3
        X = lambda v: fx0 + (math.log(v) - lo) / (hi - lo) * (fx1 - fx0)  # noqa: E731
        ax.plot([X(1), X(1)], [T - 1.8, T - 1.3], color=INK, lw=0.9, ls=(0, (3, 3)))
        for k, (hh, ll, uu) in enumerate((p["hr"], p["hr_s7"])):
            yy = T - 1.42 - k * 0.24
            col = PRIMARY if k == 0 else S7
            ax.plot([X(ll), X(uu)], [yy, yy], color=col, lw=2.2, solid_capstyle="round")
            ax.plot([X(hh)], [yy], marker="o" if k == 0 else "s", ms=6.5, color=col, mec="white", mew=1.1)
        for v in (0.8, 1.0, 1.2):
            ax.text(X(v), T - 1.93, f"{v:.1f}", fontsize=8, color=FAINT, ha="center", va="center")
        ax.text(x + 0.15, T - 2.25, f"Batch-stratified HR {p['hr_strat']:.3f} (point estimate)", fontsize=9, va="center")
        bx0, bx1 = x + 0.3, x + cw - 0.3
        BX = lambda v: bx0 + (v + 0.3) / 0.55 * (bx1 - bx0)  # noqa: E731
        yb = T - 2.55
        ax.plot([BX(0), BX(0)], [yb - 0.16, yb + 0.16], color=INK, lw=0.9, ls=(0, (3, 3)))
        ax.add_patch(Rectangle((BX(p["lobo_min"]), yb - 0.08), max(0.03, BX(p["lobo_max"]) - BX(p["lobo_min"])), 0.16,
                               fc="#dbe8f8", ec=PRIMARY, lw=1))
        ax.plot([BX(p["coef_all"])], [yb], "o", ms=5, color=PRIMARY, mec="white", mew=1)
        sc = p["lobo_sign_changes"]
        ax.text(x + 0.15, T - 2.83, f"LOBO log-HR {p['lobo_min']:+.3f} to {p['lobo_max']:+.3f}; {sc} sign change{'s' if sc != 1 else ''}",
                fontsize=8.5, color=MUTED, va="center")
        ax.text(x + 0.15, T - 3.13, f"MMR association ε² {p['mmr_e2']:.4f}, BH q {fmt_q(p['mmr_q'])}", fontsize=9, va="center")
        iy = T - 3.47
        for lab, txt in interp[p["label"]]:
            tag(ax, x + 0.15, iy, lab, fs=7)
            ax.text(x + 0.15, iy - 0.15, wrap(txt, cw - 0.3, 8.6), fontsize=8.6, va="top", linespacing=1.15)
            iy -= 0.56
    footer(ax, f"Hazard ratios per SD, log scale. Primary: ridge Cox, median representation (λ {V['lambda']['primary']:.3f}). S-7: representative-IQR summary of the same arrays (sensitivity, not replication). "
               f"Wald intervals approximate under ridge. Stratified HRs are point estimates; LOBO omits batches of {V['lobo_sizes'][0]}–{V['lobo_sizes'][1]} patients, so spread mixes batch and sample size. "
               "MMR: exploratory, outcome-free Kruskal–Wallis, n = 410, BH q; an association, not a mechanism. Coherence is not activation. Not a biomarker.")
    save(fig, "D3_pathway_interpretation")


# ---------------------------------------------------------------- D4
def d4(V):
    fig, ax = canvas("What a bulk pathway score may contain",
                     "Possible contributors to one number per patient. None was separated or quantified in this project.")
    inputs = [
        ("Tumour-cell expression programme", "Not tested", "The intended biological signal; cannot be separated from the others in bulk tissue"),
        ("Stromal, immune and endothelial cells", "Not tested", "Cell composition was not measured or deconvolved; an alternative explanation, especially for EMT"),
        ("Molecular subtype", "Observed", "All four scores associate with MMR status (exploratory, outcome-free, n = 410)"),
        ("Processing batch", "Observed", "Characterised, not resolved; residual batch structure after ComBat"),
    ]
    bx, bw, bh = 0.45, 4.6, 1.05
    ys = [5.0, 3.8, 2.6, 1.4]
    sx, sy, sw, sh = 6.15, 2.55, 2.2, 1.6
    for (title, lab, txt), y in zip(inputs, ys):
        box(ax, bx, y, bw, bh, fc="white", ec=TAG[lab] if lab != "Not tested" else FAINT, lw=1.3,
            ls="--" if lab == "Not tested" else "-")
        ax.text(bx + 0.15, y + bh - 0.25, title, fontsize=11.5, fontweight="bold", va="center")
        tag(ax, bx + bw - 1.25, y + bh - 0.25, lab, fs=7)
        ax.text(bx + 0.15, y + bh - 0.5, wrap(txt, bw - 0.3, 9.3), fontsize=9.3, color=MUTED, va="top", linespacing=1.2)
        ax.plot([bx + bw, sx], [y + bh / 2, sy + sh / 2], color=FAINT, lw=1.2)
    box(ax, sx, sy, sw, sh, fc=INK, ec=INK)
    ax.text(sx + sw / 2, sy + sh / 2 + 0.2, "Bulk pathway\nscore", fontsize=15, fontweight="bold", color="white", ha="center", va="center")
    ax.text(sx + sw / 2, sy + 0.3, "one number per patient", fontsize=9, color="#d9e2ec", ha="center", va="center")
    tag(ax, sx + 0.35, sy + sh + 0.3, "Conceptual", fs=7.5)
    ax.text(sx + 1.5, sy + sh + 0.3, "lines are not weighted", fontsize=8.5, color=MUTED, va="center")
    # observed panel
    px, pw = 8.75, W - 0.45 - 8.75
    box(ax, px, 1.15, pw, 4.95, fc=PANEL, ec=PANEL)
    ax.text(px + 0.2, 5.8, "What was observed", fontsize=13, fontweight="bold", va="center")
    m = {p["label"]: p for p in V["pathways"]}
    rows = [
        "MMR association ε² (BH q):",
        f"  EMT {m['EMT']['mmr_e2']:.4f} (1.3×10⁻⁴)",
        f"  E2F {m['E2F targets']['mmr_e2']:.4f} ({m['E2F targets']['mmr_q']:.4f})",
        f"  Wnt {m['Wnt/β-catenin']['mmr_e2']:.4f} ({m['Wnt/β-catenin']['mmr_q']:.4f})",
        f"  p53 {m['p53 pathway']['mmr_e2']:.4f} (4.7×10⁻¹⁶)",
        "Batch/extraction location effects on scores did not survive a test stratified by combat_mod (0 of 10).",
        "combat_mod matches recorded MMR in all 410 samples, but its provenance is unresolved: it is not an MMR variable.",
        "Leaving out one batch still flipped the Wnt and p53 coefficients.",
    ]
    yy = 5.4
    for r in rows:
        t = wrap(r, pw - 0.4, 9.6)
        ax.text(px + 0.2, yy, t, fontsize=9.6, va="top", linespacing=1.25)
        yy -= 0.27 * (t.count("\n") + 1) + 0.1
    footer(ax, "Conceptual schematic: no contribution was estimated, and the connecting lines carry no weights or causal direction. "
               "Sources: reports/phase8c/score_vs_annotation.csv (mmr.status); docs/PHASE8C_FROZEN_PATHWAY_SCORE_ANALYSIS.md §5.3; "
               "docs/PHASE8B_PCA_BATCH_INVESTIGATION.md (combat_mod ↔ MMR, 0 discordant of 410); docs/PHASE9_CORRECTED_FINAL_REPORT.md §12.")
    save(fig, "D4_heterogeneity_limitations")


# ---------------------------------------------------------------- D5
def d5(V):
    fig, ax = canvas("Clinical vs four-set model: the gain after optimism correction",
                     f"{V['cohort']['n']} patients, {V['cohort']['events']} recurrences; cause-specific recurrence hazard; internal validation only")
    # (a) design
    box(ax, 0.45, 3.85, 3.7, 2.35, fc=PANEL, ec=PANEL)
    ax.text(0.65, 5.95, "Models compared", fontsize=12.5, fontweight="bold", va="center")
    ax.text(0.65, 5.6, wrap("M1 clinical: age (spline), sex, stage II/III, tumour location, adjuvant chemotherapy (unpenalised)", 3.3, 9.6),
            fontsize=9.6, va="top", linespacing=1.25)
    ax.text(0.65, 4.85, wrap("M2 primary: M1 + EMT, E2F, Wnt, p53 scores (ridge-penalised)", 3.3, 9.6), fontsize=9.6, va="top", color=PRIMARY, linespacing=1.25)
    ax.text(0.65, 4.35, wrap("Five-set (+PI3K): sensitivity only, never promoted", 3.3, 9.3), fontsize=9.3, va="top", color=FIVE, linespacing=1.25)
    # (b) C-index
    a = fig.add_axes([0.37, 0.53, 0.27, 0.30])
    rows = [("M1 clinical", "M1", INK), ("M2 four-set", "M2_primary", PRIMARY), ("Five-set (sens.)", "M2_sensitivity", FIVE)]
    for i, (lab, k, c) in enumerate(rows):
        y = 2 - i
        ap, co = V["C"][k]["apparent"], V["C"][k]["corrected"]
        a.plot([co, ap], [y, y], color=GRID, lw=3)
        a.plot([ap], [y], "o", mfc="white", mec=APPARENT, ms=8, mew=1.6)
        a.plot([co], [y], "o", color=c, ms=8)
        a.text(ap + 0.0015, y + 0.18, f"{ap:.4f}", fontsize=8, color=MUTED)
        a.text(co - 0.0015, y + 0.18, f"{co:.4f}", fontsize=8, color=c, ha="right")
    a.set_yticks([2, 1, 0], [r[0] for r in rows], fontsize=9.5)
    a.set_xlim(0.592, 0.65)
    a.set_ylim(-0.6, 2.6)
    a.set_xlabel("Harrell's C (open = apparent, filled = optimism-corrected)", fontsize=8.8)
    for s in ("top", "right", "left"):
        a.spines[s].set_visible(False)
    a.tick_params(axis="y", length=0)
    a.grid(axis="x", color=GRID, lw=0.8)
    a.set_axisbelow(True)
    # (c) intervals
    b = fig.add_axes([0.70, 0.53, 0.27, 0.30])
    for r, (l, u) in enumerate(V["runs"]):
        b.plot([l, u], [4 - r, 4 - r], color=PRIMARY, lw=5, solid_capstyle="round")
    b.axvline(0, color=INK, lw=1, ls=(0, (3, 3)))
    b.axvline(V["theta"], color=INK, lw=1.8)
    b.text(V["theta"] + 0.0012, 4.55, r"$\hat{\theta}$" + f" = +{V['theta']:.6f}", fontsize=8.5)
    b.set_yticks(range(5), [f"Run {r}" for r in (4, 3, 2, 1, 0)], fontsize=9)
    b.set_xlim(-0.02, 0.06)
    b.set_ylim(-0.6, 5.0)
    b.set_xlabel("Corrected ΔC, 95 % percentile intervals", fontsize=8.8)
    for s in ("top", "right", "left"):
        b.spines[s].set_visible(False)
    b.tick_params(axis="y", length=0)
    b.grid(axis="x", color=GRID, lw=0.8)
    b.set_axisbelow(True)
    # numbers strip
    d_app = V["C"]["M2_primary"]["apparent"] - V["C"]["M1"]["apparent"]
    kx, ky = 0.45, 1.25
    box(ax, kx, ky, 7.6, 2.2, fc="white", ec=GRID)
    facts = [(f"+{d_app:.4f}", "apparent ΔC (fit and tested on the same data)", APPARENT),
             (f"+{V['theta']:.4f}", "optimism-corrected ΔC", PRIMARY),
             ("5 / 5", "corrected intervals include zero (≈ −0.009 to +0.05)", INK),
             (f"{V['slopes']['M1']:.3f} / {V['slopes']['M2_primary']:.3f}", "corrected calibration slopes, M1 / M2: residual overfitting", INK)]
    for i, (n_, t, c) in enumerate(facts):
        yy = ky + 1.85 - i * 0.48
        ax.text(kx + 0.2, yy, n_, fontsize=15, fontweight="bold", color=c, va="center")
        ax.text(kx + 2.05, yy, t, fontsize=10, va="center")
    # event bar
    ex0, ew = 8.4, W - 0.45 - 8.4
    box(ax, ex0, ky, ew, 2.2, fc=PANEL, ec=PANEL)
    ax.text(ex0 + 0.15, ky + 1.95, "Why validation stayed internal", fontsize=11, fontweight="bold", va="center")
    bars = [("GSE39582 (this cohort)", V["cohort"]["events"], PRIMARY), ("GSE14333 (exploratory)", 48, APPARENT),
            ("GSE33113 (exploratory)", 18, APPARENT), ("needed, indicative", 338, FAINT)]
    scale = (ew - 2.3) / 338
    for i, (lab, n_, c) in enumerate(bars):
        yy = ky + 1.55 - i * 0.37
        ax.text(ex0 + 0.15, yy, lab, fontsize=8.6, va="center")
        ax.add_patch(Rectangle((ex0 + 1.9, yy - 0.11), n_ * scale, 0.22, fc=c if i < 3 else "white", ec=c, lw=1, ls="-" if i < 3 else "--"))
        ax.text(ex0 + 1.95 + n_ * scale, yy, f"{'≈' if i == 3 else ''}{n_}", fontsize=8.6, va="center")
    footer(ax, "The five runs differ only in random seed: Monte Carlo replicates of one analysis, not independent cohorts or experiments "
               f"({V['B']['out']} outer × {V['B']['in']} inner resamples each; frozen verifier PASS {V['verification']['status_counts']['PASS']}/60). "
               "Intervals are conditional on 136 events; including zero does not establish no effect. The CIT partition is internal, not external validation. "
               "Required events: indicative simulation approximation at assumed C = 0.65 (Phase 5); event counts are not combined; GSE17536 excluded for patient overlap. "
               "Sources: figure_values.json (fig7, fig9); corrected_interval_results.json; docs/PHASE5_VALIDATION_PRECISION_ASSESSMENT.md §7–§8.")
    save(fig, "D5_model_comparison_uncertainty")


if __name__ == "__main__":
    V = values.load()
    for f in (d1, d2, d3, d4, d5):
        f(V)
    print("wrote", sorted(p.name for p in OUT.iterdir()))

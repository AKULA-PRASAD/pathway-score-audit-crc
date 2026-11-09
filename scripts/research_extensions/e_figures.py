"""Figures for research extensions E1–E3. Reads results/research_extensions/* only. Writes figures/research_extensions/FE*.{svg,png}."""
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import rcommon as C  # noqa: E402

INK, MUTED, GRID = "#1f2933", "#52606d", "#dde3ea"
plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 11, "svg.fonttype": "none", "svg.hashsalt": "gse39582-e",
                     "axes.edgecolor": "#8a94a3", "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": MUTED, "text.color": INK})
FOOT = "Post-hoc registered extension, outcome-free (docs/RESEARCH_EXTENSION_REGISTER.md)."


def save(fig, name):
    C.FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(C.FIG / f"{name}.svg", metadata={"Date": None})
    fig.savefig(C.FIG / f"{name}.png", dpi=200, metadata={"Software": None})
    plt.close(fig)


def fe1():
    a = pd.read_csv(C.RES / "e1_module_associations.csv")
    mods = pd.read_csv(C.RES / "e1_modules.csv")
    jac = pd.read_csv(C.RES / "e1_clusterwise_jaccard.csv").set_index("module")
    a = a[a.representation == "median"]
    cols = [("tp53.mutation", "TP53 mut"), ("braf.mutation", "BRAF mut"), ("cimp.status", "CIMP+"), ("cin.status", "CIN+"), ("mmr.status", "dMMR"),
            ("StromalScore", "Stromal est."), ("ImmuneScore", "Immune est."), ("TumorPurity", "Purity est.")]
    rows = list(mods.module)
    M = np.array([[a[(a.module == m) & (a.annotation == c)].effect.iloc[0] for c, _ in cols] for m in rows])
    fig = plt.figure(figsize=(12.2, 5.4))
    ax = fig.add_axes([0.27, 0.16, 0.55, 0.74])
    im = ax.imshow(M, cmap="RdBu_r", vmin=-0.9, vmax=0.9, aspect="auto")
    for i in range(len(rows)):
        for j in range(len(cols)):
            ax.text(j, i, f"{M[i, j]:+.2f}", ha="center", va="center", fontsize=9.5, color="white" if abs(M[i, j]) > 0.55 else INK)
    ax.set_xticks(range(len(cols)), [c[1] for c in cols], rotation=30, ha="right")
    lab = []
    for m in rows:
        r = mods[mods.module == m].iloc[0]
        lab.append(f"{m} ({r.n_genes} genes; Jaccard {jac.loc[m, 'mean_jaccard']:.2f})\n{', '.join(str(r.top_genes).split(', ')[:4])}")
    ax.set_yticks(range(len(rows)), lab, fontsize=9)
    ax.axvline(4.5, color="white", lw=3)
    ax.set_title("Gene modules within each Hallmark set: rank-biserial (genotype) | Spearman ρ (estimates)", loc="left", fontsize=11)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    fig.text(0.02, 0.02, FOOT + " GSE39582 n = 452. Registered ARI stability met only for Wnt; module-level claims are exploratory.", fontsize=9, color=MUTED)
    save(fig, "FE1_gene_modules")


def fe2():
    ws = pd.read_csv(C.RES / "e2_within_subtype_emt_stroma.csv")
    gt = pd.read_csv(C.RES / "e2_genotype_adjusted.csv")
    fig = plt.figure(figsize=(12.2, 5.0))
    a = fig.add_axes([0.12, 0.17, 0.31, 0.7])
    f = ws[ws.composition == "MCP_Fibroblasts"].reset_index(drop=True)
    y = np.arange(len(f))[::-1]
    a.errorbar(f.rho, y, xerr=[f.rho - f.lo, f.hi - f.rho], fmt="o", color="#c2366c", ms=7, capsize=3)
    a.axvline(0.7, color=MUTED, ls="--", lw=1)
    a.set_yticks(y, [f"{c} (n = {n})" for c, n in zip(f.cit, f.n)])
    a.set_xlim(0.5, 1.0)
    a.set_xlabel("Spearman ρ, EMT score vs MCP fibroblast estimate (95 % CI)")
    a.set_title("EMT tracks estimated fibroblasts in every CIT subtype", loc="left", fontsize=11)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    b = fig.add_axes([0.57, 0.24, 0.4, 0.63])
    sel = [("tp53.mutation", "p53", "TP53 → p53"), ("braf.mutation", "p53", "BRAF → p53"), ("braf.mutation", "Wnt", "BRAF → Wnt"),
           ("cimp.status", "p53", "CIMP → p53"), ("cimp.status", "Wnt", "CIMP → Wnt"), ("cin.status", "p53", "CIN → p53"), ("tp53.mutation", "E2F", "TP53 → E2F")]
    for i, (an, sc, lab) in enumerate(sel):
        r = gt[(gt.annotation == an) & (gt.score == sc)].iloc[0]
        yy = len(sel) - 1 - i
        b.plot([r.unadj_lo, r.unadj_hi], [yy + 0.15] * 2, color="#9aa5b1", lw=2)
        b.plot(r.diff_unadj, yy + 0.15, "o", color="#9aa5b1", label="unadjusted" if i == 0 else None)
        b.plot([r.adj_lo, r.adj_hi], [yy - 0.15] * 2, color="#2a78d6", lw=2)
        b.plot(r.diff_adj, yy - 0.15, "o", color="#2a78d6", label="adjusted (stromal, immune, ±dMMR)" if i == 0 else None)
    b.axvline(0, color=INK, lw=0.8)
    b.set_yticks(range(len(sel))[::-1], [s[2] for s in sel])
    b.set_xlabel("difference in score (SD units), bootstrap 95 % CI")
    b.set_title("Genotype effects before and after adjusting for estimates", loc="left", fontsize=11)
    b.legend(frameon=False, fontsize=9.5, loc="upper left", bbox_to_anchor=(0.0, -0.13), ncol=2)
    for s in ("top", "right"):
        b.spines[s].set_visible(False)
    fig.text(0.02, 0.01, FOOT + " n = 452. Estimates come from the same arrays; CIT subtypes derived from these data.", fontsize=9, color=MUTED)
    save(fig, "FE2_conditioned_structure")


def fe3():
    r = pd.read_csv(C.RES / "e3_replication_tests.csv")
    tc_path = C.RES / "e3_tcga_replication_tests.csv"
    if tc_path.exists():
        r = pd.concat([r, pd.read_csv(tc_path)], ignore_index=True)
    qc = pd.read_csv(C.RES / "e3_positive_control_qc.csv").set_index("cohort")
    qc.loc["TCGA-COAD", "passes_qc_median_r_ge_0.5"] = True
    s = json.loads((C.REP / "e3_summary.json").read_text())
    st = C.REP / "e3_tcga_summary.json"
    v_all = json.loads(st.read_text())["verdicts_qc_passing_plus_tcga"] if st.exists() else s["verdicts_excluding_qc_failed"]
    ref = {"R1": (-0.296, "TP53 mut → p53 score lower"), "R2": (0.188, "TP53 mut → E2F higher"), "R3": (-0.455, "BRAF mut → Wnt lower"),
           "R4": (0.482, "BRAF mut → p53 higher"), "R5": (None, "MSI-high → p53 higher"), "R6": (0.950, "EMT ~ fibroblast estimate (ρ)")}
    col = {"GSE41258": "#2a78d6", "GSE143985": "#1baf7a", "TCGA-COAD": "#eb6834", "GSE103479": "#9aa5b1"}
    fig = plt.figure(figsize=(12.2, 5.2))
    a = fig.add_axes([0.25, 0.16, 0.45, 0.76])
    hyps = ["R1", "R2", "R3", "R4", "R5", "R6"]
    for i, h in enumerate(hyps):
        y = len(hyps) - 1 - i
        if ref[h][0] is not None:
            a.plot(ref[h][0], y + 0.3, marker="D", color=INK, ms=7, label="GSE39582" if i == 0 else None)
        for k, (acc, c) in enumerate(col.items()):
            row = r[(r.hypothesis == h) & (r.cohort == acc) & (r.eligible == True)]  # noqa: E712
            if len(row):
                rr = row.iloc[0]
                yy = y + 0.12 - k * 0.16
                failed = not qc.loc[acc, "passes_qc_median_r_ge_0.5"]
                a.plot([rr.lo, rr.hi], [yy, yy], color=c, lw=2, ls="--" if failed else "-")
                a.plot(rr.effect, yy, "o", color=c, mfc="white" if failed else c, label=(acc + (" (failed QC)" if failed else "")) if h == "R1" else None)
    a.axvline(0, color=INK, lw=0.8)
    a.set_yticks(range(len(hyps))[::-1], [f"{h}: {ref[h][1]}" for h in hyps])
    a.set_xlim(-1, 1)
    a.set_xlabel("rank-biserial (R1–R5) or Spearman ρ (R6), bootstrap 95 % CI")
    a.legend(frameon=False, fontsize=9.5, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
    v = v_all
    fig.text(0.72, 0.5, "Pre-registered verdict\n(QC-passing cohorts incl. TCGA-COAD)\n\n" + "\n".join(f"{h}: {v[h]['verdict']}" for h in hyps), fontsize=10.5, va="top")
    fig.text(0.02, 0.01, FOOT + " Cohorts never pooled; TCGA-COAD all stages. GSE103479 failed positive-control QC (median r 0.07).", fontsize=9, color=MUTED)
    save(fig, "FE3_replication")


if __name__ == "__main__":
    fe1()
    fe2()
    fe3()
    print(sorted(p.name for p in C.FIG.glob("*.png")))

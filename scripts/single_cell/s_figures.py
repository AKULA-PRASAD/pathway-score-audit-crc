"""Figure FS1 for the single-cell extension S1-S5. Reads results/single_cell/* only; writes figures/single_cell/FS1_*.{png,svg}."""
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import sccommon as C  # noqa: E402

INK, MUTED, GRID = "#1f2933", "#52606d", "#dde3ea"
COL = {"epithelial": "#c2366c", "stromal": "#a86b12", "immune": "#5d6db8", "endothelial": "#00907e"}
COH = {"GSE132465": "#1f2933", "GSE144735": "#8a94a3", "GSE178341": "#2a78d6"}
NAME = {"GSE132465": "GSE132465 (discovery)", "GSE144735": "GSE144735", "GSE178341": "GSE178341"}
plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 10, "svg.fonttype": "none", "svg.hashsalt": "fs1",
                     "axes.edgecolor": "#8a94a3", "axes.spines.top": False, "axes.spines.right": False})


def main():
    enr = pd.read_csv(C.RES / "s1_s2_s5_enrichment.csv", dtype={"patient": str})
    coh = pd.read_csv(C.RES / "s3_coherence.csv", dtype={"patient": str})
    pb = pd.read_csv(C.RES / "s4_pseudobulk.csv")
    det = json.loads((C.RES / "s4_pseudobulk_detail.json").read_text())
    fig, ax = plt.subplots(1, 4, figsize=(16, 4.3), gridspec_kw={"width_ratios": [1.25, 1, 1, 1.05]})
    # A: S1
    a = ax[0]
    e = enr[enr["set"] == "EMT"]
    for i, acc in enumerate(C.DATASETS):
        for j, comp in enumerate(C.COMPARTMENTS):
            v = e[(e["accession"] == acc) & (e["compartment"] == comp)]["enrichment"].dropna()
            x = i * 5 + j
            a.scatter(np.full(len(v), x) + np.random.default_rng(1).uniform(-0.2, 0.2, len(v)), v, s=10, color=COL[comp], alpha=0.7, lw=0)
            a.plot([x - 0.35, x + 0.35], [v.median()] * 2, color=INK, lw=1.6)
    a.axhline(1, color=MUTED, lw=0.8, ls="--")
    a.set_yscale("log")
    a.set_xticks([1.5, 6.5, 11.5], ["GSE132465\n(discovery)", "GSE144735", "GSE178341"])
    a.set_ylabel("EMT-set enrichment vs matched random sets")
    a.set_title("A · S1: which cells supply EMT-set transcripts", loc="left", fontsize=10.5, fontweight="bold")
    for comp in C.COMPARTMENTS:
        a.scatter([], [], color=COL[comp], s=18, label=comp)
    a.set_ylim(top=60)
    a.legend(frameon=False, fontsize=8.5, loc="upper left", ncol=4, handletextpad=0.1, columnspacing=0.8)
    # B: S3
    b = ax[1]
    v = coh[coh["variant"] == "all genes"]
    for acc in C.DATASETS:
        s = v[v["accession"] == acc]
        b.scatter(s["random_p95"], s["coherence"], s=16, color=COH[acc], label=NAME[acc], lw=0)
    lim = [0, max(v["random_p95"].max(), v["coherence"].max()) * 1.1]
    b.plot(lim, lim, color=MUTED, lw=0.8, ls="--")
    b.set_xlim(lim)
    b.set_ylim(lim)
    b.set_xlabel("random-set 95th percentile (per patient)")
    b.set_ylabel("EMT-set coherence in malignant epithelium")
    b.set_title("B · S3: EMT-like programme in tumour cells?", loc="left", fontsize=10.5, fontweight="bold")
    b.text(lim[1] * 0.97, lim[1] * 0.06, "below the line = no more coherent\nthan expression-matched random genes", ha="right", fontsize=8, color=MUTED)
    b.legend(frameon=False, fontsize=8, loc="upper left")
    # C: S4
    c = ax[2]
    for acc in ("GSE132465", "GSE178341"):
        sc = pd.Series(det[acc]["score_EMT"])
        fr = pd.Series(det[acc]["stromal_fraction"]).reindex(sc.index)
        r = pb[(pb["accession"] == acc) & (pb["set"] == "EMT")].iloc[0]
        c.scatter(fr * 100, sc, s=18, color=COH[acc], lw=0, label=f"{NAME[acc]}: ρ {r['rho']:.2f} (random ≤ {r['random_p97_5']:.2f})")
    c.set_xlabel("counted stromal cells (% of tumour cells)")
    c.set_ylabel("pseudo-bulk EMT score (z)")
    c.set_title("C · S4: pseudo-bulk EMT score vs stroma", loc="left", fontsize=10.5, fontweight="bold")
    c.legend(frameon=False, fontsize=8, loc="lower right")
    # D: S5
    d = ax[3]
    m = enr[enr["accession"] == "GSE132465"].groupby(["set", "compartment"])["enrichment"].median().unstack()[C.COMPARTMENTS]
    m = m.loc[["EMT", "E2F", "Wnt", "p53"]]
    im = d.imshow(np.log2(m.to_numpy()), cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            lv = abs(np.log2(m.iloc[i, j]))
            d.text(j, i, f"{m.iloc[i, j]:.2f}", ha="center", va="center", fontsize=9, color="white" if lv > 1.5 else INK)
    d.set_xticks(range(4), C.COMPARTMENTS, rotation=20)
    d.set_yticks(range(4), m.index)
    d.set_title("D · S5: median enrichment (discovery)", loc="left", fontsize=10.5, fontweight="bold")
    cb = fig.colorbar(im, ax=d, fraction=0.046, pad=0.04)
    cb.set_label("log2 enrichment", fontsize=8.5)
    fig.text(0.01, -0.02, "Public scRNA-seq, authors' cell annotations; tumour samples of QC-passing patients (18, 4 and 39). Outcome-free; registered plan "
             "docs/SINGLE_CELL_EXTENSION_PLAN.md. mRNA only: not protein or pathway activity.", fontsize=8.5, color=MUTED)
    fig.tight_layout()
    C.FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(C.FIG / "FS1_single_cell_sources.png", dpi=200, bbox_inches="tight", metadata={"Software": None})
    fig.savefig(C.FIG / "FS1_single_cell_sources.svg", bbox_inches="tight", metadata={"Date": None})
    print("wrote", C.FIG / "FS1_single_cell_sources.png")


if __name__ == "__main__":
    main()

"""Figures for the registered exploratory analyses (X1–X8). Reads only reports/exploratory/* outputs and frozen inputs.

Writes reports/exploratory/figures/XF*.{svg,png}. Every panel states n and the analysis label.
"""
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import xcommon as X  # noqa: E402

OUTF = X.OUT / "figures"
INK, MUTED, GRID, PRIMARY = "#1f2933", "#52606d", "#dde3ea", "#2a78d6"
PW = {"EMT": "#c2366c", "E2F": "#6650c4", "Wnt": "#00907e", "p53": "#a86b12", "EMT_ex": "#e08aa8"}
plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 11, "svg.fonttype": "none", "svg.hashsalt": "gse39582-x",
                     "axes.edgecolor": "#8a94a3", "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": MUTED, "text.color": INK})


def save(fig, name):
    OUTF.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTF / f"{name}.svg", metadata={"Date": None})
    fig.savefig(OUTF / f"{name}.png", dpi=200, metadata={"Software": None})
    plt.close(fig)


def clean(a):
    for s in ("top", "right"):
        a.spines[s].set_visible(False)


def jitter_box(a, groups, labels, color):
    rng = np.random.default_rng(1)
    for i, g in enumerate(groups):
        a.boxplot([g], positions=[i], widths=0.5, showfliers=False, medianprops=dict(color=INK, lw=2),
                  boxprops=dict(color=MUTED), whiskerprops=dict(color=MUTED), capprops=dict(color=MUTED))
        a.scatter(i + rng.uniform(-0.18, 0.18, len(g)), g, s=9, color=color, alpha=0.45, lw=0)
    a.set_xticks(range(len(groups)), [f"{l}\nn = {len(g)}" for l, g in zip(labels, groups)])


def xf1():
    ctx = X.context()
    z, _ = X.cohort_scores(ctx)
    cov = ctx.cov.loc[ctx.ids]
    tp, mmr = cov["tp53.mutation"].astype(str), cov["mmr.status"].astype(str)
    res = pd.read_csv(X.OUT / "x1_molecular_anchoring" / "anchoring_tests.csv")
    sec = json.loads((X.OUT / "x1_molecular_anchoring" / "x1_summary.json").read_text())["secondary_tp53_within_pMMR"]
    fig = plt.figure(figsize=(12.2, 5.2))
    a = fig.add_axes([0.06, 0.16, 0.2, 0.68])
    jitter_box(a, [z.p53[tp == "WT"], z.p53[tp == "M"]], ["TP53 WT", "TP53 mutant"], PW["p53"])
    r = res[(res.annotation == "tp53.mutation") & (res.pathway == "p53")].iloc[0]
    a.set_title(f"All tumours: Δmedian {r.median_diff_sd:+.2f} SD\nq = {r.q_BH:.1e}", fontsize=11, loc="left")
    a.set_ylabel("p53 pathway score (SD)")
    clean(a)
    b = fig.add_axes([0.31, 0.16, 0.2, 0.68], sharey=a)
    sel = mmr == "pMMR"
    jitter_box(b, [z.p53[sel & (tp == "WT")], z.p53[sel & (tp == "M")]], ["TP53 WT", "TP53 mutant"], PW["p53"])
    b.set_title(f"pMMR only: Δmean {sec['mean_diff_sd_mutant_minus_wt']:+.2f} SD\nbatch-stratified perm. p = {sec['batch_stratified_permutation_p']:.4f}", fontsize=11, loc="left")
    clean(b)
    c = fig.add_axes([0.6, 0.16, 0.32, 0.68])
    ann = ["tp53.mutation", "kras.mutation", "braf.mutation", "cimp.status", "cin.status"]
    lab = ["TP53 mut", "KRAS mut", "BRAF mut", "CIMP+", "CIN+"]
    paths = ["EMT", "E2F", "Wnt", "p53"]
    M = np.array([[res[(res.annotation == an) & (res.pathway == p)].effect.iloc[0] for p in paths] for an in ann])
    Q = np.array([[res[(res.annotation == an) & (res.pathway == p)].q_BH.iloc[0] for p in paths] for an in ann])
    im = c.imshow(M, cmap="RdBu_r", vmin=-0.5, vmax=0.5, aspect="auto")
    for i in range(len(ann)):
        for j in range(4):
            c.text(j, i, f"{M[i, j]:+.2f}{'*' if Q[i, j] < 0.05 else ''}", ha="center", va="center", fontsize=10.5,
                   color="white" if abs(M[i, j]) > 0.3 else INK)
    c.set_xticks(range(4), paths)
    c.set_yticks(range(len(ann)), lab)
    c.set_title("Rank-biserial effect (feature vs absent); * BH q < 0.05", fontsize=11, loc="left")
    fig.colorbar(im, ax=c, fraction=0.04, pad=0.02)
    fig.text(0.06, 0.02, "Post-hoc, same data as primary. Exploratory, outcome-free (registered X1). n = 452; mutation status for 286–431. Association is not evidence of pathway function.",
             fontsize=9.5, color=MUTED)
    save(fig, "XF1_molecular_anchoring")


def xf2():
    cv = pd.read_csv(X.OUT / "x2_competing_risks" / "cumulative_incidence_curves.csv")
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.8), sharey=True)
    for a, grp in zip(axes, ["stage II", "stage III"]):
        for cause, col in [("recurrence", PRIMARY), ("death without recurrence", "#9aa5b1")]:
            s = cv[(cv.group == grp) & (cv.cause == cause)]
            a.step(s.time, s.cif, where="post", color=col, lw=2, label=cause)
            a.fill_between(s.time, s.lo, s.hi, step="post", color=col, alpha=0.15, lw=0)
        a.set_title(grp, loc="left")
        a.set_xlabel("time since surgery (unit provisional: months)")
        a.set_xlim(0, 120)
        a.grid(color=GRID)
        clean(a)
    axes[0].set_ylabel("cumulative incidence (Aalen–Johansen, 95 % CI)")
    axes[0].legend(frameon=False, loc="upper left")
    fig.text(0.06, 0.005, "Post-hoc, same data as primary. Exploratory competing-risk analysis (registered X2). n = 452: 136 recurrences, 46 deaths without recurrence. Absolute values assume the undocumented unit.",
             fontsize=9.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, "XF2_competing_risks")


def xf3():
    ctx = X.context()
    z, _ = X.cohort_scores(ctx)
    comp = pd.read_csv(X.OUT / "x5_composition" / "composition_estimates.csv", index_col=0).loc[ctx.ids]
    ex = pd.read_csv(X.OUT / "x6_scores_vs_composition" / "emt_ex_score.csv", index_col=0).loc[ctx.ids, "EMT_ex"]
    cor = pd.read_csv(X.OUT / "x6_scores_vs_composition" / "score_composition_spearman.csv")
    fig = plt.figure(figsize=(12.2, 5.0))
    a = fig.add_axes([0.06, 0.17, 0.25, 0.7])
    a.scatter(comp["MCP_Fibroblasts"], z.EMT, s=10, color=PW["EMT"], alpha=0.5, lw=0)
    r = cor[(cor.score == "EMT") & (cor.composition == "MCP_Fibroblasts")].iloc[0]
    a.set_title(f"EMT score vs fibroblast estimate\nρ = {r.rho:.2f} ({r.lo:.2f}–{r.hi:.2f})", loc="left", fontsize=11)
    a.set_xlabel("MCP-counter fibroblast score (log2, estimate)")
    a.set_ylabel("EMT score (SD)")
    clean(a)
    b = fig.add_axes([0.37, 0.17, 0.25, 0.7])
    b.scatter(comp["StromalScore"], ex, s=10, color=PW["EMT_ex"], alpha=0.6, lw=0)
    r2 = cor[(cor.score == "EMT_ex") & (cor.composition == "StromalScore")].iloc[0]
    b.set_title(f"EMT-ex (29 shared genes removed) vs\nESTIMATE stromal: ρ = {r2.rho:.2f} ({r2.lo:.2f}–{r2.hi:.2f})", loc="left", fontsize=11)
    b.set_xlabel("ESTIMATE StromalScore (estimate)")
    b.set_ylabel("EMT-ex score (SD)")
    clean(b)
    c = fig.add_axes([0.72, 0.17, 0.26, 0.7])
    rows = ["EMT", "E2F", "Wnt", "p53"]
    cols = ["StromalScore", "MCP_Fibroblasts", "MCP_Endothelial cells", "ImmuneScore", "MCP_Cytotoxic lymphocytes", "TumorPurity"]
    cl = ["Stromal", "Fibroblasts", "Endothelial", "Immune", "Cytotoxic", "Purity"]
    M = np.array([[cor[(cor.score == s) & (cor.composition == k)].rho.iloc[0] for k in cols] for s in rows])
    c.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    for i in range(4):
        for j in range(len(cols)):
            c.text(j, i, f"{M[i, j]:+.2f}", ha="center", va="center", fontsize=9.5, color="white" if abs(M[i, j]) > 0.6 else INK)
    c.set_xticks(range(len(cols)), cl, rotation=40, ha="right")
    c.set_yticks(range(4), rows)
    c.set_title("Spearman ρ, scores vs estimates", loc="left", fontsize=11)
    fig.text(0.06, 0.01, "Post-hoc, same data as primary. Exploratory (registered X5–X6). n = 452. ESTIMATE / MCP-counter: computational estimates from the same arrays, not measured fractions.",
             fontsize=9.5, color=MUTED)
    save(fig, "XF3_emt_vs_composition")


def xf4():
    s = json.loads((X.OUT / "x3_batch" / "x3_summary.json").read_text())
    lb = pd.read_csv(X.OUT / "x3_batch" / "lobo_diagnostics.csv")
    fig = plt.figure(figsize=(12.2, 4.8))
    a = fig.add_axes([0.06, 0.17, 0.4, 0.7])
    for lab in ["EMT", "E2F", "Wnt", "p53"]:
        d = lb[lb.pathway == lab]
        a.scatter(d.events_omitted, d.delta_in_se, s=36, color=PW[lab], label=lab, alpha=0.85)
    a.axhspan(-1.96, 1.96, color=GRID, alpha=0.6, lw=0)
    a.axhline(0, color=INK, lw=0.8)
    a.set_xlabel("recurrences in the omitted batch")
    a.set_ylabel("change in log-HR when batch omitted\n(multiples of full-cohort SE)")
    a.set_title("Leave-one-batch-out: 34 of 36 changes within ±1.96 SE (shaded)", loc="left", fontsize=11)
    a.legend(frameon=False, ncol=4, fontsize=10, loc="upper right")
    a.set_ylim(-3.0, 3.2)
    clean(a)
    b = fig.add_axes([0.56, 0.17, 0.4, 0.7])
    g = s["frailty"]["gamma"]
    gx, gy = np.array(g["profile_grid"]), np.array(g["profile_loglik"]) - g["marginal_loglik"]
    b.plot(gx[1:], gy[1:], marker="o", color=INK, lw=1.5)
    b.axhline(-1.92, color=MUTED, ls="--", lw=1)
    b.text(0.0055, -1.8, "95 % profile cut-off", fontsize=9.5, color=MUTED)
    b.text(0.0055, -4.9, f"θ → 0 (no frailty): {gy[0]:.1f} (not shown)", fontsize=9.5, color=MUTED)
    b.set_ylim(-5.2, 0.4)
    b.set_xscale("log")
    b.set_xlabel("gamma frailty variance θ (batch)")
    b.set_ylabel("profile log-likelihood (relative)")
    b.set_title(rf"Batch frailty: $\hat{{\theta}}$ = {g['theta_hat']:.2f}, profile 95 % ≈ {g['theta_95_profile_interval_on_grid'][0]:.2f}–{g['theta_95_profile_interval_on_grid'][1]:.2f}" + "\n" +
                f"pathway HRs unchanged (EMT {s['frailty']['gamma']['coefficients']['EMT']['HR']:.2f}, E2F {s['frailty']['gamma']['coefficients']['E2F']['HR']:.2f})",
                loc="left", fontsize=11)
    clean(b)
    fig.text(0.06, 0.005, "Post-hoc, same data as primary. Sensitivity analysis (registered X3). n = 452, 136 recurrences, 9 batches (14–103 patients, 4–52 recurrences). Frailty fits unpenalised; LOBO = frozen ridge fits.",
             fontsize=9.5, color=MUTED)
    save(fig, "XF4_batch")


def xf5():
    est = pd.read_csv(X.OUT / "x7_composition_adjusted_cox" / "adjusted_estimates.csv")
    sm = json.loads((X.OUT / "x7_composition_adjusted_cox" / "x7_summary.json").read_text())
    fig = plt.figure(figsize=(12.2, 4.8))
    a = fig.add_axes([0.1, 0.17, 0.55, 0.72])
    models = [("frozen M2 (ridge, primary)", None, "#1f2933"), ("R0 unpenalised", "R0", PRIMARY), ("R1 + ESTIMATE stromal, immune", "R1", "#eb6834"),
              ("R2 + MCP fibroblasts, endothelial, cytotoxic", "R2", "#1baf7a")]
    for i, lab in enumerate(["EMT", "E2F", "Wnt", "p53"]):
        for k, (name, m, col) in enumerate(models):
            y = 3 - i + 0.27 - k * 0.18
            if m is None:
                f = sm["frozen_ridge_M2_HR"][lab]
                h, lo, hi = f["HR"], f["lo"], f["hi"]
            else:
                r = est[(est.model == m) & (est.term == lab)].iloc[0]
                h, lo, hi = r.HR, r.lo, r.hi
            a.plot([lo, hi], [y, y], color=col, lw=2)
            a.plot([h], [y], "o", color=col, ms=6, label=name if i == 0 else None)
    a.axvline(1, color=INK, ls=(0, (3, 3)), lw=1)
    a.set_xscale("log")
    a.set_xticks([0.6, 0.8, 1, 1.25, 1.6, 2.5, 3.4], ["0.6", "0.8", "1", "1.25", "1.6", "2.5", "3.4"])
    a.minorticks_off()
    a.set_yticks([3, 2, 1, 0], ["EMT", "E2F", "Wnt", "p53"])
    a.set_xlabel("hazard ratio per SD (log scale), Wald 95 % CI")
    a.legend(frameon=False, fontsize=9.5, loc="upper right", bbox_to_anchor=(1.62, 1.0))
    clean(a)
    d = sm["discrimination"]
    txt = "Exploratory corrected C\n(single-level bootstrap, B = 200)\n" + "\n".join(
        f"{k}: {v['corrected_C']:.3f}" for k, v in d.items()) + f"\n\nVIF in R1: EMT {sm['vif']['R1']['EMT']:.1f}, stromal {sm['vif']['R1']['StromalScore']:.1f}"
    fig.text(0.73, 0.2, txt, fontsize=10.5, va="bottom")
    fig.text(0.06, 0.01, "Post-hoc, same data as primary. Exploratory (registered X7). n = 452, 136 recurrences. Not comparable with the frozen ΔC estimand. EMT and stromal terms collinear.",
             fontsize=9.5, color=MUTED)
    save(fig, "XF5_adjusted_cox")


def xf6():
    t = pd.read_csv(X.OUT / "x8_external_audit" / "eligibility_table.csv")
    t = t[t.events_II_III.notna() & t.stage_II_III_n.notna()].sort_values("events_II_III")
    fig = plt.figure(figsize=(12.2, 4.8))
    a = fig.add_axes([0.13, 0.17, 0.55, 0.72])
    y = np.arange(len(t))
    cols = ["#2a78d6" if c3 else "#9aa5b1" for c3 in t.c3_all_five_M1_predictors]
    a.barh(y, t.events_II_III, color=cols, height=0.6)
    a.axvline(136, color=INK, lw=1)
    a.text(138, len(t) - 0.6, "GSE39582 (136)", fontsize=9.5)
    a.axvline(337, color="#b42318", lw=1, ls="--")
    a.text(339, 0, "≈ 337 events for slope\nCI width 0.4 (Riley relaxed)", fontsize=9.5, color="#b42318")
    a.set_yticks(y, t.cohort)
    a.set_xlabel("stage II–III recurrence events with outcome (deposited)")
    a.set_xlim(0, 420)
    clean(a)
    fig.text(0.72, 0.62, "Blue: all five M1 predictors deposited\nGrey: one or more missing\n\nNo cohort meets all four compatibility\ncriteria (independence, compatible\nrecurrence endpoint, all predictors,\nexact score reproduction), so no\nfrozen-model validation was run.",
             fontsize=10.5, va="top")
    fig.text(0.06, 0.01, "External-cohort eligibility audit (registered X8): GEO search 2026-11-08, 554 hits, 85 screened, 14 assessed. Counts from deposited characteristics only.",
             fontsize=9.5, color=MUTED)
    save(fig, "XF6_external_audit")


if __name__ == "__main__":
    for f in (xf1, xf2, xf3, xf4, xf5, xf6):
        f()
    print("wrote", sorted(p.name for p in OUTF.glob("*.png")))

"""Stage G -- figures generated ONLY from results/s1_interval_results.json and the (read-only)
checkpoints. Writes figures/*.png, figures/*.svg and figures/FIGURE_CAPTIONS.md.

Palette (dataviz reference, validated light-mode all-pairs for slots 1-3):
  series-1 #2a78d6 (percentile), series-2 #eb6834 (basic), series-3 #1baf7a (BC z0; <3:1 contrast,
  so every use is direct-labelled and T1 carries the values); status-critical #d03b3b (EP-1 FAIL)
Usage: .venv/bin/python -B -u reports/phase9/final_project/analysis/s2_figures.py
"""
from __future__ import annotations

import json, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import FP, N_RUNS, load_run  # noqa: E402

S1, S2, S3, CRIT = "#2a78d6", "#eb6834", "#1baf7a", "#d03b3b"
INK, INK2, GRID, SURF, NEUTRAL = "#0b0b0b", "#52514e", "#e6e5e0", "#fcfcfb", "#8a8984"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
                     "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False, "font.size": 9,
                     "axes.titlesize": 10, "axes.titleweight": "bold", "lines.linewidth": 2})
PROV = "PROVISIONAL — as-run implementation; CV-seed protocol deviation not yet adjudicated"
SRC = "Source: reports/phase9/final_project/results/s1_interval_results.json"
RES = json.loads((FP / "results/s1_interval_results.json").read_text())
TH = RES["theta_hat"]["value"]
OUT = FP / "figures"
captions = []


def save(fig, name, title, caption):
    fig.text(0.01, 0.005, PROV, fontsize=7, color=INK2, ha="left", va="bottom")
    for ext in ("png", "svg"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=200 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    captions.append(f"### {name} — {title}\n\n{caption}\n\n*{SRC}. {PROV}.*\n")


# ---------------- F1 distribution of dC_corr(k), small multiples
fig, axes = plt.subplots(N_RUNS, 1, figsize=(7.2, 7.6), sharex=True)
bins = np.linspace(-0.02, 0.10, 61)
for r, ax in enumerate(axes):
    v = np.array([x["dC_corrected"] for x in load_run(r)]); R = RES["runs"][r]
    ax.hist(v, bins=bins, color=NEUTRAL, edgecolor=SURF, linewidth=0.6)
    p = R["percentile_primary"]
    for q in (p["lower"], p["upper"]):
        ax.axvline(q, color=S1, lw=2)
    ax.axvline(TH, color=INK, lw=1.5); ax.axvline(0, color=INK2, lw=1, ls=(0, (4, 3)))
    ax.set_ylabel(f"run {r}\ncount", fontsize=8)
    ax.text(0.99, 0.85, f"median {R['distribution']['dC_corrected']['median']:+.4f} · "
            f"{100*R['distribution']['dC_corrected']['frac_negative']:.1f}% < 0 · z₀ {R['z0_diagnostic']['z0']:+.2f}",
            transform=ax.transAxes, ha="right", fontsize=7.5, color=INK2)
axes[0].text(TH + 0.0008, axes[0].get_ylim()[1] * 1.02, "θ̂ = +0.0043", fontsize=7.5, va="bottom")
axes[0].text(-0.0008, axes[0].get_ylim()[1] * 1.02, "0", fontsize=7.5, va="bottom", ha="right", color=INK2)
axes[0].text(RES["runs"][0]["percentile_primary"]["upper"], axes[0].get_ylim()[1] * 0.6,
             " 95% percentile\n endpoint", fontsize=7.5, color=INK, va="top")
axes[-1].set_xlabel("ΔC_corr(k) = C_corr(M2 four-set) − C_corr(M1)   (difference in Harrell's C, unitless)")
fig.suptitle("F1. Outer-replicate distribution of the optimism-corrected ΔC, per run (B_out = 500 each)",
             x=0.01, ha="left", fontsize=10, fontweight="bold")
fig.tight_layout(rect=(0, 0.02, 1, 0.97))
save(fig, "F1_outer_distributions", "Outer-replicate distributions",
     "Histograms (bin width 0.002) of the 500 outer-replicate values ΔC_corr(k) in each run. Blue lines: 95 % "
     "percentile interval endpoints (type-7 quantiles at 0.025/0.975). Black line: θ̂ = +0.0043 "
     "(Phase 9 primary). Dashed: 0. Annotations: run median, percentage of negative replicates, z₀. "
     "Conditional on 136 events / 316 censored in the 452-patient A-19 population.")

# ---------------- F2 interval comparison
fig, ax = plt.subplots(figsize=(7.2, 4.0))
meth = (("percentile_primary", "Percentile (primary)", S1, 0.22),
        ("basic_diagnostic", "Basic (diagnostic)", S2, 0.0),
        ("bc_z0_diagnostic", "BC z₀-only (diagnostic)", S3, -0.22))
for r in range(N_RUNS):
    R = RES["runs"][r]
    for key, lab, col, off in meth:
        lo, hi = R[key]["lower"], R[key]["upper"]
        ax.plot([lo, hi], [r + off] * 2, color=col, lw=2.2, solid_capstyle="round", label=lab if r == 0 else None)
        ax.plot([lo, hi], [r + off] * 2, "|", color=col, ms=8, mew=2)
ax.axvline(TH, color=INK, lw=1.5); ax.axvline(0, color=INK2, lw=1, ls=(0, (4, 3)))
ax.text(TH + 0.0006, -0.42, "θ̂ = +0.0043", ha="left", fontsize=7.5); ax.text(-0.0006, -0.42, "0", ha="right", fontsize=7.5, color=INK2)
R0 = RES["runs"][0]
for key, lab, col, off in meth:
    ax.text(R0[key]["upper"] + 0.0012, 0 + off, lab.split(" (")[0], fontsize=7.5, va="center", color=INK)
ax.set_yticks(range(N_RUNS)); ax.set_yticklabels([f"run {r}" for r in range(N_RUNS)]); ax.set_ylim(-0.6, 4.5)
ax.invert_yaxis(); ax.set_xlim(right=0.066)
ax.set_xlabel("ΔC_corr (difference in Harrell's C, unitless)")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False, fontsize=8)
fig.suptitle("F2. 95 % intervals for ΔC_corr by method and run", x=0.01, ha="left", fontsize=10, fontweight="bold")
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
save(fig, "F2_interval_comparison", "Interval methods compared",
     "Per-run 95 % intervals: percentile (primary, plan §7), basic = [2θ̂ − q.975, 2θ̂ − q.025] and "
     "bias-corrected z₀-only (acceleration 0; diagnostic). The basic interval is the reflection of the "
     "percentile interval about θ̂ = +0.0043; their divergence (≈ 52–56 % of the percentile width at each "
     "endpoint) breaches the W6/A7 warning in every run. Values: tables/T1_intervals_per_run.csv.")

# ---------------- F3 endpoint variability + MC-SE
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))
for ax, e in zip(axes, ("lower", "upper")):
    vals = RES["between_run"][f"percentile_{e}"]
    m, sd = vals["mean"], vals["sd_ddof1"]
    ax.axhspan(m - sd, m + sd, color=GRID, zorder=0)
    ax.axhline(m, color=INK2, lw=1)
    for r in range(N_RUNS):
        k = RES["runs"][r]["kernel_MCSE_2_2"]["percentile_endpoints"][e]
        q = k["q"]; ses = [k[f"bw_x{b}"]["SE"] for b in (0.5, 1.0, 2.0)]
        ax.plot([r, r], [q - max(ses), q + max(ses)], color=S1, lw=1, alpha=0.6)
        ax.plot([r, r], [q - ses[1], q + ses[1]], color=S1, lw=2.5, solid_capstyle="round")
        ax.plot(r, q, "o", color=S1, ms=7, mec=SURF, mew=2)
    ax.set_xticks(range(N_RUNS)); ax.set_xticklabels([f"run {r}" for r in range(N_RUNS)])
    ax.set_title(f"{'Lower (2.5 %)' if e == 'lower' else 'Upper (97.5 %)'} endpoint\n"
                 f"band: mean ± between-run SD {m:+.5f} ± {sd:.5f}", loc="left", fontsize=7.5, fontweight="normal")
    ax.set_ylabel("ΔC_corr endpoint (unitless)")
fig.suptitle("F3. Endpoint Monte-Carlo variability: within-run kernel MC-SE vs between-run SD",
             x=0.01, ha="left", fontsize=10, fontweight="bold")
fig.tight_layout(rect=(0, 0.03, 1, 0.93))
save(fig, "F3_endpoint_variability", "Endpoint Monte-Carlo variability",
     "Dots: per-run percentile endpoint. Thick bar: ±1 kernel MC-SE at the Silverman bandwidth (plan §2.2); "
     "thin bar: ±1 the largest kernel MC-SE over bandwidth × {0.5, 1, 2}. Grey band: mean ± between-run SD "
     "(ddof = 1, df = 4, relative SE ≈ 35 %). These are Monte-Carlo quantities; they are not sampling "
     "uncertainty and are never summed with it (MC-4). Values: tables/T2_kernel_MCSE.csv, T3_between_run.csv.")

# ---------------- F4 MC-1 running trace
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharex=True)
for ax, (i, e) in zip(axes, enumerate(("lower", "upper"))):
    ax.axvspan(300, 500, color=GRID, zorder=0); ends = []
    for r in range(N_RUNS):
        tr = RES["runs"][r]["MC1_running_trace"]["fine_trace_step10"]
        B = [int(b) for b in tr]; y = [tr[b][i] for b in tr]
        breach = RES["runs"][r]["warnings"]["W5_breach"] and e == "upper"
        ax.plot(B, y, color=S2 if breach else NEUTRAL, lw=2 if breach else 1.2, zorder=3 if breach else 2)
        ends.append((y[-1], r, breach))
    ends.sort(); lo_, hi_ = ax.get_ylim(); gap = (hi_ - lo_) * 0.055; placed = []
    for yv, r, br in ends:
        yy = max(yv, placed[-1] + gap) if placed else yv; placed.append(yy)
        ax.text(505, yy, f"run {r}", fontsize=6.5, va="center", color=INK if br else INK2)
    ax.set_title(f"{'Lower' if e == 'lower' else 'Upper'} endpoint, first B replicates", loc="left")
    ax.set_xlabel("B (outer replicates, k order)"); ax.set_ylabel("percentile endpoint (unitless)")
    ax.set_xlim(50, 545)
axes[1].text(305, axes[1].get_ylim()[0], " final 40 % (B ≥ 300)", fontsize=7, va="bottom", color=INK2)
w1 = RES["runs"][1]["warnings"]
axes[1].set_title(f"Upper endpoint, first B replicates\nrun 1 final-40 % spread {w1['W5_spread']['upper']:.5f} "
                  f"> W5 bound {w1['W5_bound']:.5f}", loc="left", fontsize=8.5)
fig.suptitle("F4. MC-1 running trace of the percentile endpoints (W5)", x=0.01, ha="left", fontsize=10, fontweight="bold")
fig.tight_layout(rect=(0, 0.03, 1, 0.93))
save(fig, "F4_MC1_running_trace", "MC-1 running trace",
     "Percentile endpoints computed from the first B outer replicates (k order), B = 50, 60, …, 500. W5 is "
     "evaluated on the Phase 6B MC-1 grid B ∈ {300, 400, 500} (final 40 %): spread > 0.05 × width is a "
     "warning. Highlighted: run 1 upper endpoint, the only W5 breach. Values: tables/T5_warnings_hardfailures.csv.")

# ---------------- F5 EP-1 ratio grid
fig, axes = plt.subplots(1, 4, figsize=(7.2, 3.2), sharey=True)
ests = ("kernel_bw_x0.5", "kernel_bw_x1.0", "kernel_bw_x2.0", "between_run_SD")
labs = ("Kernel, h × 0.5", "Kernel, h × 1 (Silverman)", "Kernel, h × 2", "Between-run SD")
rows = [(r, e) for r in range(N_RUNS) for e in ("lower", "upper")]
for ax, est, lab in zip(axes, ests, labs):
    ax.axvline(0.05, color=INK, lw=1.2)
    for j, (r, e) in enumerate(rows):
        g = next(g for g in RES["EP1"]["grid"] if g["run"] == r and g["endpoint"] == e and g["estimator"] == est)
        col = S1 if g["pass"] else CRIT
        ax.plot(g["ratio"], j, "o" if g["pass"] else "X", color=col, ms=7 if g["pass"] else 9, mec=SURF, mew=1.5)
        if not g["pass"]:
            ax.text(g["ratio"], j - 0.45, f"FAIL {g['ratio']:.4f}", fontsize=7, ha="center", color=INK)
    ax.set_title(lab, fontsize=8, loc="left"); ax.set_xlim(0, 0.065)
    ax.set_xlabel("MC-SE / width")
axes[0].set_yticks(range(len(rows))); axes[0].set_yticklabels([f"run {r} {e}" for r, e in rows], fontsize=7.5)
axes[0].invert_yaxis()
fig.suptitle("F5. EP-1: endpoint MC-SE as a fraction of percentile-interval width (criterion ≤ 0.05)",
             x=0.01, ha="left", fontsize=10, fontweight="bold")
fig.tight_layout(rect=(0, 0.03, 1, 0.92))
save(fig, "F5_EP1_grid", "EP-1 evaluation grid",
     "Each point is one EP-1 evaluation: MC-SE(endpoint) / (upper − lower) for one run, endpoint and "
     "MC-SE estimator. 39/40 evaluations are ≤ 0.05; the exception is run 1 upper at bandwidth × 0.5 "
     "(ratio 0.0535). All 10 pass at the Silverman bandwidth and all 10 pass with the between-run SD. The "
     "plan does not state which estimator governs EP-1; the outcome is reported as UNRESOLVED. "
     "Values: tables/T4_EP1_grid.csv.")

# ---------------- F6 observed decomposition
allr = [x for r in range(N_RUNS) for x in load_run(r)]
bd = RES["bias_decomposition_observed"]
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9))
for ax, key, dval, lab in ((axes[0], "dC_app", bd["dC_app_D"], "ΔC_app(k)"),
                           (axes[1], "dOpt", bd["dOpt_D"], "ΔOpt(k)"),
                           (axes[2], "dC_corrected", bd["theta_hat"], "ΔC_corr(k)")):
    v = np.array([x[key] for x in allr])
    ax.hist(v, bins=50, color=NEUTRAL, edgecolor=SURF, linewidth=0.5)
    ax.axvline(dval, color=INK, lw=1.5); ax.axvline(v.mean(), color=S1, lw=2)
    ax.set_title(f"{lab}", loc="left", fontsize=9)
    ax.set_title(f"on D {dval:+.4f}\nmean(k) {v.mean():+.4f}", loc="right", fontsize=7, color=INK2, fontweight="normal")
    ax.set_xlabel("unitless"); ax.tick_params(labelsize=7)
axes[0].set_ylabel("count (2500 outer replicates)")
from matplotlib.lines import Line2D
fig.legend([Line2D([], [], color=INK, lw=1.5), Line2D([], [], color=S1, lw=2)],
           ["original-data value (D)", "outer-replicate mean"], loc="upper left", bbox_to_anchor=(0.01, 0.92),
           ncol=2, frameon=False, fontsize=8)
fig.suptitle("F6. Outer replicates vs original-data values: apparent ΔC, ΔOpt and ΔC_corr",
             x=0.01, ha="left", fontsize=10, fontweight="bold")
fig.tight_layout(rect=(0, 0.03, 1, 0.84))
save(fig, "F6_shift_decomposition", "Observed shift decomposition",
     "Pooled histograms over 2 500 outer replicates. Black line: original-data value (Phase 9 primary: "
     "ΔC_app(D), ΔOpt(D) = mean over 5 seeds, θ̂). Blue line: outer-replicate mean. ΔC_corr(k) = ΔC_app(k) − "
     "ΔOpt(k); the outer mean exceeds θ̂ by +0.0097, of which +0.0083 comes from a higher apparent ΔC on D*ₖ "
     "and +0.0014 from a lower ΔOpt. Descriptive only: the original-data values were produced by Phase 9 "
     "primary, under a different seeding convention (deviation record §8).")

(OUT / "FIGURE_CAPTIONS.md").write_text("# Figure captions (generated by analysis/s2_figures.py)\n\n" + "\n".join(captions))
print("figures:", sorted(p.name for p in OUT.glob("*.png")))

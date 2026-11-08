#!/usr/bin/env python3
"""Scientific results figures for the GSE39582 recurrence study (presentation only).

Draws the approved figure plan (docs/SCIENTIFIC_RESULTS_FIGURE_AUDIT.md section 4) from
EXISTING frozen outputs only. No model is fitted, no statistic is re-estimated, and no
source file is modified. The age partial-effect curve (Fig. 5b) is a deterministic display of
the stored M1 coefficients and knots through the frozen `rcs_basis`; the CIT counts (Fig. 1)
are descriptive counts from the frozen analysis bundle.

Usage (repository root):
    .venv/bin/python scripts/49_make_results_figures.py [output_dir]
Default output_dir: reports/figures/final. The directory must be absent or empty; every file
is created exclusively, so no existing figure is ever overwritten.

Outputs: Fig1..Fig9 and FigS1..FigS2 as .svg and .png, FIGURE_CAPTIONS.md,
figure_values.json (every plotted number) and figure_sources.json (input hashes, script hash,
software versions, generation parameters).
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import platform
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib  # noqa: E402
import matplotlib.ticker  # noqa: E402,F401
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from src.phase9_execution import rcs_basis  # noqa: E402  (frozen; used only for display)

# ---------------------------------------------------------------- sources
SRC = {
    "corrected_results": "reports/phase9/corrected_rerun_report/corrected_interval_results.json",
    "verification": "reports/phase9/corrected_rerun/verification/L00750.report.json",
    "execution": "reports/phase9/execution/phase9_execution_results.json",
    "execution_s7": "reports/phase9/execution_s7/phase9_s7_results.json",
    "cohort_counts": "reports/phase1/cohort_counts.csv",
    "mapping_audit": "reports/phase2/phase2_mapping_audit.json",
    "sample_qc": "reports/phase7c/sample_qc_metrics.csv",
    "batch_composition": "reports/phase7c/batch_composition.csv",
    "sample_distribution": "reports/phase8a/sample_distribution_summary.csv",
    "pca_scores": "reports/phase8a/pca_scores_gene_level.csv",
    "phase8a_summary": "reports/phase8a/phase8a_summary.json",
    "score_distributions": "reports/phase8c/score_distributions.csv",
    "score_corr_5set": "reports/phase8c/score_correlation_pearson_5set_sensitivity.csv",
    "representation": "reports/phase8c/representation_sensitivity.csv",
    "qc_report_doc": "docs/PHASE7C_QC_REPORT.md",
    "endpoint_doc": "docs/PHASE6A_ENDPOINT_RESOLUTION.md",
    "a19_freeze": "config/analysis/A19_FREEZE_RECORD.json",
    "bundle_membership": "data/processed/analysis/cohort_membership.csv",
    "bundle_covariates": "data/processed/analysis/covariates.csv",
    "bundle_endpoint": "data/processed/analysis/endpoint.csv",
    "bundle_partition": "data/processed/analysis/partition_and_batch.csv",
}
CHECKPOINTS = [f"reports/phase9/corrected_rerun/runs/run{r}.jsonl" for r in range(5)]
RUNS = (0, 1, 2, 3, 4)
THETA_EXPECTED = 0.004316342804030704

# Values transcribed from frozen documents (each string is asserted to occur in its source
# document before use, so a transcription error stops the script).
DOC_FACTS = {
    "qc_report_doc": {
        "arrays_total": (585, "the matrix holds 585"),
        "nontumoural": (19, "19 non-tumoural"),
        "stage_other": (97, "97 stage 0/I/IV or missing-stage"),
        "no_rfs": (8, "**8** stage II–III lacking documented RFS"),
        "s_s1_fence": ((4.6907, 4.8532), "[4.6907, 4.8532]"),
    },
    "endpoint_doc": {
        "died_no_recurrence": ((46, 75.0, 32.6, 23.9),
                               "| **died without recurrence** (censored) | **46** | **75.0** | **32.6%** | **23.9%** |"),
        "alive_no_recurrence": ((271, 68.0, 38.7, 42.1),
                                "| alive, no recurrence (true censoring) | 271 | 68.0 | 38.7% | 42.1% |"),
        "recurrence": ((136, 68.5, 58.8, 55.1), "| recurrence | 136 | 68.5 | 58.8% | 55.1% |"),
    },
}

PATHWAYS = ["HALLMARK_E2F_TARGETS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
            "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_P53_PATHWAY"]
PI3K = "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
SHORT = {"HALLMARK_E2F_TARGETS": "E2F targets",
         "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "EMT",
         "HALLMARK_WNT_BETA_CATENIN_SIGNALING": "Wnt/β-catenin",
         "HALLMARK_P53_PATHWAY": "p53 pathway",
         PI3K: "PI3K/AKT/mTOR"}
M1_TERMS = [("sex=Male", "Male vs female"), ("tnm.stage=3.0", "Stage III vs II"),
            ("tumor.location=proximal", "Proximal vs distal"),
            ("chemotherapy.adjuvant=Y", "Adjuvant chemotherapy (yes vs no)")]

# palette (validated: dataviz validator, light surface, all-pairs) and ink
C_PRIMARY, C_S7, C_FIVE = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, SURFACE, NEUTRAL = "#0b0b0b", "#52514e", "#e4e3df", "#ffffff", "#9a9893"
DPI = 200

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE, "svg.hashsalt": "coad-phase9-figures", "svg.fonttype": "path",
    "lines.linewidth": 2.0, "legend.frameon": False, "legend.fontsize": 8,
})


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(rel):
    with open(ROOT / rel, newline="") as f:
        return list(csv.DictReader(f))


def load_json(rel):
    return json.loads((ROOT / rel).read_bytes())


def doc_fact(doc_key, fact):
    value, needle = DOC_FACTS[doc_key][fact]
    text = (ROOT / SRC[doc_key]).read_text()
    if needle not in text:
        raise SystemExit(f"transcribed fact not found in {SRC[doc_key]}: {needle!r}")
    return value


def panel_label(ax, s):
    ax.text(-0.02, 1.04, s, transform=ax.transAxes, fontsize=11, fontweight="bold",
            ha="right", va="bottom")


# ---------------------------------------------------------------- data assembly
def assemble():
    V = {}
    cir = load_json(SRC["corrected_results"])
    ex = load_json(SRC["execution"])
    s7 = load_json(SRC["execution_s7"])
    ver = load_json(SRC["verification"])
    theta = cir["theta_hat"]
    if theta != THETA_EXPECTED:
        raise SystemExit(f"theta-hat {theta!r} != expected {THETA_EXPECTED!r}")
    if not (ver["verdict"] == "PASS" and ver["reporting_permitted"] and ver["level"] == 750):
        raise SystemExit("verification report does not permit reporting at B_out = 750")

    # ---- Fig 1: cohort flow and CIT descriptive count (frozen bundle)
    cc = {r["scenario"]: r for r in read_csv(SRC["cohort_counts"])}
    mem = read_csv(SRC["bundle_membership"])
    s3 = [r["sample_id"] for r in mem if r["in_S3"] == "True"]
    cov = {r["sample_id"]: r for r in read_csv(SRC["bundle_covariates"])}
    m1_vars = ["age.at.diagnosis_(year)", "sex", "tnm.stage", "tumor.location", "chemotherapy.adjuvant"]
    a19 = [s for s in s3 if all(cov[s][v] != "" for v in m1_vars)]
    excluded = [s for s in s3 if s not in set(a19)]
    ep = {r["sample_id"]: float(r["rfs_event_full"]) for r in read_csv(SRC["bundle_endpoint"])}
    part = {r["sample_id"]: r["dataset"] for r in read_csv(SRC["bundle_partition"])}
    freeze_excl = [o["sample_id"] for o in load_json(SRC["a19_freeze"])["excluded_observations"]]
    if excluded != freeze_excl:
        raise SystemExit(f"A-19 exclusion {excluded} != frozen record {freeze_excl}")

    def split(ids):
        out = {}
        for p in ("discovery", "validation"):
            sel = [s for s in ids if part[s] == p]
            ev = int(sum(ep[s] for s in sel))
            out[p] = {"n": len(sel), "events": ev, "censored": len(sel) - ev}
        return out
    V["fig1"] = {
        "arrays_total": doc_fact("qc_report_doc", "arrays_total"),
        "excluded_nontumoural": doc_fact("qc_report_doc", "nontumoural"),
        "excluded_stage_other": doc_fact("qc_report_doc", "stage_other"),
        "excluded_no_rfs": doc_fact("qc_report_doc", "no_rfs"),
        "S0": {"n": int(cc["S0"]["n"]), "events": int(cc["S0"]["events_full"])},
        "S3": {"n": int(cc["S3"]["n"]), "events": int(cc["S3"]["events_full"]),
               "censored": int(cc["S3"]["censored_full"])},
        "A19": {"n": len(a19), "events": int(sum(ep[s] for s in a19)),
                "censored": len(a19) - int(sum(ep[s] for s in a19))},
        "A19_excluded": [{"sample_id": s, "partition": part[s], "event": int(ep[s])} for s in excluded],
        "S3_by_partition": split(s3), "A19_by_partition": split(a19),
        "deaths_censored": doc_fact("endpoint_doc", "died_no_recurrence")[0],
        "execution_population": ex["population"],
    }
    if (V["fig1"]["A19"]["n"], V["fig1"]["A19"]["events"]) != (ex["population"]["n"], ex["population"]["events"]):
        raise SystemExit("A-19 counts do not match the execution record")

    # ---- Fig 2: preprocessing and QC
    mp = load_json(SRC["mapping_audit"])["mapping"]
    qc = read_csv(SRC["sample_qc"])
    sd = read_csv(SRC["sample_distribution"])
    batches = sorted({r["combat_batch"] for r in sd})
    comp = read_csv(SRC["batch_composition"])
    V["fig2"] = {
        "probes": {k: mp[k] for k in ("n_probes_total", "n_unique", "n_unmapped", "n_ambiguous",
                                      "n_genes_unique", "n_genes_single_probe", "n_genes_multi_probe",
                                      "max_probes_per_gene")},
        "qc_median_by_batch": {b: [float(r["median"]) for r in qc if r["combat_batch"] == b] for b in batches},
        "s_s1_fence": list(doc_fact("qc_report_doc", "s_s1_fence")),
        "qc_flags": {k: sum(1 for r in qc if r[k] == "True") for k in ("flag_S_S1", "flag_S_S2", "flag_S_S3")},
        "qc_flag_union": sum(1 for r in qc if "True" in (r["flag_S_S1"], r["flag_S_S2"], r["flag_S_S3"])),
        "qc_n": len(qc), "qc_excluded": 0,
        "batch_composition": {r["combat_batch"]: {"discovery": int(r["discovery"]),
                                                  "validation": int(r["validation"])} for r in comp},
    }

    # ---- Fig 3: exploratory structure
    pca = read_csv(SRC["pca_scores"])
    ve = load_json(SRC["phase8a_summary"])["pca"]["variance_explained_pc1_20"]
    V["fig3"] = {
        "variance_explained_pc1_20": ve,
        "pc12": [[float(r["PC1"]), float(r["PC2"])] for r in pca],
        "batch_centroids": {b: [float(np.mean([float(r["PC1"]) for r in pca if r["combat_batch"] == b])),
                                float(np.mean([float(r["PC2"]) for r in pca if r["combat_batch"] == b])),
                                sum(1 for r in pca if r["combat_batch"] == b)] for b in batches},
        "pc3_by_batch": {b: [float(r["PC3"]) for r in pca if r["combat_batch"] == b] for b in batches},
    }

    # ---- Fig 4: pathway scores
    sdist = {r["pathway"]: r for r in read_csv(SRC["score_distributions"])}
    corr_rows = read_csv(SRC["score_corr_5set"])
    corr_names = [k for k in corr_rows[0].keys() if k]
    rep = {r["pathway"]: r for r in read_csv(SRC["representation"])}
    order5 = PATHWAYS + [PI3K]
    V["fig4"] = {
        "quantiles": {p: {q: float(sdist[p][q]) for q in ("p5", "p25", "p50", "p75", "p95")} | {
            "n_genes_used": int(sdist[p]["n_genes_used"])} for p in order5},
        "corr_order": order5,
        "corr": [[float(next(r for r in corr_rows if r[""] == a)[b]) for b in order5] for a in order5],
        "representation_pearson": {p: float(rep[p]["pearson_median_vs_repiqr"]) for p in order5},
    }
    if set(corr_names) != set(order5):
        raise SystemExit("unexpected pathways in the correlation matrix")

    # ---- Fig 5: M1
    m1 = ex["apparent"]["M1"]
    V["fig5"] = {
        "terms": {t: {"hr": m1["hazard_ratios"][t], "lo": m1["hr_ci_lower"][t], "hi": m1["hr_ci_upper"][t]}
                  for t, _ in M1_TERMS},
        "age_knots": m1["age_knots"],
        "age_coef": [m1["coefficients"]["age_rcs1"], m1["coefficients"]["age_rcs2"]],
        "age_se": [m1["se"]["age_rcs1"], m1["se"]["age_rcs2"]],
    }
    k = m1["age_knots"]
    grid = np.linspace(k[0], k[2], 67)
    ref = float(k[1])
    b = np.asarray(V["fig5"]["age_coef"])
    lp = rcs_basis(grid, k) @ b - (rcs_basis(np.array([ref]), k) @ b)[0]
    V["fig5"]["age_curve"] = {"age": [float(x) for x in grid], "log_hr_vs_ref": [float(x) for x in lp],
                              "reference_age": ref}

    # ---- Fig 6: molecular HRs
    def hrs(app, names):
        return {p: {"hr": app["hazard_ratios"][p], "lo": app["hr_ci_lower"][p], "hi": app["hr_ci_upper"][p],
                    "coef": app["coefficients"][p], "se": app["se"][p]} for p in names}
    V["fig6"] = {
        "primary_median": {"role": "primary", "model": "M2 four-set, median representation",
                           "lambda": ex["apparent"]["M2_primary"]["lambda"],
                           "hr": hrs(ex["apparent"]["M2_primary"], PATHWAYS)},
        "s7_repiqr": {"role": "sensitivity", "model": "M2 four-set, representative-IQR (S-7)",
                      "lambda": s7["apparent"]["M2_primary"]["lambda"],
                      "hr": hrs(s7["apparent"]["M2_primary"], PATHWAYS)},
        "five_set": {"role": "sensitivity", "model": "M2 five-set, median representation",
                     "lambda": ex["apparent"]["M2_sensitivity"]["lambda"],
                     "hr": hrs(ex["apparent"]["M2_sensitivity"], PATHWAYS + [PI3K])},
    }

    # ---- Fig 7: discrimination
    models = [("M1", "M1 clinical", "reference"), ("M2_primary", "M2 four-set", "primary"),
              ("M2_sensitivity", "M2 five-set", "sensitivity")]
    V["fig7"] = {m: {"label": lab, "role": role, "apparent_C": ex["apparent"][m]["C"],
                     "corrected_C": ex["monte_carlo"][m]["corrected_C_mean"],
                     "mc_se": ex["monte_carlo"][m]["monte_carlo_SE"]} for m, lab, role in models}
    V["fig7"]["delta_corrected"] = theta

    # ---- Fig 8: corrected dC (verified results + checkpoints read-only)
    runs = {}
    for r in RUNS:
        x = cir["runs"][str(r)]
        runs[str(r)] = {"percentile": [x["percentile_primary"]["lower"], x["percentile_primary"]["upper"]],
                        "width": x["percentile_primary"]["width"],
                        "basic": [x["basic_diagnostic"]["lower"], x["basic_diagnostic"]["upper"]],
                        "bc": [x["bc_z0_diagnostic"]["lower"], x["bc_z0_diagnostic"]["upper"]],
                        "median": x["distribution_dC_corrected"]["median"],
                        "share_below_zero": x["distribution_dC_corrected"]["share_below_zero"],
                        "z0": x["z0"], "n": x["n_usable"]}
    outer = {}
    for r, rel in zip(RUNS, CHECKPOINTS):
        vals = []
        with open(ROOT / rel, "rb") as f:
            head = json.loads(f.readline())
            if head.get("kind") != "header" or head.get("run_id") != r:
                raise SystemExit(f"unexpected checkpoint header in {rel}")
            for line in f:
                vals.append(json.loads(line)["record"]["dC_corrected"])
        if len(vals) != 750:
            raise SystemExit(f"{rel}: {len(vals)} records, expected 750")
        outer[str(r)] = vals
        lo, hi = np.quantile(vals, [0.025, 0.975], method="linear")
        if abs(lo - runs[str(r)]["percentile"][0]) > 1e-15 or abs(hi - runs[str(r)]["percentile"][1]) > 1e-15:
            raise SystemExit(f"checkpoint {rel} does not reproduce the verified interval")
    V["fig8"] = {"theta_hat": theta, "B_out": cir["level_B_out"], "B_in": cir["B_in"], "runs": runs,
                 "between_run": cir["between_run"], "pooled_descriptive": cir["pooled_descriptive"]["percentile"],
                 "verification": {"verdict": ver["verdict"], "status_counts": ver["status_counts"]}}
    V["_outer"] = outer   # not written to figure_values (750 x 5 values); hashed instead

    # ---- Fig 9: robustness
    lobo = ex["batch_sensitivity"]["leave_one_batch_out"]
    strat = ex["batch_sensitivity"]["stratified_fit"]
    terms9 = PATHWAYS + ["tnm.stage=3.0"]
    V["fig9"] = {
        "coef_median": {p: [ex["apparent"]["M2_primary"]["coefficients"][p], ex["apparent"]["M2_primary"]["se"][p]] for p in PATHWAYS},
        "coef_s7": {p: [s7["apparent"]["M2_primary"]["coefficients"][p], s7["apparent"]["M2_primary"]["se"][p]] for p in PATHWAYS},
        "coef_stratified": {p: [strat["coefficients"][p], strat["se"][p]] for p in PATHWAYS},
        "lobo": [{"omitted_batch": x["omitted_batch"], "n_omitted": x["n_omitted"], "C_apparent": x["C_apparent"],
                  "coef": {t: x["coefficients"][t] for t in terms9}} for x in lobo],
        "all_data_coef": {t: ex["apparent"]["M2_primary"]["coefficients"][t] for t in terms9},
        "all_data_C": ex["apparent"]["M2_primary"]["C"],
        "calibration": {m: {"apparent": ex["apparent"][m]["calibration_slope"],
                            "corrected": float(np.mean([s["corrected_slope"] for s in ex["monte_carlo"][m]["per_seed"]]))}
                        for m in ("M1", "M2_primary", "M2_sensitivity")},
        "calibration_s7_M2": {"apparent": s7["apparent"]["M2_primary"]["calibration_slope"],
                              "corrected": float(np.mean([s["corrected_slope"] for s in s7["monte_carlo"]["M2_primary"]["per_seed"]]))},
    }

    # ---- S1: Monte Carlo precision at 750
    ep1 = cir["sizing_decision_at_level"]["ep1"]
    V["figS1"] = {"rows": [{"run": rw["run"], "endpoint": rw["endpoint"], "SE_by_bw": rw["SE_by_bw"],
                            "governing_SE": rw["governing_SE"], "ratio": rw["ratio"], "width": rw["width"]}
                           for rw in ep1["rows"]],
                  "rho_max": ep1["rho_max"], "threshold": ep1["fraction"], "action": cir["sizing_decision_at_level"]["action"]}

    # ---- S2: informative censoring profile (frozen Phase 6A table)
    V["figS2"] = {g: dict(zip(("n", "median_age", "stage_III_pct", "chemo_pct"), doc_fact("endpoint_doc", g)))
                  for g in ("died_no_recurrence", "alive_no_recurrence", "recurrence")}
    return V


# ---------------------------------------------------------------- figures
def fig1(V):
    d = V["fig1"]
    fig = plt.figure(figsize=(10, 5.6))
    ax = fig.add_axes([0.02, 0.05, 0.56, 0.9]); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    panel_label(ax, "a")

    def box(x, y, w, h, txt, fc="#f4f3f0", ec=INK2, bold=False):
        ax.add_patch(plt.Rectangle((x, y), w, h, fc=fc, ec=ec, lw=1.0))
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=8.5,
                fontweight="bold" if bold else "normal")

    def arrow(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.0))
    box(1.0, 8.6, 4.4, 1.0, f"GSE39582 arrays (GPL570)\nn = {d['arrays_total']}")
    box(6.1, 8.35, 3.8, 1.5, f"Excluded: {d['excluded_nontumoural']} non-tumoural;\n{d['excluded_stage_other']} stage 0/I/IV or missing;\n"
        f"{d['excluded_no_rfs']} stage II–III without RFS", fc=SURFACE)
    arrow(3.2, 8.6, 3.2, 7.6); arrow(5.4, 9.1, 6.1, 9.1)
    box(1.0, 6.6, 4.4, 1.0, f"S0: stage II–III, documented RFS\nn = {d['S0']['n']} ({d['S0']['events']} events)")
    arrow(3.2, 6.6, 3.2, 5.6)
    box(6.1, 6.6, 3.8, 1.0, f"Endpoint criteria C1 + C2\n(Phase 1): {d['S0']['n'] - d['S3']['n']} removed", fc=SURFACE)
    arrow(5.4, 7.1, 6.1, 7.1)
    box(1.0, 4.6, 4.4, 1.0, f"S3 primary cohort\nn = {d['S3']['n']} ({d['S3']['events']} events, {d['S3']['censored']} censored)")
    arrow(3.2, 4.6, 3.2, 3.6)
    ex0 = d["A19_excluded"][0]
    box(6.1, 4.6, 3.8, 1.0, f"Complete-case (A-19): 1 excluded\n{ex0['sample_id']} ({ex0['partition']}, censored;\nchemotherapy missing)", fc=SURFACE)
    arrow(5.4, 5.1, 6.1, 5.1)
    box(1.0, 2.4, 4.4, 1.2, f"A-19 analysis population\nn = {d['A19']['n']}: {d['A19']['events']} recurrences,\n{d['A19']['censored']} censored",
        fc="#dbe8f8", ec=C_PRIMARY, bold=True)
    pd_ = d["A19_by_partition"]
    arrow(2.2, 2.4, 1.6, 1.5); arrow(4.2, 2.4, 4.8, 1.5)
    box(0.0, 0.3, 3.2, 1.2, f"Discovery partition\nn = {pd_['discovery']['n']}: {pd_['discovery']['events']} events,\n{pd_['discovery']['censored']} censored", fc=SURFACE)
    box(3.6, 0.3, 3.2, 1.2, f"CIT validation partition\nn = {pd_['validation']['n']}: {pd_['validation']['events']} events,\n{pd_['validation']['censored']} censored", fc=SURFACE)
    ax.text(7.0, 0.9, "Partitions shown for description only;\nboth are analysed together.\nThe CIT partition is internal,\nnot an external test set.", fontsize=7.5, color=INK2, va="center")

    ax2 = fig.add_axes([0.60, 0.05, 0.39, 0.9]); ax2.axis("off"); panel_label(ax2, "b")
    txt = ("Endpoint and censoring\n\n"
           "• Event: first recurrence (cause-specific),\n  full follow-up (rfs_event_full / rfs_time_full).\n"
           f"• {d['deaths_censored']} deaths without recurrence (S3) are\n  censored in the deposited data; this\n  censoring is informative (Fig. S2).\n"
           "• Models estimate the cause-specific hazard\n  of recurrence, not the cumulative incidence\n  of recurrence and not absolute risk.\n"
           "• The rfs.delay time unit is undocumented:\n  only rank-based quantities are reported.\n\n"
           "Analysis\n\n"
           "• M1: clinical Cox model (6 coefficients).\n• M2: M1 + four Hallmark scores (ridge).\n"
           "• Estimand: optimism-corrected ΔC = C(M2) − C(M1).\n"
           f"• Interval: nested bootstrap, 5 runs differing\n  only in seed, {V['fig8']['B_out']} outer × {V['fig8']['B_in']} inner each.")
    ax2.text(0.0, 0.98, txt, va="top", fontsize=8.3, linespacing=1.35)
    return fig


def fig2(V):
    d = V["fig2"]
    fig, axs = plt.subplots(1, 4, figsize=(13, 3.6), gridspec_kw={"width_ratios": [1.1, 1.6, 0.9, 1.4]})
    p = d["probes"]
    ax = axs[0]; panel_label(ax, "a")
    labels = ["Unique gene", "Unmapped", "Ambiguous"]
    vals = [p["n_unique"], p["n_unmapped"], p["n_ambiguous"]]
    ax.barh(labels[::-1], vals[::-1], color=[NEUTRAL, NEUTRAL, C_PRIMARY], height=0.6)
    for i, v in enumerate(vals[::-1]):
        ax.text(v + 600, i, f"{v:,}", va="center", fontsize=8)
    ax.set_xlim(0, 56000); ax.set_xlabel("probe sets"); ax.grid(axis="y", visible=False)
    ax.set_title(f"{p['n_probes_total']:,} probes → {p['n_genes_unique']:,} genes", fontsize=9, loc="right")
    ax.text(0.98, 0.02, f"{p['n_genes_single_probe']:,} single-probe genes\n{p['n_genes_multi_probe']:,} multi-probe genes\n(median aggregated)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5, color=INK2)
    ax = axs[1]; panel_label(ax, "b")
    bs = list(d["qc_median_by_batch"])
    ax.boxplot([d["qc_median_by_batch"][b] for b in bs], showfliers=False, widths=0.55,
               medianprops=dict(color=C_PRIMARY, lw=2), boxprops=dict(color=INK2), whiskerprops=dict(color=INK2),
               capprops=dict(color=INK2))
    lo, hi = d["s_s1_fence"]
    ax.axhspan(lo, hi, color=GRID, alpha=0.7, lw=0, zorder=0)
    ax.set_xticks(range(1, len(bs) + 1), [b.replace("batch", "b") for b in bs], fontsize=7.5)
    ax.set_ylabel("per-array median (QC metric S-S1)"); ax.set_title("QC median by ComBat batch (S3, n = 453)", fontsize=9)
    ax = axs[2]; panel_label(ax, "c")
    f = d["qc_flags"]
    names = ["S-S1", "S-S2", "S-S3", "any", "excl."]
    vv = [f["flag_S_S1"], f["flag_S_S2"], f["flag_S_S3"], d["qc_flag_union"], d["qc_excluded"]]
    ax.bar(range(5), vv, color=[NEUTRAL] * 4 + [C_PRIMARY], width=0.6)
    for i, v in enumerate(vv):
        ax.text(i, v + 0.5, str(v), ha="center", fontsize=8)
    ax.set_xticks(range(5), names, fontsize=7.5); ax.set_ylabel("arrays"); ax.set_ylim(0, 34)
    ax.set_title("Soft QC flags (n = 453)", fontsize=9); ax.grid(axis="x", visible=False)
    ax = axs[3]; panel_label(ax, "d")
    comp = d["batch_composition"]; bb = list(comp)
    disc = [comp[b]["discovery"] for b in bb]; val = [comp[b]["validation"] for b in bb]
    x = np.arange(len(bb))
    ax.bar(x, disc, color=C_PRIMARY, width=0.62, label="discovery", edgecolor=SURFACE, linewidth=1)
    ax.bar(x, val, bottom=disc, color=C_S7, width=0.62, label="CIT validation", edgecolor=SURFACE, linewidth=1)
    ax.set_xticks(x, [b.replace("batch", "b") for b in bb], fontsize=7.5); ax.set_ylabel("arrays (S3)")
    ax.set_title("Batch × partition (confounded)", fontsize=9); ax.legend(loc="upper right"); ax.grid(axis="x", visible=False)
    fig.tight_layout()
    return fig


def fig3(V):
    d = V["fig3"]
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.7), gridspec_kw={"width_ratios": [1, 1.25, 1.4]})
    ax = axs[0]; panel_label(ax, "a")
    ve = np.asarray(d["variance_explained_pc1_20"]) * 100
    ax.bar(np.arange(1, 21), ve, color=NEUTRAL, width=0.7)
    ax.set_xlabel("principal component"); ax.set_ylabel("% variance"); ax.set_title(f"Scree (PC1–20: {ve.sum():.1f} %)", fontsize=9)
    ax.grid(axis="x", visible=False)
    ax = axs[1]; panel_label(ax, "b")
    pc = np.asarray(d["pc12"])
    hb = ax.hexbin(pc[:, 0], pc[:, 1], gridsize=22, cmap="Greys", mincnt=1, linewidths=0.2)
    for b, (cx, cy, n) in d["batch_centroids"].items():
        ax.plot(cx, cy, "o", ms=7, mfc=C_PRIMARY, mec=SURFACE, mew=1.5)
    ax.set_xlabel(f"PC1 ({ve[0]:.1f} %)"); ax.set_ylabel(f"PC2 ({ve[1]:.1f} %)")
    ax.set_title("Array density; blue = 9 batch centroids", fontsize=9)
    fig.colorbar(hb, ax=ax, label="arrays per cell", shrink=0.8)
    ax = axs[2]; panel_label(ax, "c")
    bs = list(d["pc3_by_batch"])
    ax.boxplot([d["pc3_by_batch"][b] for b in bs], showfliers=False, widths=0.55,
               medianprops=dict(color=C_PRIMARY, lw=2), boxprops=dict(color=INK2), whiskerprops=dict(color=INK2), capprops=dict(color=INK2))
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(range(1, len(bs) + 1), [b.replace("batch", "b") for b in bs], fontsize=7.5)
    ax.set_ylabel(f"PC3 ({ve[2]:.1f} %)"); ax.set_title("PC3 by ComBat batch (exploratory)", fontsize=9)
    fig.tight_layout()
    return fig


def fig4(V):
    d = V["fig4"]
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.0), gridspec_kw={"width_ratios": [1.3, 1.2, 1.0]})
    order = d["corr_order"]
    ax = axs[0]; panel_label(ax, "a")
    for i, p in enumerate(order):
        q = d["quantiles"][p]; col = C_FIVE if p == PI3K else C_PRIMARY
        ax.plot([q["p5"], q["p95"]], [i, i], color=col, lw=1.2)
        ax.add_patch(plt.Rectangle((q["p25"], i - 0.22), q["p75"] - q["p25"], 0.44, fc=col, ec=col, alpha=0.35))
        ax.plot([q["p50"], q["p50"]], [i - 0.22, i + 0.22], color=col, lw=2.2)
    ax.set_yticks(range(len(order)), [f"{SHORT[p]} ({d['quantiles'][p]['n_genes_used']} genes)" + (" †" if p == PI3K else "") for p in order])
    ax.invert_yaxis(); ax.axvline(0, color=INK2, lw=0.8)
    ax.set_xlabel("pathway score (mean standardised expression)"); ax.set_title("Score distributions (S3)", fontsize=9)
    ax = axs[1]; panel_label(ax, "b")
    m = np.asarray(d["corr"])
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#f0efec", "#e34948"])
    im = ax.imshow(m, cmap=cmap, vmin=-1, vmax=1)
    for i in range(len(order)):
        for j in range(len(order)):
            ax.text(j, i, f"{m[i, j]:.2f}", ha="center", va="center", fontsize=7.5)
    lab = [SHORT[p] + (" †" if p == PI3K else "") for p in order]
    ax.set_xticks(range(len(order)), lab, rotation=35, ha="right", fontsize=7.5); ax.set_yticks(range(len(order)), lab, fontsize=7.5)
    ax.grid(False); ax.set_title("Pearson correlation of scores", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.8)
    ax = axs[2]; panel_label(ax, "c")
    r = [d["representation_pearson"][p] for p in order]
    ax.barh(range(len(order)), r, color=[C_FIVE if p == PI3K else C_PRIMARY for p in order], height=0.55)
    for i, v in enumerate(r):
        ax.text(v - 0.005, i, f"{v:.3f}", va="center", ha="right", fontsize=7.5, color=SURFACE)
    ax.set_yticks(range(len(order)), lab, fontsize=7.5); ax.invert_yaxis(); ax.set_xlim(0.8, 1.0)
    ax.set_xlabel("r, median vs representative-IQR"); ax.set_title("Representation agreement", fontsize=9); ax.grid(axis="y", visible=False)
    fig.tight_layout()
    return fig


def fig5(V):
    d = V["fig5"]
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.6), gridspec_kw={"width_ratios": [1.4, 1]})
    ax = axs[0]; panel_label(ax, "a")
    for i, (t, lab) in enumerate(M1_TERMS):
        h = d["terms"][t]
        ax.plot([h["lo"], h["hi"]], [i, i], color=INK2, lw=2)
        ax.plot(h["hr"], i, "s", ms=8, mfc=INK, mec=SURFACE, mew=1.2)
        ax.text(3.9, i, f"{h['hr']:.2f} ({h['lo']:.2f}–{h['hi']:.2f})", va="center", fontsize=8)
    ax.axvline(1, color=INK2, lw=1, ls="--")
    ax.set_xscale("log"); ax.set_xlim(0.45, 6.5); ax.set_xticks([0.5, 1, 2, 3], ["0.5", "1", "2", "3"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_yticks(range(len(M1_TERMS)), [lab for _, lab in M1_TERMS]); ax.invert_yaxis()
    ax.set_xlabel("hazard ratio (95 % Wald CI), log scale"); ax.set_title("M1 clinical model: cause-specific HRs (apparent fit)", fontsize=9)
    ax = axs[1]; panel_label(ax, "b")
    c = d["age_curve"]
    ax.plot(c["age"], np.exp(c["log_hr_vs_ref"]), color=INK, lw=2)
    for kn in d["age_knots"]:
        ax.axvline(kn, color=GRID, lw=1.2, zorder=0)
    ax.axhline(1, color=INK2, lw=1, ls="--")
    ax.set_xlabel("age at diagnosis (years)"); ax.set_ylabel(f"HR vs age {c['reference_age']:.1f}")
    ax.set_title("Age (restricted cubic spline; no CI band)", fontsize=9)
    fig.tight_layout()
    return fig


def fig6(V):
    d = V["fig6"]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.0), gridspec_kw={"width_ratios": [1.15, 1]}, sharex=True)
    ax = axs[0]; panel_label(ax, "a")
    for j, (key, col, mk, off, lab) in enumerate([("primary_median", C_PRIMARY, "o", -0.14, "Primary: four-set, median"),
                                                  ("s7_repiqr", C_S7, "s", 0.14, "Sensitivity: four-set, S-7 representation")]):
        for i, p in enumerate(PATHWAYS):
            h = d[key]["hr"][p]
            ax.plot([h["lo"], h["hi"]], [i + off, i + off], color=col, lw=2)
            ax.plot(h["hr"], i + off, mk, ms=7, mfc=col if j == 0 else SURFACE, mec=col, mew=1.8, label=lab if i == 0 else None)
    for i, p in enumerate(PATHWAYS):
        h = d["primary_median"]["hr"][p]
        ax.text(1.33, i - 0.14, f"{h['hr']:.3f} ({h['lo']:.3f}–{h['hi']:.3f})", va="center", fontsize=7.5, color=C_PRIMARY)
    ax.axvline(1, color=INK2, lw=1, ls="--")
    ax.set_yticks(range(len(PATHWAYS)), [SHORT[p] for p in PATHWAYS]); ax.invert_yaxis()
    ax.set_xscale("log"); ax.set_xlim(0.72, 1.55); ax.set_xticks([0.8, 0.9, 1.0, 1.1, 1.2], ["0.8", "0.9", "1.0", "1.1", "1.2"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("HR per training-fold SD of score (95 % Wald CI)")
    ax.set_title("PRIMARY model M2 (four-set) and S-7 representation", fontsize=9)
    ax.legend(loc="lower left", bbox_to_anchor=(0, -0.42), ncol=1)
    ax = axs[1]; panel_label(ax, "b")
    names5 = PATHWAYS + [PI3K]
    for i, p in enumerate(names5):
        h = d["five_set"]["hr"][p]
        ax.plot([h["lo"], h["hi"]], [i, i], color=C_FIVE, lw=2)
        ax.plot(h["hr"], i, "^", ms=8, mfc=SURFACE, mec=C_FIVE, mew=1.8)
        ax.text(1.33, i, f"{h['hr']:.3f} ({h['lo']:.3f}–{h['hi']:.3f})", va="center", fontsize=7.5, color=INK2)
    ax.axvline(1, color=INK2, lw=1, ls="--")
    ax.set_yticks(range(len(names5)), [SHORT[p] for p in names5]); ax.invert_yaxis()
    ax.set_xlabel("HR per training-fold SD of score (95 % Wald CI)")
    ax.set_title("SENSITIVITY model: five-set (adds PI3K; correlated scores)", fontsize=9)
    ax.text(0.0, -0.3, "Not promoted: E2F upper limit 0.996 arises with collinearity\n(p53↔PI3K r = 0.64; SEs inflated); see Fig. 4b.",
            transform=ax.transAxes, fontsize=7.5, color=INK2, va="top")
    fig.tight_layout()
    return fig


def fig7(V):
    d = V["fig7"]
    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    keys = ["M1", "M2_primary", "M2_sensitivity"]
    cols = {"M1": INK, "M2_primary": C_PRIMARY, "M2_sensitivity": C_FIVE}
    for i, k in enumerate(keys):
        a, c = d[k]["apparent_C"], d[k]["corrected_C"]
        ax.annotate("", xy=(c, i), xytext=(a, i), arrowprops=dict(arrowstyle="-|>", color=NEUTRAL, lw=1.4))
        ax.plot(a, i, "o", ms=8, mfc=SURFACE, mec=cols[k], mew=1.8)
        ax.plot(c, i, "o", ms=9, mfc=cols[k], mec=SURFACE, mew=1.2)
        ax.text(a + 0.0012, i, f"apparent {a:.4f}", fontsize=7.5, color=INK2, va="center")
        ax.text(c - 0.0012, i, f"corrected {c:.4f}", fontsize=8, ha="right", va="center", fontweight="bold")
        ax.text((a + c) / 2, i - 0.2, f"optimism {a - c:.4f}", fontsize=7.2, color=INK2, ha="center")
    ax.set_yticks(range(3), [d[k]["label"] + (" (primary)" if d[k]["role"] == "primary" else
                                             " (sensitivity)" if d[k]["role"] == "sensitivity" else "") for k in keys])
    ax.set_ylim(2.5, -0.6); ax.set_xlim(0.583, 0.664); ax.set_xlabel("Harrell's C (A-19, n = 452, 136 events)")
    ax.set_title(f"Apparent (open) vs optimism-corrected (filled) C; corrected ΔC(M2 − M1) = {d['delta_corrected']:+.4f}", fontsize=9)
    fig.tight_layout()
    return fig


def fig8(V):
    d = V["fig8"]; outer = V["_outer"]; theta = d["theta_hat"]
    fig = plt.figure(figsize=(12.5, 5.2))
    gs = fig.add_gridspec(5, 2, width_ratios=[1.05, 1], hspace=0.25, wspace=0.18)
    ax = fig.add_subplot(gs[:, 0]); panel_label(ax, "a")
    for i, r in enumerate(RUNS):
        x = d["runs"][str(r)]
        lo, hi = x["percentile"]
        ax.plot([lo, hi], [i, i], color=C_PRIMARY, lw=3.2, solid_capstyle="butt")
        ax.plot([lo, lo], [i - 0.15, i + 0.15], color=C_PRIMARY, lw=2); ax.plot([hi, hi], [i - 0.15, i + 0.15], color=C_PRIMARY, lw=2)
        ax.plot(x["basic"], [i + 0.3, i + 0.3], color=NEUTRAL, lw=1.1)
        ax.plot(x["bc"], [i + 0.42, i + 0.42], color=NEUTRAL, lw=1.1, ls=(0, (3, 2)))
        ax.text(0.061, i, f"[{lo:+.4f}, {hi:+.4f}]", va="center", fontsize=8)
    ax.axvline(0, color=INK2, lw=1, ls="--"); ax.axvline(theta, color=INK, lw=1.4)
    ax.text(theta + 0.0006, -0.65, f"θ̂ = {theta:+.6f}", fontsize=8)
    ax.text(0.0006, 4.75, "0", fontsize=8, color=INK2)
    ax.set_yticks(range(5), [f"run {r}" for r in RUNS]); ax.set_ylim(4.8, -0.9)
    ax.set_xlim(-0.048, 0.085); ax.set_xlabel("optimism-corrected ΔC = C(M2) − C(M1)")
    ax.set_title("95 % percentile intervals (primary, blue); basic and BC diagnostics (grey)", fontsize=9)
    edges = np.arange(-0.016, 0.1061, 0.002)
    for i, r in enumerate(RUNS):
        a2 = fig.add_subplot(gs[i, 1])
        if i == 0:
            panel_label(a2, "b")
        x = d["runs"][str(r)]
        a2.hist(outer[str(r)], bins=edges, color=C_PRIMARY, alpha=0.55, edgecolor=SURFACE, linewidth=0.3)
        for e in x["percentile"]:
            a2.axvline(e, color=C_PRIMARY, lw=1.2)
        a2.axvline(0, color=INK2, lw=0.9, ls="--"); a2.axvline(theta, color=INK, lw=1.1)
        a2.text(0.98, 0.8, f"run {r}: median {x['median']:+.4f}; {100 * x['share_below_zero']:.1f} % < 0", transform=a2.transAxes,
                ha="right", fontsize=7.2, color=INK2)
        a2.set_xlim(-0.02, 0.085); a2.set_yticks([]); a2.grid(axis="y", visible=False)
        if i < 4:
            a2.set_xticklabels([])
        else:
            a2.set_xlabel(f"ΔC_corrected(k), {d['B_out']} outer replicates per run")
    return fig


def fig9(V):
    d = V["fig9"]
    fig, axs = plt.subplots(2, 2, figsize=(12, 7.2))
    ax = axs[0, 0]; panel_label(ax, "a")
    specs = [("coef_median", C_PRIMARY, "o", -0.2, "median (primary)"), ("coef_s7", C_S7, "s", 0.0, "S-7 representation"),
             ("coef_stratified", INK2, "D", 0.2, "batch-stratified")]
    for key, col, mk, off, lab in specs:
        for i, p in enumerate(PATHWAYS):
            c, se = d[key][p]
            ax.plot([c - 1.96 * se, c + 1.96 * se], [i + off, i + off], color=col, lw=1.8)
            ax.plot(c, i + off, mk, ms=6.5, mfc=col if key == "coef_median" else SURFACE, mec=col, mew=1.6, label=lab if i == 0 else None)
    ax.axvline(0, color=INK2, lw=1, ls="--")
    ax.set_yticks(range(4), [SHORT[p] for p in PATHWAYS]); ax.invert_yaxis()
    ax.set_xlabel("log-HR per SD (±1.96 SE)"); ax.set_title("Four-set coefficients: representation and batch stratification", fontsize=9)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3)
    ax = axs[0, 1]; panel_label(ax, "b")
    terms = PATHWAYS + ["tnm.stage=3.0"]
    nmax = max(x["n_omitted"] for x in d["lobo"])
    for i, t in enumerate(terms):
        for x in d["lobo"]:
            ax.scatter(x["coef"][t], i, s=12 + 70 * x["n_omitted"] / nmax, color=C_PRIMARY, alpha=0.55, edgecolor=SURFACE, linewidth=0.6)
        ax.plot(d["all_data_coef"][t], i, "|", ms=16, mew=2.4, color=INK)
    ax.axvline(0, color=INK2, lw=1, ls="--")
    ax.set_yticks(range(len(terms)), [SHORT.get(t, "Stage III (clinical)") for t in terms]); ax.invert_yaxis()
    ax.set_xlabel("log-HR; dot = refit omitting one batch (size ∝ omitted n); bar = all data")
    ax.set_title("Leave-one-batch-out coefficients (M2 four-set)", fontsize=9)
    ax = axs[1, 0]; panel_label(ax, "c")
    lob = sorted(d["lobo"], key=lambda x: x["omitted_batch"])
    ax.bar(range(len(lob)), [x["C_apparent"] for x in lob], color=C_PRIMARY, width=0.6)
    ax.axhline(d["all_data_C"], color=INK, lw=1.4)
    ax.text(len(lob) - 0.5, d["all_data_C"] + 0.0015, f"all data {d['all_data_C']:.3f}", ha="right", fontsize=8)
    ax.set_xticks(range(len(lob)), [f"−{x['omitted_batch'].replace('batch', 'b')}\n(n={x['n_omitted']})" for x in lob], fontsize=7)
    ax.set_ylim(0.60, 0.67); ax.set_ylabel("apparent C (M2)"); ax.set_title("Leave-one-batch-out apparent C (not optimism-corrected)", fontsize=9)
    ax.grid(axis="x", visible=False)
    ax = axs[1, 1]; panel_label(ax, "d")
    rows = [("M1", d["calibration"]["M1"], INK), ("M2 four-set (primary)", d["calibration"]["M2_primary"], C_PRIMARY),
            ("M2 four-set, S-7", d["calibration_s7_M2"], C_S7), ("M2 five-set (sensitivity)", d["calibration"]["M2_sensitivity"], C_FIVE)]
    for i, (lab, v, col) in enumerate(rows):
        ax.annotate("", xy=(v["corrected"], i), xytext=(v["apparent"], i), arrowprops=dict(arrowstyle="-|>", color=NEUTRAL, lw=1.3))
        ax.plot(v["apparent"], i, "o", ms=7, mfc=SURFACE, mec=col, mew=1.6)
        ax.plot(v["corrected"], i, "o", ms=8, mfc=col, mec=SURFACE, mew=1.0)
        ax.text(v["corrected"] - 0.006, i - 0.22, f"{v['corrected']:.3f}", ha="right", fontsize=8)
    ax.axvline(1, color=INK2, lw=1, ls="--")
    ax.set_yticks(range(4), [r[0] for r in rows]); ax.invert_yaxis(); ax.set_xlim(0.78, 1.12)
    ax.set_xlabel("calibration slope: apparent (open) → optimism-corrected (filled)"); ax.set_title("Calibration slope (log-hazard scale)", fontsize=9)
    fig.tight_layout()
    return fig


def figS1(V):
    d = V["figS1"]
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.6))
    ax = axs[0]; panel_label(ax, "a")
    for j, ep in enumerate(("lower", "upper")):
        for rw in [x for x in d["rows"] if x["endpoint"] == ep]:
            for bw, mk in (("0.5", "v"), ("1.0", "o"), ("2.0", "^")):
                ax.plot(rw["run"] + (j - 0.5) * 0.3, rw["SE_by_bw"][bw] / rw["width"], mk, ms=6,
                        mfc=C_PRIMARY if ep == "upper" else SURFACE, mec=C_PRIMARY, mew=1.3)
    ax.axhline(d["threshold"], color=INK, lw=1.2, ls="--"); ax.text(4.4, d["threshold"] + 0.001, "EP-1 limit 0.05", ha="right", fontsize=8)
    ax.set_xticks(range(5), [f"run {r}" for r in RUNS]); ax.set_ylim(0, 0.058)
    ax.set_ylabel("kernel MC-SE / interval width"); ax.set_title(f"EP-1 ratio at B_out = 750 (ρ_max = {d['rho_max']:.4f})", fontsize=9)
    ax = axs[1]; panel_label(ax, "b")
    br = V["fig8"]["between_run"]
    for j, ep in enumerate(("lower", "upper")):
        vals = br[ep]["values"]
        ax.plot(range(5), vals, "o-", color=C_PRIMARY if ep == "upper" else INK2, lw=1.2, ms=6, label=f"{ep}: SD {br[ep]['sd_ddof1']:.5f}")
    ax.set_xticks(range(5), [f"run {r}" for r in RUNS]); ax.set_ylabel("percentile endpoint"); ax.axhline(0, color=INK2, lw=0.8, ls="--")
    ax.set_title("Percentile endpoints across the five runs", fontsize=9); ax.legend(loc="center right")
    fig.tight_layout()
    return fig


def figS2(V):
    d = V["figS2"]
    groups = [("died_no_recurrence", "Died without\nrecurrence\n(censored)", C_S7), ("alive_no_recurrence", "Alive without\nrecurrence\n(censored)", NEUTRAL),
              ("recurrence", "Recurrence\n(event)", C_PRIMARY)]
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.3))
    for ax, (key, title, lim) in zip(axs, [("median_age", "Median age (years)", 85), ("stage_III_pct", "Stage III (%)", 70), ("chemo_pct", "Adjuvant chemotherapy (%)", 70)]):
        vals = [d[g][key] for g, _, _ in groups]
        ax.bar(range(3), vals, color=[c for _, _, c in groups], width=0.6)
        for i, v in enumerate(vals):
            ax.text(i, v + lim * 0.02, f"{v:g}", ha="center", fontsize=8)
        ax.set_xticks(range(3), [f"{lab}\nn = {d[g]['n']}" for g, lab, _ in groups], fontsize=7.5)
        ax.set_ylim(0, lim); ax.set_title(title, fontsize=9); ax.grid(axis="x", visible=False)
    panel_label(axs[0], "a"); panel_label(axs[1], "b"); panel_label(axs[2], "c")
    fig.tight_layout()
    return fig


FIGS = [("Fig1_study_design_cohort", fig1), ("Fig2_preprocessing_qc", fig2), ("Fig3_exploratory_structure", fig3),
        ("Fig4_pathway_scores", fig4), ("Fig5_clinical_model", fig5), ("Fig6_molecular_hazard_ratios", fig6),
        ("Fig7_discrimination", fig7), ("Fig8_corrected_deltaC", fig8), ("Fig9_robustness", fig9),
        ("FigS1_monte_carlo_precision", figS1), ("FigS2_informative_censoring", figS2)]


# ---------------------------------------------------------------- captions
def captions(V):
    f1, f6, f7, f8, f9 = V["fig1"], V["fig6"], V["fig7"], V["fig8"], V["fig9"]
    pd_ = f1["A19_by_partition"]; runs = f8["runs"]
    ivs = "; ".join(f"run {r} [{runs[str(r)]['percentile'][0]:+.4f}, {runs[str(r)]['percentile'][1]:+.4f}]" for r in RUNS)
    h = f6["primary_median"]["hr"]; h5 = f6["five_set"]["hr"]
    cal = f9["calibration"]
    C = {}
    C["Fig1_study_design_cohort"] = (
        f"**Figure 1. Study design, cohort and endpoint.** (a) Selection from the {f1['arrays_total']} GSE39582 arrays to the S0 cohort "
        f"(n = {f1['S0']['n']}), the S3 primary cohort (n = {f1['S3']['n']}; {f1['S3']['events']} recurrences, {f1['S3']['censored']} censored) and the "
        f"A-19 complete-case analysis population (n = {f1['A19']['n']}; {f1['A19']['events']} recurrences, {f1['A19']['censored']} censored), which excludes one "
        f"censored discovery-partition patient with missing adjuvant-chemotherapy status ({f1['A19_excluded'][0]['sample_id']}). Descriptive partition counts "
        f"within A-19: discovery {pd_['discovery']['n']} ({pd_['discovery']['events']} events, {pd_['discovery']['censored']} censored); CIT validation "
        f"{pd_['validation']['n']} ({pd_['validation']['events']} events, {pd_['validation']['censored']} censored). The frozen partition labels were used "
        "unchanged; the partitions are analysed together, and the CIT partition is internal, not an external test set. (b) Endpoint, censoring and analysis "
        "design. The models estimate the cause-specific hazard of recurrence, not absolute recurrence probability or cumulative incidence. "
        "Sources: reports/phase1/cohort_counts.csv; docs/PHASE7C_QC_REPORT.md (D1); frozen analysis bundle; config/analysis/A19_FREEZE_RECORD.json.")
    p = V["fig2"]["probes"]
    C["Fig2_preprocessing_qc"] = (
        f"**Figure 2. Expression preprocessing and quality control (descriptive).** (a) Of {p['n_probes_total']:,} GPL570 probe sets, {p['n_unique']:,} map "
        f"to one gene symbol and are used; {p['n_unmapped']:,} unmapped and {p['n_ambiguous']:,} ambiguous probe sets are excluded. Genes with several "
        f"probe sets ({p['n_genes_multi_probe']:,} of {p['n_genes_unique']:,}) are summarised by the median. (b) Per-array median (the S-S1 QC metric) by ComBat batch "
        "(S3, n = 453; boxes 25–75 %, whiskers 1.5 IQR, no individual points); grey band = the S-S1 fence (arrays outside it are flagged). "
        f"(c) Soft QC flags (S-S1 median, S-S2 IQR, S-S3 summed 1 − Spearman ρ): {V['fig2']['qc_flag_union']} arrays flagged by at least one criterion; none excluded, because no independent technical "
        "evidence of invalidity was found. (d) Batch composition by partition: processing batch is confounded with the discovery/CIT partition. "
        "Sources: reports/phase2/phase2_mapping_audit.json; reports/phase8a/sample_distribution_summary.csv; reports/phase7c/sample_qc_metrics.csv, batch_composition.csv.")
    ve = V["fig3"]["variance_explained_pc1_20"]
    C["Fig3_exploratory_structure"] = (
        f"**Figure 3. Exploratory expression structure (exploratory; no outcome used).** (a) Variance explained by PC1–20 of the gene-level matrix "
        f"(PC1 {100 * ve[0]:.1f} %, PC2 {100 * ve[1]:.1f} %, PC3 {100 * ve[2]:.1f} %). (b) Density of the 453 arrays in PC1–PC2; the nine ComBat-batch centroids (blue) lie close together near the origin. "
        "(c) PC3 by batch. Batch is associated with residual structure (Phase 8B: ComBat batch with PC3, ε² = 0.061); this is an association, not a "
        "demonstrated technical effect on outcome. Sources: reports/phase8a/pca_scores_gene_level.csv, phase8a_summary.json.")
    rp = V["fig4"]["representation_pearson"]
    C["Fig4_pathway_scores"] = (
        "**Figure 4. Pathway-score characterisation.** Scores are the mean of standardised member-gene mRNA levels (S3, n = 453, full-cohort "
        "standardisation for display; the models re-estimate standardisation inside every fit). They are expression summaries, not measures of pathway "
        "activation, protein activity or mutation status. In (a), boxes span 25–75 % and lines 5–95 %. † PI3K/AKT/mTOR belongs to the five-set sensitivity panel only. (a) Score distributions. "
        f"(b) Pearson correlations: p53 ↔ PI3K r = {V['fig4']['corr'][3][4]:.2f}, the collinearity that destabilises the five-set model. "
        f"(c) Agreement between median and representative-IQR summaries (r = {min(rp.values()):.3f}–{max(rp.values()):.3f}; lowest for Wnt). "
        "Sources: reports/phase8c/score_distributions.csv, score_correlation_pearson_5set_sensitivity.csv, representation_sensitivity.csv.")
    t = V["fig5"]["terms"]
    C["Fig5_clinical_model"] = (
        "**Figure 5. Clinical model M1 (apparent, unpenalised Cox fit; A-19, n = 452, 136 recurrences).** (a) Adjusted cause-specific hazard ratios "
        f"with 95 % Wald confidence intervals: stage III vs II {t['tnm.stage=3.0']['hr']:.2f} ({t['tnm.stage=3.0']['lo']:.2f}–{t['tnm.stage=3.0']['hi']:.2f}). "
        "The adjuvant-chemotherapy HR is confounded by indication and has no causal interpretation. Hazard ratios describe association, not added "
        "predictive discrimination. (b) Age enters as a restricted cubic spline (knots at the 10th/50th/90th training percentiles, grey lines) and has no "
        "single HR; the curve is the stored spline evaluated relative to the median knot. No confidence band is drawn, because the coefficient "
        "covariance was not stored. Source: reports/phase9/execution/phase9_execution_results.json (apparent.M1).")
    C["Fig6_molecular_hazard_ratios"] = (
        "**Figure 6. Molecular hazard ratios (historical single-layer Phase 9 fits; A-19, n = 452).** HR per one training-fold SD of each score, from "
        "ridge-penalised Cox models; Wald intervals from a penalised fit are approximate. (a) **Primary** four-set model M2 (median representation, "
        f"λ = {f6['primary_median']['lambda']:.3f}; filled circles): E2F {h['HALLMARK_E2F_TARGETS']['hr']:.3f} "
        f"({h['HALLMARK_E2F_TARGETS']['lo']:.3f}–{h['HALLMARK_E2F_TARGETS']['hi']:.3f}), EMT {h['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']['hr']:.3f} "
        f"({h['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']['lo']:.3f}–{h['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']['hi']:.3f}); every interval includes 1. "
        "Open squares: the S-7 representation **sensitivity** analysis (same patients and arrays, a different probe summary; not independent replication). "
        f"(b) Five-set **sensitivity** model (λ = {f6['five_set']['lambda']:.3f}): E2F {h5['HALLMARK_E2F_TARGETS']['hr']:.3f} "
        f"({h5['HALLMARK_E2F_TARGETS']['lo']:.3f}–{h5['HALLMARK_E2F_TARGETS']['hi']:.3f}) narrowly excludes 1. It coincides with collinearity introduced "
        "by PI3K and is not promoted to a primary finding. Sources: execution and S-7 results JSON.")
    C["Fig7_discrimination"] = (
        "**Figure 7. Apparent versus optimism-corrected discrimination (A-19, n = 452, 136 recurrences).** Open circles: apparent Harrell's C on the "
        "fitting data. Filled circles: C after bootstrap optimism correction (5 seeds × 500 event-stratified resamples). Corrected C: M1 "
        f"{f7['M1']['corrected_C']:.4f}, M2 four-set (primary) {f7['M2_primary']['corrected_C']:.4f}, five-set (sensitivity) "
        f"{f7['M2_sensitivity']['corrected_C']:.4f}. M2's larger optimism removes most of its apparent advantage; the corrected difference is "
        f"{f7['delta_corrected']:+.4f} (uncertainty in Fig. 8). Monte Carlo SEs (≤ {max(f7[k]['mc_se'] for k in ('M1', 'M2_primary', 'M2_sensitivity')):.5f}) "
        "exclude full-data cross-validation randomness. Source: reports/phase9/execution/phase9_execution_results.json.")
    C["Fig8_corrected_deltaC"] = (
        "**Figure 8. Optimism-corrected ΔC and its uncertainty (central result).** Corrected nested bootstrap, five runs that differ only in random seed "
        f"(not independent cohorts or experiments), each with {f8['B_out']} event-stratified outer resamples × {f8['B_in']} inner resamples; frozen "
        f"verification PASS ({f8['verification']['status_counts']['PASS']}/{sum(f8['verification']['status_counts'].values())} checks) and reproduced by an independent "
        f"re-implementation within this project from the same checkpoints (not external replication). (a) Blue: 95 % percentile intervals (primary): {ivs}. All five include zero (dashed line). Black line: the pre-specified point "
        f"estimate θ̂ = {f8['theta_hat']:+.6f}. Grey: basic (solid) and z₀-only bias-corrected (dotted) intervals, diagnostics only. Their divergence from "
        "the percentile interval reflects the bootstrap distribution lying above θ̂ (warnings A7/A8). (b) Distributions of ΔC_corrected(k) per run with "
        "percentile endpoints. These are bootstrap uncertainty intervals for the corrected ΔC. They are not Wald confidence intervals, and not the "
        "historical paired out-of-original ΔC intervals, which belong to a different quantity. The intervals are conditional on the observed 136 events "
        "and 316 censored observations, and a zero-including interval does not establish equivalence. Sources: "
        "reports/phase9/corrected_rerun_report/corrected_interval_results.json; checkpoints run0–4 (read-only; percentiles re-derived and matched).")
    C["Fig9_robustness"] = (
        "**Figure 9. Robustness and reliability (historical single-layer fits; descriptive, no intervals implied beyond those shown).** (a) Four-set "
        "coefficients under the primary median representation, the S-7 representation and a batch-stratified fit (nine ComBat strata). (b) "
        "Leave-one-batch-out refits: coefficients vary several-fold, and Wnt and p53 change sign. Omitting the largest batches removes about 23 % of the "
        "cohort, so the spread conflates batch identity with sample size. (c) Leave-one-batch-out apparent C (not optimism-corrected). (d) Calibration "
        f"slopes: corrected slopes M1 {cal['M1']['corrected']:.3f}, M2 {cal['M2_primary']['corrected']:.3f}, five-set {cal['M2_sensitivity']['corrected']:.3f} "
        f"(S-7 M2 {V['fig9']['calibration_s7_M2']['corrected']:.3f}). Values below 1 indicate residual overfitting. Batch effects are characterised, not "
        "resolved. Sources: execution JSON (batch_sensitivity, apparent, monte_carlo); S-7 JSON.")
    s1 = V["figS1"]
    C["FigS1_monte_carlo_precision"] = (
        f"**Figure S1. Monte Carlo precision of the corrected intervals at B_out = {f8['B_out']}.** (a) Kernel quantile MC-SE of each percentile endpoint "
        f"divided by interval width, at bandwidths 0.5h (▼), h (●) and 2h (▲) (Silverman); open = lower, filled = upper endpoint. The governing maximum is ρ_max = {s1['rho_max']:.4f} ≤ {s1['threshold']}, so "
        f"EP-1 is met (sizing action: {s1['action']}). (b) Percentile endpoints across the five seeds (between-run SD, ddof = 1, df = 4). These are "
        "simulation-error quantities, not sampling uncertainty. Source: corrected_interval_results.json.")
    s2 = V["figS2"]
    C["FigS2_informative_censoring"] = (
        f"**Figure S2. Informative censoring (S3 cohort, n = 453).** The {s2['died_no_recurrence']['n']} patients who died without recurrence (censored in the "
        f"deposited data) were older (median {s2['died_no_recurrence']['median_age']:g} vs {s2['alive_no_recurrence']['median_age']:g} years), less often "
        f"stage III ({s2['died_no_recurrence']['stage_III_pct']:g} % vs {s2['alive_no_recurrence']['stage_III_pct']:g} %) and less often treated with adjuvant "
        f"chemotherapy ({s2['died_no_recurrence']['chemo_pct']:g} % vs {s2['alive_no_recurrence']['chemo_pct']:g} %) than patients alive without recurrence. "
        "This is the profile of competing non-cancer mortality. The cause-specific hazard remains estimable, but it cannot be converted into an absolute "
        "recurrence probability; no competing-risk (Fine–Gray) or composite-endpoint analysis was performed. Source: docs/PHASE6A_ENDPOINT_RESOLUTION.md §2.3.")
    return C


# ---------------------------------------------------------------- output
def render(fig, fmt):
    buf = io.BytesIO()
    meta = {"Date": None} if fmt == "svg" else {"Software": None}
    fig.savefig(buf, format=fmt, dpi=DPI, metadata=meta, bbox_inches="tight", pad_inches=0.08)
    return buf.getvalue()


def write_new(path: Path, data: bytes):
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        os.write(fd, data)
    finally:
        os.close(fd)


def main():
    out = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / "reports" / "figures" / "final"
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"refusing: {out} exists and is not empty (no figure is ever overwritten)")
    V = assemble()
    out.mkdir(parents=True, exist_ok=True)
    written = {}
    for name, fn in FIGS:
        fig = fn(V)
        for fmt in ("svg", "png"):
            data = render(fig, fmt)
            write_new(out / f"{name}.{fmt}", data)
            written[f"{name}.{fmt}"] = hashlib.sha256(data).hexdigest()
        plt.close(fig)
    caps = captions(V)
    cap_md = "# Figure captions\n\nGenerated by `scripts/49_make_results_figures.py` from frozen outputs only (no model refitted). " \
             "Primary result: Figure 8. See `docs/SCIENTIFIC_RESULTS_AND_FIGURES.md`.\n\n" + \
             "\n\n".join(f"### {n}\n\n{caps[n]}" for n, _ in FIGS) + "\n"
    write_new(out / "FIGURE_CAPTIONS.md", cap_md.encode())
    outer_hash = hashlib.sha256(json.dumps(V["_outer"], sort_keys=True).encode()).hexdigest()
    vals = {k: v for k, v in V.items() if k != "_outer"}
    vals["fig8"]["outer_values_sha256"] = outer_hash
    write_new(out / "figure_values.json", (json.dumps(vals, indent=1, sort_keys=True, ensure_ascii=False) + "\n").encode())
    sources = {"script": {"path": "scripts/49_make_results_figures.py", "sha256": sha(Path(__file__).resolve())},
               "inputs": {rel: sha(ROOT / rel) for rel in list(SRC.values()) + CHECKPOINTS},
               "outputs": written,
               "parameters": {"dpi": DPI, "svg.hashsalt": plt.rcParams["svg.hashsalt"], "svg.fonttype": plt.rcParams["svg.fonttype"],
                              "font": plt.rcParams["font.family"], "palette": [C_PRIMARY, C_S7, C_FIVE],
                              "theta_expected": THETA_EXPECTED},
               "software": {"python": platform.python_version(), "matplotlib": matplotlib.__version__, "numpy": np.__version__},
               "models_fitted": 0}
    write_new(out / "figure_sources.json", (json.dumps(sources, indent=1, sort_keys=True) + "\n").encode())
    print(json.dumps({"output_dir": str(out), "files": len(written) + 3}, indent=1))


if __name__ == "__main__":
    main()

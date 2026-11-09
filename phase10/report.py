"""Build docs/PHASE10_RESULTS.md and figures from reports/phase10/*. No number is typed by hand.

Usage:  python -m phase10.report
"""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import precision_recall_curve, roc_curve  # noqa: E402

from phase10 import common as C  # noqa: E402

FIG = C.ROOT / "figures" / "phase10"
DOC = C.ROOT / "docs" / "PHASE10_RESULTS.md"
NAMES = {"enet_logistic": "Elastic-net logistic", "random_forest": "Random forest", "mlh1_only": "MLH1 only",
         "clinical_only": "Clinical only", "M1": "M1 clinical Cox", "M1_breslow": "M1 (Breslow, sksurv)",
         "M2_strict": "M2 clinical + 4 pathways", "M3_enet_cox": "M3 elastic-net Cox (genes)"}
COL = {"enet_logistic": "#2f5d8a", "random_forest": "#6b93bf", "mlh1_only": "#a7a7a7", "clinical_only": "#d0c7b8"}


def f(x, d=3):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"


def load(p):
    p = C.REPORTS / p
    return json.loads(p.read_text()) if p.exists() else None


def track_b(lines):
    s = load("mmr/summary.json")
    if not s:
        return
    lines += ["## Track B: MMR-associated classification in GSE39582 (exploratory, likely biased upward)", "",
              f"> {s['caveat']}", "",
              f"Cohort: n = {s['cohort']['n']} ({s['cohort']['dMMR']} dMMR, {s['cohort']['pMMR']} pMMR); "
              f"{s['cohort']['shared_genes']:,} shared genes; {s['design']['outer']}, inner {s['design']['inner']}. "
              f"Values are mean over repeats (SD across repeats = partition variability, not a confidence interval).", "",
              "| Model | AUROC | AUPRC (prevalence 0.144) | Brier | Scaled Brier | Calibration intercept | Calibration slope |",
              "|---|---|---|---|---|---|---|"]
    for m in ("enet_logistic", "random_forest", "mlh1_only", "clinical_only"):
        r = s["models"][m]
        lines.append(f"| {NAMES[m]} | {f(r['auroc_within_fold_pooled']['mean'])} (SD {f(r['auroc_within_fold_pooled']['sd_across_repeats'])}) | "
                     f"{f(r['auprc']['mean'])} (SD {f(r['auprc']['sd_across_repeats'])}) | {f(r['brier']['mean'])} | {f(r['scaled_brier']['mean'])} | "
                     f"{f(r['calibration_intercept']['mean'], 2)} | {f(r['calibration_slope']['mean'], 2)} |")
    lines += ["", "Paired ΔAUROC on identical positive–negative pairs (mean over repeats; SD; MC SE):", ""]
    for k, v in s["paired_delta_auroc"].items():
        lines.append(f"- {k}: {f(v['mean'], 4)} (SD {f(v['sd_across_repeats'], 4)}; MC SE {f(v['mc_se'], 4)})")
    lines += ["", f"Failed outer fits: {len(s['failed_folds'])}. Runtime {s['runtime_seconds']:.0f} s.", ""]
    sel = pd.read_csv(C.REPORTS / "mmr" / "enet_selection_stability.csv")
    if len(sel):
        lines += ["Elastic-net logistic: genes selected in the most outer fits (descriptive; not biomarkers):", "",
                  "| Gene | Selected in | Positive-sign fraction | Median coefficient |", "|---|---|---|---|"]
        for _, r in sel.head(15).iterrows():
            lines.append(f"| {r.gene} | {int(r.selected_in)}/{int(r.of_fits)} | {r.positive_sign_fraction:.2f} | {r.median_coef:+.3f} |")
        lines.append("")
    # figure: ROC + PR from repeat 0 OOF (one prediction per patient)
    oof = pd.read_csv(C.REPORTS / "mmr" / "oof_predictions.csv")
    o = oof[oof.repeat == 0].dropna(subset=["p"])
    fig, ax = plt.subplots(1, 2, figsize=(9, 4))
    for m in ("enet_logistic", "random_forest", "mlh1_only", "clinical_only"):
        g = o[o.model == m]
        fpr, tpr, _ = roc_curve(g.y, g.p)
        pr, rc, _ = precision_recall_curve(g.y, g.p)
        ax[0].plot(fpr, tpr, color=COL[m], label=NAMES[m])
        ax[1].plot(rc, pr, color=COL[m], label=NAMES[m])
    ax[0].plot([0, 1], [0, 1], ls=":", c="k", lw=0.8)
    ax[1].axhline(o[o.model == "mlh1_only"].y.mean(), ls=":", c="k", lw=0.8)
    ax[0].set(xlabel="False-positive rate", ylabel="True-positive rate", title="ROC (repeat 0, out-of-fold)")
    ax[1].set(xlabel="Recall", ylabel="Precision", title="Precision–recall (prevalence dotted)")
    ax[0].legend(fontsize=8, frameon=False)
    fig.suptitle("GSE39582 dMMR: exploratory, likely biased upward (ComBat protected MMR)", fontsize=9)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "F10_2_mmr_roc_pr.png", dpi=160)
    fig.savefig(FIG / "F10_2_mmr_roc_pr.svg")
    plt.close(fig)
    lines += ["Figure: `figures/phase10/F10_2_mmr_roc_pr.png` (repeat 0 only, so every patient appears once).", ""]


def track_b_tcga(lines):
    man = load("tcga/freeze_manifest.json")
    ev = load("tcga/evaluation.json")
    prep = load("tcga/prepare_report.json")
    lines += ["## Track B: independent-cohort, cross-platform evaluation in TCGA-COAD", ""]
    if prep:
        lines.append(f"Label-free preparation: {prep['primary_tumour_samples']} primary tumours; primary population "
                     f"(stage II–III, complete clinical) n = {prep['primary_population_n']}.")
    if man:
        lines.append(f"Freeze manifest: `reports/phase10/tcga/freeze_manifest.json` (SHA-256 `{C.sha256_file(C.REPORTS / 'tcga/freeze_manifest.json')}`), "
                     f"created {man['created_utc']}. Final elastic-net model: "
                     f"{len(man['models']['enet_logistic'].get('selected', {}))} non-zero genes, params {man['models']['enet_logistic']['best_params']}.")
    if not ev:
        lines += ["", "**Not yet evaluated.** TCGA MSI labels remain sealed until the freeze manifest is committed and the "
                  "owner approves unblinding. No TCGA performance number exists yet.", ""]
        return
    lines += ["", f"n = {ev['n']} ({ev['msi_h']} MSI-H). {ev['labels']}. Normalisation: {ev['normalisation']}. "
              f"Bootstrap: {ev['bootstrap']['resamples']} stratified patient resamples.", "",
              "| Model | AUROC (95% CI) | AUPRC (95% CI) | Brier | Calibration intercept | Calibration slope |", "|---|---|---|---|---|---|"]
    for m, r in ev["models"].items():
        lines.append(f"| {NAMES[m]} | {f(r['auroc'])} ({f(r['ci95']['auroc'][0])}–{f(r['ci95']['auroc'][1])}) | "
                     f"{f(r['auprc'])} ({f(r['ci95']['auprc'][0])}–{f(r['ci95']['auprc'][1])}) | {f(r['brier'])} | "
                     f"{f(r['calibration_intercept'], 2)} | {f(r['calibration_slope'], 2)} |")
    d = ev["enet_minus_mlh1_auroc"]
    lines += ["", f"Elastic-net − MLH1-only ΔAUROC: {f(d['point'], 3)} (paired bootstrap 95% CI {f(d['ci95_paired_bootstrap'][0], 3)} to {f(d['ci95_paired_bootstrap'][1], 3)}).", ""]


def track_a(lines):
    s = load("recurrence/summary.json")
    if not s:
        lines += ["## Track A: survival ML for recurrence", "", "**Not run.**", ""]
        return
    lines += ["## Track A: survival ML for cause-specific recurrence (exploratory)", "",
              f"Cohort n = {s['cohort']['n']}, {s['cohort']['events']} events. {s['design']['outer']}; inner {s['design']['inner']}; strict nesting. "
              "Endpoint caveats: " + "; ".join(s["endpoint_caveats"]) + ".", "",
              "| Model | Harrell C, within-fold pooled (mean over repeats) | SD across repeats | MC SE |", "|---|---|---|---|"]
    for m in ("M1", "M1_breslow", "M2_strict", "M3_enet_cox"):
        if m in s["models"]:
            r = s["models"][m]
            lines.append(f"| {NAMES[m]} | {f(r['mean'], 4)} | {f(r['sd_across_repeats'], 4)} | {f(r['mc_se'], 4)} |")
    lines += ["", "Paired ΔC on identical comparable pairs:", ""]
    for k, v in s["paired_delta_c"].items():
        lines.append(f"- {k}: {f(v['mean'], 4)} (SD {f(v['sd_across_repeats'], 4)}; MC SE {f(v['mc_se'], 4)}; {v['n_repeats']} repeats)")
    t = s["m3_tuning"]
    lines += ["", f"M3 tuning: l1_ratio chosen {t['l1_ratio_counts']}; α at a grid edge in {t['alpha_at_grid_edge']} fits; "
              f"genes selected per fit mean {f(t['genes_selected_per_fit']['mean'], 1)} (range {f(t['genes_selected_per_fit']['min'], 0)}–{f(t['genes_selected_per_fit']['max'], 0)}); "
              f"inner failures {t['inner_failures']}. M2 λ choices: {s['m2_lambda_counts']}.",
              f"N1 null (training-fold (time, event) permutation, M3): {s['n1_null']['replicates']} replicates, "
              f"mean C {f((s['n1_null']['summary'] or {}).get('mean'), 4)}; alarm = {s['n1_null']['alarm']} (rule: {s['n1_null']['alarm_rule']}).",
              f"Failed outer fits: {len(s['failed_outer_fits'])}.", "",
              f"Phase 9 reference (different estimand, not comparable): {s['phase9_reference_not_comparable']['estimand']}: "
              f"M1 {s['phase9_reference_not_comparable']['M1']:.6f}, M2 {s['phase9_reference_not_comparable']['M2_four_set']:.6f}.", ""]
    per = pd.read_csv(C.REPORTS / "recurrence" / "per_repeat.csv")
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    order = [m for m in ("M1", "M1_breslow", "M2_strict", "M3_enet_cox") if m in set(per.model)]
    for i, m in enumerate(order):
        v = per[per.model == m].c_within_fold_pooled
        ax.scatter(v, np.full(len(v), i) + np.linspace(-0.15, 0.15, len(v)), s=14, color="#2f5d8a", alpha=0.7)
        ax.plot([v.mean()] * 2, [i - 0.3, i + 0.3], color="k", lw=1.5)
    nul = C.REPORTS / "recurrence" / "n1_null_replicates.csv"
    if nul.exists() and nul.stat().st_size > 5:
        nv = pd.read_csv(nul).c_within_fold_pooled
        if len(nv):
            ax.axvspan(nv.min(), nv.max(), color="#d0c7b8", alpha=0.5, label="N1 null range (M3)")
            ax.legend(fontsize=7, frameon=False, loc="lower right")
    ax.axvline(0.5, ls=":", c="k", lw=0.8)
    ax.set_yticks(range(len(order)), [NAMES[m] for m in order], fontsize=8)
    ax.set_xlabel("Harrell C, within-fold pooled (one dot per repeat)")
    ax.set_title("Recurrence: 5-fold CV repeats (partition variability, not CIs)", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "F10_1_recurrence_models.png", dpi=160)
    fig.savefig(FIG / "F10_1_recurrence_models.svg")
    plt.close(fig)
    lines += ["Figure: `figures/phase10/F10_1_recurrence_models.png`.", ""]
    sel = C.REPORTS / "recurrence" / "m3_selection_stability.csv"
    if sel.exists() and sel.stat().st_size > 5:
        st = pd.read_csv(sel)
        if len(st):
            lines += ["M3 genes selected most often (descriptive; not biomarkers):", "", "| Gene | Selected in | Positive-sign fraction |", "|---|---|---|"]
            for _, r in st.head(10).iterrows():
                lines.append(f"| {r.gene} | {int(r.selected_in)}/{int(r.of_fits)} | {r.positive_sign_fraction:.2f} |")
            lines.append("")


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    lines = ["# Phase 10 results", "",
             "Generated by `python -m phase10.report` from `reports/phase10/`. Methods, cohorts, label definitions and limitations: "
             "`docs/PHASE10_METHODS_AND_MODEL_CARD.md`. All results are exploratory; none is a clinical claim.", ""]
    integ = C.phase9_integrity()
    lines += [f"Phase 9 integrity at report time: {integ['frozen_code_files'] - len(integ['frozen_code_mismatches'])}/{integ['frozen_code_files']} frozen files, "
              f"{integ['archive_files'] - len(integ['archive_mismatches'])}/{integ['archive_files']} archive files match.", ""]
    track_b(lines)
    track_b_tcga(lines)
    track_a(lines)
    DOC.write_text("\n".join(lines) + "\n")
    print(f"wrote {DOC.relative_to(C.ROOT)} and figures in {FIG.relative_to(C.ROOT)}")


if __name__ == "__main__":
    main()

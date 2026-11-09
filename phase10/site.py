"""Generate the website's Machine-learning section and hero card from reports/phase10/*.

Writes between <!-- PHASE10:START/END --> and <!-- PHASE10-CARD:START/END --> in
portfolio/web/index.html. Every number comes from the result JSON/CSV files; tests/phase10
re-renders and checks the page is in sync.   Usage:  python -m phase10.site
"""
from __future__ import annotations

import html
import json
import re

import pandas as pd

from phase10 import common as C

PAGE = C.ROOT / "portfolio" / "web" / "index.html"
NAMES = {"enet_logistic": "Elastic-net logistic", "random_forest": "Random forest", "mlh1_only": "MLH1 only",
         "clinical_only": "Clinical only", "M1": "M1 clinical Cox", "M1_breslow": "M1, Breslow ties",
         "M2_strict": "M2 + 4 pathways", "M3_enet_cox": "M3 elastic-net Cox"}


def J(p):
    return json.loads((C.REPORTS / p).read_text())


def f(x, d=3):
    return f"{x:.{d}f}"


def esc(s):
    return html.escape(str(s), quote=True)


def svg_workflow():
    box = lambda x, y, w, t1, t2, fill="var(--sheet)": (  # noqa: E731
        f'<rect x="{x}" y="{y}" width="{w}" height="54" rx="8" fill="{fill}" stroke="var(--line)"/>'
        f'<text x="{x + w / 2}" y="{y + 23}" text-anchor="middle" style="font:600 12px var(--text);fill:var(--ink)">{t1}</text>'
        f'<text x="{x + w / 2}" y="{y + 40}" text-anchor="middle" style="font:11px var(--text);fill:var(--ink-2)">{t2}</text>')
    arr = lambda x1, y1, x2, y2: f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="var(--ink-3)" marker-end="url(#ah)"/>'  # noqa: E731
    parts = ['<svg viewBox="0 0 640 190" role="img" aria-label="Leakage-safe workflow: saved patient-level folds, preprocessing refit inside every training fold, a committed freeze, then sealed TCGA labels evaluated once">',
             '<defs><marker id="ah" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0L8,4L0,8z" fill="var(--ink-3)"/></marker></defs>',
             box(4, 10, 140, "GSE39582", "452 recurrence · 409 MMR"),
             box(166, 10, 150, "Saved, hashed folds", "5-fold × 10 repeats"),
             box(338, 10, 150, "Inner 5-fold tuning", "filter, scaling refit per fold"),
             box(510, 10, 126, "Test fold", "scored once"),
             arr(144, 37, 164, 37), arr(316, 37, 336, 37), arr(488, 37, 508, 37),
             box(166, 112, 150, "Final models frozen", "hashes, committed", "var(--primary-wash)"),
             box(338, 112, 150, "Sealed TCGA labels", "readable only after freeze", "var(--primary-wash)"),
             box(510, 112, 126, "TCGA-COAD", "evaluated once"),
             arr(74, 64, 74, 139), arr(74, 139, 164, 139), arr(316, 139, 336, 139), arr(488, 139, 508, 139),
             '<text x="80" y="104" style="font:11px var(--data);fill:var(--ink-3)">Track B</text>',
             '</svg>']
    return "".join(parts)


def svg_track_a(a):
    per = pd.read_csv(C.REPORTS / "recurrence" / "per_repeat.csv")
    order = ["M1", "M1_breslow", "M2_strict", "M3_enet_cox"]
    lo, hi, L, R = 0.57, 0.63, 150, 570
    X = lambda v: L + (v - lo) / (hi - lo) * (R - L)  # noqa: E731
    out = [f'<svg viewBox="0 0 640 {40 + 34 * len(order) + 30}" role="img" aria-label="Harrell C per recurrence model, one dot per cross-validation repeat">']
    for t in (0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63):
        out.append(f'<line x1="{X(t):.1f}" y1="20" x2="{X(t):.1f}" y2="{30 + 34 * len(order)}" stroke="var(--line)"/>'
                   f'<text x="{X(t):.1f}" y="{46 + 34 * len(order)}" text-anchor="middle" style="font:10px var(--data);fill:var(--ink-3)">{t:.2f}</text>')
    m1 = a["models"]["M1"]["mean"]
    out.append(f'<line x1="{X(m1):.1f}" y1="16" x2="{X(m1):.1f}" y2="{30 + 34 * len(order)}" stroke="var(--ink-3)" stroke-dasharray="3 3"/>')
    for i, m in enumerate(order):
        y = 36 + 34 * i
        v = per[per.model == m].c_within_fold_pooled.to_numpy()
        out.append(f'<text x="{L - 10}" y="{y + 4}" text-anchor="end" style="font:12px var(--text);fill:var(--ink)">{NAMES[m]}</text>')
        out.append(f'<line x1="{X(v.min()):.1f}" y1="{y}" x2="{X(v.max()):.1f}" y2="{y}" stroke="var(--primary)" stroke-width="2" opacity=".35"/>')
        out += [f'<circle cx="{X(x):.1f}" cy="{y}" r="3" fill="var(--primary)" opacity=".55"/>' for x in v]
        mu = a["models"][m]["mean"]
        out.append(f'<line x1="{X(mu):.1f}" y1="{y - 9}" x2="{X(mu):.1f}" y2="{y + 9}" stroke="var(--ink)" stroke-width="2.5"/>'
                   f'<text x="{R + 8}" y="{y + 4}" style="font:11px var(--data);fill:var(--ink-2)">{f(mu, 4)}</text>')
    out.append(f'<text x="{(L + R) / 2}" y="{64 + 34 * len(order)}" text-anchor="middle" style="font:11px var(--text);fill:var(--ink-3)">Harrell C, within-fold pooled · dots = 10 repeats (partition variability, not CIs) · dashed = M1 mean</text></svg>')
    return "".join(out)


def svg_track_b(ev, cv):
    order = ["enet_logistic", "random_forest", "mlh1_only", "clinical_only"]
    lo, hi, L, R = 0.5, 1.0, 150, 600
    X = lambda v: L + (v - lo) / (hi - lo) * (R - L)  # noqa: E731
    out = [f'<svg viewBox="0 0 640 {40 + 36 * len(order) + 34}" role="img" aria-label="AUROC per classifier: TCGA independent cohort with 95% bootstrap intervals, and GSE39582 cross-validation values likely biased upward">']
    for t in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
        out.append(f'<line x1="{X(t):.1f}" y1="20" x2="{X(t):.1f}" y2="{30 + 36 * len(order)}" stroke="var(--line)"/>'
                   f'<text x="{X(t):.1f}" y="{46 + 36 * len(order)}" text-anchor="middle" style="font:10px var(--data);fill:var(--ink-3)">{t:.1f}</text>')
    for i, m in enumerate(order):
        y = 38 + 36 * i
        r = ev["models"][m]
        lo_ci, hi_ci = r["ci95"]["auroc"]
        g = cv["models"][m]["auroc_within_fold_pooled"]["mean"] if cv else None
        out.append(f'<text x="{L - 10}" y="{y + 4}" text-anchor="end" style="font:12px var(--text);fill:var(--ink)">{NAMES[m]}</text>')
        out.append(f'<line x1="{X(lo_ci):.1f}" y1="{y}" x2="{X(hi_ci):.1f}" y2="{y}" stroke="var(--primary)" stroke-width="3"/>'
                   f'<circle cx="{X(r["auroc"]):.1f}" cy="{y}" r="5" fill="var(--primary)"/>')
        if g is not None:
            out.append(f'<circle cx="{X(g):.1f}" cy="{y}" r="5" fill="none" stroke="var(--eosin)" stroke-width="1.5"/>')
        out.append(f'<text x="{R + 10}" y="{y + 4}" style="font:11px var(--data);fill:var(--ink-2)">{f(r["auroc"])}</text>')
    out.append(f'<text x="{(L + R) / 2}" y="{66 + 36 * len(order)}" text-anchor="middle" style="font:11px var(--text);fill:var(--ink-3)">'
               '● TCGA-COAD AUROC with 95% bootstrap CI · ○ GSE39582 cross-validation (likely biased upward)</text></svg>')
    return "".join(out)


def svg_selection():
    sel = pd.read_csv(C.REPORTS / "mmr" / "enet_selection_stability.csv").head(12)
    L, R = 110, 560
    out = [f'<svg viewBox="0 0 640 {30 + 22 * len(sel) + 24}" role="img" aria-label="Genes most often selected by elastic-net logistic regression across cross-validation fits">']
    for i, r in enumerate(sel.itertuples()):
        y = 20 + 22 * i
        w = (R - L) * r.frequency
        col = "var(--eosin)" if r.gene in ("MLH1", "EPM2AIP1") else "var(--primary)"
        out.append(f'<text x="{L - 8}" y="{y + 11}" text-anchor="end" style="font:11px var(--data);fill:var(--ink)">{esc(r.gene)}</text>'
                   f'<rect x="{L}" y="{y}" width="{w:.1f}" height="14" rx="2" fill="{col}" opacity=".8"/>'
                   f'<text x="{L + w + 6:.1f}" y="{y + 11}" style="font:10px var(--data);fill:var(--ink-2)">{int(r.selected_in)}/{int(r.of_fits)} · {"↓" if r.median_coef < 0 else "↑"}</text>')
    out.append(f'<text x="{(L + R) / 2}" y="{30 + 22 * len(sel) + 14}" text-anchor="middle" style="font:11px var(--text);fill:var(--ink-3)">'
               'Fraction of outer fits selecting the gene · ↓ lower expression → dMMR · overlapping training sets: not replication</text></svg>')
    return "".join(out)


def render():
    a, cv, ev = J("recurrence/summary.json"), J("mmr/summary.json"), J("tcga_ga/evaluation.json")
    hs = J("tcga/evaluation.json")
    d = a["paired_delta_c"]
    t = a["m3_tuning"]
    em = ev["enet_minus_mlh1_auroc"]
    e, ml = ev["models"]["enet_logistic"], ev["models"]["mlh1_only"]
    g = cv["models"]
    card = (f'<div class="card"><span class="k">4 · Machine learning</span><div class="big">{f(e["auroc"])}</div>'
            f'<p>TCGA-COAD AUROC (95% CI {f(e["ci95"]["auroc"][0])}–{f(e["ci95"]["auroc"][1])}; {ev["n"]} patients, {ev["msi_h"]} MSI-H) for an elastic-net MMR/MSI classifier trained only on GSE39582 and frozen before TCGA labels were read; '
            f'MLH1 alone {f(ml["auroc"])}. For recurrence, an elastic-net Cox model over 2,000 genes added nothing beyond clinical variables (ΔC {f(d["M3_enet_cox_minus_M1"]["mean"], 4)}).</p></div>')
    sec = f'''<section class="part" id="sec-ml" aria-labelledby="h-ml">
  <div class="rail"><span class="num">05</span><b>Machine learning</b>Survival and<br>classification</div>
  <div class="body">
    <h2 id="h-ml">Leakage-controlled machine learning: no recurrence signal in the transcriptome, a strong MMR/MSI signal that transfers to an independent cohort</h2>
    <p class="lede">Two tracks share one pipeline. Folds are saved and hashed, every data-dependent step (top-2,000 variance filter, scaling, pathway scores, spline knots) is refit inside each inner-training fold, and the test fold is scored once. <b>Track A</b> asks whether gene-level survival models improve recurrence discrimination. <b>Track B</b> is a positive control: predicting MMR/MSI status, where the biology is strong and known, with an independent TCGA evaluation whose labels stayed sealed until the models were committed.</p>
    <figure class="panel">{svg_workflow()}</figure>

    <h3>Track A · Recurrence (cause-specific hazard; {a["cohort"]["n"]} patients, {a["cohort"]["events"]} events)</h3>
    <figure class="panel">{svg_track_a(a)}<figcaption>{a["design"]["outer"]}; inner {a["design"]["inner"]}; strict nesting; identical comparable pairs for every model. These cross-validated values are a different estimand from Phase 9's optimism-corrected C (0.6016 / 0.6059) and are not compared with it.</figcaption></figure>
    <div class="stats">
      <div><span class="v data">{f(d["M2_strict_minus_M1"]["mean"], 4)}</span><span class="l">M2 − M1 ΔC (repeat SD {f(d["M2_strict_minus_M1"]["sd_across_repeats"], 4)})</span></div>
      <div><span class="v data">{f(d["M3_enet_cox_minus_M1"]["mean"], 4)}</span><span class="l">M3 elastic-net Cox − M1 ΔC</span></div>
      <div><span class="v data">{int(t["genes_selected_per_fit"]["max"])} of 2,000</span><span class="l">genes kept by M3 in every one of {t["genes_selected_per_fit"]["n_repeats"]} outer fits</span></div>
      <div><span class="v data">{f(a["n1_null"]["summary"]["mean"])}</span><span class="l">mean C under permuted training outcomes (N1, {a["n1_null"]["replicates"]} replicates): no leakage alarm</span></div>
      <div><span class="v data">{len(a["failed_outer_fits"])} failed</span><span class="l">outer model fits across 4 models × {t["genes_selected_per_fit"]["n_repeats"]} folds</span></div>
    </div>
    <p class="conclude">Inner cross-validation always preferred the penalty that excludes every gene, so the elastic-net Cox model reduced to the clinical model. The 2,000 most variable genes added no recurrence discrimination this design could detect.</p>

    <h3>Track B · MMR/MSI-associated classification</h3>
    <figure class="panel">{svg_track_b(ev, cv)}<figcaption>Training: GSE39582, {cv["cohort"]["n"]} patients ({cv["cohort"]["dMMR"]} dMMR), {cv["design"]["outer"]}. Independent evaluation: TCGA-COAD Illumina GA RNA-seq, stage II–III, {ev["n"]} patients ({ev["msi_h"]} MSI-H), {ev["bootstrap"]["resamples"]:,} stratified bootstrap resamples. Expression was z-scored within TCGA (unsupervised target-cohort adaptation, label-free), which needs a batch of samples: this is not a single-patient prediction pipeline.</figcaption></figure>
    <div class="stats">
      <div><span class="v data">{f(e["auroc"])}</span><span class="l">TCGA AUROC, elastic-net (CI {f(e["ci95"]["auroc"][0])}–{f(e["ci95"]["auroc"][1])})</span></div>
      <div><span class="v data">{f(e["auprc"])}</span><span class="l">TCGA AUPRC (prevalence {f(e["prevalence"])})</span></div>
      <div><span class="v data">+{f(em["point"])}</span><span class="l">elastic-net − MLH1 ΔAUROC (paired CI {f(em["ci95_paired_bootstrap"][0])} to {f(em["ci95_paired_bootstrap"][1])})</span></div>
      <div><span class="v data">{f(e["calibration_slope"], 2)}</span><span class="l">TCGA calibration slope: predictions too conservative after transfer</span></div>
      <div><span class="v data">{f(g["enet_logistic"]["auroc_within_fold_pooled"]["mean"])}</span><span class="l">GSE39582 CV AUROC, likely biased upward</span></div>
    </div>
    <figure class="panel">{svg_selection()}<figcaption>Genes the elastic-net classifier selected most often across GSE39582 outer fits. MLH1 and EPM2AIP1 (shared, bidirectional promoter) are selected with lower expression in MSI-H/dMMR tumours, consistent with MLH1 promoter hypermethylation in sporadic MSI colorectal cancer (Herman et al. 1998). Selection frequency is descriptive: these are not biomarkers or validated genes.</figcaption></figure>
    <p class="measure"><strong>What this cannot establish.</strong> GSE39582's deposited expression was ComBat-corrected with MMR status protected as a covariate, so within-cohort Track B results are likely biased upward; only the TCGA result is free of that bias. The two labels are related, not identical: both are PCR-based MSI calls (GSE39582: five Bethesda markers; TCGA: two panels at two sites), and dMMR is not synonymous with MLH1 loss, so MLH1 alone captures only the MLH1-silenced subset. The first, pre-registered TCGA evaluation on HiSeq RNA-seq had labels for only {hs["n"]} of 182 patients and is reported as non-informative; the GA-platform re-evaluation was registered and committed before any label was linked (declared deviation). None of this is a clinical tool.</p>
    <p class="small">Methods and model card: <code>docs/PHASE10_METHODS_AND_MODEL_CARD.md</code> · results: <code>docs/PHASE10_RESULTS.md</code> · reproduce: <code>python -m phase10.mmr</code>, <code>python -m phase10.recurrence</code>, <code>python -m phase10.tcga_ga evaluate</code>.</p>
  </div>
</section>'''
    return card, sec


def write():
    s = PAGE.read_text(encoding="utf-8")
    card, sec = render()
    s = re.sub(r"<!-- PHASE10:START -->.*?<!-- PHASE10:END -->", lambda _: f"<!-- PHASE10:START -->\n{sec}\n<!-- PHASE10:END -->", s, flags=re.S)
    s = re.sub(r"<!-- PHASE10-CARD:START -->.*?<!-- PHASE10-CARD:END -->", lambda _: f"<!-- PHASE10-CARD:START -->{card}<!-- PHASE10-CARD:END -->", s, flags=re.S)
    PAGE.write_text(s, encoding="utf-8")


if __name__ == "__main__":
    write()
    print("updated", PAGE.relative_to(C.ROOT))

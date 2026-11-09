"""Write audit results: JSON (machine-readable), CSV tables and a Markdown report."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items() if not isinstance(v, pd.DataFrame)}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, float) and np.isnan(o):
        return None
    return o


def write(res: dict, out, title: str = "Pathway-score audit") -> Path:
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "audit.json").write_text(json.dumps(_clean(res), indent=2))
    res["scores"].to_csv(out / "scores.csv", index_label="sample")
    if "composition" in res:
        res["composition"].to_csv(out / "composition_estimates.csv", index_label="sample")
    L = [f"# {title}", "", f"Samples: {res['n_samples']} · genes: {res['n_genes']}", ""]
    q = res["qc"]
    L += ["## 1. Data-quality gate (positive-control co-expression)", "",
          f"**{'PASS' if q['passed'] else 'FAIL'}**: median r = {q['median_r']:.2f} over {q['n_pairs']} control pairs (threshold {q['threshold']}).", ""]
    if not q["passed"]:
        L += ["> Gene-level results below are **not interpretable**: control genes that co-express in virtually every bulk tumour dataset "
              "do not co-express here, which indicates a processing or probe-labelling defect.", ""]
    L += ["## 2. Verdict per gene set", "", "| Gene set | Genes present | Verdict |", "|---|---|---|"]
    for name, s in res["gene_sets"].items():
        c = s["coverage"]
        L.append(f"| {name} | {c['n_present']}/{c['n_genes']} | {'; '.join(s['verdict'])} |")
    L.append("")
    for name, s in res["gene_sets"].items():
        if "top_genes" not in s:
            continue
        L += [f"### {name}", ""]
        L.append("Top genes (item–rest r): " + ", ".join(f"{g['gene']} {g['item_rest_r']:.2f}" for g in s["top_genes"][:8]))
        if s.get("modules"):
            L.append("Modules: " + "; ".join(f"M{m['module']} ({m['n_genes']} genes, Jaccard {m['mean_jaccard']:.2f}, {m['stability']})" for m in s["modules"]))
        if "composition_rho" in s:
            L.append("Composition ρ: " + ", ".join(f"{k} {v:+.2f}" for k, v in s["composition_rho"].items()))
        if s.get("gene_disjoint_fibroblast_test"):
            d = s["gene_disjoint_fibroblast_test"]
            L.append(f"Gene-disjoint fibroblast test: ρ {d['rho'][0]:.2f} ({d['rho'][1]:.2f}–{d['rho'][2]:.2f}); {d['n_shared_removed']} shared genes removed")
        if s.get("genotype_effects"):
            L += ["", "| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |", "|---|---|---|---|---|"]
            for lab, e in s["genotype_effects"].items():
                u = e.get("unadjusted")
                a = e.get("adjusted")
                fmt = (lambda t: f"{t[0]:+.2f} ({t[1]:+.2f} to {t[2]:+.2f})" if t else "—")
                L.append(f"| {lab} | {e['n_pos']}/{e['n_neg']} | {fmt(u)} | {fmt(a)} | {e['verdict']} |")
        L.append("")
    if res.get("pairwise"):
        L += ["## 3. Relationships between scores (partial Spearman given stromal + immune estimates)", "",
              "| Pair | Spearman | Partial | Random-set null 95 % | Exceeds null |", "|---|---|---|---|---|"]
        for p in res["pairwise"]:
            L.append(f"| {p['pair']} | {p['spearman']:+.2f} | {p['partial_given_composition']:+.2f} | "
                     f"{p['null_95_for_first_set'][0]:+.2f} to {p['null_95_for_first_set'][1]:+.2f} | {'yes' if p['exceeds_null'] else 'no'} |")
        L.append("")
    if res.get("survival_screen"):
        s = res["survival_screen"]
        L += ["## 4. Optional survival screen", "", f"Apparent ΔC {s['apparent_delta_C']:+.4f}; optimism {s['optimism_delta']:+.4f}; "
              f"**corrected ΔC {s['corrected_delta_C']:+.4f}** (n {s['n']}, events {s['events']}, B {s['B']}; {s['note']}).", ""]
    L += ["## Interpretation limits", "",
          "- Scores are averaged bulk mRNA; they do not measure pathway activity, protein function or cell-type-specific expression.",
          "- Composition values are computational estimates (ESTIMATE, MCP-counter) from the same expression data, not measured fractions.",
          "- Adjusted genotype effects may be over-adjusted if composition lies on the biological path.", ""]
    (out / "REPORT.md").write_text("\n".join(L))
    return out / "REPORT.md"

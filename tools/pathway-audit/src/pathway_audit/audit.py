"""Run the full audit for a set of gene sets and return a structured result with plain-language verdicts.

Decision rules (documented in docs/METHOD.md):
* QC: positive-control median r >= 0.5, else gene-level interpretation is flagged unsafe.
* Coherence: mean inter-gene r; < 0.05 flagged "incoherent" (score is an average of weakly related genes).
* Gene concentration: min leave-one-gene-out r < 0.95 flagged "driven by few genes".
* Composition: max |Spearman rho| of the score with any composition estimate:
  >= 0.7 "composition-dominated", 0.4–0.7 "composition-associated", < 0.4 "largely composition-independent";
  for "composition-dominated" a gene-disjoint marker test must also give rho >= 0.7 (otherwise "composition-dominated
  (shared genes)").
* Genotype effects: see association.genotype_effect.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from . import association as A
from . import composition as C
from .genes import gene_contributions, modules
from .qc import coexpression_qc
from .scoring import score_gene_set, zscore_rows


def audit(expr: pd.DataFrame, gene_sets: dict, annotations: pd.DataFrame | None = None, binary: dict | None = None,
          adjust_for: list | None = None, reference_cache=None, survival: dict | None = None,
          composition: pd.DataFrame | None = None, fibroblast_markers: list | None = None,
          n_boot: int = 1000, n_null: int = 500, module_boot: int = 200, seed: int = 0) -> dict:
    """expr: genes × samples (log scale). binary: {name: (column, positive_value, negative_value)}.
    adjust_for: annotation columns (binary 'positive=value' strings or numeric) added to composition in adjusted models.
    survival: {"time": col, "event": col, "clinical": [numeric cols]} in `annotations` (optional ΔC screen)."""
    rng = np.random.default_rng(seed)
    res = {"n_samples": int(expr.shape[1]), "n_genes": int(expr.shape[0])}
    res["qc"] = coexpression_qc(expr)
    z = zscore_rows(expr)
    comp = None
    fib = list(fibroblast_markers or [])
    if composition is not None:
        comp = composition.reindex(expr.columns).astype(float)
        res["composition_summary"] = {"estimators": list(comp.columns), "source": "user-supplied"}
    elif reference_cache is not None:
        comp = C.estimate_scores(expr, reference_cache)[["StromalScore", "ImmuneScore"]].join(C.mcp_scores(expr, reference_cache))
        res["composition_summary"] = {"estimators": list(comp.columns), "n_common_genes_estimate": comp.attrs.get("n_common_genes"),
                                      "source": "ESTIMATE + MCP-counter"}
        fib = fib or C.mcp_markers(reference_cache).get("Fibroblasts", [])
    ann = annotations.reindex(expr.columns) if annotations is not None else None
    adj = None
    adj_cols = [c for c in ("StromalScore", "ImmuneScore") if comp is not None and c in comp.columns] or (list(comp.columns) if comp is not None else [])
    if comp is not None:
        adj = comp[adj_cols].copy()
    if adjust_for and ann is not None:
        extra = {}
        for spec in adjust_for:
            col, _, pos = spec.partition("=")
            v = ann[col]
            extra[spec] = np.where(v == "", np.nan, (v == pos).astype(float)) if pos else pd.to_numeric(v.replace("", np.nan), errors="coerce").to_numpy()
        extra = pd.DataFrame(extra, index=expr.columns)
        adj = extra if adj is None else adj.join(extra)
    scores, sets_out = {}, {}
    for name, genes in gene_sets.items():
        s, cov = score_gene_set(z, genes)
        out = {"coverage": cov}
        if s is None:
            out["verdict"] = ["not scored: coverage rule not met"]
            sets_out[name] = out
            continue
        scores[name] = s
        gc = gene_contributions(z, genes)
        out["mean_inter_gene_r"] = gc.attrs["mean_inter_gene_r"]
        out["min_logo_r"] = float(gc.logo_r_with_full.min())
        out["top_genes"] = gc.head(10).to_dict("records")
        mods = modules(z, genes, min_size=10 if len(genes) >= 60 else 5, n_boot=module_boot, seed=seed)
        out["modules"] = mods.to_dict("records")
        verdict = []
        verdict.append("incoherent (mean inter-gene r < 0.05)" if out["mean_inter_gene_r"] < 0.05 else
                       f"coherent (mean inter-gene r {out['mean_inter_gene_r']:.2f})")
        if out["min_logo_r"] < 0.95:
            verdict.append(f"driven by few genes (min leave-one-out r {out['min_logo_r']:.2f})")
        if comp is not None:
            rhos = {c: float(stats.spearmanr(s, comp[c])[0]) for c in comp.columns}
            out["composition_rho"] = rhos
            top = max(rhos, key=lambda k: abs(rhos[k]))
            mx = abs(rhos[top])
            if mx >= 0.7:
                dj = A.gene_disjoint_marker_rho(z, expr, genes, fib, rng, n_boot) if fib else None
                out["gene_disjoint_fibroblast_test"] = dj
                ok = dj is not None and dj["rho"][0] >= 0.7
                verdict.append(f"composition-dominated (|rho| {mx:.2f} with {top}" + ("; holds without shared genes)" if ok else "; may reflect shared genes)"))
            elif mx >= 0.4:
                verdict.append(f"composition-associated (|rho| {mx:.2f} with {top})")
            else:
                verdict.append(f"largely composition-independent (max |rho| {mx:.2f})")
        if binary and ann is not None:
            eff = {}
            for lab, (col, pos, neg) in binary.items():
                label = ann[col].map({pos: 1.0, neg: 0.0})
                eff[lab] = A.genotype_effect(s, label, adj, rng, n_boot)
            out["genotype_effects"] = eff
        out["verdict"] = verdict
        sets_out[name] = out
    res["gene_sets"] = sets_out
    # pairwise structure with random-set null
    if comp is not None and len(scores) >= 2:
        names = list(scores)
        W = comp[adj_cols].to_numpy()
        pairs = []
        for i, a in enumerate(names):
            null = A.random_set_null(z, scores[a], int(np.median([len(g) for g in gene_sets.values()])),
                                     lambda x, y: A.partial_spearman(x, y, W), rng, n_null)
            lo, hi = np.percentile(null, [2.5, 97.5])
            for b in names[i + 1:]:
                pr = A.partial_spearman(scores[a].to_numpy(), scores[b].to_numpy(), W)
                pairs.append({"pair": f"{a} | {b}", "spearman": float(stats.spearmanr(scores[a], scores[b])[0]), "partial_given_composition": pr,
                              "null_95_for_first_set": [float(lo), float(hi)], "exceeds_null": bool(pr > hi or pr < lo)})
        res["pairwise"] = pairs
    if survival and ann is not None and scores:
        from .predictive import delta_c
        d = ann[[survival["time"], survival["event"]] + survival["clinical"]].apply(pd.to_numeric, errors="coerce")
        ok = d.notna().all(axis=1)
        sc = pd.DataFrame(scores).loc[ok.values]
        res["survival_screen"] = delta_c(d.loc[ok, survival["clinical"]], sc, d.loc[ok, survival["time"]], d.loc[ok, survival["event"]], seed=seed)
    res["scores"] = pd.DataFrame(scores)
    if comp is not None:
        res["composition"] = comp
    return res

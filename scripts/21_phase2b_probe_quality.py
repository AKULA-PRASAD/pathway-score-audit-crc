"""Phase 2B part 2: probe quality, neutral aggregation comparison,
transportability and the mapping reconciliation.

Read-only. Outcome-blind. No normalization, no batch correction.

Usage:
    .venv/bin/python scripts/21_phase2b_probe_quality.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                   # noqa: E402
from src.annotation import find_table_header            # noqa: E402
from src.cache import load_probe_matrix                 # noqa: E402
from src.gene_mapping import build_mapping              # noqa: E402

OUT = paths.ROOT / "reports" / "phase2b"


def section(t: str) -> None:
    print("\n" + "=" * 74); print(t); print("=" * 74)


def affy_suffix(probe: str) -> str:
    """Affymetrix probe-set suffix encodes hybridisation specificity.

    _at    -- probes match a single known transcript
    _s_at  -- all probes cross-hybridise to a documented transcript FAMILY
    _x_at  -- some probes cross-hybridise unpredictably (least specific)
    _a_at  -- recognises alternative transcripts of one gene
    """
    for suf in ("_x_at", "_s_at", "_a_at", "_i_at", "_g_at", "_r_at", "_b_at"):
        if probe.endswith(suf):
            return suf
    return "_at" if probe.endswith("_at") else "other"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    findings: dict = {}

    expr = load_probe_matrix()
    mapping = build_mapping()
    um = mapping.unique_map
    per_gene = um.groupby(um.values).size()
    multi = sorted(per_gene[per_gene > 1].index)
    probes_multi = um[um.isin(multi)].index

    vals = expr.to_numpy()
    q1, q3 = np.percentile(vals, [25, 75], axis=1)
    pstat = pd.DataFrame({"mean": vals.mean(axis=1),
                          "sd": vals.std(axis=1, ddof=0),
                          "var": vals.var(axis=1),
                          "iqr": q3 - q1}, index=expr.index)

    rep = pd.read_csv(OUT / "representative_probe_by_selection_set.csv",
                      index_col=0)
    sel585 = set(rep.loc[multi, "all_585"])
    pairs = pd.read_csv(paths.ROOT / "reports" / "phase2"
                        / "within_gene_probe_correlations.csv")
    wg = pd.read_csv(paths.ROOT / "reports" / "phase2"
                     / "within_gene_correlation_summary.csv",
                     index_col=0)["median_spearman_rho"]

    # ---------------------------------------------------------------
    section("3. PROBE QUALITY: SELECTED vs NOT-SELECTED (multi-probe genes)")
    # ---------------------------------------------------------------
    pm = pstat.loc[probes_multi].copy()
    pm["selected"] = pm.index.isin(sel585)
    pm["suffix"] = [affy_suffix(p) for p in pm.index]

    # mean within-gene rho of each probe against its siblings
    long = pd.concat([
        pairs[["probe_a", "spearman_rho"]].rename(columns={"probe_a": "probe"}),
        pairs[["probe_b", "spearman_rho"]].rename(columns={"probe_b": "probe"}),
    ])
    pm["mean_sibling_rho"] = long.groupby("probe")["spearman_rho"].mean(
    ).reindex(pm.index)

    print("  continuous measures (median [IQR]):")
    rows = []
    for col in ["mean", "iqr", "var", "mean_sibling_rho"]:
        s = pm[pm.selected][col]; n = pm[~pm.selected][col]
        rows.append({
            "measure": col,
            "selected_median": round(float(s.median()), 4),
            "not_selected_median": round(float(n.median()), 4),
            "selected_p25_p75": f"{s.quantile(.25):.3f}-{s.quantile(.75):.3f}",
            "not_sel_p25_p75": f"{n.quantile(.25):.3f}-{n.quantile(.75):.3f}",
        })
    qt = pd.DataFrame(rows)
    print(qt.to_string(index=False))
    qt.to_csv(OUT / "probe_quality_comparison.csv", index=False)

    print("\n  low-expression frequency:")
    for thr in (4.0, 5.0, 6.0):
        a = float((pm[pm.selected]["mean"] < thr).mean() * 100)
        b = float((pm[~pm.selected]["mean"] < thr).mean() * 100)
        print(f"    mean log2 < {thr:.1f}:  selected {a:5.1f}%   "
              f"not-selected {b:5.1f}%")

    print("\n  annotation specificity (Affymetrix probe-set suffix):")
    ct = pd.crosstab(pm["suffix"], pm["selected"], normalize="columns") * 100
    ct.columns = ["not_selected_%", "selected_%"]
    ct["n_total"] = pm["suffix"].value_counts()
    print(ct.round(2).sort_values("n_total", ascending=False).to_string())
    print("\n    _at    = probes match one known transcript (most specific)")
    print("    _s_at  = cross-hybridise to a documented transcript family")
    print("    _x_at  = some probes cross-hybridise unpredictably (least specific)")
    ct.to_csv(OUT / "probe_suffix_distribution.csv")

    x_sel = float((pm[pm.selected]["suffix"] == "_x_at").mean() * 100)
    x_not = float((pm[~pm.selected]["suffix"] == "_x_at").mean() * 100)
    print(f"\n    _x_at share: selected {x_sel:.2f}% vs not-selected {x_not:.2f}%")
    print("    NOTE: variability is NOT biological quality. A cross-hybridising")
    print("    probe can be highly variable for purely technical reasons, and")
    print("    an IQR criterion cannot distinguish that from real biology.")
    findings["probe_quality"] = {
        "rows": qt.to_dict(orient="records"),
        "x_at_selected_pct": x_sel, "x_at_not_selected_pct": x_not,
        "low_expr_selected_pct_lt5": float(
            (pm[pm.selected]["mean"] < 5).mean() * 100),
        "low_expr_not_selected_pct_lt5": float(
            (pm[~pm.selected]["mean"] < 5).mean() * 100),
    }

    # ---------------------------------------------------------------
    section("4. PROBE AGREEMENT (DESCRIPTIVE) -- Cronbach's alpha")
    # ---------------------------------------------------------------
    print("  SCOPE NOTE (corrected). Cronbach's alpha is an INTERNAL-CONSISTENCY")
    print("  statistic. It summarises how much a set of probes co-vary across")
    print("  samples. It does NOT measure accuracy, does NOT establish that an")
    print("  average is closer to true transcript abundance than a single probe,")
    print("  and is NOT a validated rule for choosing an aggregation strategy.")
    print("  It is reported here as DESCRIPTIVE evidence of probe agreement only.\n")
    print("  Assumptions it requires, and how they fare for microarray probes:")
    print("    (a) items measure ONE latent construct      -- unverifiable here;")
    print("        probes may target different transcripts/UTRs of a locus")
    print("    (b) measurement errors are UNCORRELATED     -- VIOLATED: probes")
    print("        share array, sample prep, hybridisation, scanning, and the")
    print("        depositors' joint RMA + ComBat, all of which induce")
    print("        correlated error and INFLATE alpha")
    print("    (c) items are at least tau-equivalent       -- not established")
    print("    (d) a reference measurement exists          -- ABSENT: no qPCR or")
    print("        RNA-seq on these samples, so no accuracy claim is testable\n")

    alphas, ks, genes_eval = [], [], []
    sub = expr.loc[probes_multi]
    g_of = um.loc[sub.index]
    for gene, grp in g_of.groupby(g_of.values):
        ids = list(grp.index)
        k = len(ids)
        X = sub.loc[ids].to_numpy()
        item_var = X.var(axis=1, ddof=1).sum()
        tot_var = X.sum(axis=0).var(ddof=1)
        if tot_var <= 0:
            continue
        alphas.append((k / (k - 1)) * (1 - item_var / tot_var))
        ks.append(k); genes_eval.append(gene)

    al = pd.DataFrame({"gene_symbol": genes_eval, "alpha": alphas,
                       "k_probes": ks})
    # Mean inter-item correlation implied by alpha and k (Spearman-Brown
    # rearranged). Reported ONLY to expose alpha's mechanical k-dependence --
    # it is an algebraic restatement, not an independent measurement.
    a_, k_ = al["alpha"].to_numpy(), al["k_probes"].to_numpy()
    al["implied_mean_inter_probe_r"] = a_ / (k_ - a_ * (k_ - 1))
    al["within_gene_rho"] = wg.reindex(al.gene_symbol).values
    al.to_csv(OUT / "gene_probe_agreement_alpha.csv", index=False)

    print(f"  genes evaluated: {len(al)}")
    print(al["alpha"].describe(percentiles=[.1, .25, .5, .75, .9]).round(3).to_string())
    print("\n  distribution of alpha (DESCRIPTIVE BANDS, not decision thresholds;")
    print("  the conventional 0.7 cut comes from psychometrics and has no")
    print("  established biological meaning for microarray probe sets):")
    for lo, hi, lab in [(-np.inf, 0.0, "alpha <= 0   probes unrelated/anti-correlated"),
                        (0.0, 0.5, "0.0 - 0.5    low agreement"),
                        (0.5, 0.7, "0.5 - 0.7    moderate agreement"),
                        (0.7, 0.8, "0.7 - 0.8    high agreement"),
                        (0.8, 1.01, "alpha >= 0.8 very high agreement")]:
        n = int(((al.alpha > lo) & (al.alpha <= hi)).sum())
        print(f"    {lab:48s} {n:6d}  ({n/len(al)*100:5.1f}%)")

    print("\n  -- alpha's dependence on probe count k --")
    tbl = al.groupby("k_probes").agg(
        n=("alpha", "size"), median_alpha=("alpha", "median"),
        median_implied_r=("implied_mean_inter_probe_r", "median")).round(3)
    print(tbl.head(8).to_string())
    print("\n    alpha RISES with k while the implied mean inter-probe")
    print("    correlation FALLS. The rise is the Spearman-Brown formula, not")
    print("    evidence that genes with more probes are better measured.")
    print("    This is precisely why alpha cannot serve as an aggregation rule.")

    print("\n  WHAT THIS SECTION DOES AND DOES NOT SUPPORT:")
    print("    DOES    : probe sets vary widely in how much they co-vary, from")
    print("              anti-correlated to near-identical.")
    print("    DOES NOT: show averaging is more accurate for any gene;")
    print("              show a single probe is more accurate for any gene;")
    print("              license an alpha threshold as a selection rule.")
    print("    A real answer needs an independent reference measurement on")
    print("    these samples, which does not exist.")

    findings["probe_agreement_alpha"] = {
        "interpretation": "DESCRIPTIVE internal consistency only; not an "
                          "accuracy measure and not an aggregation decision rule",
        "assumption_violations": [
            "correlated measurement error across probes (shared array, "
            "sample prep, hybridisation, joint RMA + ComBat)",
            "single-latent-construct assumption unverifiable",
            "no reference measurement available to assess accuracy",
        ],
        "n_genes": int(len(al)),
        "median_alpha": float(al.alpha.median()),
        "alpha_rises_with_k_while_implied_r_falls": True,
        "reliability_gain_claim": "WITHDRAWN -- was an algebraic restatement "
                                  "of alpha and k, not an independent measurement",
    }

    # ---------------------------------------------------------------
    section("5. MAPPING-COUNT RECONCILIATION (8,893 vs 8,894)")
    # ---------------------------------------------------------------
    hdr = find_table_header(paths.GPL570_ANNOTATION)
    raw = pd.read_csv(paths.GPL570_ANNOTATION, sep="\t", skiprows=hdr,
                      dtype=str, low_memory=False)
    ctrl = raw[raw["ID"].notna() & raw["ID"].str.startswith("!")]
    acct = pd.DataFrame([
        {"category": "unmapped (0 gene symbols)",
         "n": int((mapping.table.mapping_class == "unmapped").sum())},
        {"category": "unique (exactly 1 symbol)",
         "n": int((mapping.table.mapping_class == "unique").sum())},
        {"category": "ambiguous (2+ symbols)",
         "n": int((mapping.table.mapping_class == "ambiguous").sum())},
    ])
    acct.loc[len(acct)] = {"category": "TOTAL probes", "n": int(acct.n.sum())}
    print(acct.to_string(index=False))
    print(f"\n  raw annotation table rows          : {len(raw)}")
    print(f"  SOFT control records (ID starts '!'): {len(ctrl)} "
          f"-> {ctrl['ID'].tolist()}")
    print(f"  real probes                        : {len(raw) - len(ctrl)}")
    print(f"  expression matrix probes           : {expr.shape[0]}")
    print(f"  annotation probe set == expression : "
          f"{set(mapping.table['ID']) == set(expr.index)}")
    print("\n  EXPLANATION: script 11 uses dropna(subset=['ID']), which keeps")
    print("  '!platform_table_end' because that string is not NaN. It therefore")
    print("  counts 54,676 'probes' and 8,894 'unmapped'. src/annotation.py")
    print("  filters IDs beginning with '!', giving 54,675 and 8,893.")
    print("  8,893 is correct; 8,894 counts a SOFT terminator as a probe.")
    print("  Script 11 is UNMODIFIED; this is a reconciliation, not a fix.")
    acct.to_csv(OUT / "probe_accounting.csv", index=False)
    findings["reconciliation"] = {
        "raw_rows": int(len(raw)), "control_records": int(len(ctrl)),
        "real_probes": int(len(raw) - len(ctrl)),
        "unmapped_correct": 8893, "unmapped_script11": 8894,
        "cause": "'!platform_table_end' retained by dropna(subset=['ID'])",
    }

    # ---------------------------------------------------------------
    section("6. TRANSPORTABILITY OF EACH STRATEGY")
    # ---------------------------------------------------------------
    print("  What must travel to an independent cohort for the gene-level")
    print("  representation to be reproduced EXACTLY:\n")
    tt = pd.DataFrame([
        {"strategy": "representative_iqr / _mean",
         "frozen object": "list of 21,815 (gene -> probe_id) pairs",
         "depends on target cohort data?": "NO",
         "depends on annotation version?": "only to name genes",
         "re-derivation needed?": "none"},
        {"strategy": "median / mean",
         "frozen object": "gene -> {probe_id,...} groupings (43,163 probes)",
         "depends on target cohort data?": "NO (if all probes present)",
         "depends on annotation version?": "YES - groupings change with it",
         "re-derivation needed?": "recompute per cohort"},
    ])
    print(tt.to_string(index=False))
    print("\n  Both can be frozen on discovery. The difference is what breaks")
    print("  if the target cohort differs:")
    print("   * representative: needs 21,815 specific probes present.")
    print("   * median/mean   : needs ALL probes of each gene present, else the")
    print("     aggregate is computed over a different probe set than in")
    print("     discovery -- the same gene name denoting a different quantity.")
    print(f"\n  probes required: representative 21,815 | median/mean 43,163")
    print("\n  LIMITATIONS, both strategies:")
    print("   * Annotation version: GPL570 symbols have been revised since")
    print("     deposit. Freezing the (gene -> probe) list freezes the symbol")
    print("     vocabulary too; a target cohort annotated with a newer build")
    print("     may disagree about which gene a probe belongs to.")
    print("   * Probe availability: a same-platform GPL570 cohort has all")
    print("     54,675 probes, so neither strategy loses probes. A DIFFERENT")
    print("     Affymetrix platform carries only a subset -- this must be")
    print("     verified per cohort in Phase 4, not assumed.")
    print("   * Neither strategy addresses the deposited RMA+ComBat scale.")
    print("     Freezing a probe list does NOT make values comparable across")
    print("     cohorts. That is the Phase 4 problem and is untouched here.")
    tt.to_csv(OUT / "transportability_comparison.csv", index=False)
    findings["transportability"] = {
        "representative_probes_required": 21815,
        "aggregate_probes_required": 43163,
        "representative_frozen_object": "gene->probe_id list",
        "aggregate_frozen_object": "gene->probe-set groupings",
    }

    (OUT / "phase2b_part2.json").write_text(
        json.dumps(findings, indent=2, default=str) + "\n")
    section("DONE")
    print(f"  evidence -> {OUT.relative_to(paths.ROOT)}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

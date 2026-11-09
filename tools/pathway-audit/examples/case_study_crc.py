"""Case study: run pathway-audit on five colorectal-cancer cohorts (GSE39582, GSE41258, GSE143985, GSE103479, TCGA-COAD).

Run from the repository root with the project environment:
    .venv/bin/python tools/pathway-audit/examples/case_study_crc.py
Requires the downloaded inputs under data/ (see docs) and a reference cache (pathway-audit fetch-references).
Outputs: reports/pathway_audit_case_study/<cohort>/REPORT.md, audit.json, and SUMMARY.md.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "tools" / "pathway-audit" / "src"))
sys.path.insert(0, str(REPO / "scripts" / "research_extensions"))
sys.path.insert(0, str(REPO / "scripts" / "exploratory"))

from pathway_audit import report  # noqa: E402
from pathway_audit.audit import audit  # noqa: E402
from pathway_audit.report_html import write_html  # noqa: E402

REFS = REPO / "data" / "external" / "pathway_audit_refs"
OUT = REPO / "reports" / "pathway_audit_case_study"
NAMES = {"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "EMT", "HALLMARK_E2F_TARGETS": "E2F",
         "HALLMARK_WNT_BETA_CATENIN_SIGNALING": "Wnt", "HALLMARK_P53_PATHWAY": "p53"}


def gene_sets():
    p = json.loads((REPO / "config/pathways/m2_panel_4set.json").read_text())
    return {NAMES[k]: v["members_resolved"] for k, v in p["gene_sets"].items()}


def cohort_gse39582():
    import xcommon as X
    ctx = X.context()
    cov = ctx.cov.loc[ctx.ids]
    ann = pd.DataFrame({"tp53": cov["tp53.mutation"], "braf": cov["braf.mutation"], "kras": cov["kras.mutation"],
                        "mmr": cov["mmr.status"]}, index=ctx.ids)
    ep = X.endpoint(ctx)
    ann["time"] = ep.rfs_time_full.astype(str).values
    ann["event"] = ep.rfs_event_full.astype(str).values
    ann["age"] = cov["age.at.diagnosis_(year)"].values
    ann["male"] = (cov["sex"] == "Male").astype(int).astype(str).values
    ann["stage3"] = (cov["tnm.stage"] == "3.0").astype(int).astype(str).values
    ann["proximal"] = (cov["tumor.location"] == "proximal").astype(int).astype(str).values
    ann["chemo"] = (cov["chemotherapy.adjuvant"] == "Y").astype(int).astype(str).values
    binary = {"TP53 mutant": ("tp53", "M", "WT"), "BRAF mutant": ("braf", "M", "WT"), "KRAS mutant": ("kras", "M", "WT"), "dMMR": ("mmr", "dMMR", "pMMR")}
    return ctx.X, ann, binary, ["mmr=dMMR"], {"time": "time", "event": "event", "clinical": ["age", "male", "stage3", "proximal", "chemo"]}


def cohort_geo(acc):
    from e3_replication import COHORTS, gene_matrix, labels, platform_map, read_matrix
    a, expr = read_matrix(acc)
    G = gene_matrix(expr, platform_map(COHORTS[acc]["gpl"]), COHORTS[acc]["linear"])
    ids, L = labels(acc, a)
    if acc == "GSE41258":
        pid = a.loc[ids, "patient id"]
        ids = list(pid[~pid.duplicated()].index)
    G = G[[i for i in ids if i in G.columns]]
    ann = pd.DataFrame({k: v.reindex(G.columns).map({True: "pos", False: "neg"}).fillna("") for k, v in L.items()}, index=G.columns)
    binary = {f"{k} positive": (k, "pos", "neg") for k in L}
    return G, ann, binary, None, None


def cohort_tcga():
    from e3c_tcga import T, load
    G, L = load()
    ann = pd.DataFrame({k: v.map({True: "pos", False: "neg"}).fillna("") for k, v in L.items()}, index=G.columns)
    mut = pd.read_csv(T / "COAD_mc3_gene_level.txt.gz", sep="\t", index_col=0)
    burden = mut.sum(axis=0).reindex(G.columns)
    ann["log_mutation_burden"] = np.log1p(burden).map(lambda v: "" if pd.isna(v) else f"{v:.6f}")
    # MSI annotation is nearly absent in this RNA-seq subset, so hypermutation is represented by log mutation burden
    return G, ann, {f"{k} positive": (k, "pos", "neg") for k in ("TP53", "BRAF")}, ["log_mutation_burden"], None


def main():
    sets = gene_sets()
    cohorts = {"GSE39582": cohort_gse39582, "GSE41258": lambda: cohort_geo("GSE41258"), "GSE143985": lambda: cohort_geo("GSE143985"),
               "GSE103479": lambda: cohort_geo("GSE103479"), "TCGA-COAD": cohort_tcga}
    summary = []
    for name, fn in cohorts.items():
        G, ann, binary, adj, surv = fn()
        res = audit(G, sets, ann, binary, adj, REFS, surv, n_boot=1000, n_null=300, module_boot=100, seed=20261108)
        report.write(res, OUT / name, f"Pathway-score audit: {name}")
        write_html(res, OUT / name, f"Pathway-score audit: {name}")
        row = {"cohort": name, "samples": res["n_samples"], "qc_median_r": round(res["qc"]["median_r"], 2), "qc": "PASS" if res["qc"]["passed"] else "FAIL"}
        for s, r in res["gene_sets"].items():
            row[s] = "; ".join(r["verdict"])
        if "survival_screen" in res:
            row["survival_screen_corrected_dC"] = round(res["survival_screen"]["corrected_delta_C"], 4)
        summary.append(row)
        print(name, row["qc"], row.get("EMT", "")[:80])
    pd.DataFrame(summary).to_csv(OUT / "summary.csv", index=False)
    L = ["# pathway-audit case study: five colorectal cancer cohorts", "",
         "Gene sets: the four frozen Hallmark sets of the GSE39582 study. Composition: ESTIMATE + MCP-counter (computational estimates).", "",
         "| Cohort | Samples | QC (median r) | EMT | p53 |", "|---|---|---|---|---|"]
    for r in summary:
        L.append(f"| {r['cohort']} | {r['samples']} | {r['qc']} ({r['qc_median_r']}) | {r.get('EMT', '')} | {r.get('p53', '')} |")
    L += ["", "Per-cohort reports: `<cohort>/REPORT.md`."]
    (OUT / "SUMMARY.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()

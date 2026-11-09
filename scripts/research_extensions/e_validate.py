"""Validation checks for research extensions E1–E3. Writes reports/research_extensions/validation_checks.json; exits non-zero on failure."""
import hashlib
import json
import sys

import numpy as np
import pandas as pd

import rcommon as C
from rcommon import X

checks = []


def check(name, expected, observed, ok):
    checks.append({"check": name, "expected": expected, "observed": observed, "result": "PASS" if ok else "FAIL"})


def main():
    reg = json.loads((C.REP / "registration.json").read_text())
    cur = hashlib.sha256((X.ROOT / reg["plan"]).read_bytes()).hexdigest()
    check("extension plan unchanged since registration", reg["sha256"], cur, cur == reg["sha256"])
    prov = json.loads((C.EXT / "PROVENANCE.json").read_text())
    bad = [f["file"] for f in prov["files"] if hashlib.sha256((C.EXT / f["file"]).read_bytes()).hexdigest() != f["sha256"]]
    check("E3 downloads match recorded checksums", "all match", bad or "all match", not bad)
    ctx = X.context()
    check("GSE39582 cohort", "452", len(ctx.ids), len(ctx.ids) == 452)
    # E1: module scores derived from frozen gene sets only
    g = pd.read_csv(C.RES / "e1_gene_contributions.csv")
    ok = all(set(g.gene[g.pathway == lab]) <= set(ctx.panel["gene_sets"][pid]["members_resolved"]) for pid, lab in C.SETS)
    check("E1 genes are members of the frozen gene sets", True, ok, ok)
    check("E1 gene counts per set", "198/195/41/194", "/".join(str((g.pathway == lab).sum()) for _, lab in C.SETS),
          [int((g.pathway == lab).sum()) for _, lab in C.SETS] == [198, 195, 41, 194])
    s1 = json.loads((C.REP / "e1_summary.json").read_text())
    check("E1 H-E1 evaluability follows registered ARI rule (p53 median ARI < 0.7 → not evaluable)", "consistent",
          f"ARI {s1['stability']['p53']['median_ARI']:.3f}, evaluable {s1['H_E1']['evaluable']}",
          (s1["stability"]["p53"]["median_ARI"] >= 0.7) == s1["H_E1"]["evaluable"])
    # E1 full-set scores equal frozen cohort scores
    z, _ = X.cohort_scores(ctx)
    Z = C.zrows(ctx.X)
    d = max(float(np.max(np.abs(C.score(Z, ctx.panel["gene_sets"][pid]["members_resolved"]).to_numpy() - z[lab].to_numpy()))) for pid, lab in C.SETS)
    check("E1 scoring helper reproduces frozen cohort scores", "< 1e-10", d, d < 1e-10)
    # E2
    ws = pd.read_csv(C.RES / "e2_within_subtype_emt_stroma.csv")
    check("E2 within-subtype n sums to 452", 452, int(ws[ws.composition == "MCP_Fibroblasts"].n.sum()), int(ws[ws.composition == "MCP_Fibroblasts"].n.sum()) == 452)
    # E3
    s3 = json.loads((C.REP / "e3_summary.json").read_text())
    ci = s3["cohorts"]
    check("E3 GSE41258 primary tumours only", 186, ci["GSE41258"]["n_tumours"], ci["GSE41258"]["n_tumours"] == 186)
    check("E3 GSE143985 / GSE103479 tumour counts", "91 / 156", f"{ci['GSE143985']['n_tumours']} / {ci['GSE103479']['n_tumours']}",
          ci["GSE143985"]["n_tumours"] == 91 and ci["GSE103479"]["n_tumours"] == 156)
    cov = pd.read_csv(C.RES / "e3_gene_coverage.csv")
    check("E3 frozen coverage rule met for all scores", "all", int(cov.coverage_rule_met.sum()), bool(cov.coverage_rule_met.all()))
    r = pd.read_csv(C.RES / "e3_replication_tests.csv")
    check("E3 six hypotheses × three cohorts recorded", 18, len(r), len(r) == 18)
    check("E3 QC verdict set present", "present", "verdicts_excluding_qc_failed" in s3, "verdicts_excluding_qc_failed" in s3)
    # independence (no shared GSM with GSE39582)
    gsm39582 = set(X.soft_characteristics().index)
    shared = [acc for acc in ("GSE41258", "GSE143985", "GSE103479")
              if set(pd.read_csv(C.RES / f"e3_scores_{acc}.csv").sample_id) & gsm39582]
    check("E3 cohorts share no samples with GSE39582", "none", shared or "none", not shared)
    # frozen primary untouched
    iv = json.loads((X.ROOT / "reports/phase9/corrected_rerun_report/corrected_interval_results.json").read_text())
    check("frozen θ̂ unchanged", 0.004316, round(iv["theta_hat"], 6), round(iv["theta_hat"], 6) == 0.004316)
    check("three extension figures (PNG + SVG)", 3, len(list(C.FIG.glob("FE*.png"))), len(list(C.FIG.glob("FE*.png"))) == 3 and len(list(C.FIG.glob("FE*.svg"))) == 3)
    ch = json.loads((C.REP / "challenge_tests.json").read_text())
    check("challenge: gene-disjoint R6 rho >= 0.9 in all three QC-passing cohorts", ">= 0.9",
          {k: round(ch["C2_R6_gene_disjoint"][k][0], 3) for k in ("GSE39582", "GSE41258", "GSE143985")},
          all(ch["C2_R6_gene_disjoint"][k][0] >= 0.9 for k in ("GSE39582", "GSE41258", "GSE143985")))
    for d in ("RESEARCH_QUESTION_REASSESSMENT", "BIOLOGICAL_CONTRIBUTION_STRENGTHENING_PLAN", "RESEARCH_EXTENSION_REGISTER",
              "FINAL_RESEARCH_CONTRIBUTION", "INDEPENDENT_SCIENTIFIC_CHALLENGE"):
        f = X.ROOT / "docs" / f"{d}.md"
        check(f"deliverable docs/{d}.md exists", True, f.is_file(), f.is_file())
    rt = json.loads((C.REP / "registration_tcga.json").read_text())
    cur_t = hashlib.sha256((X.ROOT / rt["plan"]).read_bytes()).hexdigest()
    check("TCGA addendum unchanged since registration", rt["sha256"], cur_t, cur_t == rt["sha256"])
    tp = json.loads((X.ROOT / "data/external/tcga_coad_xena/PROVENANCE.json").read_text())
    badt = [f["file"] for f in tp["files"] if hashlib.sha256((X.ROOT / "data/external/tcga_coad_xena" / f["file"]).read_bytes()).hexdigest() != f["sha256"]]
    check("TCGA downloads match recorded checksums", "all match", badt or "all match", not badt)
    ts = json.loads((C.REP / "e3_tcga_summary.json").read_text())
    check("TCGA passes positive-control QC", ">= 0.5", round(ts["qc"]["median_r"], 3), ts["qc"]["median_r"] >= 0.5)
    check("TCGA primary tumours", 286, ts["n_primary_tumours"], ts["n_primary_tumours"] == 286)
    tch = json.loads((C.REP / "tcga_challenge_tests.json").read_text())
    base = tch["C6_random_baseline_R6"]
    ok_b = all(ch["C2_R6_gene_disjoint"][k][0] > base[k]["p97_5"] for k in ("GSE39582", "GSE41258", "GSE143985")) and ts["R6_gene_disjoint"][0] > base["TCGA-COAD"]["p97_5"]
    check("R6 exceeds each cohort's random-gene-set 97.5th percentile", True, {k: round(v["p97_5"], 2) for k, v in base.items()}, ok_b)
    r1 = tch["C5_tcga_hypermutation"]["R1"]["excluding_hypermutated"]["rank_biserial"]
    check("TCGA R1 robust to excluding hypermutated tumours (CI < 0)", "< 0", round(r1[2], 3), r1[2] < 0)
    X.write_json(C.REP / "validation_checks.json", {"n": len(checks), "n_fail": sum(c["result"] == "FAIL" for c in checks), "checks": checks})
    for c in checks:
        print(c["result"], c["check"], "|", c["observed"])
    if any(c["result"] == "FAIL" for c in checks):
        sys.exit(1)


if __name__ == "__main__":
    main()

"""Tests on synthetic data with known structure (no network, no references needed)."""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pathway_audit import association as A
from pathway_audit.audit import audit
from pathway_audit.cli import main as cli_main
from pathway_audit.genes import gene_contributions, modules
from pathway_audit.qc import coexpression_qc
from pathway_audit.scoring import score_gene_set, zscore_rows


def synthetic(n=300, seed=1):
    """Bulk-like matrix: a 'stromal' factor drives collagens and a STROMA set; a proliferation factor drives a PROLIF set;
    a binary 'mutation' lowers a TARGET set directly; noise genes elsewhere."""
    rng = np.random.default_rng(seed)
    stroma = rng.normal(size=n)
    prolif = rng.normal(size=n)
    mut = rng.random(n) < 0.4
    rows = {}
    for g in ["COL1A1", "COL1A2", "COL3A1", "SPARC", "COL5A1"] + [f"STR{i}" for i in range(40)]:
        rows[g] = 8 + 1.0 * stroma + 0.4 * rng.normal(size=n)
    for g in ["MKI67", "TOP2A", "BUB1", "CCNB1"] + [f"PRO{i}" for i in range(40)]:
        rows[g] = 7 + 1.0 * prolif + 0.5 * rng.normal(size=n)
    tcell = rng.normal(size=n)
    for g in ["CD2", "CD3E"]:
        rows[g] = 5 + tcell + 0.4 * rng.normal(size=n)
    epi = rng.normal(size=n)
    for g in ["CDX2", "VIL1"]:
        rows[g] = 9 + epi + 0.4 * rng.normal(size=n)
    for i in range(30):
        rows[f"TGT{i}"] = 6 - 1.0 * mut + 0.8 * rng.normal(size=n)
    for i in range(300):
        rows[f"NOISE{i}"] = 5 + rng.normal(size=n)
    expr = pd.DataFrame(rows).T
    expr.columns = [f"S{i}" for i in range(n)]
    ann = pd.DataFrame({"mut": np.where(mut, "M", "WT"), "stroma_true": stroma}, index=expr.columns)
    return expr, ann


def test_qc_passes_clean_and_fails_misaligned():
    expr, _ = synthetic()
    assert coexpression_qc(expr)["passed"]
    rng = np.random.default_rng(0)
    bad = expr.copy()
    bad.index = rng.permutation(bad.index.to_numpy())   # values no longer match their gene labels
    q = coexpression_qc(bad)
    assert not q["passed"] and q["median_r"] < 0.3


def test_coverage_rule():
    expr, _ = synthetic()
    z = zscore_rows(expr)
    s, cov = score_gene_set(z, [f"STR{i}" for i in range(10)] + [f"MISSING{i}" for i in range(10)])
    assert s is None and not cov["coverage_ok"]
    s, cov = score_gene_set(z, [f"STR{i}" for i in range(20)])
    assert s is not None and abs(s.std(ddof=1) - 1) < 1e-9


def test_gene_contributions_and_modules():
    expr, _ = synthetic()
    z = zscore_rows(expr)
    genes = [f"STR{i}" for i in range(20)] + [f"PRO{i}" for i in range(20)]
    gc = gene_contributions(z, genes)
    assert gc.logo_r_with_full.min() > 0.95
    m = modules(z, genes, n_boot=30)
    assert len(m) == 2 and set(m.n_genes) == {20} and (m.stability == "stable").all()


def test_genotype_effect_explained_vs_robust():
    rng = np.random.default_rng(3)
    n = 400
    conf = rng.normal(size=n)
    lab = pd.Series((conf + rng.normal(size=n) > 0.5).astype(float))
    y_conf = pd.Series(conf + 0.5 * rng.normal(size=n))                     # driven only by the confounder
    y_direct = pd.Series(-1.0 * lab.to_numpy() + 0.5 * rng.normal(size=n))  # driven by the label itself
    adj = pd.DataFrame({"conf": conf})
    e1 = A.genotype_effect(y_conf, lab, adj, rng, B=300)
    e2 = A.genotype_effect(y_direct, lab, adj, rng, B=300)
    assert e1["verdict"] == "explained by adjustment covariates"
    assert e2["verdict"] == "robust to adjustment"


def test_gene_disjoint_and_partial_null():
    expr, ann = synthetic()
    z = zscore_rows(expr)
    rng = np.random.default_rng(0)
    set_genes = [f"STR{i}" for i in range(30)] + ["COL1A1", "COL3A1"]
    d = A.gene_disjoint_marker_rho(z, expr, set_genes, ["COL1A1", "COL1A2", "COL3A1", "SPARC"], rng, B=200)
    assert d["n_shared_removed"] == 2 and d["rho"][0] > 0.9
    W = ann[["stroma_true"]].to_numpy()
    target = z.loc[[f"PRO{i}" for i in range(20)]].mean(axis=0)
    null = A.random_set_null(z, target, 30, lambda a, b: A.partial_spearman(a, b, W), rng, n=50, pool=[f"NOISE{i}" for i in range(300)])
    assert abs(np.median(null)) < 0.1


def test_audit_end_to_end_without_references(tmp_path):
    expr, ann = synthetic()
    sets = {"STROMA": [f"STR{i}" for i in range(30)], "TARGET": [f"TGT{i}" for i in range(30)], "TOO_SMALL": ["STR0", "X1", "X2"]}
    res = audit(expr, sets, ann, binary={"MUT": ("mut", "M", "WT")}, n_boot=200, module_boot=20)
    assert res["qc"]["passed"]
    t = res["gene_sets"]["TARGET"]["genotype_effects"]["MUT"]
    assert t["unadjusted"][2] < 0 and t["verdict"] == "associated"
    assert res["gene_sets"]["TOO_SMALL"]["verdict"][0].startswith("not scored")
    from pathway_audit.report import write
    p = write(res, tmp_path / "rep")
    assert "Data-quality gate" in p.read_text()


def test_cli(tmp_path):
    expr, ann = synthetic(n=120)
    expr.to_csv(tmp_path / "expr.tsv", sep="\t")
    ann.to_csv(tmp_path / "ann.csv")
    (tmp_path / "sets.gmt").write_text("TARGET\tna\t" + "\t".join(f"TGT{i}" for i in range(30)) + "\n")
    rc = cli_main(["run", "--expr", str(tmp_path / "expr.tsv"), "--gene-sets", str(tmp_path / "sets.gmt"), "--out", str(tmp_path / "out"),
                   "--annotations", str(tmp_path / "ann.csv"), "--binary", "MUT=mut:M:WT", "--log2", "no"])
    assert rc == 0 and (tmp_path / "out" / "audit.json").exists()


@pytest.mark.skipif(not os.environ.get("PATHWAY_AUDIT_REFS"), reason="set PATHWAY_AUDIT_REFS to a fetch-references cache")
def test_estimate_reimplementation_matches_package():
    from pathway_audit.composition import validate_estimate
    assert validate_estimate(Path(os.environ["PATHWAY_AUDIT_REFS"])) < 1e-6


def test_delta_c_noise_scores_near_zero():
    pytest.importorskip("lifelines")
    from pathway_audit.predictive import delta_c
    rng = np.random.default_rng(5)
    n = 300
    x = rng.normal(size=n)
    t = rng.exponential(np.exp(-0.7 * x))
    e = (rng.random(n) < 0.6).astype(float)
    r = delta_c(pd.DataFrame({"x": x}), pd.DataFrame({"noise1": rng.normal(size=n), "noise2": rng.normal(size=n)}), t, e, B=30)
    assert r["apparent_delta_C"] >= -0.01 and r["corrected_delta_C"] < r["apparent_delta_C"]


def test_demo_recovers_planted_truth(tmp_path):
    from pathway_audit.demo import run_demo
    res, html_path = run_demo(tmp_path / "demo")
    v = {k: " ".join(s["verdict"]) for k, s in res["gene_sets"].items()}
    assert "composition-dominated" in v["STROMAL_PROGRAMME"] and "holds without shared genes" in v["STROMAL_PROGRAMME"]
    assert "incoherent" in v["RANDOM_SET"]
    ge = res["gene_sets"]
    assert ge["TARGET_PROGRAMME"]["genotype_effects"]["Target mutation"]["verdict"] == "robust to adjustment"
    assert ge["CONFOUNDED_PROGRAMME"]["genotype_effects"]["Driver mutation"]["verdict"] == "explained by adjustment covariates"
    t = html_path.read_text()
    assert t.startswith("<!doctype html>") and "<svg" in t and "http" not in t.split("<main>")[1]


def test_cli_demo(tmp_path):
    assert cli_main(["demo", "--out", str(tmp_path / "d")]) == 0
    assert (tmp_path / "d" / "report.html").exists() and (tmp_path / "d" / "REPORT.md").exists()

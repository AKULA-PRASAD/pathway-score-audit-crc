"""Results figures (scripts/49_make_results_figures.py): traceability, labelling,
determinism and non-overwrite.

Every check re-reads the frozen sources independently of the script. The script is run twice,
into fresh temporary directories, and must produce byte-identical files. Requires the locally
held frozen outputs, analysis bundle and checkpoints (skipped when they are absent).
"""
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "49_make_results_figures.py"
CIR = ROOT / "reports/phase9/corrected_rerun_report/corrected_interval_results.json"
EX = ROOT / "reports/phase9/execution/phase9_execution_results.json"
S7 = ROOT / "reports/phase9/execution_s7/phase9_s7_results.json"
ASRUN = ROOT / "reports/phase9/final_project/results/s1_interval_results.json"
BUNDLE = ROOT / "data/processed/analysis"
CKPT = ROOT / "reports/phase9/corrected_rerun/runs/run0.jsonl"
THETA = 0.004316342804030704
FIGS = ["Fig1_study_design_cohort", "Fig2_preprocessing_qc", "Fig3_exploratory_structure", "Fig4_pathway_scores",
        "Fig5_clinical_model", "Fig6_molecular_hazard_ratios", "Fig7_discrimination", "Fig8_corrected_deltaC",
        "Fig9_robustness", "FigS1_monte_carlo_precision", "FigS2_informative_censoring"]
PATHWAYS = ["HALLMARK_E2F_TARGETS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
            "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_P53_PATHWAY"]

pytestmark = pytest.mark.skipif(not (BUNDLE.is_dir() and CKPT.is_file() and CIR.is_file()),
                                reason="frozen outputs, bundle or checkpoints not available locally")


def tree_hashes(d: Path) -> dict:
    return {p.relative_to(d).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(d.rglob("*")) if p.is_file()}


def protected_hashes() -> dict:
    roots = [ROOT / "reports/phase8a/figures", ROOT / "reports/phase8b/figures", ROOT / "reports/phase8c/figures",
             ROOT / "reports/phase9/final_project/figures"]
    out = {}
    for r in roots:
        out.update({f"{r.name}/{k}": v for k, v in tree_hashes(r).items()})
    for f in (CIR, EX, S7):
        out[f.name] = hashlib.sha256(f.read_bytes()).hexdigest()
    return out


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    before = protected_hashes()
    outs = []
    for i in (1, 2):
        d = tmp_path_factory.mktemp(f"figs{i}") / "out"
        r = subprocess.run([sys.executable, "-B", str(SCRIPT), str(d)], capture_output=True, text=True, cwd=ROOT, timeout=600)
        assert r.returncode == 0, r.stderr[-2000:]
        outs.append(d)
    after = protected_hashes()
    return outs, before, after


@pytest.fixture(scope="module")
def values(runs):
    return json.loads((runs[0][0] / "figure_values.json").read_text())


def test_all_outputs_present_and_deterministic(runs):
    a, b = runs[0]
    expected = {f"{n}.{x}" for n in FIGS for x in ("svg", "png")} | {"FIGURE_CAPTIONS.md", "figure_values.json", "figure_sources.json"}
    assert set(tree_hashes(a)) == expected
    assert tree_hashes(a) == tree_hashes(b), "repeated execution must produce byte-identical outputs"


def test_no_existing_figure_or_frozen_result_modified(runs):
    _, before, after = runs
    assert before == after


def test_refuses_non_empty_output_directory(runs, tmp_path):
    d = tmp_path / "occupied"
    d.mkdir()
    (d / "keep.txt").write_text("existing")
    r = subprocess.run([sys.executable, "-B", str(SCRIPT), str(d)], capture_output=True, text=True, cwd=ROOT, timeout=120)
    assert r.returncode != 0 and "refusing" in (r.stderr + r.stdout)
    assert (d / "keep.txt").read_text() == "existing" and len(list(d.iterdir())) == 1


def test_corrected_intervals_and_theta_match_verified_results(values):
    cir = json.loads(CIR.read_text())
    assert values["fig8"]["theta_hat"] == THETA == cir["theta_hat"]
    for r in range(5):
        x = cir["runs"][str(r)]
        v = values["fig8"]["runs"][str(r)]
        assert v["percentile"] == [x["percentile_primary"]["lower"], x["percentile_primary"]["upper"]]
        assert v["basic"] == [x["basic_diagnostic"]["lower"], x["basic_diagnostic"]["upper"]]
        assert v["bc"] == [x["bc_z0_diagnostic"]["lower"], x["bc_z0_diagnostic"]["upper"]]
    assert values["fig8"]["B_out"] == 750 and values["fig8"]["verification"]["verdict"] == "PASS"


def test_superseded_as_run_results_not_plotted(values):
    asrun = json.loads(ASRUN.read_text())
    old_lower = set(asrun["between_run"]["percentile_lower"]["values"])
    old_upper = set(asrun["between_run"]["percentile_upper"]["values"])
    for r in range(5):
        lo, hi = values["fig8"]["runs"][str(r)]["percentile"]
        assert lo not in old_lower and hi not in old_upper


def test_outer_distributions_reproduce_intervals(values):
    vals = []
    with open(ROOT / "reports/phase9/corrected_rerun/runs/run2.jsonl", "rb") as f:
        f.readline()
        vals = [json.loads(l)["record"]["dC_corrected"] for l in f]
    lo, hi = np.quantile(vals, [0.025, 0.975], method="linear")
    assert [lo, hi] == pytest.approx(values["fig8"]["runs"]["2"]["percentile"], abs=1e-15)


def test_hazard_ratios_match_sources_and_roles(values):
    ex, s7 = json.loads(EX.read_text()), json.loads(S7.read_text())
    f6 = values["fig6"]
    assert f6["primary_median"]["role"] == "primary"
    assert f6["s7_repiqr"]["role"] == "sensitivity" and f6["five_set"]["role"] == "sensitivity"
    for key, src, model in (("primary_median", ex, "M2_primary"), ("s7_repiqr", s7, "M2_primary"), ("five_set", ex, "M2_sensitivity")):
        for p, h in f6[key]["hr"].items():
            a = src["apparent"][model]
            assert (h["hr"], h["lo"], h["hi"]) == (a["hazard_ratios"][p], a["hr_ci_lower"][p], a["hr_ci_upper"][p])
    assert set(f6["primary_median"]["hr"]) == set(PATHWAYS)
    assert "HALLMARK_PI3K_AKT_MTOR_SIGNALING" in f6["five_set"]["hr"]


def test_figure_text_labels_primary_and_sensitivity(runs):
    caps = (runs[0][0] / "FIGURE_CAPTIONS.md").read_text()
    for needle in ("**Primary** four-set model", "Five-set **sensitivity** model", "not independent replication",
                   "differ only in random seed", "not Wald confidence intervals", "paired out-of-original",
                   "not absolute recurrence probability", "not added predictive discrimination",
                   "conditional on the observed 136 events"):
        assert needle in caps, needle
    assert "seed 0" not in (runs[0][0] / "Fig8_corrected_deltaC.svg").read_text()


def test_discrimination_and_calibration_match_execution(values):
    ex = json.loads(EX.read_text())
    for m in ("M1", "M2_primary", "M2_sensitivity"):
        assert values["fig7"][m]["apparent_C"] == ex["apparent"][m]["C"]
        assert values["fig7"][m]["corrected_C"] == ex["monte_carlo"][m]["corrected_C_mean"]
        cs = np.mean([s["corrected_slope"] for s in ex["monte_carlo"][m]["per_seed"]])
        assert values["fig9"]["calibration"][m]["corrected"] == pytest.approx(cs, abs=0)
    lobo = ex["batch_sensitivity"]["leave_one_batch_out"]
    assert [x["C_apparent"] for x in values["fig9"]["lobo"]] == [x["C_apparent"] for x in lobo]


def test_m1_terms_and_age_curve(values):
    ex = json.loads(EX.read_text())
    m1 = ex["apparent"]["M1"]
    for t, h in values["fig5"]["terms"].items():
        assert (h["hr"], h["lo"], h["hi"]) == (m1["hazard_ratios"][t], m1["hr_ci_lower"][t], m1["hr_ci_upper"][t])
    c = values["fig5"]["age_curve"]
    k = m1["age_knots"]
    assert c["reference_age"] == k[1]
    i = int(np.argmin(np.abs(np.asarray(c["age"]) - k[1])))
    assert abs(c["log_hr_vs_ref"][i]) < 0.02   # curve passes through 0 at the reference age


def test_cit_counts_reconcile_with_frozen_partition(values):
    def rows(name):
        with open(BUNDLE / name, newline="") as f:
            return list(csv.DictReader(f))
    s3 = [r["sample_id"] for r in rows("cohort_membership.csv") if r["in_S3"] == "True"]
    cov = {r["sample_id"]: r for r in rows("covariates.csv")}
    ep = {r["sample_id"]: float(r["rfs_event_full"]) for r in rows("endpoint.csv")}
    part = {r["sample_id"]: r["dataset"] for r in rows("partition_and_batch.csv")}
    a19 = [s for s in s3 if cov[s]["chemotherapy.adjuvant"] != "" and all(
        cov[s][v] != "" for v in ("age.at.diagnosis_(year)", "sex", "tnm.stage", "tumor.location"))]
    for p in ("discovery", "validation"):
        sel = [s for s in a19 if part[s] == p]
        ev = int(sum(ep[s] for s in sel))
        assert values["fig1"]["A19_by_partition"][p] == {"n": len(sel), "events": ev, "censored": len(sel) - ev}
    v = values["fig1"]
    assert v["A19_by_partition"]["discovery"] == {"n": 353, "events": 107, "censored": 246}
    assert v["A19_by_partition"]["validation"] == {"n": 99, "events": 29, "censored": 70}
    assert v["A19_excluded"] == [{"sample_id": "GSM972262", "partition": "discovery", "event": 0}]
    assert (v["A19"]["n"], v["A19"]["events"]) == (452, 136)


def test_sources_record_hashes(runs):
    src = json.loads((runs[0][0] / "figure_sources.json").read_text())
    assert src["models_fitted"] == 0
    for rel, h in src["inputs"].items():
        assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == h, rel
    assert src["script"]["sha256"] == hashlib.sha256(SCRIPT.read_bytes()).hexdigest()

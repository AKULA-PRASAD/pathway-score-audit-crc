"""Report / presentation / table / figure consistency tests (read-only)."""
import csv, hashlib, json, re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[5]
FP = ROOT / "reports/phase9/final_project"
DOCS = [FP / "documentation/PHASE9_FINAL_REPORT.md", FP / "presentation/PHASE9_TECHNICAL_PRESENTATION.html",
        FP / "PROJECT_COMPLETION_REPORT.md"]
FIGS = ["F1_outer_distributions", "F2_interval_comparison", "F3_endpoint_variability",
        "F4_MC1_running_trace", "F5_EP1_grid", "F6_shift_decomposition"]


@pytest.fixture(scope="module")
def res():
    return json.loads((FP / "results/s1_interval_results.json").read_text())


@pytest.mark.parametrize("p", DOCS, ids=lambda p: p.name)
def test_rendered_no_placeholders_and_labelled(p):
    t = p.read_text()
    assert "{{" not in t                      # no unrendered template placeholder
    assert "PROVISIONAL" in t


def test_report_contains_every_per_run_percentile_interval(res):
    t = (FP / "documentation/PHASE9_FINAL_REPORT.md").read_text()
    for x in res["runs"]:
        p = x["percentile_primary"]
        assert f"[{p['lower']:+.5f}, {p['upper']:+.5f}]" in t


def test_claim_all_intervals_include_zero_is_true(res):
    for x in res["runs"]:
        for k in ("percentile_primary", "basic_diagnostic", "bc_z0_diagnostic"):
            assert x[k]["lower"] < 0 < x[k]["upper"]


def test_claims_about_warnings_are_true(res):
    w = {x["run_id"]: x["warnings"] for x in res["runs"]}
    assert all(v["W6_breach"] and v["W7_breach"] for v in w.values())          # "all runs"
    assert [r for r, v in w.items() if v["W5_breach"]] == [1]                  # "run 1"
    assert not any(v["W3_breach"] for v in w.values())
    assert not any(g["breach"] for g in res["W4"])
    assert all(x["z0_diagnostic"]["z0"] < 0 for x in res["runs"])             # "above theta"
    assert res["bias_decomposition_observed"]["shift_dC_corrected"] > res["theta_hat"]["value"]  # "larger than theta"


def test_claims_about_ep1_are_true(res):
    e = res["EP1"]
    assert e["by_estimator"]["kernel_bw_x1.0"]["all_pass"] and e["by_estimator"]["between_run_SD"]["all_pass"]
    assert [(g["run"], g["endpoint"], g["estimator"]) for g in e["failing_cases"]] == [(1, "upper", "kernel_bw_x0.5")]


def test_figures_and_captions_exist():
    cap = (FP / "figures/FIGURE_CAPTIONS.md").read_text()
    for f in FIGS:
        assert (FP / f"figures/{f}.png").stat().st_size > 10_000 and (FP / f"figures/{f}.svg").exists()
        assert f"### {f}" in cap


def test_report_references_only_existing_files():
    t = (FP / "documentation/PHASE9_FINAL_REPORT.md").read_text()
    for m in re.findall(r"\]\((\.\./[^)]+)\)", t):
        assert (FP / "documentation" / m).resolve().exists(), m


def test_csv_tables_reproduce_json(res):
    with open(FP / "tables/T1_intervals_per_run.csv") as fh:
        for row in csv.DictReader(fh):
            x = res["runs"][int(row["run"])]
            assert float(row["pct_lower"]) == x["percentile_primary"]["lower"]
            assert float(row["bc_upper"]) == x["bc_z0_diagnostic"]["upper"]


def test_deviation_records():
    v1 = ROOT / "reports/phase9/h5_spotcheck/CV_SEED_DEVIATION_DRAFT.md"
    assert hashlib.sha256(v1.read_bytes()).hexdigest().startswith("2866571cecf471af")   # retained unchanged
    v2 = (FP / "documentation/CV_SEED_DEVIATION_DRAFT_v2.md").read_text()
    assert v2.startswith("# DRAFT — NOT YET ADJUDICATED")
    for needle in ("line 364", "lines 117", "142", "Possible future investigations", "Disclosure considerations"):
        assert needle in v2

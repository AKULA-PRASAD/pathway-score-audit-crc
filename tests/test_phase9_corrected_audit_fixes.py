"""Regression tests for the defects found by the independent engineering audit."""
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pytest

from src.phase9_corrected import gate as gate_mod
from src.phase9_corrected import reliability as rel
from src.phase9_corrected import statistics as st
from src.phase9_corrected.config import (PLAN_BASE, ConfigurationError, CorrectedConfig,
                                         canonical_json)
from src.phase9_corrected.data import LAM_GRID
from src.phase9_corrected.engine import compute_outer_replicate
from src.phase9_corrected.orchestration import run_corrected
from src.phase9_corrected.rng import PlanStreams
from src.phase9_corrected.verification import (v1_self_consistency, v2_record_assignment,
                                               v5_spot_check)
from tests.p9c_support import FakePipeline, approved_test_config, synthetic_context

ROOT = Path(__file__).resolve().parents[1]
S = PlanStreams(PLAN_BASE, 0, "independent", 7)


def rec(k=0, B_in=6, fake=None, ctx=None):
    return json.loads(canonical_json(compute_outer_replicate(
        ctx or synthetic_context(), S, fake or FakePipeline(), 0, k, B_in, LAM_GRID)))


@pytest.fixture(scope="module")
def sc():
    spec = importlib.util.spec_from_file_location(
        "fp_stats_core_ro2", ROOT / "reports/phase9/final_project/analysis/stats_core.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ---- M2: plan-fixed sizes
@pytest.mark.parametrize("field,value", [("B_in", 10), ("B_in", 200), ("B_out_initial", 600)])
def test_plan_sizes_enforced(field, value):
    c = approved_test_config(ROOT, "scratch_p9c_test_out")
    c[field] = value
    if field == "B_in":
        c["B_in_max"] = max(value, 100)
    with pytest.raises(ConfigurationError, match=field):
        CorrectedConfig(c, ROOT).validate_for_execution()


# ---- minor: resume not blocked by worker count; .git/.venv protected; non-empty dir refused
def test_worker_count_not_in_config_hash():
    a = approved_test_config(ROOT, "x_out", n_workers=8)
    b = approved_test_config(ROOT, "x_out", n_workers=3)
    assert CorrectedConfig(a, ROOT).sha256() == CorrectedConfig(b, ROOT).sha256()
    b["ep1_option"] = "A"
    assert CorrectedConfig(a, ROOT).sha256() != CorrectedConfig(b, ROOT).sha256()


@pytest.mark.parametrize("d", [".git/p9c", ".venv/p9c"])
def test_git_and_venv_protected(d):
    with pytest.raises(ConfigurationError, match="protected"):
        CorrectedConfig(approved_test_config(ROOT, d), ROOT).validate_for_execution()


def test_existing_foreign_output_dir_refused_before_any_write():
    cfg = CorrectedConfig(approved_test_config(ROOT, "reports/phase1"), ROOT)
    with pytest.raises(ConfigurationError, match="holds no p9c manifest"):
        run_corrected(cfg, 0, 500)


# ---- major: V1 fit-consumption check is no longer tautological
def test_v1_detects_fit_that_consumed_a_different_sample():
    r = rec()
    assert v1_self_consistency([r], synthetic_context().e, S)["status"] == "PASS"
    t = json.loads(json.dumps(r))
    t["inner"][2]["fits"]["M1"]["consumed_multiplicity_sha256"] = "0" * 64
    t["apparent"]["M2"]["consumed_multiplicity_sha256"] = "1" * 64
    res = v1_self_consistency([t], synthetic_context().e, S)
    assert res["status"] == "FAIL" and res["n_mismatches"] == 2


def test_engine_consumed_hash_is_computed_from_fit_inputs():
    """A fit fed the wrong ids must show up: the consumed hash comes from the id list."""
    class Spy(FakePipeline):
        seen = []

        def pipeline_fit(self, X, cov, panel, train_ids, *a, **kw):
            Spy.seen.append(list(train_ids))
            return super().pipeline_fit(X, cov, panel, train_ids, *a, **kw)
    r = rec(fake=Spy())
    from src.phase9_corrected.resampling import INNER_MULT_TAG, tagged_sha256
    ctx = synthetic_context()
    pos = {s: i for i, s in enumerate(ctx.ids)}
    # fits: M1 apparent, M2 apparent, then (M1, M2) per b
    for b in range(6):
        ids_b = Spy.seen[2 + 2 * b]
        h = tagged_sha256(INNER_MULT_TAG, np.bincount([pos[s] for s in ids_b], minlength=len(ctx.ids)))
        assert h == r["inner"][b]["fits"]["M1"]["consumed_multiplicity_sha256"] == r["inner"][b]["resample"]["multiplicity_sha256"]


# ---- minor: V2 when M2 legitimately skips CV
def test_v2_accepts_m2_without_cv_and_flags_inconsistent_consumption():
    class NoCV(FakePipeline):
        def pipeline_fit(self, X, cov, panel, train_ids, apply_ids, t, e, grid, rng):
            m, a, b, d = super().pipeline_fit(X, cov, None, train_ids, apply_ids, t, e, grid, rng)
            return m, a, b, d                    # all pathways "dropped": no CV, no inner_cv_c
    r = rec(fake=NoCV())
    assert r["apparent"]["M2"]["cv_rng_consumed"] is False
    assert v2_record_assignment([r], S)["status"] == "PASS"
    t = json.loads(json.dumps(r))
    t["inner"][0]["fits"]["M2"]["inner_cv_c"] = 0.6        # claims CV ran but nothing drawn
    assert v2_record_assignment([t], S)["status"] == "FAIL"


def test_v5_blocked_when_selected_replicates_absent():
    assert v5_spot_check({0: {}}, lambda k: {}, selected_k=(0, 25))["status"] == "BLOCKED"


# ---- minor: statistics edge cases
def test_silverman_degenerate_sample_raises():
    with pytest.raises(ValueError):
        st.silverman_bandwidth(np.zeros(50))


def test_z0_tie_counts_reported():
    c = st.z0_counts(np.array([1.0, 2.0, 2.0, 3.0]), 2.0)
    assert c == {"n_less": 1, "n_tie": 2, "B": 4}


# ---- major: post-run reliability layer
def to_hist_shape(r):
    pairs = [row for row in r["inner"] if row["paired_complete"]]
    return {"inner_C": {m: [[row["fits"][m]["C_fit"], row["fits"][m]["C_eval"]] for row in pairs]
                        for m in ("M1", "M2")},
            "dC_corrected": r["dC_corrected"], "dC_app": r["dC_app"]}


def test_bounding_matches_final_project_convention(sc):
    calls = {"n": 0}

    def fail_when(ids, panel):
        if panel is not None:
            calls["n"] += 1
            return "x" if calls["n"] == 4 else None
    r = rec(B_in=30, fake=FakePipeline(fail_when))
    assert r["B_in_paired_complete"] == 29
    mine, theirs = rel.bounding_record(r), sc.bounding_record(to_hist_shape(r), 30)
    assert mine["n_missing"] == theirs["n_missing"] == 1
    for f in ("dC_min", "dC_max", "width"):
        assert math.isclose(mine[f], theirs[f], rel_tol=0, abs_tol=1e-15)
    assert mine["dC_min"] <= r["dC_corrected"] <= mine["dC_max"]


def test_w3_matches_final_project_convention(sc):
    recs = [rec(k) for k in range(8)]
    mine = rel.w3_inner_noise(recs, list(range(8)))
    theirs = sc.inner_noise([to_hist_shape(r) for r in recs])
    assert math.isclose(mine["noise_var"], theirs["mean_inner_noise_var"], rel_tol=1e-12)
    assert math.isclose(mine["SD_inflation"], theirs["SD_inflation"], rel_tol=1e-12)


def test_6_5_h1_h2_require_explicit_mcse_and_classify():
    recs = [rec(k) for k in range(6)]
    with pytest.raises(ValueError):
        rel.section_6_5(recs, None)
    s = rel.section_6_5(recs, 1e-3)
    assert s["n_excluded"] == 0 and s["H1"] == "not triggered" and s["usable_k"] == list(range(6))
    bad = json.loads(json.dumps(recs))
    bad[0]["status"] = "unreliable_apparent_failed"
    assert rel.section_6_5(bad, 1e-3)["H1"] == "UNRELIABLE"           # 1/6 > 5 %
    with pytest.raises(ValueError):
        rel.h2(recs, s["usable_k"], None)
    h = rel.h2(recs, s["usable_k"], {"lower": 1e-3, "upper": 1e-3})
    assert h["vacuous"] and h["H2"] == "not triggered"
    assert rel.w4(0.003, 0.001)["W4"] and not rel.w4(0.0015, 0.001)["W4"]


# ---- major: pre-execution gate actually runs
def test_pre_execution_gate_passes_on_small_declared_domain(tmp_path):
    c = approved_test_config(ROOT, "scratch_p9c_test_out",
                             ladder={"rungs": [750], "beyond_ladder_rule": "stop_and_label_imprecise"},
                             run_ids=[0, 1])
    cfg = CorrectedConfig(c, ROOT)
    g = gate_mod.run_pre_execution_gate(cfg, "m" * 64, out_path=tmp_path / "g.json")
    assert g["all_pass"], {k: (v if isinstance(v, str) else v.get("status")) for k, v in g.items()}
    assert g["V2_domain"]["n_states"] == 2 * 750 * (2 + 3 * 100) + 5
    assert g["V1_outer_historical"]["n_compared"] == 1000
    assert not (ROOT / "scratch_p9c_test_out").exists()


def test_pre_execution_gate_fails_on_numpy_mismatch(tmp_path, monkeypatch):
    monkeypatch.setattr(gate_mod, "PINNED_NUMPY", "0.0.0")
    monkeypatch.setattr(gate_mod, "enumerate_stream_states",
                        lambda *a, **k: {"pass": True, "collisions": [], "n_states": 0})
    monkeypatch.setattr(gate_mod, "v1_outer_historical_all", lambda *a, **k: {"status": "PASS"})
    cfg = CorrectedConfig(approved_test_config(ROOT, "scratch_p9c_test_out"), ROOT)
    g = gate_mod.run_pre_execution_gate(cfg, "m" * 64, out_path=tmp_path / "g.json")
    assert g["V2_numpy"]["status"] == "FAIL" and g["V2_structural"]["status"] == "FAIL"
    assert not g["all_pass"]

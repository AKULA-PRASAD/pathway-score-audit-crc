"""Synthetic fixtures and a fast stand-in pipeline for the corrected Phase 9 tests.

NOT scientific data. The fake pipeline mimics the frozen pipeline's interface and its
RNG behaviour (M1 draws nothing; M2 draws a CV fold permutation and picks lambda
from it), so stream assignment, pairing and aggregation can be tested quickly.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from src.phase9_corrected.data import AGE_COL, AnalysisContext
from src.phase9_execution import c_index as real_c_index


class FakeFailure(RuntimeError):
    def __init__(self, reason):
        super().__init__(reason)
        self.reason = reason


class FakePipeline:
    """fail_when(train_ids, panel) -> reason or None, for failure injection."""
    name = "fake"

    def __init__(self, fail_when=None):
        self.fail_when = fail_when
        self.failure_types = (FakeFailure,)
        self.c_index = real_c_index

    def pipeline_fit(self, X, cov, panel, train_ids, apply_ids, t_train, e_train, lam_grid, rng):
        if self.fail_when is not None:
            why = self.fail_when(train_ids, panel)
            if why:
                raise FakeFailure(why)
        x_tr = cov.loc[list(train_ids), "x"].to_numpy(float)
        x_ap = cov.loc[list(apply_ids), "x"].to_numpy(float)
        if panel is None:
            lam, cols_tr, cols_ap = 0.0, x_tr[:, None], x_ap[:, None]
            beta = np.array([1.0 + 0.01 * float(np.mean(x_tr))])
        else:
            perm = rng.permutation(len(train_ids))            # consumes the CV generator
            lam = float(lam_grid[int(perm[0]) % len(lam_grid)])
            z_tr = cov.loc[list(train_ids), "z"].to_numpy(float)
            z_ap = cov.loc[list(apply_ids), "z"].to_numpy(float)
            cols_tr, cols_ap = np.c_[x_tr, z_tr], np.c_[x_ap, z_ap]
            beta = np.array([1.0, 0.5 / (1.0 + lam)])
        diag = {"lambda": lam, "age_knots": [float(np.min(x_tr)), float(np.median(x_tr)),
                                              float(np.max(x_tr))],
                "age_linear_fallback": False}
        if panel is not None:
            diag["inner_cv_c"] = 0.6
        return {"beta": beta}, cols_tr, cols_ap, diag

    @staticmethod
    def linear_predictor(m, d):
        return np.asarray(d) @ m["beta"]


def synthetic_context(n=60, seed=12345, events_frac=0.35) -> AnalysisContext:
    g = np.random.default_rng(seed)                      # fixture data only
    ids = [f"S{i:03d}" for i in range(n)]
    x = g.normal(size=n)
    z = g.normal(size=n)
    t = np.round(g.exponential(scale=np.exp(-0.7 * x - 0.3 * z)) * 100, 1) + 1.0
    e = (g.random(n) < events_frac).astype(float)
    e[:3] = 1.0
    cov = pd.DataFrame({"x": x, "z": z}, index=pd.Index(ids, name="sample_id"))
    cov[AGE_COL] = np.round(60 + 10 * x, 1).astype(str)
    cov["sex"] = np.where(g.random(n) < 0.5, "Male", "Female")
    cov["tnm.stage"] = np.where(g.random(n) < 0.5, "3.0", "2.0")
    cov["tumor.location"] = np.where(g.random(n) < 0.5, "proximal", "distal")
    cov["chemotherapy.adjuvant"] = np.where(g.random(n) < 0.5, "Y", "N")
    return AnalysisContext(ids, cov, t.astype(float), e, None, {"panel": "synthetic"})


def synthetic_real_pipeline_context(n=80, seed=777) -> AnalysisContext:
    """Synthetic context usable by the REAL frozen pipeline (clinical covariates as strings,
    a two-pathway panel and a small expression matrix)."""
    g = np.random.default_rng(seed)
    ids = [f"R{i:03d}" for i in range(n)]
    genes = [f"G{j}" for j in range(50)]          # >= 20 genes per set (frozen coverage rule)
    X = pd.DataFrame(g.normal(size=(50, n)), index=genes, columns=ids)
    risk = 0.6 * X.loc[genes[:25]].mean().to_numpy() + 0.3 * g.normal(size=n)
    t = np.round(g.exponential(scale=np.exp(-risk)) * 100, 1) + 1.0
    e = (g.random(n) < 0.45).astype(float)
    age = 65 + 8 * g.normal(size=n)
    cov = pd.DataFrame(index=pd.Index(ids, name="sample_id"))
    cov[AGE_COL] = np.round(age, 1).astype(str)
    cov["sex"] = np.where(g.random(n) < 0.5, "Male", "Female")
    cov["tnm.stage"] = np.where(g.random(n) < 0.5, "3.0", "2.0")
    cov["tumor.location"] = np.where(g.random(n) < 0.5, "proximal", "distal")
    cov["chemotherapy.adjuvant"] = np.where(g.random(n) < 0.5, "Y", "N")
    panel = {"scoring": {"sigma_ddof": 1},
             "gene_sets": {"P1": {"members_resolved": genes[:25], "n_source": 25},
                           "P2": {"members_resolved": genes[25:], "n_source": 25}}}
    return AnalysisContext(ids, cov, t.astype(float), e, X, panel)


def approved_test_config(root, output_dir, **over) -> dict:
    """A fully-populated config for TESTS ONLY. Every value marked approved here is a
    test fixture, not an owner decision."""
    from src.phase9_corrected.config import template_config
    c = template_config()
    c.update({"cv_seed_deviation_adjudicated": True,
              "apparent_cv": {"convention": "independent", "T": 7},
              "B_in_max": 100, "as_run_variant": False, "ep1_option": "B",
              "w5": {"action": "W5-B", "corrective_rule": None,
                     "definition": "adjudication_D5_as_proposed",
                     "exclusion_rule": "usable_within_k_positions_v1"},
              "ladder": {"rungs": [750, 1000, 1500, 2000, 3000],
                         "beyond_ladder_rule": "stop_and_label_imprecise"},
              "output_dir": str(output_dir), "verification_gate": "V1-V5_H4_as_implemented_p9c-1",
              "runtime_budget_hours_per_run": 1.0})
    c["approvals"] = {f: {"approved": True, "date": "TEST", "reference": "TEST FIXTURE"}
                      for f in ("cv_seed_deviation_adjudicated", "apparent_cv", "B_in_max",
                                "as_run_variant", "ep1_option", "w5", "ladder", "output_dir",
                                "verification_gate", "runtime_budget_hours_per_run")}
    c["execution_authorization"] = {"authorized": True, "date": "TEST", "reference": "TEST FIXTURE"}
    c.update(over)
    return c


def synthetic_gate_body(config, msha, fail_check=None) -> dict:
    """A SYNTHETIC six-check gate result with the structure gate.py writes. NOT the gate:
    no check is evaluated. `fail_check` names one check to report as FAIL."""
    import numpy as _np
    from src.phase9_corrected import orchestration as _orch
    k_max = max([config.B_out_initial] + list((config.get("ladder") or {}).get("rungs") or []))
    n = _orch.expected_v2_domain_states(config.run_ids, k_max, config.B_in_max, config.convention)
    body = {"key": {"config_sha256": config.sha256(), "manifest_sha256": msha, "numpy": _np.__version__},
            "H4": "PASS", "historical_preservation": "PASS",
            "V2_numpy": {"status": "PASS", "detail": {"all_pass": True}},
            "V2_structural": {"status": "PASS", "detail": {"holds": True}},
            "V2_domain": {"status": "PASS", "pass": True, "collisions": [], "n_states": n, "n_distinct": n,
                          "domain": {"run_ids": list(config.run_ids), "k_max": k_max,
                                     "b_max": config.B_in_max, "convention": config.convention,
                                     "T": config.T}},
            "V1_outer_historical": {"status": "PASS", "n_compared": 500 * len(config.run_ids),
                                    "mismatches": []},
            "all_pass": True}
    if fail_check in ("H4", "historical_preservation"):
        body[fail_check] = "FAIL"
    elif fail_check:
        body[fail_check]["status"] = "FAIL"
    if fail_check:
        body["all_pass"] = False
    return body


def make_synthetic_gate(fail_check=None, calls=None):
    """Test double for orchestration.run_pre_execution_gate: writes synthetic evidence to
    `out_path` exactly as the real gate writes its result file, and returns it."""
    from src.phase9_corrected.config import canonical_json as _cj

    def gate(config, msha, progress=None, out_path=None):
        from pathlib import Path as _P
        body = synthetic_gate_body(config, msha, fail_check)
        p = _P(out_path) if out_path else _P(config.output_dir) / "gate" / "pre_execution_gate.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(_cj(body))
        if calls is not None:
            calls.append(str(p))
        return body
    return gate


synthetic_gate = make_synthetic_gate()

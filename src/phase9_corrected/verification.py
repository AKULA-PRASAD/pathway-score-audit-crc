"""Verification gate V1-V5 (H4 lives in integrity.py, as a separate hard check).

Every check returns {"check", "status", ...} with status in
PASS | FAIL | NOT_APPLICABLE | BLOCKED | NOT_COMPARABLE.

Historical limitation (stated, not hidden): the historical checkpoints store only a
16-hex OUTER multiplicity hash. For inner resamples there is no stored historical
evidence, so V1 against history can establish (a) outer identity from stored hashes and
(b) inner CODE-PATH identity (new driver vs frozen driver function regenerating from the
same streams) -- never what the 2026-11-03/04 processes actually drew. V3's M1 inner C
match is corroborating, not proving, evidence (C is not injective in the resample).
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

from .resampling import draw_inner, draw_outer
from .rng import HistoricalAsRunStreams, state_fingerprint

FROZEN_DRIVER = "scripts/40_phase9_dc_interval.py"


def load_frozen_driver(root: Path):
    """Load the frozen driver as a module WITHOUT calling main() (as the H5 spot-check did)."""
    spec = importlib.util.spec_from_file_location("p9_dc_driver_frozen_ro", Path(root) / FROZEN_DRIVER)
    drv = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = drv
    spec.loader.exec_module(drv)
    return drv


def load_historical_run(root: Path, run_id: int) -> dict:
    p = Path(root) / "reports" / "phase9" / "dc_interval" / f"run{run_id}_outer.jsonl"
    rows = [json.loads(l) for l in p.read_text().splitlines() if l.strip()]
    return {r["k"]: r for r in rows}


def _hex(x):
    return None if x is None else float(x).hex()


# ------------------------------------------------------------------ V1
def v1_self_consistency(records, e: np.ndarray, streams) -> dict:
    """Regenerate every outer and inner resample from the recorded streams and require
    the stored hashes / counts to match exactly (corrected run's own evidence)."""
    n = len(e)
    bad = []
    n_inner = 0
    for rec in records:
        k = rec["k"]
        if rec["streams"]["outer"]["state_sha256"] != state_fingerprint(streams.outer(k)):
            bad.append((k, None, "outer stream fingerprint"))
        o = draw_outer(e, streams.generator(streams.outer(k)))
        ev = o.evidence()
        for f in ("multiplicity_sha256", "legacy_multiplicity_sha16", "n_distinct", "events", "censored"):
            if rec["outer"][f] != ev[f]:
                bad.append((k, None, f"outer {f}"))
        if rec["outer"]["multiplicity"] != ev["multiplicity"]:
            bad.append((k, None, "outer multiplicity vector"))
        eO = e[o.indices]
        for row in rec.get("inner", []):
            b = row["b"]
            n_inner += 1
            J = draw_inner(eO, o.indices, n, streams.generator(streams.inner(k, b)))
            if row["resample"] != J.evidence():
                bad.append((k, b, "inner resample evidence"))
            for m, fd in row["fits"].items():
                if fd["resample_index_sha256"] != J.index_sha256:
                    bad.append((k, b, f"fit {m} labelled with a different resample"))
                if fd.get("consumed_multiplicity_sha256") != J.multiplicity_sha256:
                    bad.append((k, b, f"fit {m} consumed a different resample"))
        for m, fd in (rec.get("apparent") or {}).items():
            if fd.get("consumed_multiplicity_sha256") != rec["outer"]["multiplicity_sha256"]:
                bad.append((k, None, f"apparent {m} consumed a different sample"))
    return {"check": "V1-self", "status": "PASS" if not bad else "FAIL",
            "n_outer": len(records), "n_inner": n_inner, "mismatches": bad[:50],
            "n_mismatches": len(bad)}


def v1_outer_vs_historical(records, historical: dict) -> dict:
    """Outer identity against stored historical hashes (k present in both)."""
    common = [r for r in records if r["k"] in historical]
    bad = [r["k"] for r in common
           if r["outer"]["legacy_multiplicity_sha16"] != historical[r["k"]]["multiplicity_sha"]]
    st = "NOT_APPLICABLE" if not common else ("PASS" if not bad else "FAIL")
    return {"check": "V1-outer-historical", "status": st, "n_compared": len(common),
            "mismatched_k": bad, "evidence": "stored historical multiplicity_sha (16 hex)"}


def v1_inner_vs_frozen_path(records, e: np.ndarray, frozen_driver, base: int, run_id: int,
                            B_in_hist: int) -> dict:
    """Inner CODE-PATH identity: regenerate the as-run inner indices with the FROZEN
    driver's own strat_boot and its as-run stream construction, compare index-for-index
    with the corrected run's stored inner evidence. Labelled as code-path evidence."""
    from .resampling import INNER_INDEX_TAG, tagged_sha256
    hs = HistoricalAsRunStreams(base, run_id, B_in_hist)
    bad, n = [], 0
    for rec in records:
        k = rec["k"]
        I = frozen_driver.strat_boot(e, np.random.default_rng(hs.outer(k)))
        eO = e[I]
        for row in rec.get("inner", []):
            b = row["b"]
            if b >= B_in_hist:
                continue
            J = frozen_driver.strat_boot(eO, np.random.default_rng(hs.inner(k, b)))
            n += 1
            if tagged_sha256(INNER_INDEX_TAG, J) != row["resample"]["index_sha256"]:
                bad.append((k, b))
    st = "NOT_APPLICABLE" if n == 0 else ("PASS" if not bad else "FAIL")
    return {"check": "V1-inner-frozen-code-path", "status": st, "n_compared": n,
            "mismatches": bad[:50],
            "scope": "establishes identical inner draws between corrected and frozen code paths; "
                     "does NOT establish what the historical processes drew (no stored evidence)"}


# ------------------------------------------------------------------ V2 (record-level)
def _expected_consumption(model, fd):
    """M1 never draws. M2 draws iff the ridge CV ran (it is skipped when every pathway
    is dropped, src/phase9_execution.py:195-202); inner_cv_c is recorded only then."""
    return model == "M2" and fd.get("inner_cv_c") is not None


def v2_record_assignment(records, streams) -> dict:
    """Recorded stream fingerprints equal the prescribed construction; the CV stream of
    every fit differs from its resample stream; only M2 consumes CV draws."""
    bad = []
    for rec in records:
        k = rec["k"]
        if rec["streams"]["apparent_cv"]["state_sha256"] != state_fingerprint(streams.apparent_cv(k)):
            bad.append((k, None, "apparent_cv fingerprint"))
        if (getattr(streams, "convention", None) == "independent"
                and rec["streams"]["apparent_cv"]["state_sha256"] == rec["streams"]["outer"]["state_sha256"]):
            bad.append((k, None, "apparent CV state equals outer resample state"))
        for m, fd in (rec.get("apparent") or {}).items():
            if fd["ok"] and fd["cv_rng_consumed"] != _expected_consumption(m, fd):
                bad.append((k, None, f"apparent {m} CV consumption {fd['cv_rng_consumed']}"))
        for row in rec.get("inner", []):
            b = row["b"]
            s = row["streams"]
            if s["inner"]["state_sha256"] != state_fingerprint(streams.inner(k, b)):
                bad.append((k, b, "inner fingerprint"))
            if s["inner_cv"]["state_sha256"] != state_fingerprint(streams.inner_cv(k, b)):
                bad.append((k, b, "inner_cv fingerprint"))
            if s["inner_cv"]["state_sha256"] == s["inner"]["state_sha256"]:
                bad.append((k, b, "inner CV state equals inner resample state"))
            for m, fd in row["fits"].items():
                if fd["ok"] and fd["cv_rng_consumed"] != _expected_consumption(m, fd):
                    bad.append((k, b, f"inner {m} CV consumption {fd['cv_rng_consumed']}"))
    return {"check": "V2-record-assignment", "status": "PASS" if not bad else "FAIL",
            "n_mismatches": len(bad), "mismatches": bad[:50]}


# ------------------------------------------------------------------ V3
def v3_m1_invariance(records, historical: dict) -> dict:
    """M1 quantities that must be invariant to the CV-seed correction (M1 draws no CV
    numbers): apparent C, lambda (= 0), knots, per-b inner C at the stored 8-decimal
    precision, and opt_M1_mean when the paired set is identical."""
    bad, nc, n = [], [], 0
    for rec in records:
        k = rec["k"]
        h = historical.get(k)
        if h is None:
            continue
        n += 1
        if rec["status"] != "ok" or h["status"] != "ok":
            nc.append((k, "status not ok"))
            continue
        if _hex(rec["C_app"]["M1"]) != _hex(h["C_app"]["M1"]):
            bad.append((k, "C_app.M1"))
        if rec["apparent"]["M1"]["lam"] != h["lambda_apparent"]["M1"]:
            bad.append((k, "lambda_apparent.M1"))
        if rec["apparent"]["M1"]["age_knots"] != h["age_knots"]["M1"]:
            bad.append((k, "age_knots.M1"))
        # historical inner_C lists hold paired-complete b only, in b order, with no b label;
        # positional mapping is valid only where every historical b was paired-complete.
        if h["B_in_paired_complete"] != h["B_in_attempted"]:
            nc.append((k, "historical paired set incomplete: positional b mapping unavailable"))
            continue
        hist_m1 = h["inner_C"]["M1"]
        for row in rec["inner"]:
            b = row["b"]
            if b >= len(hist_m1):
                continue
            fd = row["fits"]["M1"]
            if not fd["ok"]:
                bad.append((k, f"inner b={b} M1 failed (historical succeeded)"))
                continue
            new = np.round(np.array([fd["C_fit"], fd["C_eval"]]), 8).tolist()
            if new != hist_m1[b]:
                bad.append((k, f"inner_C.M1 b={b}"))
        same_set = (rec["B_in_attempted"] == h["B_in_attempted"]
                    and all(r["paired_complete"] for r in rec["inner"]))
        if same_set:
            if _hex(rec["opt_mean"]["M1"]) != _hex(h["opt_M1_mean"]):
                bad.append((k, "opt_M1_mean"))
        else:
            nc.append((k, "paired set differs: opt_M1_mean not comparable"))
    st = "NOT_APPLICABLE" if n == 0 else ("PASS" if not bad else "FAIL")
    return {"check": "V3", "status": st, "n_compared": n, "failures": bad[:50],
            "n_failures": len(bad), "not_comparable": nc[:50], "n_not_comparable": len(nc)}


# ------------------------------------------------------------------ V4
def v4_as_run_apparent(records, historical: dict, config_convention: str, as_run_variant) -> dict:
    """Applies ONLY to a computation that uses the as-run apparent convention: the optional
    paired variant, or the primary if the owner chose the as_run convention. Never compares
    the independent-stream primary against the as-run values."""
    if config_convention == "as_run":
        field_app, field_dc = "apparent", "dC_app"
    elif as_run_variant is True:
        field_app, field_dc = "apparent_as_run", "dC_app_as_run"
    elif as_run_variant is None:
        return {"check": "V4", "status": "BLOCKED", "reason": "as_run_variant decision not supplied"}
    else:
        return {"check": "V4", "status": "NOT_APPLICABLE",
                "reason": "independent primary without the as-run variant"}
    bad, n = [], 0
    for rec in records:
        h = historical.get(rec["k"])
        if h is None or h["status"] != "ok" or field_dc not in rec:
            continue
        n += 1
        if _hex(rec[field_dc]) != _hex(h["dC_app"]):
            bad.append((rec["k"], field_dc))
        for m in ("M1", "M2"):
            fd = rec[field_app][m]
            if _hex(fd["C_fit"]) != _hex(h["C_app"][m]) or fd["lam"] != h["lambda_apparent"][m]:
                bad.append((rec["k"], f"{field_app}.{m}"))
    st = "NOT_APPLICABLE" if n == 0 else ("PASS" if not bad else "FAIL")
    return {"check": "V4", "status": st, "field": field_app, "n_compared": n, "failures": bad[:50]}


# ------------------------------------------------------------------ V5
V5_SELECTED_K = tuple(range(0, 500, 25))       # plan line 383: 0, 25, ..., 475 of run 0


def v5_spot_check(stored: dict, recompute, selected_k=V5_SELECTED_K) -> dict:
    """Recompute selected replicates in a fresh call and require the FULL record (after the
    same JSON round-trip as the checkpoint) and dC_corrected (float.hex) to be identical."""
    from .config import canonical_json
    missing = [k for k in selected_k if k not in stored]
    if missing:
        return {"check": "V5", "status": "BLOCKED", "selected_k": list(selected_k),
                "reason": f"selected replicates absent: {missing[:10]}"}
    rows = []
    for k in selected_k:
        if k not in stored:
            rows.append({"k": k, "present": False, "identical": False})
            continue
        new = json.loads(canonical_json(recompute(k)))
        old = stored[k]
        diff = sorted(x for x in set(old) | set(new) if old.get(x) != new.get(x))
        rows.append({"k": k, "present": True, "identical": not diff, "differing_fields": diff,
                     "dC_hex_equal": _hex(old.get("dC_corrected")) == _hex(new.get("dC_corrected"))})
    ok = all(r["identical"] and r.get("dC_hex_equal", False) for r in rows)
    return {"check": "V5", "status": "PASS" if ok else "FAIL", "selected_k": list(selected_k),
            "rows": rows}

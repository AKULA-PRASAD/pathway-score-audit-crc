"""Pre-execution gate, run by `run_corrected` before any replicate is computed.

  V2-numpy      NumPy version and SeedSequence behaviour (empirical)
  V2-structural input-level disjointness argument (holds only if V2-numpy passes)
  V2-domain     exhaustive initial-state enumeration over k < ladder top, b < B_in_max
  V1-outer      corrected outer resamples reproduce all stored historical hashes
  H4            frozen inputs, and historical preservation
The result is written once to <output_dir>/gate/pre_execution_gate.json and reused only
when the configuration hash, manifest hash and NumPy version are identical.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .config import canonical_json
from .data import load_endpoint_only
from .integrity import h4_check, historical_preservation_check
from .resampling import draw_outer
from .rng import (PlanStreams, enumerate_stream_states, structural_domain_report,
                  verify_numpy_seed_behaviour)
from .verification import load_historical_run

PINNED_NUMPY = "2.5.4"            # requirements.txt:69; version used for the historical runs


def v1_outer_historical_all(root, base, run_ids, B_hist=500) -> dict:
    ids, e = load_endpoint_only(root)
    bad, n = [], 0
    for r in run_ids:
        hist = load_historical_run(root, r)
        s = PlanStreams(base, r, "as_run")         # outer stream is convention-independent
        for k in range(B_hist):
            n += 1
            if draw_outer(e, s.generator(s.outer(k))).legacy_multiplicity_sha16 != hist[k]["multiplicity_sha"]:
                bad.append((r, k))
    return {"check": "V1-outer-historical", "status": "PASS" if not bad else "FAIL",
            "n_compared": n, "mismatches": bad[:50]}


def run_pre_execution_gate(config, manifest_sha256: str, progress=None, out_path=None) -> dict:
    root = config.root
    out = Path(out_path) if out_path else Path(config.output_dir) / "gate" / "pre_execution_gate.json"
    key = {"config_sha256": config.sha256(), "manifest_sha256": manifest_sha256,
           "numpy": np.__version__}
    if out.exists():
        prev = json.loads(out.read_text())
        if prev.get("key") == key and prev.get("all_pass"):
            return prev
    res = {"key": key}
    res["H4"] = h4_check(root, "gate", raise_on_fail=False)["status"]
    res["historical_preservation"] = historical_preservation_check(root)["status"]
    nb = verify_numpy_seed_behaviour(PINNED_NUMPY, config.base, config.T)
    res["V2_numpy"] = {"status": "PASS" if nb["all_pass"] else "FAIL", "detail": nb}
    sd = structural_domain_report(config.base, config.run_ids, config.T, config.convention,
                                  nb["all_pass"])
    res["V2_structural"] = {"status": "PASS" if sd["holds"] else "FAIL", "detail": sd}
    k_max = max([config.B_out_initial] + list(config.get("ladder")["rungs"]))
    dom = enumerate_stream_states(PlanStreams(config.base, config.run_ids[0], config.convention,
                                              config.T), config.run_ids, k_max, config.B_in_max,
                                  progress=progress)
    res["V2_domain"] = {"status": "PASS" if dom["pass"] else "FAIL",
                        **{k: v for k, v in dom.items() if k != "collisions"},
                        "collisions": dom["collisions"][:50]}
    res["V1_outer_historical"] = v1_outer_historical_all(root, config.base, config.run_ids)
    res["all_pass"] = all(x == "PASS" if isinstance(x, str) else x["status"] == "PASS"
                          for kk, x in res.items() if kk not in ("key",))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(canonical_json(res))
    return res

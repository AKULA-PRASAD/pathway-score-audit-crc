"""Evidence for the CV-seed protocol deviation (documentation only; no model fitting).

For real replicate seeds of run 0, compares the PCG64 starting state and the first raw 64-bit
outputs of:
  (a) the resampling generator           default_rng(ss)              (as executed)
  (b) the CV generator as executed       default_rng(ss)  (fresh)     (driver lines 117 / 142)
  (c) the CV generator as planned        default_rng(ss.spawn(2)[0])  (plan section 8.2 line 364)
at the apparent level (ss = ss_outer) and the inner level (ss = ss_inner[b]).
Also confirms spawn(k+1)[k] == spawn(500)[k] (outer derivation equivalence).
Writes ONLY validation/cv_seed_state_evidence.json.
"""
import json
from pathlib import Path

import numpy as np

FP = Path(__file__).resolve().parents[1]
BASE, RUN = 2026110300, 0
out = {"numpy": np.__version__, "cases": []}
for k in (0, 25, 475):
    root = np.random.SeedSequence(entropy=[BASE, RUN])
    ss_outer = root.spawn(k + 1)[k]
    alt = np.random.SeedSequence(entropy=[BASE, RUN]).spawn(500)[k]
    equiv = ss_outer.spawn_key == alt.spawn_key and \
        ss_outer.generate_state(8).tolist() == alt.generate_state(8).tolist()
    levels = [("apparent", ss_outer)]
    ss_inner = np.random.SeedSequence(entropy=[BASE, RUN]).spawn(k + 1)[k].spawn(100)
    levels += [(f"inner_b{b}", ss_inner[b]) for b in (0, 99)]
    for name, ss in levels:
        a = np.random.default_rng(ss); b_ = np.random.default_rng(ss)
        c = np.random.default_rng(np.random.SeedSequence(entropy=ss.entropy, spawn_key=ss.spawn_key).spawn(2)[0])
        sa, sb, sc_ = a.bit_generator.state, b_.bit_generator.state, c.bit_generator.state
        ra = a.bit_generator.random_raw(4).tolist(); rb = b_.bit_generator.random_raw(4).tolist()
        rc = c.bit_generator.random_raw(4).tolist()
        out["cases"].append({"k": k, "level": name, "spawn_key": list(ss.spawn_key),
                             "outer_spawn_equivalent": equiv,
                             "executed_cv_state_equals_resample_state": sa == sb,
                             "planned_cv_state_equals_resample_state": sa == sc_,
                             "first4_raw_resample": [str(v) for v in ra],
                             "first4_raw_cv_executed": [str(v) for v in rb],
                             "first4_raw_cv_planned": [str(v) for v in rc]})
out["summary"] = {
    "n_cases": len(out["cases"]),
    "all_outer_spawn_equivalent": all(c["outer_spawn_equivalent"] for c in out["cases"]),
    "all_executed_cv_equals_resample": all(c["executed_cv_state_equals_resample_state"] for c in out["cases"]),
    "any_planned_cv_equals_resample": any(c["planned_cv_state_equals_resample_state"] for c in out["cases"])}
(FP / "validation" / "cv_seed_state_evidence.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out["summary"]))

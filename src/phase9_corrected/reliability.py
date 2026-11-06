"""Post-run reliability layer: plan sections 6.4, 6.5 and H1, H2, W3, W4.

The "endpoint MC-SE" against which 6.4/6.5/H2 are judged is NOT fixed by the plan once
three bandwidths exist (it depends on the EP-1 reading, Decision 4). Every function that
needs it takes it as a REQUIRED argument; nothing is defaulted.
"""
from __future__ import annotations

import math

import numpy as np

from .statistics import percentile_ci


def paired_differences(rec) -> np.ndarray:
    """d_b = opt_M2,b - opt_M1,b over paired-complete b, in b order."""
    d = []
    for row in rec["inner"]:
        if row["paired_complete"]:
            f1, f2 = row["fits"]["M1"], row["fits"]["M2"]
            d.append((f2["C_fit"] - f2["C_eval"]) - (f1["C_fit"] - f1["C_eval"]))
    return np.asarray(d, dtype=np.float64)


def bounding_record(rec) -> dict:
    """Plan 6.4 under convention 13: missing paired differences imputed at the
    within-replicate min and max. Vacuous (width 0) when no b is missing."""
    B_in = rec["B_in_attempted"]
    d = paired_differences(rec)
    n_miss = B_in - d.size
    if n_miss == 0:
        v = rec["dC_corrected"]
        return {"n_missing": 0, "dC_min": v, "dC_max": v, "width": 0.0}
    lo = (d.sum() + n_miss * d.min()) / B_in
    hi = (d.sum() + n_miss * d.max()) / B_in
    return {"n_missing": int(n_miss), "dC_min": rec["dC_app"] - hi, "dC_max": rec["dC_app"] - lo,
            "width": float(hi - lo)}


def section_6_5(records, endpoint_mcse: float) -> dict:
    """Outer replicates excluded as unreliable: apparent failure, paired < 95 % of B_in,
    or bounding width > endpoint MC-SE. H1: > 5 % excluded -> UNRELIABLE."""
    if endpoint_mcse is None:
        raise ValueError("endpoint MC-SE for 6.5 must be supplied explicitly")
    ex = {"apparent_failed": [], "inner_shortfall": [], "bounding_width": []}
    for r in records:
        if r["status"] == "unreliable_apparent_failed":
            ex["apparent_failed"].append(r["k"])
        elif r["status"] == "unreliable_inner_shortfall":
            ex["inner_shortfall"].append(r["k"])
        elif bounding_record(r)["width"] > endpoint_mcse:
            ex["bounding_width"].append(r["k"])
    n_ex = sum(len(v) for v in ex.values())
    frac = n_ex / len(records)
    return {"excluded": ex, "n_excluded": n_ex, "fraction": frac,
            "H1": "UNRELIABLE" if frac > 0.05 else "not triggered",
            "usable_k": sorted(set(r["k"] for r in records) - {k for v in ex.values() for k in v})}


def h2(records, usable_k, endpoint_mcse: dict) -> dict:
    """Plan 6.4 decision rule / H2: interval recomputed under the min and max bounds; if
    either endpoint moves by more than its MC-SE the interval is UNRELIABLE."""
    if not endpoint_mcse or set(endpoint_mcse) != {"lower", "upper"}:
        raise ValueError("endpoint MC-SE for both endpoints must be supplied explicitly")
    byk = {r["k"]: r for r in records}
    v = np.array([byk[k]["dC_corrected"] for k in usable_k])
    b = [bounding_record(byk[k]) for k in usable_k]
    p0 = percentile_ci(v)
    pmin = percentile_ci([x["dC_min"] for x in b])
    pmax = percentile_ci([x["dC_max"] for x in b])
    move = {"lower": max(abs(pmin[0] - p0[0]), abs(pmax[0] - p0[0])),
            "upper": max(abs(pmin[1] - p0[1]), abs(pmax[1] - p0[1]))}
    trig = any(move[e] > endpoint_mcse[e] for e in move)
    return {"move": move, "endpoint_mcse": endpoint_mcse,
            "H2": "UNRELIABLE" if trig else "not triggered",
            "vacuous": all(x["n_missing"] == 0 for x in b)}


def w3_inner_noise(records, usable_k) -> dict:
    """Convention 14: noise var = mean_k var(d_b, ddof=1)/n_b; SD inflation =
    sqrt(Var_obs / (Var_obs - noise)) - 1; warning if > 5 %."""
    byk = {r["k"]: r for r in records}
    nv = []
    for k in usable_k:
        d = paired_differences(byk[k])
        nv.append(float(np.var(d, ddof=1)) / d.size)
    noise = float(np.mean(nv))
    obs = float(np.var([byk[k]["dC_corrected"] for k in usable_k], ddof=1))
    infl = math.sqrt(obs / (obs - noise)) - 1 if obs > noise else math.inf
    return {"noise_var": noise, "observed_var": obs, "SD_inflation": infl, "W3": infl > 0.05}


def w4(kernel_se: float, between_sd: float) -> dict:
    r = kernel_se / between_sd if between_sd > 0 else math.inf
    return {"ratio": r, "W4": r > 2 or r < 0.5}

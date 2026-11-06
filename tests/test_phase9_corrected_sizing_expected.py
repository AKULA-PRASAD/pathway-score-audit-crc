"""D11-T1 (owner decision 2026-11-06): independent expected-value tests for the Decision 7B
sizing computation (`orchestration.compute_sizing_decision`).

Independence. Expected values come from the pure-Python reference helper below (`math`
only; no NumPy, no project import). `reports/phase9/final_project/analysis/stats_core.py`
is a separate implementation but shares the same formulas and NumPy primitives, so it is
NOT used as the oracle here. The helper re-implements the published formulas (type-7
quantiles, sample SD ddof=1, Silverman 1986 eq. 3.31, Gaussian KDE, kernel quantile SE per
plan §2.2) independently in code and arithmetic -- not independently in specification.
Selection logic (populations, exclusions, MIN/MAX choices, H1, W5 positions, ladder,
actions) is derived by hand in the comments and recomputed explicitly in `expected()`.

Floating point. Helper and production sum in different orders (Python left-to-right vs
NumPy pairwise), so statistical values are compared with a relative tolerance of 1e-10.
Exact `float.hex()` identity is required only where both sides evaluate the same
expression (B_req from the production rho_max) and for every discrete outcome.

Scope. The 95 % paired-completeness rule is decided in the ENGINE (engine.py, status
`unreliable_inner_shortfall`; tested in test_phase9_corrected_engine.py) -- sizing only
consumes the resulting statuses, so it is not re-derived here.
"""
import math
from pathlib import Path

import pytest

from src.phase9_corrected.config import CorrectedConfig
from src.phase9_corrected.orchestration import compute_sizing_decision
from tests.p9c_support import approved_test_config

ROOT = Path(__file__).resolve().parents[1]
RTOL = 1e-10
FRACTION = 0.05                                   # EP-1 / W5 bound fraction (plan line 84)
BW = (0.5, 1.0, 2.0)                              # plan §2.2 bandwidth multipliers
LADDER = [750, 1000, 1500, 2000, 3000]            # Decision 7A
TRACE_GRID = (50, 100, 200, 300, 400, 500, 750, 1000, 1500, 2000, 3000)   # W5 grid (D5 / 6B)
W5A = {"action": "W5-A", "corrective_rule": "next_rung_above_current",
       "definition": "adjudication_D5_as_proposed", "exclusion_rule": "usable_within_k_positions_v1"}


# ================================================================ pure-Python reference helper
def q7(x, p):
    """Type-7 quantile: h = (n-1)p, linear interpolation between order statistics."""
    s = sorted(x)
    n = len(s)
    h = (n - 1) * p
    j = math.floor(h)
    g = h - j
    return s[j] if j + 1 >= n else s[j] + g * (s[j + 1] - s[j])


def sd(x):
    m = sum(x) / len(x)
    return math.sqrt(sum((v - m) ** 2 for v in x) / (len(x) - 1))


def silverman(x):
    return 0.9 * min(sd(x), (q7(x, 0.75) - q7(x, 0.25)) / 1.34) * len(x) ** -0.2


def kde(x, t, h):
    return sum(math.exp(-0.5 * ((t - v) / h) ** 2) for v in x) / (len(x) * h * math.sqrt(2 * math.pi))


def endpoint_se(x):
    """{endpoint: {bw multiplier: kernel SE}} and the percentile width, plan §2.2."""
    lo, hi = q7(x, 0.025), q7(x, 0.975)
    h = silverman(x)
    se = {e: {m: math.sqrt(p * (1 - p) / len(x)) / kde(x, q, m * h) for m in BW}
          for e, p, q in (("lower", 0.025, lo), ("upper", 0.975, hi))}
    return se, hi - lo


def bounding_width(rec):
    """Plan §6.4 convention 13: missing paired diffs imputed at the within-replicate min/max;
    width = n_missing * (max d - min d) / B_in (0 when nothing is missing)."""
    B_in = rec["B_in_attempted"]
    d = [(row["fits"]["M2"]["C_fit"] - row["fits"]["M2"]["C_eval"])
         - (row["fits"]["M1"]["C_fit"] - row["fits"]["M1"]["C_eval"])
         for row in rec["inner"] if row["paired_complete"]]
    n_miss = B_in - len(d)
    return 0.0 if n_miss == 0 else n_miss * (max(d) - min(d)) / B_in


def expected(by, B, beyond="top_up_to_max_rung_then_reassess"):
    """Independent derivation of the Decision 5A/5B/6A/6B/7A/7B outcome."""
    per = {}
    for r in sorted(by):
        recs = [by[r][k] for k in range(B)]
        # §6.5 MC-SE population (DEC-D): every status "ok" record (each has dC_corrected)
        P = [x["dC_corrected"] for x in recs if x["status"] == "ok"]
        se_P, _ = endpoint_se(P)
        mcse = {e: min(se_P[e].values()) for e in se_P}          # MIN over bandwidths (5B)
        t65 = min(mcse.values())                                  # MIN of the two endpoints
        # exclusions: criteria 1-2 (status) and 3 (bounding width > t65, computed BEFORE)
        excl = sorted([x["k"] for x in recs if x["status"] != "ok"]
                      + [x["k"] for x in recs if x["status"] == "ok" and bounding_width(x) > t65])
        U = [recs[k]["dC_corrected"] for k in range(B) if k not in excl]
        se_U, W = endpoint_se(U)
        ratio = {e: max(se_U[e].values()) / W for e in se_U}     # EP-1 Option B: MAX (5A)
        h1 = len(excl) / B > 0.05
        pts = [g for g in TRACE_GRID if 0.6 * B <= g <= B]
        trace, n_usable = {}, {}
        for g in pts:                                             # 6B: positions k < g kept
            Ug = [recs[k]["dC_corrected"] for k in range(g) if k not in excl]
            trace[g], n_usable[g] = (q7(Ug, 0.025), q7(Ug, 0.975)), len(Ug)
        spread = {i: max(t[i] for t in trace.values()) - min(t[i] for t in trace.values()) for i in (0, 1)}
        bound = FRACTION * W
        per[r] = {"P": P, "se_P": se_P, "mcse": mcse, "t65": t65, "excl": excl, "U": U,
                  "se_U": se_U, "W": W, "ratio": ratio, "h1": h1, "pts": pts, "trace": trace,
                  "n_usable": n_usable, "bound": bound,
                  "w5_unevaluable": h1, "breach": (None if h1 else max(spread.values()) > bound)}
    rho_max = max(max(v["ratio"].values()) for v in per.values())       # 7B: all runs/endpoints
    ep1_fail = rho_max > FRACTION
    B_req = B * (rho_max / FRACTION) ** 2
    if not ep1_fail:
        ep1_next, ep1_action = None, "stop_ep1_pass"
    elif B >= LADDER[-1]:
        ep1_next, ep1_action = None, "stop_label_imprecise"
    else:
        cand = [L for L in LADDER if L >= B_req and L > B]
        if cand:
            ep1_next, ep1_action = cand[0], "top_up"
        elif beyond == "top_up_to_max_rung_then_reassess":
            ep1_next, ep1_action = LADDER[-1], "top_up"
        else:
            ep1_next, ep1_action = None, "stop_label_imprecise"
    if any(v["w5_unevaluable"] for v in per.values()):
        return {"per": per, "rho_max": rho_max, "B_req": B_req, "ep1_action": ep1_action,
                "ep1_next": ep1_next, "combined": None, "action": "stop_for_owner_review"}
    nxt = [L for L in LADDER if L > B]
    w5_next = nxt[0] if (any(v["breach"] for v in per.values()) and nxt) else None
    cands = [x for x in (ep1_next, w5_next) if x]
    combined = max(cands) if cands else None
    if combined:
        action = "top_up"
    elif ep1_action != "stop_ep1_pass":
        action = ep1_action
    elif any(v["breach"] for v in per.values()):
        action = "report_persisting_breach"
    else:
        action = "stop_ep1_pass"
    return {"per": per, "rho_max": rho_max, "B_req": B_req, "ep1_action": ep1_action,
            "ep1_next": ep1_next, "combined": combined, "action": action}


# ================================================================ deterministic synthetic fixtures
def rank(n, salt):
    """A fixed pseudo-random permutation (multiplicative hashing; no RNG)."""
    order = sorted(range(n), key=lambda k: ((k + 1) * 2654435761 + salt * 40503) % 4294967296)
    pos = [0] * n
    for i, k in enumerate(order):
        pos[k] = i
    return pos


def values(n, a, scale, shift, salt, drift=0.0, up=1.0):
    """dC_k = shift + scale * sign(t)|t|^a * (up if t > 0) + drift * k / n, with
    t = tan(pi (u - 1/2)), u = (rank_k + 1/2)/n. a controls tail heaviness (a = 1: Cauchy)."""
    pos = rank(n, salt)
    out = []
    for k in range(n):
        t = math.tan(math.pi * ((pos[k] + 0.5) / n - 0.5))
        out.append(shift + scale * math.copysign(abs(t) ** a, t) * (up if t > 0 else 1.0)
                   + drift * k / n)
    return out


def record(r, k, v, status="ok"):
    d = {"run_id": r, "k": k, "status": status, "B_in_attempted": 0, "inner": []}
    if status == "ok":
        d.update(dC_corrected=v, dC_app=v)
    return d


def give_bounding_width(rec, width):
    """3 paired-complete inner rows with differences (0, 0, D) out of B_in = 4: width = D/4."""
    def row(dd):
        return {"paired_complete": True,
                "fits": {"M1": {"C_fit": 0.6, "C_eval": 0.6}, "M2": {"C_fit": 0.6 + dd, "C_eval": 0.6}}}
    rec["B_in_attempted"] = 4
    rec["inner"] = [row(0.0), row(0.0), row(4 * width)]


def build(n, specs, bad=None):
    by = {}
    for r, sp in enumerate(specs):
        x = values(n, **sp)
        by[r] = {k: record(r, k, x[k], (bad or {}).get((r, k), "ok")) for k in range(n)}
    return by


def cfg(beyond="top_up_to_max_rung_then_reassess"):
    return CorrectedConfig(approved_test_config(ROOT, "scratch_unused", w5=W5A,
                                                ladder={"rungs": LADDER, "beyond_ladder_rule": beyond}), ROOT)


def close(a, b):
    return math.isclose(a, b, rel_tol=RTOL, abs_tol=0.0)


def light(r, n_drift=0.0, shift=0.0):
    return {"a": 0.3, "scale": 0.010, "shift": shift, "salt": r + 1, "drift": n_drift}


# ================================================================ helper sanity (hand-calculated)
def test_reference_helper_hand_values():
    assert q7([1, 2, 3, 4], 0.25) == 1.75                    # h = 0.75 -> 1 + 0.75*(2-1)
    assert q7([1, 2, 3, 4], 0.975) == 3.925                  # h = 2.925 -> 3 + 0.925
    assert math.isclose(sd([1, 2, 3, 4]), math.sqrt(5 / 3))  # var = 5/3 (ddof = 1)
    assert math.isclose(kde([0.0], 0.0, 1.0), 1 / math.sqrt(2 * math.pi))
    # Silverman for [1..4]: s = 1.2910, IQR = 3.25 - 1.75 = 1.5 -> 1.5/1.34 = 1.1194 (smaller)
    assert math.isclose(silverman([1, 2, 3, 4]), 0.9 * (1.5 / 1.34) * 4 ** -0.2)


# ================================================================ F1: clean top-up to 1000
@pytest.fixture(scope="module")
def F1():
    """Run 0: asymmetric moderately heavy tails (a = 0.43, upper tail x1.3); apparent fit
    failed at k = 10 and 400, inner shortfall at k = 250 (no dC_corrected); bounding widths:
    k = 100 gets 1.3 * t65 (EXCLUDED under MIN-MIN), k = 200 gets 0.5 * t65 (kept).
    Runs 1-4: light tails (a = 0.3). Hand derivation: rho_max ~ 0.064 at 0.5h in run 0 ->
    B_req = 500 (rho/0.05)^2 ~ 830 -> smallest rung >= B_req and > 500 is 1000; using the
    Silverman bandwidth alone (Option A) would give ~712 -> 750, so Option B matters."""
    bad = {(0, 10): "unreliable_apparent_failed", (0, 400): "unreliable_apparent_failed",
           (0, 250): "unreliable_inner_shortfall"}
    specs = [{"a": 0.43, "scale": 0.010, "shift": 0.0, "salt": 1, "up": 1.3}] + \
        [light(r, shift=0.002 * r) for r in range(1, 5)]
    by = build(500, specs, bad)
    P0 = [x["dC_corrected"] for x in by[0].values() if x["status"] == "ok"]
    se, _ = endpoint_se(P0)
    t65 = min(min(v.values()) for v in se.values())
    give_bounding_width(by[0][100], 1.3 * t65)
    give_bounding_width(by[0][200], 0.5 * t65)
    return by, se, t65


def test_F1_section_6_5_threshold_population_and_exclusions(F1):
    by, se, t65 = F1
    # fixture sanity: 1.3 * t65 discriminates the approved MIN-MIN rule from
    #   max over endpoints of the bandwidth-minimum, and from the bandwidth-MAX (EP-1) rule
    alt_endpoint_max = max(min(v.values()) for v in se.values())
    alt_bw_max = min(max(v.values()) for v in se.values())
    assert 1.3 * t65 < alt_endpoint_max and 1.3 * t65 < alt_bw_max
    exp = expected(by, 500)
    d = compute_sizing_decision(cfg(), 500, by)
    p0 = d["per_run"]["0"]
    # population: 497 "ok" records (k = 10, 250, 400 have no dC_corrected)
    assert len(exp["per"][0]["P"]) == 497
    for e in ("lower", "upper"):                               # MIN over bandwidths
        assert close(p0["endpoint_mcse_min"][e], exp["per"][0]["mcse"][e])
        assert close(p0["endpoint_mcse_min"][e], min(se[e].values()))
    assert close(min(p0["endpoint_mcse_min"].values()), t65)  # MIN of the two endpoints
    # threshold computed BEFORE the bounding exclusion: with k = 100 removed it would differ
    P_after = [x["dC_corrected"] for k, x in by[0].items() if x["status"] == "ok" and k != 100]
    se_after, _ = endpoint_se(P_after)
    assert not close(min(min(v.values()) for v in se_after.values()), t65)
    # exact excluded positions: status 10, 250, 400 + bounding 100 (200 kept)
    assert exp["per"][0]["excl"] == [10, 100, 250, 400] and p0["n_excluded_6_5"] == 4
    # n_usable by hand: k < 300 loses {10, 100, 250} -> 297; k < 400 loses the same three
    # (400 is not < 400) -> 397; all 500 lose four -> 496
    assert p0["w5"]["n_usable"] == {"300": 297, "400": 397, "500": 496}
    assert {str(g): n for g, n in exp["per"][0]["n_usable"].items()} == p0["w5"]["n_usable"]


def test_F1_ep1_option_b_rho_breq_and_rung(F1):
    by, _, _ = F1
    exp = expected(by, 500)
    d = compute_sizing_decision(cfg(), 500, by)
    for row in d["ep1"]["rows"]:                               # EP-1 MAX over bandwidths
        e_se = exp["per"][row["run"]]["se_U"][row["endpoint"]]
        for m in BW:
            assert close(row["SE_by_bw"][str(m)], e_se[m])
        assert row["governing_SE"] == max(row["SE_by_bw"].values())
        assert close(row["ratio"], exp["per"][row["run"]]["ratio"][row["endpoint"]])
    assert close(d["ep1"]["rho_max"], exp["rho_max"])          # 7B: max over 5 runs x 2 ends
    assert d["ep1"]["rho_max"] == max(r["ratio"] for r in d["ep1"]["rows"])
    B_req = d["ep1_decision"]["B_req"]
    assert B_req.hex() == (500 * (d["ep1"]["rho_max"] / FRACTION) ** 2).hex()
    assert close(B_req, exp["B_req"]) and 750 < B_req <= 1000
    # Option A (Silverman h only) would select 750: the MAX rule is what selects 1000
    rho_A = max(exp["per"][r]["se_U"][e][1.0] / exp["per"][r]["W"] for r in exp["per"] for e in ("lower", "upper"))
    assert 500 * (rho_A / FRACTION) ** 2 < 750
    assert d["ep1_decision"]["B_next"] == exp["ep1_next"] == 1000
    assert d["combined_B_next"] == exp["combined"] == 1000 and d["action"] == "top_up"
    assert [d["per_run"][str(r)]["w5"]["breach"] for r in range(5)] == \
        [exp["per"][r]["breach"] for r in range(5)]


# ================================================================ F2: B_req beyond the ladder
def test_F2_beyond_ladder_selects_3000_directly():
    by = build(500, [{"a": 1.0, "scale": 0.010, "shift": 0.0, "salt": 1}] + [light(r) for r in range(1, 5)])
    exp = expected(by, 500)
    d = compute_sizing_decision(cfg(), 500, by)
    assert exp["B_req"] > 3000 and close(d["ep1_decision"]["B_req"], exp["B_req"])
    assert d["ep1_decision"]["B_next"] == 3000 and d["combined_B_next"] == 3000 == exp["combined"]
    assert d["action"] == "top_up"


# ================================================================ F3: the top rung (B_cur = 3000)
def test_F3_ep1_failure_at_3000_is_imprecise():
    by = build(3000, [{"a": 1.0, "scale": 0.010, "shift": 0.0, "salt": 1}] + [light(r) for r in range(1, 5)])
    exp = expected(by, 3000)
    d = compute_sizing_decision(cfg(), 3000, by)
    assert exp["rho_max"] > FRACTION and close(d["ep1"]["rho_max"], exp["rho_max"])
    assert d["action"] == exp["action"] == "stop_label_imprecise" and d["combined_B_next"] is None
    assert d["per_run"]["0"]["w5"]["points"] == [2000, 3000]


def test_F3_w5_breach_at_3000_is_persisting_breach_not_imprecise():
    by = build(3000, [light(r) for r in range(4)] + [light(4, n_drift=0.012)])
    exp = expected(by, 3000)
    d = compute_sizing_decision(cfg(), 3000, by)
    assert d["ep1"]["pass"] is True and exp["rho_max"] <= FRACTION
    assert [d["per_run"][str(r)]["w5"]["breach"] for r in range(5)] == [False] * 4 + [True]
    assert d["action"] == exp["action"] == "report_persisting_breach"
    assert d["action"] != "stop_label_imprecise" and d["combined_B_next"] is None


# ================================================================ F4: H1 -> W5 UNEVALUABLE -> owner review
def test_F4_h1_makes_w5_unevaluable_and_stops_for_owner_review():
    bad = {(3, 19 * j): "unreliable_apparent_failed" for j in range(26)}   # 26/500 = 5.2 % > 5 %
    by = build(500, [{"a": 0.48, "scale": 0.010, "shift": 0.0, "salt": 1}] + [light(r) for r in range(1, 5)], bad)
    exp = expected(by, 500)
    d = compute_sizing_decision(cfg(), 500, by)
    p3 = d["per_run"]["3"]
    assert exp["per"][3]["h1"] is True and p3["H1"] == "UNRELIABLE"
    assert p3["w5"]["status"] == "UNEVALUABLE" and p3["w5_decision"]["action"] == "stop_for_owner_review"
    assert d["ep1_decision"]["action"] == exp["ep1_action"] == "top_up"     # EP-1 alone: top-up
    assert d["action"] == exp["action"] == "stop_for_owner_review" and d["combined_B_next"] is None


# ================================================================ F5: W5-A breach drives the top-up
def test_F5_w5_breach_positions_no_backfill_and_next_rung():
    bad = {(2, 350): "unreliable_apparent_failed", (2, 420): "unreliable_inner_shortfall"}
    by = build(500, [light(0), light(1), light(2, n_drift=0.015), light(3), light(4)], bad)
    exp = expected(by, 500)
    d = compute_sizing_decision(cfg(), 500, by)
    w = d["per_run"]["2"]["w5"]
    assert d["ep1"]["pass"] is True
    assert w["points"] == [300, 400, 500]
    assert w["n_usable"] == {"300": 300, "400": 399, "500": 498}       # positions kept
    for g in (300, 400, 500):
        lo, hi = exp["per"][2]["trace"][g]
        assert close(w["trace"][str(g)][0], lo) and close(w["trace"][str(g)][1], hi)
    # no back-fill: the first 399 USABLE values would include k = 400, which is not < 400
    usable_seq = [by[2][k]["dC_corrected"] for k in range(500) if by[2][k]["status"] == "ok"]
    backfilled = (q7(usable_seq[:400], 0.025), q7(usable_seq[:400], 0.975))
    assert not (close(w["trace"]["400"][0], backfilled[0]) and close(w["trace"]["400"][1], backfilled[1]))
    assert [d["per_run"][str(r)]["w5"]["breach"] for r in range(5)] == [False, False, True, False, False]
    assert exp["per"][2]["breach"] is True and close(w["bound"], exp["per"][2]["bound"])
    assert d["combined_B_next"] == exp["combined"] == 750 and d["action"] == "top_up"

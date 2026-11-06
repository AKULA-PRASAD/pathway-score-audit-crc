"""Owner Decision 6B (2026-11-06): W5 exclusion rule 'usable_within_k_positions_v1',
and the p9c-config-2 schema that requires it. Deterministic fixtures, explicit values.
The final-project stats_core.py is loaded read-only as the historical reference."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from src.phase9_corrected import statistics as st
from src.phase9_corrected.config import (CONFIG_SCHEMA, CorrectedConfig,
                                         UnresolvedDecisionError)
from tests.p9c_support import approved_test_config

ROOT = Path(__file__).resolve().parents[1]
RULE = "usable_within_k_positions_v1"
LADDER = {"rungs": [750, 1000, 1500, 2000, 3000], "beyond_ladder_rule": "stop_and_label_imprecise"}


@pytest.fixture(scope="module")
def sc():
    spec = importlib.util.spec_from_file_location(
        "fp_stats_core_ro_w5", ROOT / "reports/phase9/final_project/analysis/stats_core.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def x500(seed=11, n=500):
    return np.random.default_rng(seed).normal(0.0, 0.014, n)


def as_map(x, excluded=()):
    return {k: (None if k in excluded else float(v)) for k, v in enumerate(x)}


def width(v):
    lo, hi = st.percentile_ci(v)
    return hi - lo


# 1. zero exclusions: bit-identical to the legacy path and to stats_core
def test_zero_exclusions_bit_identical_to_legacy_and_stats_core(sc):
    for seed in (11, 12, 13):
        x = x500(seed)
        legacy = st.w5_evaluate(as_map(x), 500, width(x))
        new = st.w5_evaluate(as_map(x), 500, exclusion_rule=RULE, h1_triggered=False)
        assert new["status"] == "EVALUATED"
        for f in ("points", "trace", "spread", "bound", "breach"):
            assert new[f] == legacy[f], f
        tr = sc.running_trace(x)
        fin = [tr[B] for B in sc.MC1_FINAL]
        assert new["spread"]["lower"] == max(t[0] for t in fin) - min(t[0] for t in fin)
        assert new["spread"]["upper"] == max(t[1] for t in fin) - min(t[1] for t in fin)
        lo, hi = sc.percentile_ci(x)
        assert new["bound"] == 0.05 * (hi - lo)


# 2-4. excluded positions contribute nothing; positions kept; n_usable exact
EXCL = (10, 350, 420)


def test_excluded_positions_do_not_contribute_and_n_usable_exact():
    x = x500()
    w = st.w5_evaluate(as_map(x, EXCL), 500, exclusion_rule=RULE, h1_triggered=False)
    assert w["n_usable"] == {"300": 299, "400": 398, "500": 497}
    assert w["n_usable_full"] == 497
    for B in (300, 400, 500):
        expect = [x[k] for k in range(B) if k not in EXCL]
        assert w["trace"][str(B)] == list(st.percentile_ci(np.array(expect)))
    # the excluded values cannot matter: replace them with extreme numbers -> same result
    y = x.copy()
    y[list(EXCL)] = 1e6
    assert st.w5_evaluate(as_map(y, EXCL), 500, exclusion_rule=RULE, h1_triggered=False) == w


def test_positions_preserved_no_backfill():
    x = x500()
    base = st.w5_evaluate(as_map(x), 500, exclusion_rule=RULE, h1_triggered=False)
    w = st.w5_evaluate(as_map(x, (405, 470)), 500, exclusion_rule=RULE, h1_triggered=False)
    assert w["trace"]["300"] == base["trace"]["300"]          # k >= 400: earlier points intact
    assert w["trace"]["400"] == base["trace"]["400"]
    assert w["n_usable"] == {"300": 300, "400": 400, "500": 498}
    # an exclusion below 300 shrinks B=300 to 299 values -- it is NOT refilled with k=300
    v = st.w5_evaluate(as_map(x, (0,)), 500, exclusion_rule=RULE, h1_triggered=False)
    assert v["n_usable"]["300"] == 299
    assert v["trace"]["300"] == list(st.percentile_ci(x[1:300]))
    assert v["trace"]["300"] != list(st.percentile_ci(x[1:301]))


# 5-6. bound from the usable full-B_out population; conflicting width rejected
def test_bound_uses_usable_full_population_and_conflicting_width_rejected():
    x = x500()
    usable = np.array([x[k] for k in range(500) if k not in EXCL])
    w = st.w5_evaluate(as_map(x, EXCL), 500, exclusion_rule=RULE, h1_triggered=False)
    assert w["width"] == width(usable) and w["bound"] == 0.05 * width(usable)
    assert w["width"] != width(x)                              # differs from all-500 width
    ok = st.w5_evaluate(as_map(x, EXCL), 500, width(usable), exclusion_rule=RULE,
                        h1_triggered=False)                     # consistent width accepted
    assert ok == w
    with pytest.raises(ValueError, match="conflicts"):
        st.w5_evaluate(as_map(x, EXCL), 500, width(x), exclusion_rule=RULE, h1_triggered=False)


# 7-8. UNEVALUABLE triggers
def test_h1_triggered_is_unevaluable():
    w = st.w5_evaluate(as_map(x500(), range(0, 500, 19)), 500, exclusion_rule=RULE,
                       h1_triggered=True)
    assert w["status"] == "UNEVALUABLE" and w["breach"] is None
    assert any("H1" in r for r in w["reasons"])
    assert "spread" not in w and "trace" not in w                # no substitute calculation
    with pytest.raises(ValueError, match="h1_triggered"):
        st.w5_evaluate(as_map(x500()), 500, exclusion_rule=RULE)


def test_fewer_than_two_usable_is_unevaluable():
    x = x500()
    excl = set(range(299))                                     # B=300 keeps only k=299
    w = st.w5_evaluate(as_map(x, excl), 500, exclusion_rule=RULE, h1_triggered=False)
    assert w["n_usable"]["300"] == 1
    assert w["status"] == "UNEVALUABLE" and w["breach"] is None
    assert any("300" in r for r in w["reasons"])
    two = st.w5_evaluate(as_map(x, set(range(298))), 500, exclusion_rule=RULE,
                         h1_triggered=False)                    # exactly two: evaluable
    assert two["n_usable"]["300"] == 2 and two["status"] == "EVALUATED"


# 9-10. decisions and combined sizing
def test_unevaluable_stops_for_owner_review_under_both_actions():
    u = {"status": "UNEVALUABLE", "breach": None}
    for action, rule in (("W5-A", "next_rung_above_current"), ("W5-B", None)):
        assert st.w5_decision(u, action, rule, 500, LADDER) == \
            {"action": "stop_for_owner_review", "B_next": None}
    with pytest.raises(UnresolvedDecisionError):                # action still required
        st.w5_decision(u, None, None, 500, LADDER)


def test_combined_sizing_refuses_when_owner_review_required():
    stop = {"action": "stop_for_owner_review", "B_next": None}
    with pytest.raises(st.OwnerReviewRequired):
        st.combined_next_size({"action": "top_up", "B_next": 750}, stop)
    with pytest.raises(st.OwnerReviewRequired):
        st.combined_next_size(stop, {"action": "top_up", "B_next": 1000})
    assert st.combined_next_size({"action": "top_up", "B_next": 750},
                                 {"action": "top_up", "B_next": 1000}) == 1000


# 11. default behaviour unchanged and still gated
def test_omitting_exclusion_rule_keeps_gated_behaviour():
    x = as_map(x500(), (7,))
    with pytest.raises(UnresolvedDecisionError):
        st.w5_evaluate(x, 500, 0.05)
    with pytest.raises(ValueError, match="width is required"):
        st.w5_evaluate(as_map(x500()), 500)
    with pytest.raises(UnresolvedDecisionError, match="unknown"):
        st.w5_evaluate(x, 500, exclusion_rule="first_B_usable", h1_triggered=False)


# 12-14. configuration schema p9c-config-2
def problems(c):
    return CorrectedConfig(c, ROOT).validate_structure()


@pytest.mark.parametrize("val", [None, "", "first_B_usable", "USABLE_WITHIN_K_POSITIONS_V1"])
def test_config_rejects_missing_or_unknown_exclusion_rule(val):
    c = approved_test_config(ROOT, "scratch_p9c_test_out")
    if val is None:
        del c["w5"]["exclusion_rule"]
    else:
        c["w5"]["exclusion_rule"] = val
    assert any("w5.exclusion_rule" in p for p in problems(c))


def test_config_accepts_approved_rule_under_schema_2():
    assert CONFIG_SCHEMA == "p9c-config-2"
    c = approved_test_config(ROOT, "scratch_p9c_test_out")
    assert c["schema"] == "p9c-config-2" and c["w5"]["exclusion_rule"] == RULE
    assert problems(c) == []
    CorrectedConfig(c, ROOT).validate_for_execution()
    c["w5"] = {"action": "W5-A", "corrective_rule": "next_rung_above_current",
               "definition": "adjudication_D5_as_proposed", "exclusion_rule": RULE}
    assert problems(c) == []


@pytest.mark.parametrize("schema", ["p9c-config-1", "p9c-config-3", None])
def test_config_rejects_other_schema_versions(schema):
    c = approved_test_config(ROOT, "scratch_p9c_test_out", schema=schema)
    assert any("schema must be 'p9c-config-2'" in p for p in problems(c))


def test_blank_template_still_has_no_w5_and_is_not_executable():
    from src.phase9_corrected.config import template_config
    t = template_config()
    assert t["schema"] == "p9c-config-2" and t["w5"] is None
    assert CorrectedConfig(json.loads(json.dumps(t)), ROOT).missing_owner_items()

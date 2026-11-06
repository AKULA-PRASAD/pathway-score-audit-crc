"""RNG architecture of the corrected Phase 9 pipeline (no bootstrap is run).
T = 7 below is a TEST VALUE, not an approved domain tag."""
import numpy as np
import pytest

from src.phase9_corrected.config import PLAN_BASE, ConfigurationError, UINT32_MAX
from src.phase9_corrected.rng import (HistoricalAsRunStreams, PlanStreams, assembled_words,
                                      coerce_words, enumerate_stream_states, state_fingerprint,
                                      structural_domain_report, verify_numpy_seed_behaviour)

B = PLAN_BASE
T_TEST = 7
S = np.random.SeedSequence


def test_numpy_behaviour_verified_for_pinned_version():
    chk = verify_numpy_seed_behaviour("2.5.4", B, T_TEST)
    assert chk["numpy_version"]["pass"]
    assert chk["coercion_model"]["pass"]
    assert chk["zero_padding_hazard_present"]["pass"]       # why T must be non-zero
    assert chk["T_stream_distinct_from_outer"]["pass"]
    assert chk["explicit_spawn_key_equals_spawn_list"]["pass"]
    assert chk["spawn_is_stateful"]["pass"]
    assert chk["all_pass"]


def test_version_mismatch_is_reported_not_hidden():
    assert not verify_numpy_seed_behaviour("0.0.0", B, T_TEST)["all_pass"]


def test_zero_T_would_alias_outer_stream_and_is_rejected():
    # the hazard is real ...
    assert state_fingerprint(S([B, 3, 0], spawn_key=(11,))) == state_fingerprint(S([B, 3], spawn_key=(11,)))
    # ... so the construction refuses it
    for bad in (0, -1, UINT32_MAX + 1, True, None, 1.5):
        with pytest.raises(ConfigurationError):
            PlanStreams(B, 0, "independent", bad)


def test_index_above_uint32_aliases_and_is_rejected():
    # NumPy splits ints >= 2**32 into two words: inner (k, 2**32 + a) == reserve (k, a, 1)
    a, k = 5, 9
    assert state_fingerprint(S([B, 0], spawn_key=(k, 2 ** 32 + a))) == \
        state_fingerprint(S([B, 0], spawn_key=(k, a, 1)))
    s = PlanStreams(B, 0, "independent", T_TEST)
    with pytest.raises(ConfigurationError):
        s.inner(k, 2 ** 32 + a)
    with pytest.raises(ConfigurationError):
        s.apparent_cv(2 ** 32)


def test_plan_streams_equal_plan_spawn_tree():
    """Explicit construction == plan section 8.2 spawn list (lines 361-364)."""
    for r in (0, 4):
        s = PlanStreams(B, r, "independent", T_TEST)
        root = S([B, r])
        outer = root.spawn(30)
        for k in (0, 17, 29):
            assert state_fingerprint(s.outer(k)) == state_fingerprint(outer[k])
            inner = S([B, r], spawn_key=(k,)).spawn(100)
            for b in (0, 63, 99):
                assert state_fingerprint(s.inner(k, b)) == state_fingerprint(inner[b])
                cv = S([B, r], spawn_key=(k, b)).spawn(2)
                assert state_fingerprint(s.inner_cv(k, b)) == state_fingerprint(cv[0])
                assert state_fingerprint(s.inner_reserve(k, b)) == state_fingerprint(cv[1])


def test_historical_streams_reproduce_frozen_driver_construction_and_deviation():
    h = HistoricalAsRunStreams(B, 2, 100)
    p = PlanStreams(B, 2, "independent", T_TEST)
    for k in (0, 250):
        ss_outer = S(entropy=[B, 2]).spawn(k + 1)[k]                 # driver lines 96-97
        assert state_fingerprint(h.outer(k)) == state_fingerprint(ss_outer)
        assert state_fingerprint(h.apparent_cv(k)) == state_fingerprint(ss_outer)   # lines 116-117
        ss_inner = ss_outer.spawn(100)                                # line 130
        for b in (0, 99):
            assert state_fingerprint(h.inner(k, b)) == state_fingerprint(ss_inner[b])
            assert state_fingerprint(h.inner_cv(k, b)) == state_fingerprint(ss_inner[b])  # deviation
            # corrected: resample stream unchanged, CV stream changed
            assert state_fingerprint(p.inner(k, b)) == state_fingerprint(ss_inner[b])
            assert state_fingerprint(p.inner_cv(k, b)) != state_fingerprint(ss_inner[b])
        assert state_fingerprint(p.outer(k)) == state_fingerprint(ss_outer)
        assert state_fingerprint(p.apparent_cv(k)) != state_fingerprint(ss_outer)


def test_as_run_convention_uses_outer_state_for_apparent():
    p = PlanStreams(B, 1, "as_run")
    assert state_fingerprint(p.apparent_cv(4)) == state_fingerprint(p.outer(4))
    with pytest.raises(ConfigurationError):
        PlanStreams(B, 1, "as_run", 7)
    with pytest.raises(ConfigurationError):
        PlanStreams(B, 1, None)


def test_fresh_generators_are_independent_objects_with_identical_start():
    p = PlanStreams(B, 0, "independent", T_TEST)
    g1, g2 = p.generator(p.inner_cv(3, 4)), p.generator(p.inner_cv(3, 4))
    a = g1.permutation(50)
    assert g1 is not g2
    assert np.array_equal(a, g2.permutation(50))       # consuming g1 did not advance g2


def test_coercion_model_matches_numpy():
    from numpy.random.bit_generator import _coerce_to_uint32_array as c
    for x in (0, 1, B, UINT32_MAX, UINT32_MAX + 1, 2 ** 70 + 3):
        assert tuple(int(w) for w in c(x)) == coerce_words(x)


def test_assembled_words_model_predicts_numpy_equalities():
    """Model-predicted equal word arrays <=> equal states, on a grid of near-miss pairs."""
    cases = [(([B, 0], (1,)), ([B, 0, 0], (1,))),
             (([B, 0], ()), ([B, 0, 0, 0], ())),
             (([B, 0, T_TEST], (1,)), ([B, 0], (1,))),
             (([B, 0], (1, 2)), ([B, 0], (1, 2, 0))),
             (([B, 0], (1, 2 ** 32 + 2)), ([B, 0], (1, 2, 1)))]
    for (e1, k1), (e2, k2) in cases:
        same_model = assembled_words(e1, k1) == assembled_words(e2, k2)
        same_state = state_fingerprint(S(e1, spawn_key=k1)) == state_fingerprint(S(e2, spawn_key=k2))
        assert same_model == same_state


def test_structural_report_requires_numpy_verification_and_nonzero_T():
    assert structural_domain_report(B, range(5), T_TEST, "independent", True)["holds"]
    assert not structural_domain_report(B, range(5), T_TEST, "independent", False)["holds"]
    assert not structural_domain_report(B, range(5), 0, "independent", True)["holds"]


def test_structural_word_patterns_disjoint_on_grid():
    words = {}
    for r in range(5):
        for k in range(6):
            fam = {"outer": assembled_words([B, r], (k,)),
                   "apparent": assembled_words([B, r, T_TEST], (k,))}
            for b in range(4):
                fam[f"inner{b}"] = assembled_words([B, r], (k, b))
                fam[f"cv{b}"] = assembled_words([B, r], (k, b, 0))
                fam[f"res{b}"] = assembled_words([B, r], (k, b, 1))
            for name, w in fam.items():
                assert w not in words, (name, r, k, words.get(w))
                words[w] = (name, r, k)


def test_enumeration_requires_declared_B_in_max():
    p = PlanStreams(B, 0, "independent", T_TEST)
    with pytest.raises(ConfigurationError):
        enumerate_stream_states(p, [0], 3, None)


def test_enumeration_small_domain_no_collision_and_detects_planted_collision():
    p = PlanStreams(B, 0, "independent", T_TEST)
    res = enumerate_stream_states(p, [0, 1, 2, 3, 4], 12, 6)
    assert res["pass"] and res["n_states"] == res["n_distinct"] == 5 * 12 * (2 + 3 * 6) + 5
    assert res["claim"].startswith("no duplicate initial state within the enumerated domain")

    class Aliased(PlanStreams):                      # a deliberately broken construction
        def apparent_cv(self, k):
            return self.outer(k)
    bad = Aliased(B, 0, "independent", T_TEST)
    assert not enumerate_stream_states(bad, [0], 3, 2, include_primary_seeds=False)["pass"]


def test_primary_phase9_seeds_distinct_from_all_corrected_streams():
    p = PlanStreams(B, 0, "independent", T_TEST)
    res = enumerate_stream_states(p, [0, 1, 2, 3, 4], 3, 3, include_primary_seeds=True)
    assert res["pass"]

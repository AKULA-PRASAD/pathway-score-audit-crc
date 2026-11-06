"""Resampling module and V1 outer/inner evidence (reads frozen inputs and historical
checkpoints READ-ONLY; draws resample indices only -- no model is fitted)."""
from pathlib import Path

import numpy as np
import pytest

from src.phase9_corrected.config import PLAN_BASE
from src.phase9_corrected.data import load_endpoint_only
from src.phase9_corrected.resampling import (INNER_INDEX_TAG, OUTER_MULT_TAG, draw_inner,
                                             draw_outer, legacy_multiplicity_sha16,
                                             strat_boot, tagged_sha256)
from src.phase9_corrected.rng import HistoricalAsRunStreams, PlanStreams
from src.phase9_corrected.verification import load_frozen_driver, load_historical_run

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def frozen():
    return load_frozen_driver(ROOT)


@pytest.fixture(scope="module")
def endpoint():
    return load_endpoint_only(ROOT)


def test_strat_boot_identical_to_frozen_driver(frozen):
    e = (np.random.default_rng(1).random(97) < 0.3).astype(float)
    for s in range(20):
        a = strat_boot(e, np.random.default_rng(s))
        b = frozen.strat_boot(e, np.random.default_rng(s))
        assert np.array_equal(a, b)


def test_outer_evidence_properties(endpoint):
    ids, e = endpoint
    o = draw_outer(e, np.random.default_rng(3))
    assert o.multiplicity.sum() == len(e) and o.events == int(e.sum())
    assert o.censored == int((e == 0).sum())
    assert o.n_distinct == int((o.multiplicity > 0).sum())
    assert o.legacy_multiplicity_sha16 == legacy_multiplicity_sha16(o.indices, len(e))
    assert o.multiplicity_sha256 == tagged_sha256(OUTER_MULT_TAG, o.multiplicity)


def test_inner_evidence_properties_and_tag_separation(endpoint):
    ids, e = endpoint
    o = draw_outer(e, np.random.default_rng(3))
    eO = e[o.indices]
    J = draw_inner(eO, o.indices, len(e), np.random.default_rng(4))
    assert J.events == int(eO.sum())                         # event stratification preserved
    assert np.array_equal(J.original, o.indices[J.positions])
    assert J.index_sha256 == tagged_sha256(INNER_INDEX_TAG, J.positions)
    assert J.index_sha256 != tagged_sha256(OUTER_MULT_TAG, J.positions)   # domain-separated
    J2 = draw_inner(eO, o.indices, len(e), np.random.default_rng(4))
    assert J2.evidence() == J.evidence()
    J3 = draw_inner(eO, o.indices, len(e), np.random.default_rng(5))
    assert J3.index_sha256 != J.index_sha256


def test_dataset_shape(endpoint):
    ids, e = endpoint
    assert len(ids) == 452 and int(e.sum()) == 136        # plan line 139 / 143


@pytest.mark.parametrize("run_id", [0, 1, 2, 3, 4])
def test_V1_outer_corrected_streams_reproduce_all_historical_hashes(endpoint, run_id):
    """Every historical outer multiplicity_sha (5 x 500) is reproduced by the corrected
    outer stream: the outer resample is unchanged by the correction."""
    ids, e = endpoint
    hist = load_historical_run(ROOT, run_id)
    assert sorted(hist) == list(range(500))
    p = PlanStreams(PLAN_BASE, run_id, "as_run")
    bad = [k for k in range(500)
           if draw_outer(e, p.generator(p.outer(k))).legacy_multiplicity_sha16
           != hist[k]["multiplicity_sha"]]
    assert bad == []


def test_V1_inner_corrected_path_equals_frozen_path(endpoint, frozen):
    """Inner CODE-PATH identity (no historical inner evidence exists)."""
    ids, e = endpoint
    for run_id in (0, 3):
        p, h = PlanStreams(PLAN_BASE, run_id, "as_run"), HistoricalAsRunStreams(PLAN_BASE, run_id, 100)
        for k in (0, 25, 499):
            o = draw_outer(e, p.generator(p.outer(k)))
            I_f = frozen.strat_boot(e, np.random.default_rng(h.outer(k)))
            assert np.array_equal(o.indices, I_f)
            eO = e[o.indices]
            for b in range(100):
                J = draw_inner(eO, o.indices, len(e), p.generator(p.inner(k, b)))
                Jf = frozen.strat_boot(eO, np.random.default_rng(h.inner(k, b)))
                assert np.array_equal(J.positions, Jf)

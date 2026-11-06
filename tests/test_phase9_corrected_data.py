"""Corrected data loader equals the frozen driver's init_worker (read-only)."""
from pathlib import Path

import numpy as np

from src.phase9_corrected.data import load_context, load_endpoint_only
from src.phase9_corrected.verification import load_frozen_driver

ROOT = Path(__file__).resolve().parents[1]


def test_load_context_equals_frozen_init_worker():
    drv = load_frozen_driver(ROOT)
    drv._G.clear()
    drv.init_worker()
    ctx = load_context(ROOT)
    assert ctx.ids == drv._G["ids"]
    assert np.array_equal(ctx.t, drv._G["t"]) and np.array_equal(ctx.e, drv._G["e"])
    assert ctx.X.equals(drv._G["X"]) and ctx.cov.equals(drv._G["cov"])
    assert ctx.panel == drv._G["panel"]
    ids, e = load_endpoint_only(ROOT)
    assert ids == ctx.ids and np.array_equal(e, ctx.e)

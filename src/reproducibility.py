"""Reproducibility helpers: seeding and environment capture.

Important distinction maintained throughout this project:

  * COMPUTATIONAL reproducibility -- the same inputs and code produce
    bit-identical outputs on this machine. That is what this module and the
    Phase 0 test suite establish.

  * BIOLOGICAL / NORMALIZATION validity -- whether the deposited values are a
    sound basis for the scientific question. GSE39582 was deposited as
    "RMA normalisation + ComBat" applied jointly across all 585 arrays, which
    has consequences for cross-cohort transport. Nothing in Phase 0 addresses
    that; it is taken up in Phases 1 and 4.

A pipeline can be perfectly reproducible and still scientifically invalid.
These two properties are tracked separately and must never be conflated.
"""

from __future__ import annotations

import os
import platform
import random
import sys
from datetime import datetime, timezone
from importlib import import_module

DEFAULT_SEED = 20261102

_TRACKED_PACKAGES = (
    "numpy", "pandas", "scipy", "sklearn", "lifelines", "pyarrow", "pytest",
)


def set_seeds(seed: int = DEFAULT_SEED) -> int:
    """Seed Python and NumPy RNGs. Returns the seed for provenance logging."""
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    try:
        import numpy as np

        np.random.seed(seed)
    except ImportError:
        pass
    return seed


def package_versions() -> dict[str, str]:
    """Version of each tracked package, or 'not installed'."""
    out: dict[str, str] = {}
    for name in _TRACKED_PACKAGES:
        try:
            out[name] = getattr(import_module(name), "__version__", "unknown")
        except ImportError:
            out[name] = "not installed"
    return out


def environment_report() -> dict:
    """Machine-readable snapshot of the runtime, for embedding in manifests."""
    return {
        "captured_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version.split()[0],
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "packages": package_versions(),
    }

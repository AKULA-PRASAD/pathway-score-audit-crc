"""Event-stratified resampling and canonical resampling evidence.

`strat_boot` reproduces the frozen driver's function (scripts/40_phase9_dc_interval.py
lines 61-66) exactly; equivalence is tested against the frozen file.

Canonical serialisation for hashes (p9c-1):
  sha256( TAG + b"\\0" + numpy.asarray(x, dtype="<i8").tobytes() )
with TAG one of OUTER_MULT_TAG, INNER_INDEX_TAG, INNER_MULT_TAG. The legacy 16-hex
outer hash (driver lines 103, 106-107) is also computed, byte-for-byte as the
historical driver did, so V1 can compare with the historical checkpoints.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass

import numpy as np

OUTER_MULT_TAG = b"p9c-1/outer-multiplicity"
INNER_INDEX_TAG = b"p9c-1/inner-positions"
INNER_MULT_TAG = b"p9c-1/inner-multiplicity-original-ids"


def strat_boot(ev_mask, rng):
    """Event-stratified patient bootstrap. Returns POSITIONAL indices with repeats.
    [IMPL] identical to the frozen driver (lines 61-66)."""
    idx = np.arange(len(ev_mask))
    out = [rng.choice(idx[ev_mask == lab], size=int((ev_mask == lab).sum()), replace=True)
           for lab in (1.0, 0.0)]
    return np.sort(np.concatenate(out))


def tagged_sha256(tag: bytes, arr) -> str:
    return hashlib.sha256(tag + b"\0" + np.asarray(arr, dtype="<i8").tobytes()).hexdigest()


def legacy_multiplicity_sha16(I, n: int) -> str:
    """Byte-for-byte the historical `multiplicity_sha` (driver lines 103, 106-107)."""
    mult = np.bincount(I, minlength=n).astype(int)
    return hashlib.sha256(mult.tobytes()).hexdigest()[:16]


@dataclass
class OuterResample:
    indices: np.ndarray
    multiplicity: np.ndarray
    multiplicity_sha256: str
    legacy_multiplicity_sha16: str
    n_distinct: int
    events: int
    censored: int

    def evidence(self) -> dict:
        return {"multiplicity": [int(x) for x in self.multiplicity],
                "multiplicity_sha256": self.multiplicity_sha256,
                "legacy_multiplicity_sha16": self.legacy_multiplicity_sha16,
                "n_distinct": self.n_distinct, "events": self.events, "censored": self.censored}


def draw_outer(e: np.ndarray, rng) -> OuterResample:
    I = strat_boot(e, rng)
    n = len(e)
    mult = np.bincount(I, minlength=n)
    eO = e[I]
    return OuterResample(I, mult, tagged_sha256(OUTER_MULT_TAG, mult),
                         legacy_multiplicity_sha16(I, n), int((mult > 0).sum()),
                         int(eO.sum()), int((eO == 0).sum()))


@dataclass
class InnerResample:
    positions: np.ndarray           # J: positional indices into D*_k
    original: np.ndarray            # I[J]: indices into the original cohort
    index_sha256: str
    multiplicity_sha256: str
    n_distinct: int
    events: int
    censored: int

    def evidence(self) -> dict:
        return {"index_sha256": self.index_sha256,
                "multiplicity_sha256": self.multiplicity_sha256,
                "n_distinct": self.n_distinct, "events": self.events, "censored": self.censored}


def draw_inner(eO: np.ndarray, I: np.ndarray, n: int, rng) -> InnerResample:
    J = strat_boot(eO, rng)
    orig = I[J]
    eJ = eO[J]
    return InnerResample(J, orig, tagged_sha256(INNER_INDEX_TAG, J),
                         tagged_sha256(INNER_MULT_TAG, np.bincount(orig, minlength=n)),
                         int(len(np.unique(orig))), int(eJ.sum()), int((eJ == 0).sum()))


def binary_cell_counts(cov, ids, binaries: dict) -> dict:
    """Counts of both levels of each clinical binary in a fitting sample (raw inputs for
    the plan's `min_binary_cell`, line 270, whose exact definition the plan leaves open)."""
    sub = cov.loc[list(ids)]
    return {c: {"pos": int((sub[c].to_numpy() == pos).sum()),
                "neg": int((sub[c].to_numpy() != pos).sum())} for c, (pos, _r) in binaries.items()}

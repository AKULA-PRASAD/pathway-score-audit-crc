"""RNG stream architecture for the corrected rerun.

Stream families (all SeedSequences constructed EXPLICITLY from (entropy, spawn_key),
never through the stateful `spawn()` counter, so a stream is a pure function of its id):

  family                   entropy                spawn_key    class   consumer
  outer(k)                 [BASE, run_id]         (k,)         [PLAN 361-362] outer strat_boot
  inner(k, b)              [BASE, run_id]         (k, b)       [PLAN 363]     inner strat_boot
  inner_cv(k, b)           [BASE, run_id]         (k, b, 0)    [PLAN 364]     inner 5-fold CV
  inner_reserve(k, b)      [BASE, run_id]         (k, b, 1)    [PLAN 364]     never consumed
  apparent_cv(k)  indep.   [BASE, run_id, T]      (k,)         [PROP D2]      apparent CV on D*_k
  apparent_cv(k)  as_run   = outer(k)                          [IMPL 116-117] (verification/variant)

`root.spawn(n)[k]` has spawn_key (k,) when root is fresh, so the explicit construction
equals the plan's spawn list; this equivalence is VERIFIED empirically
(`verify_numpy_seed_behaviour`), not assumed.

A fresh Generator is created for every consuming operation and every model
([PROP]: mirrors the historical per-model pattern; only M2 consumes CV draws).

Guarantees are separated explicitly:
  * structural  -- `structural_domain_report`: input-level distinctness of the assembled
                   entropy words, given NumPy behaviour that must itself be verified;
  * empirical   -- `enumerate_stream_states`: no duplicate PCG64 initial state inside a
                   DECLARED finite domain; nothing is claimed outside that domain;
  * not provided -- distinct initial states do not prove non-overlapping output
                   sequences (NumPy: "very probably non-overlapping").
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Iterable, Optional

import numpy as np

from .config import UINT32_MAX, ConfigurationError

POOL_SIZE = 4                                       # NumPy SeedSequence default
PHASE9_PRIMARY_SEEDS = (20261102, 20261103, 20261104, 20261105, 20261106)  # 38_phase9_execute.py:32
FAMILIES = ("outer", "inner", "inner_cv", "inner_reserve", "apparent_cv")


@dataclass(frozen=True)
class StreamId:
    family: str
    run_id: int
    k: Optional[int] = None
    b: Optional[int] = None

    def canonical(self) -> str:
        return f"{self.family}/r{self.run_id}/k{self.k}/b{self.b}"


def state_fingerprint(ss: np.random.SeedSequence) -> str:
    """sha256 of the PCG64 initial (state, inc) that default_rng(ss) would start from."""
    st = np.random.PCG64(ss).state["state"]
    return hashlib.sha256(f"{st['state']}:{st['inc']}".encode()).hexdigest()


def seedseq_descriptor(ss: np.random.SeedSequence) -> dict:
    ent = ss.entropy
    return {"entropy": list(ent) if isinstance(ent, (list, tuple)) else ent,
            "spawn_key": list(ss.spawn_key)}


class PlanStreams:
    """Streams for the corrected (plan-conforming) rerun."""

    kind = "plan"

    def __init__(self, base: int, run_id: int, convention: str, T: Optional[int] = None):
        if convention not in ("independent", "as_run"):
            raise ConfigurationError("apparent-CV convention must be explicit")
        if convention == "independent" and not (isinstance(T, int) and not isinstance(T, bool)
                                                and 1 <= T <= UINT32_MAX):
            raise ConfigurationError("independent apparent-CV stream requires T in [1, 2**32-1]")
        if convention == "as_run" and T is not None:
            raise ConfigurationError("T must be None under the as_run convention")
        self.base, self.run_id, self.convention, self.T = base, run_id, convention, T

    def _ss(self, key) -> np.random.SeedSequence:
        for x in key:
            if not (0 <= x <= UINT32_MAX):
                raise ConfigurationError("stream indices must lie in [0, 2**32-1]")
        return np.random.SeedSequence(entropy=[self.base, self.run_id], spawn_key=tuple(key))

    def outer(self, k):
        return self._ss((k,))

    def inner(self, k, b):
        return self._ss((k, b))

    def inner_cv(self, k, b):
        return self._ss((k, b, 0))

    def inner_reserve(self, k, b):
        return self._ss((k, b, 1))

    def apparent_cv(self, k):
        if self.convention == "as_run":
            return self.outer(k)
        if not 0 <= k <= UINT32_MAX:
            raise ConfigurationError("k out of range")
        return np.random.SeedSequence(entropy=[self.base, self.run_id, self.T], spawn_key=(k,))

    def apparent_cv_as_run(self, k):
        """Optional paired variant (Decision 2): the as-run apparent convention."""
        return self.outer(k)

    @staticmethod
    def generator(ss) -> np.random.Generator:
        return np.random.default_rng(ss)          # always a FRESH generator


class HistoricalAsRunStreams:
    """Reconstruction of the frozen driver's streams (scripts/40_phase9_dc_interval.py
    lines 96-98, 116-117, 130-142). FOR VERIFICATION ONLY (V1 frozen-path comparison,
    engine-equivalence tests). The execution entry point refuses this class."""

    kind = "historical_as_run"

    def __init__(self, base: int, run_id: int, B_in: int):
        self.base, self.run_id, self.B_in = base, run_id, B_in

    def outer(self, k):
        return np.random.SeedSequence(entropy=[self.base, self.run_id]).spawn(k + 1)[k]

    def inner(self, k, b):
        return self.outer(k).spawn(self.B_in)[b]

    def inner_cv(self, k, b):                    # the documented deviation: CV = resample state
        return self.inner(k, b)

    def apparent_cv(self, k):
        return self.outer(k)

    apparent_cv_as_run = apparent_cv

    @staticmethod
    def generator(ss):
        return np.random.default_rng(ss)


# ------------------------------------------------------------------ structural analysis
def coerce_words(x: int) -> tuple:
    """Model of numpy.random.bit_generator._coerce_to_uint32_array for one int
    (lowest 32-bit word first; documented in the installed docstring)."""
    if x < 0:
        raise ValueError("negative entropy is rejected by NumPy")
    if x == 0:
        return (0,)
    w = []
    while x:
        w.append(x & UINT32_MAX)
        x >>= 32
    return tuple(w)


def assembled_words(entropy: Iterable[int], spawn_key: Iterable[int]) -> tuple:
    """Model of SeedSequence's assembled entropy, normalised for pool-equivalence:
      - run entropy padded with zeros to POOL_SIZE when a spawn key is present
        (upstream `get_assembled_entropy`; NOT inspectable in the installed binary);
      - arrays shorter than POOL_SIZE behave as if zero-padded (upstream `mix_entropy`).
    Both behaviours are checked empirically by `verify_numpy_seed_behaviour`."""
    run = tuple(w for x in entropy for w in coerce_words(int(x)))
    sk = tuple(w for x in spawn_key for w in coerce_words(int(x)))
    if sk and len(run) < POOL_SIZE:
        run = run + (0,) * (POOL_SIZE - len(run))
    out = run + sk
    if len(out) < POOL_SIZE:
        out = out + (0,) * (POOL_SIZE - len(out))
    return out


def structural_domain_report(base: int, run_ids, T: Optional[int], convention: str,
                             numpy_behaviour_verified: bool) -> dict:
    """Input-level disjointness argument, valid for ALL k, b < 2**32 (no finite B_in bound):

      outer      (B, r, 0, 0, k)          length 5, word[2] = 0
      apparent   (B, r, T, 0, k)          length 5, word[2] = T != 0
      inner      (B, r, 0, 0, k, b)       length 6
      inner_cv   (B, r, 0, 0, k, b, 0)    length 7, last word 0
      reserve    (B, r, 0, 0, k, b, 1)    length 7, last word 1
      P9 seeds   (s, 0, 0, 0)             length 4, s != B
    Distinct families differ in length or in a fixed word; within a family the free
    words are the stream indices. This proves distinct HASH INPUTS only. Distinct pools
    and states follow only up to hash collisions (probability negligible, not zero)."""
    cond = {
        "base_single_word": 0 <= base <= UINT32_MAX,
        "run_ids_single_word": all(0 <= r <= UINT32_MAX for r in run_ids),
        "indices_single_word": "required at construction (PlanStreams._ss rejects >= 2**32)",
        "base_not_a_primary_seed": base not in PHASE9_PRIMARY_SEEDS,
        "numpy_behaviour_verified": bool(numpy_behaviour_verified),
    }
    if convention == "independent":
        cond["T_nonzero_single_word"] = isinstance(T, int) and 1 <= T <= UINT32_MAX
    holds = all(v is True or isinstance(v, str) for v in cond.values())
    return {"guarantee": "input-level distinctness of assembled entropy words",
            "scope": "all k, b in [0, 2**32-1]; all run_ids listed",
            "conditions": cond, "holds": holds,
            "not_guaranteed": ["distinct pools / PCG64 states (hash collision possible, negligible)",
                               "non-overlapping output sequences"]}


# ------------------------------------------------------------------ empirical checks
def verify_numpy_seed_behaviour(expected_numpy: str, base: int, T: Optional[int],
                                sample_k=(0, 1, 499, 2999), sample_b=(0, 99)) -> dict:
    """Empirical verification of the NumPy behaviour the structural argument relies on."""
    checks = {}
    checks["numpy_version"] = {"observed": np.__version__, "expected": expected_numpy,
                               "pass": np.__version__ == expected_numpy}
    try:
        from numpy.random.bit_generator import _coerce_to_uint32_array as c
        cases = [base, 0, 1, UINT32_MAX, UINT32_MAX + 1, 2 ** 64 + 5] + ([T] if T else [])
        ok = all(tuple(int(w) for w in c(x)) == coerce_words(x) for x in cases)
        checks["coercion_model"] = {"pass": ok}
    except ImportError:
        checks["coercion_model"] = {"pass": False, "detail": "private coercion function unavailable"}
    fp = state_fingerprint
    S = np.random.SeedSequence
    pad, distinct, spawn_eq = True, True, True
    for r in (0, 4):
        # positive controls: the zero-padding hazard exists (why T must be non-zero)
        pad &= fp(S([base, r])) == fp(S([base, r, 0, 0]))
        for k in sample_k:
            pad &= fp(S([base, r], spawn_key=(k,))) == fp(S([base, r, 0], spawn_key=(k,)))
            if T:
                distinct &= fp(S([base, r, T], spawn_key=(k,))) != fp(S([base, r], spawn_key=(k,)))
            ch = S([base, r]).spawn(k + 1)[k]
            spawn_eq &= tuple(ch.spawn_key) == (k,) and fp(ch) == fp(S([base, r], spawn_key=(k,)))
            for b in sample_b:
                inn = S([base, r]).spawn(k + 1)[k].spawn(b + 1)[b]   # fresh parents: spawn() is stateful
                spawn_eq &= fp(inn) == fp(S([base, r], spawn_key=(k, b)))
                cv = inn.spawn(2)
                spawn_eq &= fp(cv[0]) == fp(S([base, r], spawn_key=(k, b, 0)))
                spawn_eq &= fp(cv[1]) == fp(S([base, r], spawn_key=(k, b, 1)))
    # spawn() is stateful: a second call yields different children
    p = S([base, 0])
    first, second = p.spawn(1)[0], p.spawn(1)[0]
    checks["zero_padding_hazard_present"] = {"pass": pad,
        "meaning": "[B,r] == [B,r,0] == [B,r,0,0]: a zero domain tag would alias the outer stream"}
    checks["T_stream_distinct_from_outer"] = {"pass": distinct if T else None,
                                              "evaluated": bool(T)}
    checks["explicit_spawn_key_equals_spawn_list"] = {"pass": spawn_eq}
    checks["spawn_is_stateful"] = {"pass": tuple(first.spawn_key) != tuple(second.spawn_key)}
    checks["all_pass"] = all(v.get("pass") is True for v in checks.values()
                             if isinstance(v, dict) and v.get("pass") is not None)
    return checks


def enumerate_stream_states(streams: PlanStreams, run_ids, k_max: int, b_max: Optional[int],
                            include_primary_seeds=True, progress=None) -> dict:
    """Exhaustive PCG64 initial-state enumeration over a DECLARED finite domain:
    k in [0, k_max), b in [0, b_max), all families, all run_ids, plus the Phase 9 primary
    seeds. Requires b_max explicitly (B_in,max is an owner decision)."""
    if b_max is None:
        raise ConfigurationError("B_in_max not declared: the empirical domain is undefined")
    seen, collisions, n = {}, [], 0

    def add(fpv, sid):
        nonlocal n
        n += 1
        if fpv in seen:
            collisions.append((seen[fpv], sid))
        else:
            seen[fpv] = sid

    for r in run_ids:
        s = type(streams)(streams.base, r, streams.convention, streams.T)
        for k in range(k_max):
            add(state_fingerprint(s.outer(k)), f"outer/r{r}/k{k}")
            if s.convention == "independent":
                add(state_fingerprint(s.apparent_cv(k)), f"apparent_cv/r{r}/k{k}")
            for b in range(b_max):
                add(state_fingerprint(s.inner(k, b)), f"inner/r{r}/k{k}/b{b}")
                add(state_fingerprint(s.inner_cv(k, b)), f"inner_cv/r{r}/k{k}/b{b}")
                add(state_fingerprint(s.inner_reserve(k, b)), f"inner_reserve/r{r}/k{k}/b{b}")
            if progress and k % 100 == 0:
                progress(r, k)
    if include_primary_seeds:
        for sd in PHASE9_PRIMARY_SEEDS:
            add(state_fingerprint(np.random.SeedSequence(sd)), f"phase9_primary_seed/{sd}")
    return {"kind": "empirical", "n_states": n, "n_distinct": len(seen),
            "collisions": collisions, "pass": not collisions,
            "domain": {"run_ids": list(run_ids), "k_max": k_max, "b_max": b_max,
                       "convention": streams.convention, "T": streams.T},
            "claim": "no duplicate initial state within the enumerated domain only"}

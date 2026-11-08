"""Phase 8C recommendation 4: the scorer moved to src/ WITHOUT changing.

Two independent proofs that the extraction was behaviour-preserving:

  1. SOURCE identity -- the function source is byte-identical to the pre-move
     version, pinned by sha256 in tests/fixtures/premove_scorer_reference.json.
  2. BEHAVIOURAL identity -- scores recomputed on representative inputs match
     bit-for-bit (sha256 of the float64 buffer) against outputs captured from
     the pre-move implementation, together with its diagnostics and its error
     contract.

The fixture holds DIGESTS, not a copy of the code, so there is no second
implementation to drift.
"""

from __future__ import annotations

import hashlib
import importlib.util
import inspect
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.pathway_scoring import score_panel                   # noqa: E402

FIXTURE = ROOT / "tests" / "fixtures" / "premove_scorer_reference.json"
ANALYSIS = ROOT / "data" / "processed" / "analysis"
PANELS = ROOT / "config" / "pathways"


@pytest.fixture(scope="module")
def ref():
    return json.loads(FIXTURE.read_text())


@pytest.fixture(scope="module")
def data():
    mem = pd.read_csv(ANALYSIS / "cohort_membership.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    X = pd.read_parquet(ANALYSIS / "expression_median_S0.parquet")[s3].astype(np.float64)
    panels = {n: json.loads((PANELS / f"m2_panel_{n}.json").read_text())
              for n in ("4set", "5set")}
    return X, s3, panels


def digest(df: pd.DataFrame) -> dict:
    return {"sha256_float64": hashlib.sha256(
                np.ascontiguousarray(df.values, dtype=np.float64).tobytes()).hexdigest(),
            "columns": list(df.columns), "index_head": list(df.index[:3]),
            "index_tail": list(df.index[-3:]), "shape": list(df.shape)}


# ==================== 1. source identity ====================

def test_M01_source_is_byte_identical_to_premove(ref):
    got = hashlib.sha256(inspect.getsource(score_panel).encode()).hexdigest()
    assert got == ref["source_sha256"], (
        "the scorer body changed during the move -- this was a pure refactor")


def test_M02_signature_unchanged(ref):
    assert str(inspect.signature(score_panel)) == ref["signature"]


def test_M03_single_implementation_no_second_copy():
    """scripts/32 must RE-EXPORT, not redefine."""
    src = (ROOT / "scripts" / "32_phase8a_eda.py").read_text()
    assert "def score_panel(" not in src, "a second copy of the scorer exists"
    assert "from src.pathway_scoring import score_panel" in src
    s34 = (ROOT / "scripts" / "34_phase8c_pathway_scores.py").read_text()
    assert "def score_panel(" not in s34
    assert "from src.pathway_scoring import score_panel" in s34
    # the audit oracle is deliberately separate and must survive
    assert "def score_panel_independent(" in s34


def test_M04_callers_resolve_to_the_same_object():
    for rel, attr in (("scripts/32_phase8a_eda.py", "score_panel"),
                      ("scripts/34_phase8c_pathway_scores.py", "canonical_scorer")):
        spec = importlib.util.spec_from_file_location(rel, ROOT / rel)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        obj = getattr(m, attr)
        obj = obj() if attr == "canonical_scorer" else obj
        assert obj is score_panel, rel


def test_M05_module_is_importable_without_path_loading():
    import src.pathway_scoring as mod
    assert mod.score_panel is score_panel
    assert mod.__all__ == ["score_panel"]
    assert score_panel.__module__ == "src.pathway_scoring"


# ==================== 2. behavioural identity ====================

@pytest.mark.parametrize("pn", ["4set", "5set"])
def test_M06_full_cohort_scores_bit_identical(ref, data, pn):
    X, s3, panels = data
    got, _ = score_panel(X, panels[pn], train_ids=s3)
    assert digest(got) == ref["scenarios"][pn]["full_cohort"]


@pytest.mark.parametrize("pn", ["4set", "5set"])
def test_M07_training_fold_fit_bit_identical(ref, data, pn):
    X, s3, panels = data
    got, _ = score_panel(X, panels[pn], train_ids=sorted(s3)[:300])
    assert digest(got) == ref["scenarios"][pn]["fold_fit"]


@pytest.mark.parametrize("pn", ["4set", "5set"])
def test_M08_fold_isolation_deletion_bit_identical(ref, data, pn):
    X, s3, panels = data
    train = sorted(s3)[:300]
    held = [s for s in s3 if s not in set(train)]
    got, _ = score_panel(X.drop(columns=held), panels[pn], train_ids=train)
    assert digest(got) == ref["scenarios"][pn]["fold_fit_heldout_deleted"]


@pytest.mark.parametrize("pn", ["4set", "5set"])
def test_M09_fold_isolation_poisoning_bit_identical(ref, data, pn):
    X, s3, panels = data
    train = sorted(s3)[:300]
    held = [s for s in s3 if s not in set(train)]
    pois = X.copy()
    pois.loc[:, held] += 1000.0
    got, _ = score_panel(pois, panels[pn], train_ids=train)
    assert digest(got) == ref["scenarios"][pn]["fold_fit_heldout_poisoned"]


@pytest.mark.parametrize("pn", ["4set", "5set"])
def test_M10_training_params_still_ignore_heldout(data, pn):
    """Re-derived here, not read from the fixture."""
    X, s3, panels = data
    train = sorted(s3)[:300]
    held = [s for s in s3 if s not in set(train)]
    a, _ = score_panel(X, panels[pn], train_ids=train)
    b, _ = score_panel(X.drop(columns=held), panels[pn], train_ids=train)
    pois = X.copy(); pois.loc[:, held] += 1000.0
    c, _ = score_panel(pois, panels[pn], train_ids=train)
    assert np.array_equal(a.loc[train].values, b.loc[train].values)
    assert np.array_equal(a.loc[train].values, c.loc[train].values)
    assert not np.allclose(a.loc[held].values, c.loc[held].values)


@pytest.mark.parametrize("pn", ["4set", "5set"])
def test_M11_ddof_override_bit_identical(ref, data, pn):
    X, s3, panels = data
    got, _ = score_panel(X, panels[pn], train_ids=s3, ddof=0)
    assert digest(got) == ref["scenarios"][pn]["ddof0_override"]


@pytest.mark.parametrize("pn", ["4set", "5set"])
def test_M12_diagnostics_identical(ref, data, pn):
    X, s3, panels = data
    _, dfull = score_panel(X, panels[pn], train_ids=s3)
    _, dfold = score_panel(X, panels[pn], train_ids=sorted(s3)[:300])
    assert json.loads(json.dumps(dfull, default=str)) == ref["scenarios"][pn]["diag_full"]
    assert json.loads(json.dumps(dfold, default=str)) == ref["scenarios"][pn]["diag_fold"]


@pytest.mark.parametrize("pn", ["4set", "5set"])
def test_M13_membership_and_ordering_preserved(ref, data, pn):
    X, s3, panels = data
    got, _ = score_panel(X, panels[pn], train_ids=s3)
    assert list(got.columns) == ref["scenarios"][pn]["full_cohort"]["columns"]
    assert list(got.columns) == list(panels[pn]["gene_sets"])
    assert list(got.index) == s3


@pytest.mark.parametrize("pn", ["4set", "5set"])
def test_M14_zero_variance_handling_identical(ref, data, pn):
    X, s3, panels = data
    sc = ref["scenarios"][pn]
    Xz = X.copy()
    Xz.loc[sc["victim_gene"]] = 7.0
    _, dz = score_panel(Xz, panels[pn], train_ids=s3)
    assert {k: v["n_zero_variance_dropped"] for k, v in dz.items()} == \
        sc["zero_variance_injection"]


def test_M15_error_contract_identical(ref, data):
    X, s3, panels = data
    train = sorted(s3)[:300]
    unpinned = json.loads(json.dumps(panels["4set"]))
    del unpinned["scoring"]["sigma_ddof"]
    cases = {"no_train_ids": lambda: score_panel(X, panels["4set"]),
             "train_ids_none": lambda: score_panel(X, panels["4set"], train_ids=None),
             "unpinned_ddof": lambda: score_panel(X, unpinned, train_ids=train),
             "one_train_sample": lambda: score_panel(X, panels["4set"],
                                                     train_ids=train[:1])}
    for label, fn in cases.items():
        expect = ref["error_contract"][label]
        with pytest.raises(Exception) as ei:
            fn()
        assert type(ei.value).__name__ == expect["raised"], label
        if expect["raised"] == "ValueError":
            assert str(ei.value) == expect["message"], label


def test_M16_sigma_ddof_still_read_from_spec(data):
    X, s3, panels = data
    train = sorted(s3)[:300]
    a, _ = score_panel(X, panels["4set"], train_ids=train)
    b, _ = score_panel(X, panels["4set"], train_ids=train, ddof=1)
    c, _ = score_panel(X, panels["4set"], train_ids=train, ddof=0)
    assert panels["4set"]["scoring"]["sigma_ddof"] == 1
    assert np.array_equal(a.values, b.values)
    assert not np.allclose(a.values, c.values)


def test_M17_fixture_pins_real_content(ref):
    """Guard against a fixture that would pass vacuously."""
    assert len(ref["source_sha256"]) == 64
    for pn in ("4set", "5set"):
        s = ref["scenarios"][pn]
        assert len(s["full_cohort"]["columns"]) == (4 if pn == "4set" else 5)
        assert s["full_cohort"]["shape"] == [453, 4 if pn == "4set" else 5]
        assert s["full_cohort"]["sha256_float64"] != s["fold_fit"]["sha256_float64"]
    assert set(ref["error_contract"]) == {"no_train_ids", "train_ids_none",
                                          "unpinned_ddof", "one_train_sample"}

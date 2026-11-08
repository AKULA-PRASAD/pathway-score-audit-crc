"""Amendment M2-A1: the training-fold sigma estimator is pinned to ddof=1.

Verifies the five properties the amendment is supposed to guarantee:

  * both specifications DECLARE ddof=1;
  * the scoring implementation USES the declared value and refuses to fall
    back on a library default;
  * mu/sigma are fold-specific, estimated from training samples only;
  * held-out samples cannot influence standardization;
  * membership and model roles are unchanged by the amendment.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

PANELS = ROOT / "config" / "pathways"
PRIMARY = PANELS / "m2_panel_4set.json"
SENSITIVITY = PANELS / "m2_panel_5set.json"
RECORD = PANELS / "M2_FREEZE_RECORD.json"
ANALYSIS = ROOT / "data" / "processed" / "analysis"

AMENDMENT_ID = "M2-A1"


@pytest.fixture(scope="module")
def specs():
    return (json.loads(PRIMARY.read_text()), json.loads(SENSITIVITY.read_text()))


@pytest.fixture(scope="module")
def record():
    return json.loads(RECORD.read_text())


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location(
        "phase8a_scoring", ROOT / "scripts" / "32_phase8a_eda.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def expr():
    mem = pd.read_csv(ANALYSIS / "cohort_membership.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    X = pd.read_parquet(ANALYSIS / "expression_median_S0.parquet")[s3]
    return X.astype(np.float64), s3


# ============ 1. both specifications declare ddof = 1 ============

def test_D01_both_panels_declare_sigma_ddof_1(specs):
    for d in specs:
        assert d["scoring"]["sigma_ddof"] == 1, d["panel"]
        assert isinstance(d["scoring"]["sigma_ddof"], int)


def test_D02_freeze_record_declares_sigma_ddof_1(record):
    assert record["scoring"]["sigma_ddof"] == 1


def test_D03_declaration_forbids_library_defaults(specs, record):
    for d in (*specs, record):
        rule = d["scoring"]["sigma_ddof_rule"]
        assert "n_train - 1" in rule
        assert "MUST NOT rely on a library default" in rule


def test_D04_formula_states_the_denominator(specs, record):
    for d in (*specs, record):
        assert "n_train - 1" in d["scoring"]["formula"]


# ============ 2. the implementation uses the declared value ============

def test_D05_scoring_reads_ddof_from_the_specification(mod, expr, specs):
    """ddof must come from the panel, not from a default in the signature."""
    X, s3 = expr
    primary, _ = specs
    train = sorted(s3)[:300]
    from_spec, _ = mod.score_panel(X, primary, train_ids=train)
    explicit_1, _ = mod.score_panel(X, primary, train_ids=train, ddof=1)
    explicit_0, _ = mod.score_panel(X, primary, train_ids=train, ddof=0)
    assert np.array_equal(from_spec.values, explicit_1.values)
    assert not np.allclose(from_spec.values, explicit_0.values)


def test_D06_scoring_refuses_a_panel_that_does_not_pin_ddof(mod, expr, specs):
    """The failure mode the amendment exists to prevent: a silent default."""
    X, s3 = expr
    unpinned = json.loads(json.dumps(specs[0]))
    del unpinned["scoring"]["sigma_ddof"]
    with pytest.raises(ValueError, match="sigma_ddof"):
        mod.score_panel(X, unpinned, train_ids=sorted(s3)[:300])


def test_D07_declared_value_is_actually_the_sample_sd(mod, expr, specs):
    """Recompute one gene's standardisation by hand against the declaration."""
    X, s3 = expr
    primary, _ = specs
    train = sorted(s3)[:300]
    gene = primary["gene_sets"]["HALLMARK_E2F_TARGETS"]["members_resolved"][0]
    row = X.loc[gene, train].values
    n = len(train)
    manual_sd = float(np.sqrt(((row - row.mean()) ** 2).sum() / (n - 1)))
    population_sd = float(np.sqrt(((row - row.mean()) ** 2).sum() / n))
    assert manual_sd != population_sd
    sd_used = float(X.loc[[gene], train].std(axis=1,
                    ddof=primary["scoring"]["sigma_ddof"]).iloc[0])
    assert np.isclose(sd_used, manual_sd, rtol=0, atol=1e-12)


def test_D08_ddof_choice_is_a_constant_rescaling(mod, expr, specs):
    """The amendment's stated numerical consequence, verified."""
    X, s3 = expr
    primary, _ = specs
    train = sorted(s3)[:300]
    a, _ = mod.score_panel(X, primary, train_ids=train, ddof=1)
    b, _ = mod.score_panel(X, primary, train_ids=train, ddof=0)
    ratio = b.values / a.values
    assert np.allclose(ratio, np.sqrt(len(train) / (len(train) - 1)), atol=1e-10)
    # rank order preserved -> correlations unchanged
    assert np.allclose(a.corr().values, b.corr().values, atol=1e-10)


# ============ 3. fold-specific training mu and sigma ============

def test_D09_mu_and_sigma_are_fold_specific(mod, expr, specs):
    """Two different training folds must give different standardisation."""
    X, s3 = expr
    primary, _ = specs
    fold_a, fold_b = sorted(s3)[:200], sorted(s3)[200:400]
    sa, _ = mod.score_panel(X, primary, train_ids=fold_a)
    sb, _ = mod.score_panel(X, primary, train_ids=fold_b)
    assert not np.allclose(sa.values, sb.values), (
        "scores identical across different training folds -- mu/sigma are not "
        "being re-estimated per fold")


def test_D10_training_population_is_mandatory(mod, expr, specs):
    X, _ = expr
    with pytest.raises(TypeError):
        mod.score_panel(X, specs[0])
    with pytest.raises(ValueError, match="mandatory"):
        mod.score_panel(X, specs[0], train_ids=None)


def test_D11_mu_sigma_come_from_training_samples_only(mod, expr, specs):
    """A training-fold sample must standardise to its own fold's statistics."""
    X, s3 = expr
    primary, _ = specs
    train = sorted(s3)[:250]
    scores, _ = mod.score_panel(X, primary, train_ids=train)
    genes = [g for g in primary["gene_sets"]["HALLMARK_P53_PATHWAY"]
             ["members_resolved"] if g in X.index]
    Xtr = X.loc[genes, train]
    mu, sd = Xtr.mean(axis=1), Xtr.std(axis=1, ddof=1)
    keep = sd.index[sd.values > 0]
    manual = ((X.loc[keep] .sub(mu[keep], axis=0))
              .div(sd[keep], axis=0)).mean(axis=0)
    assert np.allclose(scores["HALLMARK_P53_PATHWAY"].values, manual.values,
                       atol=1e-12)
    # training-fold scores must be centred on ~0 by construction
    assert abs(scores["HALLMARK_P53_PATHWAY"].loc[train].mean()) < 1e-10


# ============ 4. held-out samples cannot influence standardization ============

def test_D12_heldout_rows_cannot_change_training_parameters(mod, expr, specs):
    """Deleting held-out samples from the input must change nothing."""
    X, s3 = expr
    for panel in specs:
        train = sorted(s3)[:300]
        held = [s for s in s3 if s not in set(train)]
        full, _ = mod.score_panel(X, panel, train_ids=train)
        trimmed, _ = mod.score_panel(X.drop(columns=held), panel, train_ids=train)
        assert np.array_equal(full.loc[train].values, trimmed.loc[train].values), \
            panel["panel"]


def test_D13_perturbing_heldout_values_cannot_change_training_scores(mod, expr, specs):
    """Stronger than deletion: corrupt the held-out data and refit."""
    X, s3 = expr
    primary, _ = specs
    train = sorted(s3)[:300]
    held = [s for s in s3 if s not in set(train)]
    clean, _ = mod.score_panel(X, primary, train_ids=train)
    poisoned = X.copy()
    poisoned.loc[:, held] = poisoned.loc[:, held] + 1000.0
    dirty, _ = mod.score_panel(poisoned, primary, train_ids=train)
    assert np.array_equal(clean.loc[train].values, dirty.loc[train].values), (
        "held-out values leaked into the training-fold standardisation")
    # and the held-out scores MUST move, or the test proves nothing
    assert not np.allclose(clean.loc[held].values, dirty.loc[held].values)


def test_D14_full_cohort_fit_differs_from_fold_fit(mod, expr, specs):
    """If these agreed, fold isolation would be unobservable."""
    X, s3 = expr
    primary, _ = specs
    train = sorted(s3)[:300]
    held = [s for s in s3 if s not in set(train)]
    fold, _ = mod.score_panel(X, primary, train_ids=train)
    whole, _ = mod.score_panel(X, primary, train_ids=s3)
    assert np.abs(fold.loc[held].values - whole.loc[held].values).max() > 0


# ============ 5. membership and roles unchanged ============

def test_D15_panel_membership_unchanged_by_the_amendment(specs, record):
    primary, sensitivity = specs
    assert sorted(primary["gene_sets"]) == sorted(
        record["specifications"]["primary"]["gene_sets"])
    assert sorted(sensitivity["gene_sets"]) == sorted(
        record["specifications"]["sensitivity"]["gene_sets"])
    assert primary["n_sets"] == 4 and sensitivity["n_sets"] == 5
    assert primary["union_unique_genes"] == 607
    assert sensitivity["union_unique_genes"] == 702


def test_D16_resolved_members_byte_identical_to_pre_amendment(specs):
    """Per-set member lists and counts must be untouched."""
    expected = {
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": (200, 198),
        "HALLMARK_E2F_TARGETS": (200, 195),
        "HALLMARK_WNT_BETA_CATENIN_SIGNALING": (42, 41),
        "HALLMARK_P53_PATHWAY": (200, 194),
        "HALLMARK_PI3K_AKT_MTOR_SIGNALING": (105, 104),
    }
    for d in specs:
        for name, gs in d["gene_sets"].items():
            n_source, n_res = expected[name]
            assert gs["n_source"] == n_source, name
            assert len(dict.fromkeys(gs["members_resolved"])) == n_res, name


def test_D17_roles_unchanged(specs, record):
    primary, sensitivity = specs
    assert primary["role"] == "PRIMARY"
    assert sensitivity["role"].startswith("SENSITIVITY")
    assert record["approval"]["primary_panel"] == "4-set"
    assert record["approval"]["sensitivity_panel"] == "5-set"
    assert (record["approval"]["sensitivity_only_addition"]
            == "HALLMARK_PI3K_AKT_MTOR_SIGNALING")
    assert "HALLMARK_PI3K_AKT_MTOR_SIGNALING" not in primary["gene_sets"]
    assert "HALLMARK_PI3K_AKT_MTOR_SIGNALING" in sensitivity["gene_sets"]


def test_D18_immutability_clause_survives(record):
    clause = record["approval"]["immutability_clause"]
    assert "must NOT be promoted to the primary panel" in clause
    assert "must NOT be removed from the sensitivity panel" in clause


def test_D19_coverage_rules_unchanged(specs):
    for d in specs:
        mc = d["scoring"]["minimum_set_coverage"]
        assert ">=80% of source members must resolve" in mc
        assert ">=20 genes must remain" in mc
        assert d["scoring"]["missing_gene_members"].endswith("never imputed")
        assert "EXCLUDED for that fold" in d["scoring"]["zero_variance_genes"]


# ============ amendment record integrity ============

def test_D20_amendment_recorded_in_both_panels_and_the_record(specs, record):
    for d in (*specs, record):
        hist = d["amendment_history"]
        assert [a["amendment_id"] for a in hist][0] == AMENDMENT_ID
        a = hist[0]
        assert a["date"] == "2026-11-02"
        assert a["approved_by"] == "project owner"
        assert "ddof = 1" in a["approved_convention"]
        assert "AMBIGUOUS" in a["previous_state"]
        assert "ddof=0" in a["previous_state"] and "ddof=1" in a["previous_state"]
        assert a["membership_changed"] is False
        assert a["roles_changed"] is False
        assert a["outcome_data_used"] is False
        assert "sample-standard-deviation scaling" in a["expected_numerical_consequence"]
        assert "population-standard-deviation" in a["expected_numerical_consequence"]
        assert any("Phase 8A" in r for r in a["rationale"])
        assert any("reproducibility" in r for r in a["rationale"])


def test_D21_specification_revision_bumped(specs, record):
    """M2-A1 took the spec to revision 2; later amendments only increase it."""
    rev = record["specification_revision"]
    assert rev >= 2
    assert rev == 1 + len(record["amendment_history"])
    for d in specs:
        assert d["specification_revision"] == rev
    for role in ("primary", "sensitivity"):
        assert record["specifications"][role]["specification_revision"] == rev


def test_D22_original_hashes_preserved_not_overwritten(record):
    """Every superseded hash survives; later amendments append, never replace."""
    orig = record["original_freeze"]["sha256"]
    assert orig["m2_panel_4set.json"].startswith("bab99da9f139db34")
    assert orig["m2_panel_5set.json"].startswith("6058f733b2acda1c")
    for role, pre in (("primary", "bab99da9f139db34"),
                      ("sensitivity", "6058f733b2acda1c")):
        sup = record["specifications"][role]["superseded_sha256"]
        assert any(e["sha256"].startswith(pre) for e in sup), role
        # M2-A1's entry must still be there after M2-A2 appended its own
        assert AMENDMENT_ID in {e["superseded_by"] for e in sup}, role
        assert len(sup) == len(record["amendment_history"]), role


def test_D23_new_hashes_recorded_and_correct(record):
    """The recorded current hash tracks the file through every amendment."""
    for role, path in (("primary", PRIMARY), ("sensitivity", SENSITIVITY)):
        s = record["specifications"][role]
        got = hashlib.sha256(path.read_bytes()).hexdigest()
        assert s["sha256"] == got, role
        assert got not in [e["sha256"] for e in s["superseded_sha256"]]


def test_D24_panels_still_marked_frozen_and_immutable(specs):
    for d in specs:
        assert d["status"] == "FROZEN"
        assert d["immutable"] is True
        assert d["outcome_data_used"] is False


def test_D25_amendment_is_permitted_by_the_modification_policy(specs):
    """A revision informed by outcome data would reclassify M2 as exploratory."""
    for d in specs:
        a = d["amendment_history"][0]
        assert "non-outcome corrections" in a["permitted_under"]
        assert "outcome-blind" in a["outcome_blind_statement"].lower()
        assert "M3-exploratory" in a["outcome_blind_statement"]


def test_D26_earlier_phase_attestations_are_not_rewritten():
    """A-19 and the QC decisions record must still hold the ORIGINAL hashes."""
    a19 = json.loads((ROOT / "config/analysis/A19_FREEZE_RECORD.json").read_text())
    qc = json.loads((ROOT / "config/analysis/QC_DECISIONS_RECORD.json").read_text())
    assert a19["upstream_frozen_hashes"]["m2_panel_4set.json"].startswith(
        "bab99da9f139db34")
    assert qc["frozen_artifacts_unmodified"]["m2_panel_5set.json"].startswith(
        "6058f733b2acda1c")


def test_D27_resolver_accepts_amendment_but_rejects_drift():
    from src.m2_amendment import resolve_hash, sha256
    ok, why = resolve_hash(PRIMARY, sha256(PRIMARY))
    assert ok and why == "current"
    ok, why = resolve_hash(
        PRIMARY, "bab99da9f139db34ccc2d5f5f954cda9e92fc6cf33f41a076024e5146b065dfe")
    assert ok and AMENDMENT_ID in why
    ok, why = resolve_hash(PRIMARY, "f" * 64)
    assert not ok and "UNDOCUMENTED DRIFT" in why

"""Phase 10 safeguards: cohorts, leakage, folds, metrics, controls, sealed labels, Phase 9 integrity.

Run with the Phase 10 environment:  ~/.venvs/coad-phase10/bin/python -m pytest tests/phase10 -q
"""
import json

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import roc_auc_score

from phase10 import common as C

DATA = (C.ROOT / "data" / "processed" / "analysis").is_dir()
needs_data = pytest.mark.skipif(not DATA, reason="needs the local data bundle")


# ---------------------------------------------------------------- 5. metrics vs references
def test_auc_pair_counts_match_sklearn():
    rng = np.random.default_rng(1)
    for _ in range(20):
        y = rng.integers(0, 2, 60)
        s = rng.integers(0, 8, 60).astype(float)          # many tied scores
        num, den = C.auc_pair_counts(y, s)
        assert num / den == pytest.approx(roc_auc_score(y, s), abs=1e-12)


def test_harrell_pair_counts_match_lifelines():
    from lifelines.utils import concordance_index
    rng = np.random.default_rng(2)
    for _ in range(20):
        t = rng.integers(1, 15, 70).astype(float)          # heavy tied times
        e = rng.integers(0, 2, 70)
        r = rng.integers(0, 6, 70).astype(float)           # tied risks
        num, den = C.harrell_pair_counts(t, e, r)
        assert num / den == pytest.approx(concordance_index(t, -r, e), abs=1e-12)


def test_calibration_matches_logistic_reference():
    from sklearn.linear_model import LogisticRegression
    rng = np.random.default_rng(3)
    p = rng.uniform(0.02, 0.98, 400)
    y = (rng.uniform(size=400) < p ** 1.5).astype(int)
    _, slope = C.calibration(y, p)
    lo = np.log(p / (1 - p))
    ref = LogisticRegression(C=np.inf, tol=1e-12, max_iter=10000).fit(lo[:, None], y).coef_[0, 0]
    assert slope == pytest.approx(ref, rel=1e-5)
    a, _ = C.calibration((rng.uniform(size=400) < p).astype(int), p)
    assert abs(a) < 0.5                                  # well-calibrated draws -> intercept near 0


# ---------------------------------------------------------------- 4/11. folds
def test_folds_partition_stratify_and_are_deterministic():
    ids = [f"s{i}" for i in range(409)]
    y = np.r_[np.ones(59), np.zeros(350)]
    a = C.make_folds(ids, y, "mmr", 3, 5, 123)
    b = C.make_folds(ids, y, "mmr", 3, 5, 123)
    pd.testing.assert_frame_equal(a, b)
    for r, g in a.groupby("repeat"):
        assert sorted(g.sample_id) == sorted(ids)          # every patient exactly once per repeat
        pos = g.assign(y=pd.Series(y, index=ids).loc[g.sample_id].to_numpy()).groupby("fold").y.sum()
        assert pos.max() - pos.min() <= 1                  # stratified
    assert not a.equals(C.make_folds(ids, y, "mmr", 3, 5, 124))


# ---------------------------------------------------------------- synthetic CV helpers
class Spy(C.ColumnSelect):
    seen = []

    def fit(self, X, y=None):
        Spy.seen.append(set(pd.DataFrame(X).index))
        return self


def _synthetic(n=409, p=150, informative=10, signal=True, seed=0):
    rng = np.random.default_rng(seed)
    X = pd.DataFrame(rng.normal(size=(n, p)), columns=[f"g{i}" for i in range(p)],
                     index=pd.Index([f"s{i}" for i in range(n)], name="sample_id"))
    X["MLH1"] = rng.normal(size=n)
    # planted genes vary more (SD 2), as biologically driven genes do, so the outcome-blind
    # variance filter retains them; the test then checks learning, not the filter's luck.
    X.iloc[:, :informative] *= 2.0
    lin = X.iloc[:, :informative].sum(axis=1) / 2.0 * (1.2 if signal else 0.0)
    y = pd.Series((lin + rng.logistic(size=n) > np.quantile(lin + rng.logistic(size=n), 0.855)).astype(int),
                  index=X.index)
    return X, y


def _cfg(k=50):
    cfg = C.load_config()
    cfg["variance_top_k"] = k
    cfg["mmr"]["enet"]["C"] = [0.01, 0.1, 1.0]
    cfg["mmr"]["enet"]["l1_ratio"] = [0.5]
    cfg["mmr"]["rf"]["n_estimators"] = 50
    return cfg


def _run(X, y, cfg, repeats=1, only=None, spy=False):
    from phase10 import mmr as M
    genes = [c for c in X.columns]
    Xc = pd.DataFrame({"age": np.arange(len(X), dtype=float) % 7}, index=X.index)
    folds = C.make_folds(list(X.index), y.to_numpy(), "mmr", repeats, 5, cfg["seed_base"])

    def factory(cv, rs):
        ms = M.build_models(cfg, genes, cv, rs)
        if spy:
            from sklearn.pipeline import Pipeline
            from sklearn.linear_model import LogisticRegression
            from sklearn.model_selection import GridSearchCV
            ms = {"spy_enet": (GridSearchCV(Pipeline([("spy", Spy(genes)), ("var", C.VarianceFilter(cfg["variance_top_k"])),
                                                      ("sc", __import__("sklearn.preprocessing").preprocessing.StandardScaler()),
                                                      ("clf", LogisticRegression(solver="saga", l1_ratio=0.5, max_iter=2000))]),
                                             {"clf__C": [0.1, 1.0]}, cv=cv, n_jobs=1), "genes")}
        return ms
    return folds, M.run_cv(cfg, list(X.index), y, {"genes": X, "clinical": Xc}, folds, factory, only=only, log=lambda *_: None)


# ---------------------------------------------------------------- 3. no test rows in any fit
def test_no_test_rows_reach_any_preprocessing_fit():
    X, y = _synthetic(n=200, p=60)
    Spy.seen = []
    folds, _ = _run(X, y, _cfg(20), only=(0, 2), spy=True)
    test_ids = set(folds[(folds.repeat == 0) & (folds.fold == 2)].sample_id)
    assert Spy.seen, "spy never fitted"
    for s in Spy.seen:                                      # inner fits AND the refit
        assert not (s & test_ids)


@needs_data
def test_track_a_variance_filter_sees_only_training_rows(monkeypatch):
    from phase10 import recurrence as R
    from src.phase9_corrected.data import load_context
    ctx = load_context(C.ROOT)
    seen = []

    class SpyVF(C.VarianceFilter):
        def fit(self, X, y=None):
            seen.append(set(X.index))
            return super().fit(X, y)
    monkeypatch.setattr(R.C, "VarianceFilter", SpyVF)
    a, b = ctx.ids[:300], ctx.ids[300:]
    R.m3_design(ctx, ctx.X.T, a, b, 100)
    assert seen == [set(a)]


# ---------------------------------------------------------------- 4/11. identical folds, reproducibility
def test_identical_folds_across_models_and_reproducible():
    X, y = _synthetic(n=200, p=60)
    cfg = _cfg(20)
    _, (fd1, oof1, _) = _run(X, y, cfg)
    _, (fd2, oof2, _) = _run(X, y, cfg)
    pd.testing.assert_frame_equal(oof1, oof2)               # fixed seeds -> identical predictions
    sets = oof1.groupby("model").sample_id.apply(lambda s: tuple(s))
    assert len(set(sets)) == 1                              # same patients, same order, every model
    per = fd1.groupby(["repeat", "fold"])[["n_test", "n_pos_test"]].nunique()
    assert (per == 1).all().all()


# ---------------------------------------------------------------- 6/7. synthetic positive and negative controls
def test_synthetic_signal_recovered():
    from phase10 import mmr as M
    X, y = _synthetic(signal=True, seed=5)
    fd, oof, det = _run(X, y, _cfg(80))[1]
    en = fd[fd.model == "enet_logistic"]
    assert en.auc_num.sum() / en.auc_den.sum() > 0.8
    sel = M.selection_stability(det, len(en))
    top = set(sel.head(10).gene)
    assert len(top & {f"g{i}" for i in range(10)}) >= 6       # planted genes dominate selection


def test_synthetic_null_near_chance():
    X, y = _synthetic(signal=False, seed=6)
    fd = _run(X, y, _cfg(80))[1][0]
    for m in ("enet_logistic", "random_forest"):
        g = fd[fd.model == m]
        assert abs(g.auc_num.sum() / g.auc_den.sum() - 0.5) < 0.12


def test_n1_permutation_moves_time_event_pairs_together():
    from phase10.recurrence import permute_pairs
    t = np.arange(50, dtype=float)
    e = (np.arange(50) % 3 == 0).astype(float)
    t2, e2 = permute_pairs(t, e, 7)
    pairs = set(zip(t, e))
    assert set(zip(t2, e2)) == pairs and not np.array_equal(t, t2)
    assert e2.sum() == e.sum()


# ---------------------------------------------------------------- M2 definition preserved
@needs_data
def test_m2_strict_wrapper_equals_frozen_pipeline_fit_at_fixed_lambda():
    from phase10 import recurrence as R
    from src.phase9_corrected.data import LAM_GRID, load_context
    from src.phase9_execution import pipeline_fit
    ctx = load_context(C.ROOT)
    tr = ctx.ids[:360]
    idx = R._idx(ctx.ids, tr)
    lam = float(LAM_GRID[5])
    mine = R.m2_fit(ctx, tr, ctx.t[idx], ctx.e[idx], lam)
    ref, *_ = pipeline_fit(ctx.X, ctx.cov, ctx.panel, tr, tr, ctx.t[idx], ctx.e[idx], np.array([lam]),
                           np.random.default_rng(0))
    np.testing.assert_allclose(mine.params_.to_numpy(), ref.params_.to_numpy(), rtol=0, atol=1e-12)


# ---------------------------------------------------------------- 1/2. cohorts and labels
@needs_data
def test_track_b_cohort_counts_and_labels():
    from phase10 import mmr as M
    genes = (C.REPORTS / "shared_genes.txt").read_text().split() if (C.REPORTS / "shared_genes.txt").exists() else M.shared_gene_list()
    ids, y, Xg, Xc = M.load_gse(genes)
    assert len(ids) == 409 and len(set(ids)) == 409
    assert int(y.sum()) == 59 and int((1 - y).sum()) == 350
    assert y.notna().all() and Xc.notna().all().all() and np.isfinite(Xg.to_numpy()).all()
    assert len(genes) == 16322 and "MLH1" in genes


@needs_data
def test_track_a_cohort_counts():
    from src.phase9_corrected.data import load_context
    ctx = load_context(C.ROOT)
    assert len(ctx.ids) == len(set(ctx.ids)) == 452
    assert int(ctx.e.sum()) == 136 and int((1 - ctx.e).sum()) == 316
    assert np.isfinite(ctx.t).all() and set(np.unique(ctx.e)) <= {0.0, 1.0}


@needs_data
def test_tcga_population_label_free_and_primary_tumours():
    p = C.REPORTS / "tcga" / "population_primary_stage2_3.csv"
    if not p.exists():
        pytest.skip("run `python -m phase10.tcga prepare` first")
    pop = pd.read_csv(p).sample_id
    assert pop.is_unique and (pop.str[13:15] == "01").all() and pop.str[:12].is_unique
    rep = json.loads((C.REPORTS / "tcga" / "prepare_report.json").read_text())
    assert rep["label_free"] and rep["primary_population_n"] == len(pop)


# ---------------------------------------------------------------- 8/9. sealed labels, label-free scaling
@needs_data
def test_tcga_labels_blocked_without_valid_committed_freeze(tmp_path, monkeypatch):
    from phase10 import tcga as T
    cfg = C.load_config()
    if not T.MANIFEST.exists():
        with pytest.raises(T.SealedLabelError):
            T.load_sealed_labels(cfg)
        return
    monkeypatch.setattr(T, "APPROVAL", tmp_path / "approval.json")
    with pytest.raises(T.SealedLabelError, match="approval"):
        T.load_sealed_labels(cfg)                           # no approval file
    (tmp_path / "approval.json").write_text(json.dumps({"manifest_sha256": "0" * 64}))
    with pytest.raises(T.SealedLabelError, match="does not match"):
        T.load_sealed_labels(cfg)                           # wrong manifest hash
    (tmp_path / "approval.json").write_text(json.dumps({"manifest_sha256": C.sha256_file(T.MANIFEST)}))
    import subprocess
    committed = subprocess.run(["git", "-C", str(C.ROOT), "ls-files", "--error-unmatch",
                                str(T.MANIFEST.relative_to(C.ROOT))], capture_output=True).returncode == 0
    if not committed:
        with pytest.raises(T.SealedLabelError, match="committed"):
            T.load_sealed_labels(cfg)                       # approved but not committed


def test_tcga_adaptation_uses_expression_only(monkeypatch):
    from phase10 import tcga as T
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    rng = np.random.default_rng(4)
    Xs = pd.DataFrame(rng.normal(size=(80, 5)), columns=list("abcde"))
    p = Pipeline([("sel", C.ColumnSelect(list("abcde"))), ("var", C.VarianceFilter(3)), ("sc", StandardScaler()),
                  ("clf", LogisticRegression())]).fit(Xs, rng.integers(0, 2, 80))
    Xt = pd.DataFrame(rng.normal(5, 3, size=(40, 5)), columns=list("abcde"))
    real_read = pd.read_csv

    def guard(path, *a, **k):
        assert "clinicalMatrix" not in str(path), "adaptation must not read the clinical/label file"
        return real_read(path, *a, **k)
    monkeypatch.setattr(pd, "read_csv", guard)
    ad, nzero = T.adapt(p, Xt)
    cols = p.named_steps["var"].columns_
    np.testing.assert_allclose(ad.named_steps["sc"].mean_, Xt[cols].mean().to_numpy())
    np.testing.assert_allclose(ad.named_steps["sc"].scale_, Xt[cols].std(ddof=0).to_numpy())
    np.testing.assert_array_equal(ad.named_steps["clf"].coef_, p.named_steps["clf"].coef_)   # frozen
    assert nzero == 0 and p.named_steps["sc"].mean_[0] != ad.named_steps["sc"].mean_[0]      # original untouched


# ---------------------------------------------------------------- 10. Phase 9 integrity
@needs_data
def test_phase9_hashes_unchanged():
    r = C.phase9_integrity()
    assert r["frozen_code_files"] == 27 and r["frozen_code_mismatches"] == []
    assert r["archive_files"] == 52 and r["archive_mismatches"] == []


# ---------------------------------------------------------------- 12. results on the registered cohort/label
def test_reported_results_use_registered_cohorts_and_labels():
    b = C.REPORTS / "mmr" / "summary.json"
    a = C.REPORTS / "recurrence" / "summary.json"
    m = C.REPORTS / "tcga" / "freeze_manifest.json"
    if not (a.exists() or b.exists() or m.exists()):
        pytest.skip("no Phase 10 results yet")
    if b.exists():
        s = json.loads(b.read_text())
        assert s["cohort"] == {"n": 409, "dMMR": 59, "pMMR": 350, "shared_genes": 16322}
        assert "ComBat" in s["caveat"]
    if a.exists():
        s = json.loads(a.read_text())
        assert s["cohort"] == {"n": 452, "events": 136, "censored": 316}
        assert not any(k for k in json.dumps(s).lower().split('"') if k in ("brier", "auc_t", "td_auc"))
    if m.exists():
        s = json.loads(m.read_text())
        assert s["label_definition"]["column"] == "MSI_updated_Oct62011"
        assert s["label_definition"]["positive"] == ["MSI-H"] and s["label_definition"]["negative"] == ["MSS", "MSI-L"]
        assert s["training"]["n"] == 409

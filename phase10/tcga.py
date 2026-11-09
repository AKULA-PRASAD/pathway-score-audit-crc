"""Track B, independent-cohort, cross-platform evaluation of MMR/MSI-associated expression
classification in TCGA-COAD.

GSE39582 label = mmr.status (dMMR/pMMR). TCGA label = MSI_updated_Oct62011 (MSI-H vs MSS+MSI-L).
These are different, associated constructs; this is not validation of a dMMR classifier.

Stages (each a subcommand):
  prepare   label-free: expression, primary tumours, allow-listed clinical columns, populations
  freeze    fit final GSE39582 models; write freeze_manifest.json (models, features, hashes)
  evaluate  sealed: reads TCGA MSI labels only if the manifest is approved AND committed, once

Normalisation: expression features are z-scored with the TCGA evaluation population's own gene
means/SDs (unsupervised target-cohort adaptation; label-free). This needs a batch of target
samples and is NOT a single-patient prospective pipeline.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import json
import re
import subprocess
from datetime import datetime, timezone

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss

from phase10 import common as C
from phase10 import mmr as M

OUT = C.REPORTS / "tcga"
MANIFEST = OUT / "freeze_manifest.json"
APPROVAL = OUT / "UNBLINDING_APPROVAL.json"
CLIN = C.TCGA_DIR / "COAD_clinicalMatrix"
EXPR = C.TCGA_DIR / "HiSeqV2.gz"


class SealedLabelError(PermissionError):
    pass


# ---------------------------------------------------------------- label-free preparation
def stage_roman(s):
    m = re.match(r"^\s*Stage\s+(IV|III|II|I)", str(s))
    return m.group(1) if m else None


def read_allow_listed_clinical(cfg):
    allow = cfg["tcga"]["clinical_allow_list"]
    header = CLIN.open().readline().rstrip("\n").split("\t")
    assert cfg["tcga"]["label_column"] not in allow
    clin = pd.read_csv(CLIN, sep="\t", usecols=allow, dtype=str)
    assert list(clin.columns) == [c for c in header if c in allow]       # nothing else was loaded
    return clin.set_index("sampleID")


def load_tcga_expression(genes=None):
    G = pd.read_csv(EXPR, sep="\t", index_col=0)
    G = G.loc[:, [c for c in G.columns if c[13:15] == "01"]]
    G.columns.name = "sample_id"
    return (G.loc[list(genes)] if genes is not None else G).T.astype(float)


def prepare(cfg):
    genes = (C.REPORTS / "shared_genes.txt").read_text().split()
    Xt = load_tcga_expression(genes)
    assert Xt.index.str[:12].is_unique, "one primary tumour per patient expected"
    clin = read_allow_listed_clinical(cfg).reindex(Xt.index)
    t = cfg["tcga"]
    loc = clin["anatomic_neoplasm_subdivision"]
    Xc = pd.DataFrame({
        "age": pd.to_numeric(clin["age_at_initial_pathologic_diagnosis"], errors="coerce"),
        "male": clin["gender"].map({"MALE": 1.0, "FEMALE": 0.0}),
        "stage3": clin["pathologic_stage"].map(stage_roman).map({"II": 0.0, "III": 1.0}),
        "proximal": np.where(loc.isin(t["proximal_sites"]), 1.0, np.where(loc.isin(t["distal_sites"]), 0.0, np.nan)),
    }, index=Xt.index)
    stage = clin["pathologic_stage"].map(stage_roman)
    complete = Xc[["age", "male", "proximal"]].notna().all(axis=1)
    primary = Xt.index[stage.isin(t["primary_stages"]) & complete]
    allst = Xt.index[stage.notna()]
    OUT.mkdir(parents=True, exist_ok=True)
    pd.Series(sorted(primary), name="sample_id").to_csv(OUT / "population_primary_stage2_3.csv", index=False)
    pd.Series(sorted(allst), name="sample_id").to_csv(OUT / "population_all_stages.csv", index=False)
    zero_frac = (Xt.loc[primary] == 0).mean()
    rep = {"label_free": True, "primary_tumour_samples": int(Xt.shape[0]),
           "stage_counts": stage.value_counts(dropna=False).to_dict(),
           "location_unmapped": sorted(set(loc.dropna()) - set(t["proximal_sites"]) - set(t["distal_sites"])),
           "primary_population_n": int(len(primary)), "all_stage_population_n": int(len(allst)),
           "excluded_from_primary": {"stage_not_II_III_or_missing": int((~stage.isin(t["primary_stages"])).sum()),
                                     "incomplete_clinical_among_stage_II_III": int((stage.isin(t["primary_stages"]) & ~complete).sum())},
           "shared_genes": len(genes), "genes_all_zero_in_primary": int((zero_frac == 1).sum()),
           "genes_over_half_zero_in_primary": int((zero_frac > 0.5).sum()),
           "expression_value_range": [float(Xt.min().min()), float(Xt.max().max())]}
    C.write_json(OUT / "prepare_report.json", rep)
    Xc.loc[primary].to_csv(OUT / "clinical_primary_label_free.csv")
    return rep


# ---------------------------------------------------------------- target-cohort adaptation
def adapt(est, Xt_genes):
    """Copy of a fitted gene pipeline whose StandardScaler uses TCGA population moments.

    Gene selection and coefficients stay frozen; only per-gene location/scale are re-estimated
    on the target expression (no labels)."""
    p = copy.deepcopy(M._pipe(est))
    cols = p.named_steps["var"].columns_ if "var" in p.named_steps else list(p.named_steps["sel"].columns)
    sc = p.named_steps["sc"]
    mu = Xt_genes[cols].mean(axis=0).to_numpy()
    sd = Xt_genes[cols].std(axis=0, ddof=0).to_numpy()
    zero = sd == 0
    sd[zero] = 1.0
    sc.mean_, sc.scale_, sc.var_ = mu, sd, sd ** 2
    return p, int(zero.sum())


# ---------------------------------------------------------------- freeze
def freeze(cfg):
    genes = (C.REPORTS / "shared_genes.txt").read_text().split()
    ids, y, Xg, Xc = M.load_gse(genes)
    models = M.fit_final(cfg, genes, ids, y, {"genes": Xg, "clinical": Xc})
    mdir = OUT / "frozen_models"
    mdir.mkdir(parents=True, exist_ok=True)
    entries = {}
    for name, (est, which) in models.items():
        path = mdir / f"{name}.joblib"
        joblib.dump(M._pipe(est), path)
        info = M.describe_fit(name, est)
        p = M._pipe(est)
        feats = p.named_steps["var"].columns_ if "var" in p.named_steps else (
            list(p.named_steps["sel"].columns) if "sel" in p.named_steps else list(Xc.columns))
        entries[name] = {"file": str(path.relative_to(C.ROOT)), "sha256": C.sha256_file(path), "input": which,
                         "feature_list": feats, "n_features": len(feats), **info}
    # Only the files that define Track B training, preprocessing and this evaluation are sealed.
    sealed = ["common.py", "mmr.py", "tcga.py"]
    code = {f"phase10/{n}": C.sha256_file(C.ROOT / "phase10" / n) for n in sealed}
    code[str(C.CONFIG.relative_to(C.ROOT))] = C.sha256_file(C.CONFIG)
    pop = OUT / "population_primary_stage2_3.csv"
    man = {
        "schema": "phase10-tcga-freeze-1",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "statement": "Frozen before any TCGA MSI label was linked to TCGA expression. Evaluate once; no refit, retune or model selection after unblinding.",
        "training": {"cohort": "GSE39582, Phase 9 452-patient cohort with mmr.status", "n": len(ids),
                     "dMMR": int(y.sum()), "pMMR": int((1 - y).sum()),
                     "ids_sha256": C.frame_sha(pd.DataFrame({"id": ids}))},
        "evaluation_population": {"file": str(pop.relative_to(C.ROOT)), "sha256": C.sha256_file(pop),
                                  "definition": "TCGA-COAD primary tumour (-01), pathologic stage II-III, complete age/sex/location"},
        "label_definition": {"column": cfg["tcga"]["label_column"], "positive": cfg["tcga"]["positive"],
                             "negative": cfg["tcga"]["negative"], "other_values": "excluded (missing or unexpected)"},
        "preprocessing": {"gene_space": "symbols present in both GSE39582 (median-of-probes, log2 RMA+ComBat as deposited) and TCGA HiSeqV2",
                          "shared_genes_file": "reports/phase10/shared_genes.txt",
                          "shared_genes_sha256": C.sha256_file(C.REPORTS / "shared_genes.txt"),
                          "gse_training": "VarianceFilter(top 2000, training rows) -> StandardScaler(training rows)",
                          "tcga_adaptation": "per-gene z-score with TCGA evaluation-population mean/SD (label-free); clinical features use the GSE-fitted scaler",
                          "clinical_mapping": {"stage": "Stage II* -> 0, Stage III* -> 1", "proximal_sites": cfg["tcga"]["proximal_sites"],
                                               "distal_sites": cfg["tcga"]["distal_sites"]}},
        "models": entries,
        "data_sha256": {"gse_expression": C.sha256_file(C.GSE_EXPR), "tcga_expression": C.sha256_file(EXPR),
                        "tcga_clinical": C.sha256_file(CLIN)},
        "code_sha256": code, "environment": C.environment_record(), "git_head": C.git_head(),
        "bootstrap_resamples": cfg["tcga"]["bootstrap_resamples"],
    }
    C.write_json(MANIFEST, man)
    return C.sha256_file(MANIFEST), man


# ---------------------------------------------------------------- sealed label access
def verify_freeze(require_git=True):
    """Raise SealedLabelError unless the freeze is approved, intact and (optionally) committed."""
    if not MANIFEST.exists():
        raise SealedLabelError("no freeze manifest")
    if not APPROVAL.exists():
        raise SealedLabelError("no unblinding approval file")
    msha = C.sha256_file(MANIFEST)
    appr = json.loads(APPROVAL.read_text())
    if appr.get("manifest_sha256") != msha:
        raise SealedLabelError("approval does not match the current manifest")
    man = json.loads(MANIFEST.read_text())
    checks = {e["file"]: e["sha256"] for e in man["models"].values()}
    checks.update(man["code_sha256"])
    checks[man["evaluation_population"]["file"]] = man["evaluation_population"]["sha256"]
    for rel, h in checks.items():
        if C.sha256_file(C.ROOT / rel) != h:
            raise SealedLabelError(f"frozen file changed: {rel}")
    if require_git:
        rel = str(MANIFEST.relative_to(C.ROOT))
        tracked = subprocess.run(["git", "-C", str(C.ROOT), "ls-files", "--error-unmatch", rel], capture_output=True)
        clean = subprocess.run(["git", "-C", str(C.ROOT), "diff", "--quiet", "HEAD", "--", rel])
        if tracked.returncode or clean.returncode:
            raise SealedLabelError("freeze manifest is not committed unchanged")
    return man


def load_sealed_labels(cfg, require_git=True):
    man = verify_freeze(require_git)
    col = cfg["tcga"]["label_column"]
    lab = pd.read_csv(CLIN, sep="\t", usecols=["sampleID", col], dtype=str).set_index("sampleID")[col]
    pos, neg = set(cfg["tcga"]["positive"]), set(cfg["tcga"]["negative"])
    y = lab.map(lambda v: 1 if v in pos else (0 if v in neg else np.nan))
    return man, y


# ---------------------------------------------------------------- evaluation (once)
def _metrics(y, p, clip):
    num, den = C.auc_pair_counts(y, p)
    a, b = C.calibration(y, p, clip)
    prev = float(np.mean(y))
    br = brier_score_loss(y, p)
    return {"auroc": num / den, "auprc": average_precision_score(y, p), "prevalence": prev, "brier": br,
            "scaled_brier": 1 - br / (prev * (1 - prev)), "calibration_intercept": a, "calibration_slope": b}


def evaluate(cfg, require_git=True):
    if (OUT / "evaluation.json").exists():
        raise SealedLabelError("evaluation already performed once; refusing to re-run")
    man, ylab = load_sealed_labels(cfg, require_git)
    genes = (C.REPORTS / "shared_genes.txt").read_text().split()
    pop = pd.read_csv(C.ROOT / man["evaluation_population"]["file"])["sample_id"].tolist()
    Xt = load_tcga_expression(genes).loc[pop]
    Xc = pd.read_csv(OUT / "clinical_primary_label_free.csv", index_col=0).loc[pop]
    y = ylab.reindex(pop)
    keep = y.notna()
    clip = cfg["mmr"]["prob_clip"]
    preds, res, notes = {}, {}, {}
    for name, e in man["models"].items():
        est = joblib.load(C.ROOT / e["file"])
        if e["input"] == "genes":
            p_ad, nzero = adapt(est, Xt.loc[keep[keep].index])
            notes[name] = {"genes_zero_sd_in_tcga": nzero}
            preds[name] = p_ad.predict_proba(Xt.loc[keep[keep].index])[:, 1]
        else:
            preds[name] = est.predict_proba(Xc.loc[keep[keep].index][e["feature_list"]])[:, 1]
    yy = y[keep].astype(int).to_numpy()
    for name, p in preds.items():
        res[name] = _metrics(yy, p, clip)
    B = cfg["tcga"]["bootstrap_resamples"]
    rng = np.random.default_rng(C.seed(cfg["seed_base"], 3, 1))
    pos, neg = np.flatnonzero(yy == 1), np.flatnonzero(yy == 0)
    boot = {n: [] for n in preds}
    dpair = []
    for _ in range(B):
        idx = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]
        aucs = {n: C.auc_pair_counts(yy[idx], preds[n][idx]) for n in preds}
        for n in preds:
            boot[n].append({"auroc": aucs[n][0] / aucs[n][1], "auprc": average_precision_score(yy[idx], preds[n][idx]),
                            "brier": brier_score_loss(yy[idx], preds[n][idx])})
        dpair.append(aucs["enet_logistic"][0] / aucs["enet_logistic"][1] - aucs["mlh1_only"][0] / aucs["mlh1_only"][1])
    for n in preds:
        bdf = pd.DataFrame(boot[n])
        res[n]["ci95"] = {k: [float(bdf[k].quantile(.025)), float(bdf[k].quantile(.975))] for k in bdf}
    out = {"title": "Independent-cohort, cross-platform evaluation of MMR/MSI-associated expression classification",
           "labels": "GSE39582 mmr.status (training) vs TCGA MSI_updated_Oct62011, MSI-H vs MSS+MSI-L (evaluation): different, associated constructs",
           "normalisation": "unsupervised target-cohort adaptation (TCGA population gene z-scores); not a single-patient pipeline",
           "freeze_manifest_sha256": C.sha256_file(MANIFEST), "n": int(keep.sum()), "msi_h": int(yy.sum()),
           "excluded_missing_label": int((~keep).sum()), "bootstrap": {"resamples": B, "type": "stratified patient bootstrap"},
           "models": res, "enet_minus_mlh1_auroc": {"point": res["enet_logistic"]["auroc"] - res["mlh1_only"]["auroc"],
                                                     "ci95_paired_bootstrap": [float(np.quantile(dpair, .025)), float(np.quantile(dpair, .975))]},
           "adaptation_notes": notes}
    pd.DataFrame({"sample_id": keep[keep].index, "y": yy, **{n: p for n, p in preds.items()}}).to_csv(OUT / "predictions.csv", index=False)
    C.write_json(OUT / "evaluation.json", out)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["prepare", "freeze", "evaluate"])
    a = ap.parse_args(argv)
    cfg = C.load_config()
    if a.stage == "prepare":
        print(json.dumps(prepare(cfg), indent=2))
    elif a.stage == "freeze":
        sha, man = freeze(cfg)
        print("freeze manifest sha256:", sha)
        for n, e in man["models"].items():
            print(f"  {n}: {e['n_features']} features, params={e.get('best_params')}, "
                  f"selected={len(e.get('selected', {})) if 'selected' in e else '-'}")
    else:
        print(json.dumps(evaluate(cfg), indent=2, default=float))


if __name__ == "__main__":
    main()

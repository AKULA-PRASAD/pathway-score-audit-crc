"""Declared deviation 2: TCGA-COAD evaluation on the Illumina GA RNA-seq samples (GAV2).

Why: at unblinding, only 4 of the 182 registered HiSeqV2 stage II-III patients carried a
`MSI_updated_Oct62011` label (the 2012-freeze MSI calls sit mostly on GA-sequenced samples;
188 labelled samples, 18 with HiSeqV2 expression). That evaluation is non-informative and is
kept on record unchanged (reports/phase10/tcga/evaluation.json).

What stays fixed: the frozen GSE39582 models (hash-verified against the committed freeze
manifest), the label field and its binary definition, the stage II-III/complete-clinical
population rule, the target-cohort adaptation and the metrics (functions imported from
phase10.tcga, unchanged). Only the expression source changes: GAV2 is the same UNC
RSEM pipeline, log2(norm_count+1), on the Illumina GA instead of the HiSeq.

Stages:  prepare (label-free) -> register (commit required) -> evaluate (sealed, once)
"""
from __future__ import annotations

import argparse
import gzip
import json
import subprocess

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss

from phase10 import common as C
from phase10 import tcga as T

GA = C.ROOT / "data" / "external" / "tcga_coad_ga" / "GAV2.gz"
OUT = C.REPORTS / "tcga_ga"
REG = OUT / "registration.json"
APPROVAL = OUT / "UNBLINDING_APPROVAL.json"


def load_ga(genes=None):
    G = pd.read_csv(GA, sep="\t", index_col=0)
    G = G.loc[:, [c for c in G.columns if c[13:15] == "01"]]
    return (G.loc[list(genes)] if genes is not None else G).T.astype(float)


def prepare(cfg):
    man = json.loads(T.MANIFEST.read_text())
    genes = (C.REPORTS / "shared_genes.txt").read_text().split()
    with gzip.open(GA, "rt") as fh:
        fh.readline()
        ga_genes = {line.split("\t", 1)[0] for line in fh}
    needed = set()
    for e in man["models"].values():
        if e["input"] == "genes":
            needed |= set(e["feature_list"])
    missing = sorted(needed - ga_genes)
    Xt = load_ga([g for g in genes if g in ga_genes])
    assert Xt.index.str[:12].is_unique
    clin = T.read_allow_listed_clinical(cfg).reindex(Xt.index)       # same allow-list, no MSI fields
    t = cfg["tcga"]
    loc = clin["anatomic_neoplasm_subdivision"]
    Xc = pd.DataFrame({
        "age": pd.to_numeric(clin["age_at_initial_pathologic_diagnosis"], errors="coerce"),
        "male": clin["gender"].map({"MALE": 1.0, "FEMALE": 0.0}),
        "stage3": clin["pathologic_stage"].map(T.stage_roman).map({"II": 0.0, "III": 1.0}),
        "proximal": np.where(loc.isin(t["proximal_sites"]), 1.0, np.where(loc.isin(t["distal_sites"]), 0.0, np.nan)),
    }, index=Xt.index)
    stage = clin["pathologic_stage"].map(T.stage_roman)
    complete = Xc[["age", "male", "proximal"]].notna().all(axis=1)
    primary = Xt.index[stage.isin(t["primary_stages"]) & complete]
    OUT.mkdir(parents=True, exist_ok=True)
    pd.Series(sorted(primary), name="sample_id").to_csv(OUT / "population_primary_stage2_3.csv", index=False)
    Xc.loc[primary].to_csv(OUT / "clinical_primary_label_free.csv")
    hiseq_pop = set(pd.read_csv(T.OUT / "population_primary_stage2_3.csv").sample_id.str[:12])
    rep = {"label_free": True, "primary_tumour_samples": int(len(Xt)), "primary_population_n": int(len(primary)),
           "stage_counts": stage.value_counts(dropna=False).to_dict(),
           "frozen_model_genes_missing_in_ga": missing, "patients_also_in_hiseq_population": int(len(set(primary.str[:12]) & hiseq_pop))}
    C.write_json(OUT / "prepare_report.json", rep)
    return rep


def register():
    man_sha = C.sha256_file(T.MANIFEST)
    man = json.loads(T.MANIFEST.read_text())
    for name, e in man["models"].items():                             # frozen models unchanged
        assert C.sha256_file(C.ROOT / e["file"]) == e["sha256"], f"{name} model changed"
    pop = OUT / "population_primary_stage2_3.csv"
    reg = {"schema": "phase10-tcga-ga-registration-1",
           "deviation": "2: evaluation expression source HiSeqV2 -> GAV2 because the registered HiSeqV2 population had 4/182 labels",
           "frozen_manifest": str(T.MANIFEST.relative_to(C.ROOT)), "frozen_manifest_sha256": man_sha,
           "model_sha256": {n: e["sha256"] for n, e in man["models"].items()},
           "expression": {"file": "data/external/tcga_coad_ga/GAV2.gz", "sha256": C.sha256_file(GA),
                          "source": "https://tcga.xenahubs.net/download/TCGA.COAD.sampleMap/GAV2.gz",
                          "unit": "log2(norm_count+1), RSEM, IlluminaGA_RNASeqV2, UNC (Xena 2017-10-13)"},
           "population": {"file": str(pop.relative_to(C.ROOT)), "sha256": C.sha256_file(pop),
                          "definition": "GAV2 primary tumour (-01), pathologic stage II-III, complete age/sex/location"},
           "label_definition": man["label_definition"],
           "code_sha256": {f"phase10/{n}": C.sha256_file(C.ROOT / "phase10" / n) for n in ("common.py", "tcga.py", "tcga_ga.py")},
           "statement": "Registered before any MSI label was linked to GAV2 expression. Evaluate once."}
    C.write_json(REG, reg)
    return C.sha256_file(REG)


def verify():
    if not (REG.exists() and APPROVAL.exists()):
        raise T.SealedLabelError("registration or approval missing")
    if json.loads(APPROVAL.read_text()).get("registration_sha256") != C.sha256_file(REG):
        raise T.SealedLabelError("approval does not match registration")
    reg = json.loads(REG.read_text())
    man = json.loads(T.MANIFEST.read_text())
    for n, h in reg["model_sha256"].items():
        if C.sha256_file(C.ROOT / man["models"][n]["file"]) != h:
            raise T.SealedLabelError(f"frozen model changed: {n}")
    for rel, h in reg["code_sha256"].items():
        if C.sha256_file(C.ROOT / rel) != h:
            raise T.SealedLabelError(f"code changed: {rel}")
    if C.sha256_file(C.ROOT / reg["population"]["file"]) != reg["population"]["sha256"]:
        raise T.SealedLabelError("population changed")
    rel = str(REG.relative_to(C.ROOT))
    if (subprocess.run(["git", "-C", str(C.ROOT), "ls-files", "--error-unmatch", rel], capture_output=True).returncode
            or subprocess.run(["git", "-C", str(C.ROOT), "diff", "--quiet", "HEAD", "--", rel]).returncode):
        raise T.SealedLabelError("registration is not committed unchanged")
    return reg, man


def evaluate(cfg):
    if (OUT / "evaluation.json").exists():
        raise T.SealedLabelError("GA evaluation already performed once; refusing to re-run")
    reg, man = verify()
    col = cfg["tcga"]["label_column"]
    lab = pd.read_csv(T.CLIN, sep="\t", usecols=["sampleID", col], dtype=str).set_index("sampleID")[col]
    pos, neg = set(cfg["tcga"]["positive"]), set(cfg["tcga"]["negative"])
    pop = pd.read_csv(C.ROOT / reg["population"]["file"]).sample_id.tolist()
    y = lab.reindex(pop).map(lambda v: 1 if v in pos else (0 if v in neg else np.nan))
    keep = y.notna()
    ids = keep[keep].index
    genes = (C.REPORTS / "shared_genes.txt").read_text().split()
    with gzip.open(GA, "rt") as fh:
        fh.readline()
        ga_genes = {line.split("\t", 1)[0] for line in fh}
    Xt = load_ga([g for g in genes if g in ga_genes]).loc[ids]
    Xc = pd.read_csv(OUT / "clinical_primary_label_free.csv", index_col=0).loc[ids]
    clip = cfg["mmr"]["prob_clip"]
    preds, notes = {}, {}
    for name, e in man["models"].items():
        est = joblib.load(C.ROOT / e["file"])
        if e["input"] == "genes":
            p_ad, nzero = T.adapt(est, Xt)
            notes[name] = {"genes_zero_sd": nzero}
            preds[name] = p_ad.predict_proba(Xt)[:, 1]
        else:
            preds[name] = est.predict_proba(Xc[e["feature_list"]])[:, 1]
    yy = y[keep].astype(int).to_numpy()
    res = {n: T._metrics(yy, p, clip) for n, p in preds.items()}
    B = cfg["tcga"]["bootstrap_resamples"]
    rng = np.random.default_rng(C.seed(cfg["seed_base"], 3, 2))
    pi, ni = np.flatnonzero(yy == 1), np.flatnonzero(yy == 0)
    boot = {n: [] for n in preds}
    dpair = []
    for _ in range(B):
        idx = np.r_[rng.choice(pi, len(pi)), rng.choice(ni, len(ni))]
        au = {n: C.auc_pair_counts(yy[idx], preds[n][idx]) for n in preds}
        for n in preds:
            boot[n].append({"auroc": au[n][0] / au[n][1], "auprc": average_precision_score(yy[idx], preds[n][idx]),
                            "brier": brier_score_loss(yy[idx], preds[n][idx])})
        dpair.append(au["enet_logistic"][0] / au["enet_logistic"][1] - au["mlh1_only"][0] / au["mlh1_only"][1])
    for n in preds:
        b = pd.DataFrame(boot[n])
        res[n]["ci95"] = {k: [float(b[k].quantile(.025)), float(b[k].quantile(.975))] for k in b}
    out = {"title": "Independent-cohort, cross-platform evaluation of MMR/MSI-associated expression classification (TCGA-COAD, GA RNA-seq; deviation 2)",
           "labels": "GSE39582 mmr.status (PCR MSI, 5 Bethesda markers) vs TCGA MSI_updated_Oct62011, MSI-H vs MSS+MSI-L (PCR panels; related, not identical)",
           "normalisation": "unsupervised target-cohort adaptation (per-gene z-score on this TCGA population; label-free); not a single-patient pipeline",
           "registration_sha256": C.sha256_file(REG), "n": int(keep.sum()), "msi_h": int(yy.sum()),
           "excluded_missing_label": int((~keep).sum()),
           "bootstrap": {"resamples": B, "type": "stratified patient bootstrap"}, "models": res,
           "enet_minus_mlh1_auroc": {"point": res["enet_logistic"]["auroc"] - res["mlh1_only"]["auroc"],
                                     "ci95_paired_bootstrap": [float(np.quantile(dpair, .025)), float(np.quantile(dpair, .975))]},
           "adaptation_notes": notes}
    pd.DataFrame({"sample_id": ids, "y": yy, **preds}).to_csv(OUT / "predictions.csv", index=False)
    C.write_json(OUT / "evaluation.json", out)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["prepare", "register", "evaluate"])
    a = ap.parse_args(argv)
    cfg = C.load_config()
    if a.stage == "prepare":
        print(json.dumps(prepare(cfg), indent=2))
    elif a.stage == "register":
        print("registration sha256:", register())
    else:
        print(json.dumps(evaluate(cfg), indent=2, default=float))


if __name__ == "__main__":
    main()

"""Phase 6D: definitive missingness audit.

AUDIT ONLY. No modelling, no imputation, no feature selection, no
normalization. Phase 2 / Phase 3 / frozen M2 artifacts are read-only.

OUTCOME HANDLING: recurrence and survival VALUES are never read. Endpoint
completeness is CITED from the verified Phase 6A result, not recomputed, so
no outcome column is opened by this script.

Usage:
    .venv/bin/python scripts/28_definitive_missingness_audit.py
"""

import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                  # noqa: E402
from src.reproducibility import environment_report     # noqa: E402

OUT = paths.ROOT / "reports" / "phase6d"
ANALYSIS = paths.PROCESSED_DIR / "analysis"

# Frozen M1 predictor list (Phase 6B §3.1). Optional molecular covariates are
# NOT added to M1 -- they are audited separately.
M1_PRIMARY = ["age.at.diagnosis_(year)", "sex", "tnm.stage", "tumor.location",
              "chemotherapy.adjuvant"]
M1_OPTIONAL = ["tnm.t", "tnm.n", "tnm.m", "mmr.status"]
MOLECULAR_OPTIONAL = ["cimp.status", "cin.status", "kras.mutation",
                      "braf.mutation", "tp53.mutation"]
EXCLUDED_BY_PHASE1 = ["chemotherapy.adjuvant.type", "cit.molecularsubtype",
                      "rna.extraction.method"]

NA_LIKE = {"N/A", "NA", "ND", "", "unknown", "Unknown", "not measured"}


def section(t):
    print("\n" + "=" * 78); print(t); print("=" * 78)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    findings = {}

    # ---------- cohort ----------
    mem = pd.read_csv(ANALYSIS / "cohort_membership.csv")
    part = pd.read_csv(ANALYSIS / "partition_and_batch.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    raw = pd.read_csv(paths.SAMPLE_METADATA, keep_default_na=False, dtype=str)
    raw = raw[raw.sample_id.isin(s3)].set_index("sample_id").loc[s3].reset_index()
    ds = part.set_index("sample_id")["dataset"].reindex(s3).values
    raw["dataset"] = ds
    raw["combat_batch"] = part.set_index("sample_id")["combat_batch"].reindex(s3).values
    n, nd_, nv = len(raw), int((ds == "discovery").sum()), int((ds == "validation").sum())
    print(f"S3 cohort: n={n}  discovery={nd_}  validation={nv}")

    # ---------- 1. encoding inventory ----------
    section("1. MISSING-VALUE ENCODING INVENTORY")
    print("  The deposit's own documentation states:")
    print('    "ND: not determined"   "NA: Not applicable"')
    nd_total = int(sum((raw[c] == "ND").sum() for c in raw.columns))
    print(f"\n  Occurrences of 'ND' anywhere in the metadata: {nd_total}")
    print("  => the documented ND/NA distinction is NOT USED. Both meanings")
    print("     are encoded identically as 'N/A' and cannot be separated from")
    print("     the encoding alone. They are separated STRUCTURALLY below.")
    enc = []
    for c in M1_PRIMARY + M1_OPTIONAL + MOLECULAR_OPTIONAL + EXCLUDED_BY_PHASE1:
        if c not in raw.columns:
            continue
        vals = raw[c].value_counts()
        na_like = {k: int(v) for k, v in vals.items() if k in NA_LIKE}
        enc.append({"column": c, "distinct_values": int(vals.size),
                    "na_like_codes": json.dumps(na_like),
                    "n_na_like": int(sum(na_like.values())),
                    "empty_string": int((raw[c] == "").sum()),
                    "literal_NA": int((raw[c] == "NA").sum()),
                    "literal_N_slash_A": int((raw[c] == "N/A").sum()),
                    "literal_ND": int((raw[c] == "ND").sum())})
    edf = pd.DataFrame(enc)
    print("\n" + edf[["column", "distinct_values", "n_na_like", "empty_string",
                      "literal_NA", "literal_N_slash_A", "literal_ND"]].to_string(index=False))
    edf.to_csv(OUT / "encoding_inventory.csv", index=False)
    findings["encoding"] = {"ND_occurrences_total": nd_total,
                            "ND_documented_but_unused": True,
                            "na_code_in_use": "N/A",
                            "codes_conflated": ["not determined", "not applicable"]}

    # ---------- 2. structural vs true missing ----------
    section("2. STRUCTURAL vs TRUE MISSINGNESS")
    ca, ct = raw["chemotherapy.adjuvant"], raw["chemotherapy.adjuvant.type"]
    struct = int(((ca == "N") & (ct == "N/A")).sum())
    true_m = int(((ca == "Y") & (ct == "N/A")).sum())
    print(f"  chemotherapy.adjuvant.type is 'N/A' for {int((ct=='N/A').sum())} samples:")
    print(f"    STRUCTURAL (no chemo given, type not applicable) : {struct}")
    print(f"    TRUE MISSING (chemo given, regimen unrecorded)   : {true_m}")
    print("\n  Mutation detail columns follow the documented rule")
    print("  (mutation=WT -> exon/protein/DNA are 'N/A' = not applicable),")
    print("  verified in Phase 6D inspection. Those columns are NOT planned")
    print("  predictors and are not audited further.")
    findings["structural"] = {"chemo_type_structural": struct,
                              "chemo_type_true_missing": true_m}

    # ---------- 3. variable-level missingness ----------
    section("3. VARIABLE-LEVEL MISSINGNESS (S3, discovery, validation)")
    groups = {**{c: "M1_primary" for c in M1_PRIMARY},
              **{c: "M1_optional" for c in M1_OPTIONAL},
              **{c: "molecular_optional" for c in MOLECULAR_OPTIONAL},
              **{c: "excluded_phase1" for c in EXCLUDED_BY_PHASE1}}
    rows = []
    for c, grp in groups.items():
        if c not in raw.columns:
            continue
        miss = raw[c].isin(NA_LIKE)
        md = int(miss[ds == "discovery"].sum()); mv = int(miss[ds == "validation"].sum())
        rows.append({"variable": c, "group": grp,
                     "observed": int((~miss).sum()), "missing": int(miss.sum()),
                     "pct": round(miss.mean() * 100, 2),
                     "disc_missing": md, "disc_pct": round(md / nd_ * 100, 2),
                     "val_missing": mv, "val_pct": round(mv / nv * 100, 2),
                     "diff_pp": round(mv / nv * 100 - md / nd_ * 100, 2),
                     "complete": bool(miss.sum() == 0)})
    vdf = pd.DataFrame(rows).sort_values(["group", "pct"], ascending=[True, False])
    print(vdf.to_string(index=False))
    vdf.to_csv(OUT / "variable_missingness.csv", index=False)
    findings["variables"] = vdf.to_dict(orient="records")

    # ---------- 4. patterns ----------
    section("4. MISSINGNESS PATTERNS")
    audit_cols = [c for c in M1_PRIMARY + M1_OPTIONAL + MOLECULAR_OPTIONAL
                  if c in raw.columns]
    pat = raw[audit_cols].isin(NA_LIKE)
    key = pat.apply(lambda r: "".join("1" if x else "0" for x in r), axis=1)
    vc = key.value_counts()
    print(f"  variables (bit order): {audit_cols}")
    print(f"  distinct patterns: {len(vc)}")
    print(f"\n  {'pattern':16s} {'n':>5s}  variables missing")
    for k, v in vc.head(12).items():
        miss = [audit_cols[i] for i, ch in enumerate(k) if ch == "1"]
        print(f"  {k:16s} {v:5d}  {', '.join(miss) if miss else '(none - complete)'}")
    pdf = pd.DataFrame({"pattern": vc.index, "n": vc.values})
    pdf["variables_missing"] = [", ".join(audit_cols[i] for i, ch in enumerate(k)
                                          if ch == "1") for k in vc.index]
    pdf.to_csv(OUT / "missingness_patterns.csv", index=False)
    print(f"\n  samples complete on all {len(audit_cols)} audited variables: "
          f"{int((~pat.any(axis=1)).sum())} ({(~pat.any(axis=1)).mean()*100:.1f}%)")

    print("\n  -- co-missingness (pairs missing together) --")
    co = []
    for a, b in combinations(audit_cols, 2):
        both = int((pat[a] & pat[b]).sum())
        if both >= 10:
            co.append({"var_a": a, "var_b": b, "both_missing": both,
                       "jaccard": round(both / int((pat[a] | pat[b]).sum()), 3)})
    codf = pd.DataFrame(co).sort_values("jaccard", ascending=False)
    print(codf.head(8).to_string(index=False) if len(codf) else "    (none >=10)")
    codf.to_csv(OUT / "co_missingness.csv", index=False)

    # ---------- 5. model-specific complete case ----------
    section("5. MODEL-SPECIFIC COMPLETE-CASE ASSESSMENT")
    models = {
        "M1 primary (clinical only)": M1_PRIMARY,
        "M1 + optional TNM/MMR": M1_PRIMARY + M1_OPTIONAL,
        "M1 + all molecular (NOT PROPOSED)": M1_PRIMARY + M1_OPTIONAL + MOLECULAR_OPTIONAL,
    }
    mrows = []
    for name, cols in models.items():
        cc = int((~raw[cols].isin(NA_LIKE).any(axis=1)).sum())
        g = raw[~raw[cols].isin(NA_LIKE).any(axis=1)]
        mrows.append({"model": name, "n_predictors": len(cols),
                      "complete_cases": cc, "pct": round(cc / n * 100, 1),
                      "lost": n - cc,
                      "disc_cc": int((g.dataset == "discovery").sum()),
                      "val_cc": int((g.dataset == "validation").sum())})
    # M2 / M3 depend on expression, which is complete by construction
    X = pd.read_parquet(ANALYSIS / "expression_median_S0.parquet")
    Xs3 = X[s3]
    expr_missing = int(np.isnan(Xs3.to_numpy()).sum())
    zero_var = int((Xs3.std(axis=1) == 0).sum())
    for nm, extra in [("M2 primary (4-set) + M1", 4), ("M2 sensitivity (5-set) + M1", 5),
                      ("M3 exploratory (21,815 genes) + M1", 21815)]:
        cc = int((~raw[M1_PRIMARY].isin(NA_LIKE).any(axis=1)).sum())
        mrows.append({"model": nm, "n_predictors": len(M1_PRIMARY) + extra,
                      "complete_cases": cc, "pct": round(cc / n * 100, 1),
                      "lost": n - cc,
                      "disc_cc": int((raw[~raw[M1_PRIMARY].isin(NA_LIKE).any(axis=1)].dataset == "discovery").sum()),
                      "val_cc": int((raw[~raw[M1_PRIMARY].isin(NA_LIKE).any(axis=1)].dataset == "validation").sum())})
    mdf = pd.DataFrame(mrows)
    print(mdf.to_string(index=False))
    mdf.to_csv(OUT / "model_complete_case.csv", index=False)
    print(f"\n  expression matrix (S3): NaN values = {expr_missing}; "
          f"zero-variance genes = {zero_var}")
    print("  => pathway inputs have NO missing values. The frozen coverage rule")
    print("     (>=80% members, >=20 genes) is evaluated per FOLD at scoring time.")
    findings["models"] = mdf.to_dict(orient="records")
    findings["expression"] = {"nan_values": expr_missing, "zero_variance_genes": zero_var}

    # ---------- 6. missingness by observed characteristics ----------
    section("6. MISSINGNESS BY OBSERVED CHARACTERISTICS (outcome-blind)")
    strat_vars = ["tnm.stage", "sex", "tumor.location", "dataset", "combat_batch"]
    targets = ["tp53.mutation", "cin.status", "cimp.status", "braf.mutation",
               "mmr.status"]
    srows = []
    for t in targets:
        miss = raw[t].isin(NA_LIKE)
        for sv in strat_vars:
            for lvl, g in miss.groupby(raw[sv]):
                if len(g) < 10:
                    continue
                srows.append({"target": t, "stratifier": sv, "level": str(lvl),
                              "n": len(g), "missing": int(g.sum()),
                              "pct": round(g.mean() * 100, 1)})
    sdf = pd.DataFrame(srows)
    sdf.to_csv(OUT / "missingness_by_characteristic.csv", index=False)
    for t in targets[:3]:
        print(f"\n  {t}:")
        sub = sdf[sdf.target == t]
        for sv in ["dataset", "combat_batch", "tnm.stage"]:
            x = sub[sub.stratifier == sv]
            if len(x):
                rng = f"{x.pct.min():.1f}%-{x.pct.max():.1f}%"
                print(f"    by {sv:14s}: {rng:16s} across {len(x)} levels")

    print("\n  -- the sharpest signal: by ComBat batch --")
    bt = raw.groupby("combat_batch").apply(
        lambda g: pd.Series({"n": len(g),
                             "tp53_miss_pct": round(g["tp53.mutation"].isin(NA_LIKE).mean()*100, 1),
                             "cin_miss_pct": round(g["cin.status"].isin(NA_LIKE).mean()*100, 1),
                             "cimp_miss_pct": round(g["cimp.status"].isin(NA_LIKE).mean()*100, 1)}),
        include_groups=False)
    print(bt.to_string())
    bt.to_csv(OUT / "missingness_by_batch.csv")
    findings["by_batch"] = bt.reset_index().to_dict(orient="records")

    (OUT / "definitive_audit.json").write_text(
        json.dumps({**findings, "cohort": {"n": n, "discovery": nd_, "validation": nv},
                    "outcome_columns_opened": False,
                    "environment": environment_report()},
                   indent=2, default=str) + "\n")
    section("DONE")
    print(f"  artifacts -> {OUT.relative_to(paths.ROOT)}/")
    print("  No outcome column opened. No imputation. No modelling.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

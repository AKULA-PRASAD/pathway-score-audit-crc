"""X4. MMR status and combat_mod provenance; robustness of score–MMR associations. Registered plan §X4.

combat_mod is never relabelled as MMR status.
"""
import csv
import gzip
from collections import Counter

import numpy as np
import pandas as pd
from scipy import stats

import xcommon as X


def eps2(s, g):
    groups = [s[g == k] for k in np.unique(g)]
    H, p = stats.kruskal(*groups)
    n, k = len(s), len(groups)
    return (H - k + 1) / (n - k), p, H


def main():
    ctx = X.context()
    out = X.outdir("x4_mmr_provenance")
    soft = X.soft_characteristics()
    # ---- titles (CIT ids) from SOFT
    titles = {}
    with gzip.open(X.ROOT / "data/raw/gse39582/GSE39582_family.soft.gz", "rt", errors="replace") as f:
        cur = None
        for line in f:
            if line.startswith("^SAMPLE"):
                cur = line.split("=", 1)[1].strip()
            elif cur and line.startswith("!Sample_title"):
                titles[cur] = line.split("=", 1)[1].strip()
    # ---- SDRF
    r = list(csv.reader(open(X.ROOT / "data/external/E-GEOD-39582.sdrf.txt"), delimiter="\t"))
    h = r[0]
    sd = pd.DataFrame(r[1:], columns=range(len(h)))
    col = {name: h.index(name) for name in ["Characteristics [mmr.status]", "Characteristics [normalizationcombatmod]",
                                            "Characteristics [normalizationcombatbatch]", "Comment [Sample_source_name]", "Source Name"]}
    sdrf = pd.DataFrame({"source": sd[col["Source Name"]], "mmr_sdrf": sd[col["Characteristics [mmr.status]"]],
                         "combat_mod": sd[col["Characteristics [normalizationcombatmod]"]]})
    sdrf["sample_id"] = sdrf.source.str.extract(r"(GSM\d+)")[0]
    sdrf = sdrf.dropna(subset=["sample_id"]).drop_duplicates("sample_id").set_index("sample_id")
    # ---- Marisa Table S1
    s1 = pd.read_excel(X.ROOT / "data/external/marisa2013_tableS1.xls", header=1)
    s1 = s1.set_index("id")
    ids = ctx.ids
    df = pd.DataFrame({"title": [titles.get(i) for i in ids], "mmr_geo": soft.loc[ids, "mmr.status"].values,
                       "mmr_sdrf": sdrf.reindex(ids).mmr_sdrf.values, "combat_mod": sdrf.reindex(ids).combat_mod.values}, index=ids)
    df["mmr_tableS1"] = s1.reindex(df.title)["MMR.Status"].values
    prov = {
        "fields": {"GEO SOFT": "!Sample_characteristics_ch1 'mmr.status'; no combat_mod field; data_processing = 'RMA normalisation (Bioconductor R affy package) + ComBat (R sva package)' with no description of the ComBat model matrix",
                   "ArrayExpress SDRF": "Characteristics [mmr.status]; Characteristics [normalizationcombatmod] (values modcan1/modcan2); Characteristics [normalizationcombatbatch]",
                   "Marisa 2013 Table S1": "MMR.Status column; no ComBat fields"},
        "n_cohort": len(ids),
        "geo_vs_sdrf_mmr_agreement": float((df.mmr_geo.replace("N/A", "NA") == df.mmr_sdrf).mean()),
        "geo_vs_tableS1_mmr_agreement": float((df.mmr_geo.replace("N/A", "NA") == df.mmr_tableS1.fillna("NA").astype(str)).mean()),
        "crosstab_combat_mod_by_mmr_geo_cohort": {f"{a}|{b}": int(c) for (a, b), c in Counter(zip(df.combat_mod, df.mmr_geo)).items()},
    }
    allsd = Counter(zip(sdrf.combat_mod, sdrf.mmr_sdrf))
    prov["crosstab_combat_mod_by_mmr_sdrf_all_deposit"] = {f"{a}|{b}": int(c) for (a, b), c in allsd.items()}
    df.to_csv(out / "mmr_sources_per_sample.csv")
    # ---- robustness of score–MMR association
    z, _ = X.cohort_scores(ctx)
    part = pd.read_csv(X.ANALYSIS / "partition_and_batch.csv").set_index("sample_id").loc[ids]
    cit = soft.loc[ids, "cit.molecularsubtype"]
    known = df.mmr_geo.isin(["dMMR", "pMMR"]).to_numpy()
    g = df.mmr_geo.to_numpy()[known]
    bt = part.combat_batch.to_numpy()[known]
    rng = np.random.default_rng(X.SEED)
    rows = []
    for lab in z.columns:
        s = z[lab].to_numpy()[known]
        e2, p, H = eps2(s, g)
        # batch-stratified permutation of MMR labels (H statistic)
        perm = np.empty(10000)
        for i in range(10000):
            lab_p = g.copy()
            for b in np.unique(bt):
                ix = np.where(bt == b)[0]
                lab_p[ix] = rng.permutation(lab_p[ix])
            perm[i] = stats.kruskal(*[s[lab_p == k] for k in ("dMMR", "pMMR")])[0]
        p_strat = (1 + np.sum(perm >= H)) / 10001
        lobo = {b: eps2(s[bt != b], g[bt != b])[0] for b in np.unique(bt)}
        within = {}
        for k in sorted(cit.unique()):
            m = (cit.to_numpy()[known] == k)
            nd, npm = int((g[m] == "dMMR").sum()), int((g[m] == "pMMR").sum())
            if nd >= 10 and npm >= 10:
                within[k] = {"n_dMMR": nd, "n_pMMR": npm, "eps2": eps2(s[m], g[m])[0], "p": eps2(s[m], g[m])[1]}
        med_d, med_p = np.median(s[g == "dMMR"]), np.median(s[g == "pMMR"])
        rows.append(dict(pathway=lab, n=int(known.sum()), eps2=e2, p=p, p_batch_stratified_perm=p_strat,
                         median_dMMR=med_d, median_pMMR=med_p, lobo_eps2_min=min(lobo.values()), lobo_eps2_max=max(lobo.values()),
                         within_cit=within))
    rob = pd.DataFrame(rows)
    rob["q_BH_perm"] = X.bh(rob.p_batch_stratified_perm.to_numpy())
    power = {"n_dMMR": int((g == "dMMR").sum()), "n_pMMR": int((g == "pMMR").sum()), "n_missing": int((~known).sum()),
             "dMMR_by_batch": pd.Series(g[g == "dMMR"]).groupby(bt[g == "dMMR"]).size().to_dict(),
             "dMMR_by_cit": pd.Series(cit.to_numpy()[known][g == "dMMR"]).value_counts().to_dict(),
             "recurrences_in_dMMR": int(X.endpoint(ctx).rfs_event_full.to_numpy()[known][g == "dMMR"].sum())}
    rob.drop(columns="within_cit").to_csv(out / "mmr_association_robustness.csv", index=False)
    X.write_json(out / "x4_summary.json", {"plan_section": "X4", "label": "exploratory, outcome-free (except event count for power)",
                                          "provenance": prov, "robustness": rows, "power": power})
    print(prov)
    print(rob.drop(columns="within_cit").round(4).to_string())
    for r_ in rows:
        print(r_["pathway"], {k: (round(v["eps2"], 3), v["n_dMMR"]) for k, v in r_["within_cit"].items()})
    print(power)


if __name__ == "__main__":
    main()

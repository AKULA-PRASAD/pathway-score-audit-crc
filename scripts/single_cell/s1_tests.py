"""Single-cell steps 1-5: QC gate, then the registered questions S1-S5 (docs/SINGLE_CELL_EXTENSION_PLAN.md).

Outcome-free. Reads the aggregates written by s0_build.py. Decision rules are applied exactly as registered:
  S1/S2  stromal vs epithelial enrichment of set transcripts (paired Wilcoxon in discovery; direction in replication)
  S3     within-malignant-epithelium coherence vs matched random sets, with and without an ambient-RNA control
  S4     pseudo-bulk EMT score vs counted stromal fraction (Spearman, bootstrap CI, random-set baseline)
  S5     compartment sources of the E2F, Wnt and p53 sets (descriptive)
Writes results/single_cell/*.csv and summary.json.
"""
import json

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

import sccommon as C


def load(acc):
    d = C.PROC / acc
    groups = pd.read_csv(d / "groups.csv", dtype={"patient": str})
    genes = pd.read_csv(d / "genes.csv")["gene"].astype(str).to_numpy()
    sums = pd.DataFrame(np.load(d / "group_sums.npy").T, index=genes).groupby(level=0).sum()  # genes x groups, duplicates summed
    det = pd.Series(np.load(d / "tumour_detected_cells.npy"), index=genes).groupby(level=0).sum()
    summ = json.loads((d / "build_summary.json").read_text())
    return groups, sums, det, summ


def compartment_matrix(groups, sums, cls, patient=None):
    """genes x compartments UMI sums for one sample class (and optionally one patient)."""
    g = groups[(groups["cls"] == cls) & ((groups["patient"] == patient) if patient is not None else True)]
    out = pd.DataFrame(0.0, index=sums.index, columns=C.COMPARTMENTS)
    for comp in C.COMPARTMENTS:
        cols = g.index[g["compartment"] == comp]
        if len(cols):
            out[comp] = sums.iloc[:, cols].sum(axis=1).to_numpy()
    return out


def qc_gate(acc, groups, sums):
    t = groups[groups["cls"] == "tumour"]
    counts = t.pivot_table(index="patient", columns="compartment", values="n_cells", aggfunc="sum", fill_value=0)
    for comp in C.COMPARTMENTS:
        if comp not in counts:
            counts[comp] = 0
    retained = sorted(counts.index[(counts["epithelial"] >= C.MIN_EPI) & (counts["stromal"] >= C.MIN_STROMAL)])
    cm = compartment_matrix(groups, sums, "tumour")
    cpm = cm / cm.sum(axis=0) * 1e6
    markers = {}
    for comp, gl in C.QC_MARKERS.items():
        for g in gl:
            if g in cpm.index:
                markers[g] = {"expected": comp, "top": cpm.loc[g].idxmax(), "pass": bool(cpm.loc[g].idxmax() == comp)}
    comps_present = [c for c in C.COMPARTMENTS if counts[c].sum() > 0]
    each_comp_has_marker = all(any(m["expected"] == c for m in markers.values()) for c in C.COMPARTMENTS)
    passed = (len(comps_present) == 4) and each_comp_has_marker and all(m["pass"] for m in markers.values()) and len(retained) >= 3
    return {"accession": acc, "passed": bool(passed), "compartments_present": comps_present, "markers": markers,
            "n_patients_tumour": int(len(counts)), "n_patients_retained": len(retained), "retained": retained,
            "cells_per_patient": counts[C.COMPARTMENTS].astype(int).to_dict(orient="index")}


def expressed_pool(det, summ, sums, groups):
    pool = det.index[det >= 0.01 * summ["n_tumour_cells"]]
    t = groups.index[groups["cls"] == "tumour"]
    tot = sums.iloc[:, t].sum(axis=1)
    cpm = (tot / tot.sum() * 1e6).loc[pool]
    dec = pd.Series(pd.qcut(cpm.rank(method="first"), 10, labels=False), index=pool)
    return list(pool), dec


def random_sets(target, pool, dec, rng, n=C.N_RANDOM, exclude=()):
    """n random gene sets matched to the target on size and mean-expression decile; returns a sparse indicator (n x len(pool))."""
    excl = set(target) | set(exclude)
    pos = {g: i for i, g in enumerate(pool)}
    by_dec = {d: np.array([pos[g] for g in dec.index[dec == d] if g not in excl]) for d in range(10)}
    need = dec.loc[[g for g in target if g in pos]].value_counts()
    rows, cols = [], []
    for r in range(n):
        for d, k in need.items():
            pick = rng.choice(by_dec[d], int(k), replace=False)
            rows += [r] * len(pick)
            cols += list(pick)
    return sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, len(pool)))


def enrichment(acc, groups, sums, pool, dec, sets, retained, rng):
    rows = []
    for label, gl in sets.items():
        target = [g for g in gl if g in set(pool)]
        R = random_sets(target, pool, dec, rng)
        for p in retained:
            X = compartment_matrix(groups, sums, "tumour", p).loc[pool]
            obs = X.loc[target].sum(axis=0)
            obs_share = obs / obs.sum()
            rnd = R @ X.to_numpy()
            rnd_share = (rnd / rnd.sum(axis=1, keepdims=True)).mean(axis=0)
            for i, comp in enumerate(C.COMPARTMENTS):
                e = obs_share[comp] / rnd_share[i] if rnd_share[i] > 0 else np.nan
                rows.append({"accession": acc, "set": label, "n_genes": len(target), "patient": p, "compartment": comp,
                             "share": obs_share[comp], "random_share": rnd_share[i], "enrichment": e})
    return pd.DataFrame(rows)


def paired_test(df, label):
    w = df[df["set"] == label].pivot(index="patient", columns="compartment", values="enrichment")
    diff = (w["stromal"] - w["epithelial"]).dropna()
    p = stats.wilcoxon(diff, alternative="two-sided").pvalue if len(diff) >= 3 and (diff != 0).any() else np.nan
    return {"n_patients": int(len(diff)), "median_stromal_enrichment": float(w["stromal"].median()),
            "median_epithelial_enrichment": float(w["epithelial"].median()), "median_difference": float(diff.median()),
            "wilcoxon_p": float(p) if p == p else None}


def normal_stromal_epithelial_ratio(groups, sums, pool):
    cm = compartment_matrix(groups, sums, "normal").loc[pool]
    if cm["stromal"].sum() == 0 or cm["epithelial"].sum() == 0:
        return None
    cpm = cm / cm.sum(axis=0) * 1e6
    return (cpm["stromal"] + 1) / (cpm["epithelial"] + 1)


def coherence(acc, sets, pool, dec, ambient_excl, retained, rng):
    d = C.PROC / acc
    E = sp.load_npz(d / "s3_malignant_epi_counts.npz").tocsr()
    meta = pd.read_csv(d / "s3_malignant_epi_cells.csv", dtype={"patient": str})
    genes = pd.read_csv(d / "genes.csv")["gene"].astype(str).to_numpy()
    first = pd.Series(np.arange(len(genes)), index=genes)
    # duplicate symbols: sum their columns
    uniq, inv = np.unique(genes, return_inverse=True)
    agg = sp.csr_matrix((np.ones(len(genes)), (np.arange(len(genes)), inv)), shape=(len(genes), len(uniq)))
    E = (E @ agg).tocsr()
    col = pd.Series(np.arange(len(uniq)), index=uniq)
    del first
    lib = meta["library_size"].to_numpy()
    rows = []
    for variant, excl in (("all genes", set()), ("ambient-RNA control", set(ambient_excl))):
        target_all = [g for g in sets["EMT"] if g in set(pool) and g not in excl]
        vpool = [g for g in pool if g not in excl]
        R = random_sets(target_all, vpool, dec.loc[vpool], rng)
        for p, sub in meta.groupby("patient", sort=True):
            if p not in retained or len(sub) < C.MIN_EPI:
                continue
            idx = sub.index.to_numpy()
            X = E[idx][:, col.loc[vpool].to_numpy()].toarray()
            X = np.log1p(X / lib[idx, None] * 1e4)
            Rk = stats.rankdata(X, axis=0)  # average ranks across cells, per gene
            sd = Rk.std(axis=0)
            ok = sd > 0
            Z = np.zeros_like(Rk)
            Z[:, ok] = (Rk[:, ok] - Rk[:, ok].mean(axis=0)) / sd[ok]
            pos = {g: i for i, g in enumerate(vpool)}
            tgt = np.array([pos[g] for g in target_all if ok[pos[g]]])
            n = Z.shape[0]

            def mean_r(sel_cols):
                k = len(sel_cols)
                s = Z[:, sel_cols].sum(axis=1)
                return (s @ s / n - k) / (k * (k - 1))

            obs = mean_r(tgt)
            Rok = R[:, ok]
            Zok = Z[:, ok]
            S = (Rok @ Zok.T)  # n_random x cells: sum of z over each random set (constant genes dropped)
            k_r = np.asarray(Rok.sum(axis=1)).ravel()
            rnd = ((np.asarray(S) ** 2).sum(axis=1) / n - k_r) / (k_r * (k_r - 1))
            rows.append({"accession": acc, "variant": variant, "patient": p, "n_cells": int(n), "n_genes": int(len(tgt)),
                         "coherence": float(obs), "random_median": float(np.median(rnd)), "random_p95": float(np.percentile(rnd, 95)),
                         "exceeds_p95": bool(obs > np.percentile(rnd, 95))})
    return pd.DataFrame(rows)


def pseudobulk(acc, groups, sums, pool, dec, sets, retained, rng):
    t = groups[(groups["cls"] == "tumour") & groups["patient"].isin(retained)]
    P = pd.DataFrame({p: sums.iloc[:, t.index[t["patient"] == p]].sum(axis=1) for p in retained})
    L = np.log2(P / P.sum(axis=0) * 1e6 + 1).loc[pool]
    sd = L.std(axis=1)
    Zg = L.sub(L.mean(axis=1), axis=0).div(sd.replace(0, np.nan), axis=0)
    cells = t.pivot_table(index="patient", columns="compartment", values="n_cells", aggfunc="sum", fill_value=0)
    frac = (cells.get("stromal", 0) / cells.sum(axis=1)).reindex(retained)
    out = {"patients": retained, "stromal_fraction": frac.to_dict()}
    rows = []
    for label in ("EMT", "EMT-ex6"):
        target = [g for g in sets[label] if g in Zg.index and sd[g] > 0]
        score = Zg.loc[target].mean(axis=0)
        rho, p = stats.spearmanr(score, frac)
        boot = []
        for _ in range(2000):
            i = rng.integers(0, len(retained), len(retained))
            if len(set(i)) > 2:
                boot.append(stats.spearmanr(score.iloc[i], frac.iloc[i])[0])
        R = random_sets(target, pool, dec, rng)
        Zv = Zg.fillna(0).to_numpy()
        kz = np.asarray(R.sum(axis=1)).ravel()
        rs = np.asarray(R @ Zv) / kz[:, None]
        base = [stats.spearmanr(r, frac)[0] for r in rs]
        rows.append({"accession": acc, "set": label, "n_patients": len(retained), "n_genes": len(target), "rho": float(rho), "p": float(p),
                     "ci_lo": float(np.nanpercentile(boot, 2.5)), "ci_hi": float(np.nanpercentile(boot, 97.5)),
                     "random_median": float(np.nanmedian(base)), "random_p97_5": float(np.nanpercentile(base, 97.5))})
        out[f"score_{label}"] = score.to_dict()
    return pd.DataFrame(rows), out


def holm(pvals):
    items = sorted([(p, k) for k, p in pvals.items() if p is not None])
    m, adj, run = len(items), {}, 0.0
    for i, (p, k) in enumerate(items):
        run = max(run, min(1.0, (m - i) * p))
        adj[k] = run
    return adj


def main():
    sets, _ = C.gene_sets()
    rng = np.random.default_rng(C.SEED)
    C.RES.mkdir(parents=True, exist_ok=True)
    gates, enr, coh, pb, pb_detail = {}, [], [], [], {}
    for acc in C.DATASETS:
        groups, sums, det, summ = load(acc)
        gate = qc_gate(acc, groups, sums)
        gates[acc] = gate
        print(acc, "QC", "PASS" if gate["passed"] else "FAIL", gate["n_patients_retained"], "patients retained", flush=True)
        if not gate["passed"]:
            continue
        pool, dec = expressed_pool(det, summ, sums, groups)
        enr.append(enrichment(acc, groups, sums, pool, dec, sets, gate["retained"], rng))
        ratio = normal_stromal_epithelial_ratio(groups, sums, pool)
        amb = list(ratio.index[ratio >= ratio.quantile(0.9)]) if ratio is not None else []
        gates[acc]["ambient_control_genes_excluded"] = len(amb)
        gates[acc]["emt_genes_excluded_by_ambient_control"] = sorted(set(amb) & set(sets["EMT"]))
        coh.append(coherence(acc, sets, pool, dec, amb, set(gate["retained"]), rng))
        r, detail = pseudobulk(acc, groups, sums, pool, dec, sets, gate["retained"], rng)
        pb.append(r)
        pb_detail[acc] = detail
        print(acc, "done", flush=True)
    enr = pd.concat(enr, ignore_index=True)
    coh = pd.concat(coh, ignore_index=True)
    pb = pd.concat(pb, ignore_index=True)
    enr.to_csv(C.RES / "s1_s2_s5_enrichment.csv", index=False)
    coh.to_csv(C.RES / "s3_coherence.csv", index=False)
    pb.to_csv(C.RES / "s4_pseudobulk.csv", index=False)
    C.write_json(C.RES / "qc_gate.json", gates)
    C.write_json(C.RES / "s4_pseudobulk_detail.json", pb_detail)
    passing = [a for a in C.DATASETS if gates[a]["passed"]]
    disc, reps = "GSE132465", [a for a in passing if a != "GSE132465"]
    res = {"label": "registered single-cell extension S1-S5; outcome-free", "cohorts_passing_qc": passing}
    for q, label in (("S1", "EMT"), ("S2", "EMT-ex6")):
        per = {a: paired_test(enr[enr["accession"] == a], label) for a in passing}
        d = per.get(disc)
        ok = bool(d and d["wilcoxon_p"] is not None and d["wilcoxon_p"] < 0.05 and d["median_difference"] > 0
                  and any(per[a]["median_difference"] > 0 for a in reps))
        res[q] = {"per_cohort": per, "supported": ok}
    s3 = {}
    for variant in ("all genes", "ambient-RNA control"):
        v = coh[coh["variant"] == variant]
        s3[variant] = {a: {"n_patients": int((v["accession"] == a).sum()),
                           "fraction_exceeding_p95": float(v.loc[v["accession"] == a, "exceeds_p95"].mean()),
                           "median_coherence": float(v.loc[v["accession"] == a, "coherence"].median()),
                           "median_random_p95": float(v.loc[v["accession"] == a, "random_p95"].median())} for a in passing}
    detected = all(s3[v][disc]["fraction_exceeding_p95"] > 0.5 and any(s3[v][a]["fraction_exceeding_p95"] > 0.5 for a in reps)
                   for v in s3)
    nd = s3["all genes"][disc]
    s3_p = stats.binomtest(int(round(nd["fraction_exceeding_p95"] * nd["n_patients"])), nd["n_patients"], 0.05, alternative="greater").pvalue
    res["S3"] = {"per_variant": s3, "verdict": "coherent tumour-cell EMT-like programme detected" if detected else "not detected at this resolution",
                 "binomial_p_discovery_vs_5pct": float(s3_p)}
    s4d = pb[(pb["accession"] == disc) & (pb["set"] == "EMT")].iloc[0]
    res["S4"] = {"rows": pb.to_dict(orient="records"), "supported": bool(s4d["rho"] > 0 and s4d["ci_lo"] > 0)}
    s5 = enr[enr["set"].isin(["E2F", "Wnt", "p53", "EMT"])].groupby(["accession", "set", "compartment"])["enrichment"].median().unstack()
    res["S5_median_enrichment"] = {f"{a}|{s}": row.to_dict() for (a, s), row in s5.iterrows()}
    pv = {"S1": res["S1"]["per_cohort"][disc]["wilcoxon_p"], "S2": res["S2"]["per_cohort"][disc]["wilcoxon_p"], "S3": res["S3"]["binomial_p_discovery_vs_5pct"],
          "S4": float(s4d["p"])}
    res["p_unadjusted_discovery"] = pv
    res["p_holm_discovery"] = holm(pv)
    C.write_json(C.RES / "summary.json", res)
    print(json.dumps({k: res[k] for k in ("S1", "S2")}, indent=1, default=str)[:3000])
    print("S3", res["S3"]["verdict"], json.dumps(res["S3"]["per_variant"], default=str))
    print("S4", pb.to_string())


if __name__ == "__main__":
    main()

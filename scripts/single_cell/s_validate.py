"""Validation checks for the single-cell extension S1-S5. Writes results/single_cell/validation_checks.json; exits non-zero on failure.

Independent recomputation where possible: V1 re-reads the raw discovery matrix for two patients and
rebuilds their compartment sums without the streaming builder; V4 recomputes the discovery S4
correlation with separate code; the rest check the registered decision rules against the stored values.
"""
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats

import sccommon as C

checks = []


def check(name, ok, detail):
    checks.append({"check": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name, "|", detail)


summ = json.loads((C.RES / "summary.json").read_text())
gate = json.loads((C.RES / "qc_gate.json").read_text())
enr = pd.read_csv(C.RES / "s1_s2_s5_enrichment.csv", dtype={"patient": str})
coh = pd.read_csv(C.RES / "s3_coherence.csv", dtype={"patient": str})
pb = pd.read_csv(C.RES / "s4_pseudobulk.csv")

# V1: raw re-read of two discovery patients
acc = "GSE132465"
cells = C.cell_table(acc)
retained = gate[acc]["retained"]
pats = retained[:2]
want = cells[(cells["patient"].isin(pats)) & (cells["cls"] == "tumour")]
f = next((C.EXT / acc).glob("*raw_UMI*"))
raw = pd.read_csv(f, sep="\t", index_col=0, usecols=["Index"] + list(want["cell"]))
raw = raw.groupby(level=0).sum()
groups = pd.read_csv(C.PROC / acc / "groups.csv", dtype={"patient": str})
genes = pd.read_csv(C.PROC / acc / "genes.csv")["gene"].astype(str).to_numpy()
sums = pd.DataFrame(np.load(C.PROC / acc / "group_sums.npy").T, index=genes).groupby(level=0).sum()
maxdiff = 0.0
for p in pats:
    for comp in C.COMPARTMENTS:
        ids = want.loc[(want["patient"] == p) & (want["compartment"] == comp), "cell"]
        r = raw[list(ids)].sum(axis=1) if len(ids) else pd.Series(0.0, index=raw.index)
        gi = groups.index[(groups["cls"] == "tumour") & (groups["patient"] == p) & (groups["compartment"] == comp)]
        b = sums.iloc[:, gi].sum(axis=1) if len(gi) else pd.Series(0.0, index=sums.index)
        maxdiff = max(maxdiff, float((r.reindex(sums.index).fillna(0) - b).abs().max()))
check("V1 raw re-read of 2 discovery patients reproduces streamed compartment sums", maxdiff == 0, f"max |diff| = {maxdiff}")

# V2: observed EMT shares recomputed directly for the same patients
sets, _ = C.gene_sets()
det = pd.Series(np.load(C.PROC / acc / "tumour_detected_cells.npy"), index=genes).groupby(level=0).sum()
n_t = json.loads((C.PROC / acc / "build_summary.json").read_text())["n_tumour_cells"]
pool = set(det.index[det >= 0.01 * n_t])
tg = [g for g in sets["EMT"] if g in pool]
d2 = 0.0
for p in pats:
    tot = {comp: float(raw.loc[raw.index.isin(tg), list(want.loc[(want["patient"] == p) & (want["compartment"] == comp), "cell"])].to_numpy().sum())
           for comp in C.COMPARTMENTS}
    s = sum(tot.values())
    for comp in C.COMPARTMENTS:
        st = enr[(enr["accession"] == acc) & (enr["set"] == "EMT") & (enr["patient"] == p) & (enr["compartment"] == comp)]["share"].iloc[0]
        d2 = max(d2, abs(tot[comp] / s - st))
check("V2 observed EMT compartment shares recomputed from raw counts", d2 < 1e-9, f"max |diff| = {d2:.2e}")

# V3: cell counts per compartment equal the annotation table
cnt = want.groupby(["patient", "compartment"]).size()
ok3 = all(int(cnt.get((p, c), 0)) == int(gate[acc]["cells_per_patient"][p][c]) for p in pats for c in C.COMPARTMENTS)
check("V3 per-patient compartment cell counts match the annotation file", ok3, f"patients {pats}")

# V4: discovery S4 rho recomputed with separate code
det4 = json.loads((C.RES / "s4_pseudobulk_detail.json").read_text())[acc]
sc = pd.Series(det4["score_EMT"])
fr = pd.Series(det4["stromal_fraction"]).reindex(sc.index)
rx, ry = sc.rank(), fr.rank()
rho = float(np.corrcoef(rx, ry)[0, 1])
stored = float(pb[(pb["accession"] == acc) & (pb["set"] == "EMT")]["rho"].iloc[0])
check("V4 discovery S4 Spearman rho recomputed (rank-Pearson)", abs(rho - stored) < 1e-12, f"{rho:.6f} vs {stored:.6f}")

# V5: QC gate markers all pass in every analysed cohort
for a in summ["cohorts_passing_qc"]:
    m = gate[a]["markers"]
    check(f"V5 QC positive-control markers in own compartment ({a})", all(v["pass"] for v in m.values()), f"{sum(v['pass'] for v in m.values())}/{len(m)} markers")

# V6: retained patients satisfy the gate thresholds
for a in summ["cohorts_passing_qc"]:
    cp = gate[a]["cells_per_patient"]
    ok6 = all(cp[p]["epithelial"] >= C.MIN_EPI and cp[p]["stromal"] >= C.MIN_STROMAL for p in gate[a]["retained"])
    check(f"V6 retained patients meet >= {C.MIN_EPI} epithelial and >= {C.MIN_STROMAL} stromal cells ({a})", ok6, f"{len(gate[a]['retained'])} patients")

# V7: S3 analysed only retained patients
ok7 = all(set(coh.loc[coh["accession"] == a, "patient"]) <= set(gate[a]["retained"]) for a in summ["cohorts_passing_qc"])
check("V7 S3 restricted to QC-retained patients", ok7, "all cohorts")

# V8-V11: decision rules re-applied to stored values
disc, reps = "GSE132465", [a for a in summ["cohorts_passing_qc"] if a != "GSE132465"]
for q in ("S1", "S2"):
    pc = summ[q]["per_cohort"]
    rule = pc[disc]["wilcoxon_p"] < 0.05 and pc[disc]["median_difference"] > 0 and any(pc[a]["median_difference"] > 0 for a in reps)
    check(f"V8 {q} verdict follows the registered rule", rule == summ[q]["supported"], f"supported={summ[q]['supported']}")
pv = summ["S3"]["per_variant"]
rule3 = all(pv[v][disc]["fraction_exceeding_p95"] > 0.5 and any(pv[v][a]["fraction_exceeding_p95"] > 0.5 for a in reps) for v in pv)
check("V9 S3 verdict follows the registered rule", (summ["S3"]["verdict"].startswith("coherent")) == rule3, summ["S3"]["verdict"])
r4 = pb[(pb["accession"] == disc) & (pb["set"] == "EMT")].iloc[0]
check("V10 S4 verdict follows the registered rule", (r4["rho"] > 0 and r4["ci_lo"] > 0) == summ["S4"]["supported"], f"rho {r4['rho']:.3f}, CI lo {r4['ci_lo']:.3f}")
ph = summ["p_holm_discovery"]
okh = all(ph[k] >= summ["p_unadjusted_discovery"][k] - 1e-15 for k in ph)
check("V11 Holm-adjusted p-values are >= unadjusted", okh, json.dumps({k: f"{v:.2e}" for k, v in ph.items()}))

n_fail = sum(not c["pass"] for c in checks)
C.write_json(C.RES / "validation_checks.json", {"n": len(checks), "n_fail": n_fail, "checks": checks})
print(f"{len(checks) - n_fail}/{len(checks)} checks pass")
sys.exit(1 if n_fail else 0)

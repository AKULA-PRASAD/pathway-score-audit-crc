"""E3 deviation: positive-control co-expression QC (added after GSE103479 showed implausible results; declared post hoc).

Median Pearson r over canonical, tissue-independent co-expressed gene pairs in each cohort (and in GSE39582 for reference).
Verdicts in e3 are reported with and without any cohort failing QC (median r < 0.5).
"""
import numpy as np
import pandas as pd

import rcommon as C
from rcommon import X
from e3_replication import COHORTS, gene_matrix, labels, platform_map, read_matrix

PAIRS = [("COL1A1", "COL1A2"), ("COL1A1", "COL3A1"), ("COL1A2", "COL3A1"), ("MKI67", "TOP2A"), ("CD2", "CD3E"),
         ("CDX2", "VIL1"), ("SPARC", "COL5A1"), ("BUB1", "CCNB1")]


def qc(G):
    rs = {f"{a}-{b}": float(np.corrcoef(G.loc[a], G.loc[b])[0, 1]) for a, b in PAIRS if a in G.index and b in G.index}
    return rs, float(np.median(list(rs.values())))


def main():
    rows = []
    ctx = X.context()
    rs, med = qc(ctx.X)
    rows.append(dict(cohort="GSE39582 (reference)", median_r=med, **rs))
    for acc, cfg in COHORTS.items():
        ann, expr = read_matrix(acc)
        G = gene_matrix(expr, platform_map(cfg["gpl"]), cfg["linear"])
        ids, _ = labels(acc, ann)
        G = G[[i for i in ids if i in G.columns]]
        rs, med = qc(G)
        rows.append(dict(cohort=acc, median_r=med, **rs))
    t = pd.DataFrame(rows)
    t["passes_qc_median_r_ge_0.5"] = t.median_r >= 0.5
    t.to_csv(C.RES / "e3_positive_control_qc.csv", index=False)
    print(t.round(2).to_string())


if __name__ == "__main__":
    main()


def verdicts_excluding_failed():
    import json
    t = pd.read_csv(C.RES / "e3_positive_control_qc.csv")
    failed = set(t.loc[~t["passes_qc_median_r_ge_0.5"], "cohort"])
    res = pd.read_csv(C.RES / "e3_replication_tests.csv")
    res = res[~res.cohort.isin(failed)]
    out = {}
    for hid in ["R1", "R2", "R3", "R4", "R5", "R6"]:
        h = res[(res.hypothesis == hid) & (res.eligible == True)]  # noqa: E712
        n_el = len(h)
        n_rep = int(h.replicated_in_cohort.astype(bool).sum()) if n_el else 0
        n_opp = int(h.opposite_in_cohort.astype(bool).sum()) if n_el else 0
        v = "NOT TESTABLE" if n_el == 0 else "REPLICATED" if (n_rep >= 2 and n_opp == 0) else "PARTIALLY REPLICATED" if (n_rep == 1 and n_opp == 0) else "NOT REPLICATED"
        out[hid] = {"eligible_cohorts": n_el, "consistent": n_rep, "opposite": n_opp, "verdict": v}
    s = json.loads((C.REP / "e3_summary.json").read_text())
    s["qc_failed_cohorts"] = sorted(failed)
    s["verdicts_excluding_qc_failed"] = out
    X.write_json(C.REP / "e3_summary.json", s)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    verdicts_excluding_failed()

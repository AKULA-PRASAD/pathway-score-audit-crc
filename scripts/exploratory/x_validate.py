"""Validation checks for the registered exploratory analyses. Writes reports/exploratory/validation_checks.json and exits non-zero on failure.

Each check ties an exploratory output to an input file or to a frozen primary result.
"""
import json
import sys

import numpy as np
import pandas as pd

import xcommon as X

checks = []


def check(name, expected, observed, ok):
    checks.append({"check": name, "expected": expected, "observed": observed, "result": "PASS" if ok else "FAIL"})


def main():
    O = X.OUT
    ctx = X.context()
    reg = json.loads((O / "plan_registration.json").read_text())
    check("registered plan unchanged since registration", reg["sha256"], X.sha256(X.ROOT / reg["plan"]), reg["sha256"] == X.sha256(X.ROOT / reg["plan"]))
    check("cohort size and events", "452 / 136", f"{len(ctx.ids)} / {int(ctx.e.sum())}", len(ctx.ids) == 452 and int(ctx.e.sum()) == 136)
    # scores reproduce Phase 8C exactly on S3
    from src.pathway_scoring import score_panel
    mem = pd.read_csv(X.ANALYSIS / "cohort_membership.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    Xs3 = pd.read_parquet(X.ANALYSIS / "expression_median_S0.parquet")[s3].astype(float)
    raw, _ = score_panel(Xs3, ctx.panel, train_ids=s3)
    ref = pd.read_csv(X.ROOT / "reports/phase8c/scores_4set_primary_reproduced_full_cohort.csv", index_col=0)
    d = float((raw[ref.columns] - ref.loc[raw.index]).abs().max().max())
    check("frozen pathway definitions reproduce Phase 8C scores", "max |diff| < 1e-12", d, d < 1e-12)
    # X1
    x1 = pd.read_csv(O / "x1_molecular_anchoring/anchoring_tests.csv")
    check("X1 number of primary tests", 24, len(x1), len(x1) == 24)
    s1 = json.loads((O / "x1_molecular_anchoring/x1_summary.json").read_text())
    check("X1 TP53 group sizes", "156 M / 130 WT", s1["missingness"]["tp53.mutation"]["levels"], s1["missingness"]["tp53.mutation"]["levels"] == {"M": 156, "WT": 130})
    # X2
    s2 = json.loads((O / "x2_competing_risks/x2_summary.json").read_text())
    g = s2["event_coding_gate"]
    check("X2 event coding reproduces Phase 6A", "136 recurrences, 46 deaths w/o recurrence", f"{g['recurrences']}, {g['deaths_without_recurrence']}",
          g["recurrences"] == 136 and g["deaths_without_recurrence"] == 46)
    cs = pd.read_csv(O / "x2_competing_risks/cause_specific_cox.csv")
    st = cs[(cs.model == "CS recurrence: clinical") & (cs.term == "tnm.stage=3.0")].iloc[0]
    fz = json.loads((X.ROOT / "reports/phase9/execution/phase9_execution_results.json").read_text())
    m1 = fz["apparent"]["M1"]
    check("X2 unpenalised clinical model reproduces frozen M1 stage HR", round(m1["hazard_ratios"]["tnm.stage=3.0"], 4), round(st.HR, 4),
          abs(st.HR - m1["hazard_ratios"]["tnm.stage=3.0"]) < 1e-3)
    # X3
    s3j = json.loads((O / "x3_batch/x3_summary.json").read_text())
    check("X3 batch count and sizes", "9 batches, 14–103", f"{s3j['n_batches']}, {min(s3j['batch_sizes'].values())}–{max(s3j['batch_sizes'].values())}",
          s3j["n_batches"] == 9 and min(s3j["batch_sizes"].values()) == 14 and max(s3j["batch_sizes"].values()) == 103)
    check("X3 events per batch sum to 136", 136, sum(s3j["events_per_batch"].values()), sum(s3j["events_per_batch"].values()) == 136)
    # X4
    s4 = json.loads((O / "x4_mmr_provenance/x4_summary.json").read_text())
    ct = s4["provenance"]["crosstab_combat_mod_by_mmr_sdrf_all_deposit"]
    check("X4 combat_mod vs MMR in deposit: no discordant known-MMR samples", 0, ct.get("modcan1|pMMR", 0) + ct.get("modcan2|dMMR", 0),
          ct.get("modcan1|pMMR", 0) + ct.get("modcan2|dMMR", 0) == 0)
    check("X4 MMR agreement GEO vs SDRF vs Table S1", "1.0 / 1.0", f"{s4['provenance']['geo_vs_sdrf_mmr_agreement']} / {s4['provenance']['geo_vs_tableS1_mmr_agreement']}",
          s4["provenance"]["geo_vs_sdrf_mmr_agreement"] == 1.0 and s4["provenance"]["geo_vs_tableS1_mmr_agreement"] == 1.0)
    # X5
    s5 = json.loads((O / "x5_composition/x5_summary.json").read_text())
    check("X5 ESTIMATE re-implementation validation gate", "< 1e-6", s5["estimate_validation_gate"]["max_rel_err_scores"], s5["estimate_validation_gate"]["passed"])
    comp = pd.read_csv(O / "x5_composition/composition_estimates.csv", index_col=0)
    check("X5 composition estimates for all 452", 452, len(comp), len(comp) == 452 and set(comp.index) == set(ctx.ids))
    check("X5 all MCP-counter marker probe sets present", "all", s5["mcp_probe_mode_coverage"],
          all(v["n_present"] == v["n_markers"] for v in s5["mcp_probe_mode_coverage"].values()))
    # X6
    s6 = json.loads((O / "x6_scores_vs_composition/x6_summary.json").read_text())
    check("X6 EMT–stromal overlap", "29 of 198 genes", f"{s6['overlap']['EMT']['stromal_union']['n']} of {s6['overlap']['EMT']['n_members']}",
          s6["overlap"]["EMT"]["stromal_union"]["n"] == 29 and s6["overlap"]["EMT"]["n_members"] == 198)
    # X7
    s7 = json.loads((O / "x7_composition_adjusted_cox/x7_summary.json").read_text())
    c0 = s7["discrimination"]["C0"]["apparent_C"]
    check("X7 clinical model apparent C reproduces frozen M1 apparent C", round(fz["apparent"]["M1"]["C"], 4), round(c0, 4), abs(c0 - fz["apparent"]["M1"]["C"]) < 1e-4)
    check("X7 bootstrap: no failed refits", 0, sum(v["n_fail"] for v in s7["discrimination"].values()), sum(v["n_fail"] for v in s7["discrimination"].values()) == 0)
    # X8
    t8 = pd.read_csv(O / "x8_external_audit/eligibility_table.csv")
    r14 = t8[t8.cohort == "GSE14333"].iloc[0]
    check("X8 GSE14333 stage II–III events reproduce Phase 5 (48)", 48, int(r14.events_II_III), int(r14.events_II_III) == 48)
    check("X8 GSE33113 events reproduce Phase 5 (18)", 18, int(t8[t8.cohort == "GSE33113"].events_II_III.iloc[0]), int(t8[t8.cohort == "GSE33113"].events_II_III.iloc[0]) == 18)
    check("X8 no cohort eligible for frozen-model validation", 0, int(t8.frozen_model_validation_eligible.sum()), int(t8.frozen_model_validation_eligible.sum()) == 0)
    # frozen primary results untouched
    iv = json.loads((X.ROOT / "reports/phase9/corrected_rerun_report/corrected_interval_results.json").read_text())
    check("frozen θ̂ unchanged", 0.004316, round(iv["theta_hat"], 6), round(iv["theta_hat"], 6) == 0.004316)
    # figures exist
    figs = sorted(p.name for p in (O / "figures").glob("XF*.png"))
    check("six exploratory figures (PNG + SVG)", 6, len(figs), len(figs) == 6 and len(list((O / "figures").glob("XF*.svg"))) == 6)
    X.write_json(O / "validation_checks.json", {"n": len(checks), "n_fail": sum(c["result"] == "FAIL" for c in checks), "checks": checks})
    for c in checks:
        print(c["result"], c["check"], "|", c["observed"])
    if any(c["result"] == "FAIL" for c in checks):
        sys.exit(1)


if __name__ == "__main__":
    main()

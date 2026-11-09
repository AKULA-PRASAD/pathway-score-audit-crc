"""Portfolio values, read directly from frozen project outputs (read-only).

Every number used in the portfolio diagrams and slides comes from this module.
Nothing here fits a model or recomputes a statistic: values are read from frozen
JSON/CSV outputs, except the set-coherence values, which are quoted from
docs/PHASE8C_FROZEN_PATHWAY_SCORE_ANALYSIS.md §3.1 and checked against that text
by portfolio/src/qc_check.py.
"""
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

SRC = {
    "figure_values": "reports/figures/final/figure_values.json",
    "execution": "reports/phase9/execution/phase9_execution_results.json",
    "s7": "reports/phase9/execution_s7/phase9_s7_results.json",
    "intervals": "reports/phase9/corrected_rerun_report/corrected_interval_results.json",
    "panel4": "config/pathways/m2_panel_4set.json",
    "panel5": "config/pathways/m2_panel_5set.json",
    "annotation": "reports/phase8c/score_vs_annotation.csv",
    "representation": "reports/phase8c/representation_sensitivity.csv",
    "phase8c_doc": "docs/PHASE8C_FROZEN_PATHWAY_SCORE_ANALYSIS.md",
}

PATHWAYS = [
    ("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "EMT"),
    ("HALLMARK_E2F_TARGETS", "E2F targets"),
    ("HALLMARK_WNT_BETA_CATENIN_SIGNALING", "Wnt/β-catenin"),
    ("HALLMARK_P53_PATHWAY", "p53 pathway"),
]

# Quoted from PHASE8C §3.1 (implied mean inter-gene correlation); verified by qc_check.py.
COHERENCE = {"EMT": 0.241, "E2F targets": 0.304, "Wnt/β-catenin": 0.041, "p53 pathway": 0.033}


def _json(key):
    return json.loads((ROOT / SRC[key]).read_text())


def _csv(key):
    with open(ROOT / SRC[key], newline="") as f:
        return list(csv.DictReader(f))


def load():
    fv, ex, s7, iv = _json("figure_values"), _json("execution"), _json("s7"), _json("intervals")
    p4, p5 = _json("panel4"), _json("panel5")
    ann, rep = _csv("annotation"), _csv("representation")
    V = {}
    f1 = fv["fig1"]
    V["cohort"] = {
        "arrays": f1["arrays_total"], "S3_n": f1["S3"]["n"], "S3_events": f1["S3"]["events"],
        "n": f1["A19"]["n"], "events": f1["A19"]["events"], "censored": f1["A19"]["censored"],
        "discovery": f1["A19_by_partition"]["discovery"], "cit": f1["A19_by_partition"]["validation"],
        "deaths_censored": f1["deaths_censored"],
    }
    pr = fv["fig2"]["probes"]
    V["probes"] = {"total": pr["n_probes_total"], "single_gene": pr["n_unique"], "unmapped": pr["n_unmapped"],
                   "ambiguous": pr["n_ambiguous"], "genes": pr["n_genes_unique"],
                   "multi_probe_genes": pr["n_genes_multi_probe"]}
    V["panel"] = {
        "union4": p4["union_unique_genes"], "union5": p5["union_unique_genes"],
        "cond4": p4["collinearity_all461_descriptive_only"]["condition_number"],
        "cond5": p5["collinearity_all461_descriptive_only"]["condition_number"],
        "maxr4": p4["collinearity_all461_descriptive_only"]["max_abs_r"],
        "maxr5": p5["collinearity_all461_descriptive_only"]["max_abs_r"],
        "sizes": {lab: p4["gene_sets"][pid]["n_resolved_unique"] for pid, lab in PATHWAYS},
        "version": p4["database"]["version"], "frozen_utc": p4["frozen_utc"],
    }
    a1, a2, a5 = ex["apparent"]["M1"], ex["apparent"]["M2_primary"], ex["apparent"]["M2_sensitivity"]
    s2 = s7["apparent"]["M2_primary"]
    strat = ex["batch_sensitivity"]["stratified_fit"]["coefficients"]
    lobo = ex["batch_sensitivity"]["leave_one_batch_out"]
    annm = {r["pathway"]: r for r in ann if r["annotation"] == "mmr.status"}
    repm = {r["pathway"]: r for r in rep}
    V["pathways"] = []
    for pid, lab in PATHWAYS:
        lv = [x["coefficients"][pid] for x in lobo]
        all_c = a2["coefficients"][pid]
        V["pathways"].append({
            "id": pid, "label": lab, "genes": V["panel"]["sizes"][lab], "coherence": COHERENCE[lab],
            "hr": (a2["hazard_ratios"][pid], a2["hr_ci_lower"][pid], a2["hr_ci_upper"][pid]),
            "hr_s7": (s2["hazard_ratios"][pid], s2["hr_ci_lower"][pid], s2["hr_ci_upper"][pid]),
            "hr_strat": math.exp(strat[pid]),
            "coef_all": all_c, "lobo_min": min(lv), "lobo_max": max(lv),
            "lobo_sign_changes": sum(1 for x in lv if (x > 0) != (all_c > 0)),
            "mmr_e2": float(annm[pid]["location_epsilon_squared"]), "mmr_q": float(annm[pid]["location_q_BH"]),
            "mmr_n": int(annm[pid]["n"]),
            "rep_r": float(repm[pid]["pearson_median_vs_repiqr"]),
        })
    e2f = "HALLMARK_E2F_TARGETS"
    V["five_set_e2f"] = (a5["hazard_ratios"][e2f], a5["hr_ci_lower"][e2f], a5["hr_ci_upper"][e2f])
    V["lambda"] = {"primary": a2["lambda"], "s7": s2["lambda"], "five": a5["lambda"]}
    V["m1_terms"] = {k: (a1["hazard_ratios"][k], a1["hr_ci_lower"][k], a1["hr_ci_upper"][k]) for k in a1["hazard_ratios"]}
    f7, f9 = fv["fig7"], fv["fig9"]
    V["C"] = {m: {"apparent": f7[m]["apparent_C"], "corrected": f7[m]["corrected_C"]} for m in ("M1", "M2_primary", "M2_sensitivity")}
    V["slopes"] = {m: f9["calibration"][m]["corrected"] for m in ("M1", "M2_primary", "M2_sensitivity")}
    V["theta"] = iv["theta_hat"]
    V["runs"] = [(iv["runs"][str(r)]["percentile_primary"]["lower"], iv["runs"][str(r)]["percentile_primary"]["upper"]) for r in range(5)]
    V["B"] = {"out": fv["fig8"]["B_out"], "in": fv["fig8"]["B_in"]}
    V["verification"] = fv["fig8"]["verification"]
    V["lobo_C"] = (min(x["C_apparent"] for x in lobo), max(x["C_apparent"] for x in lobo))
    V["lobo_sizes"] = (min(x["n_omitted"] for x in lobo), max(x["n_omitted"] for x in lobo))
    V["lobo_batches"] = [{"batch": x["omitted_batch"], "n": x["n_omitted"],
                          "coef": {lab: x["coefficients"][pid] for pid, lab in PATHWAYS}} for x in lobo]
    V["mc"] = {"rho_max": fv["figS1"]["rho_max"], "threshold": fv["figS1"]["threshold"]}
    V["stage_hr"] = V["m1_terms"]["tnm.stage=3.0"]
    return V


if __name__ == "__main__":
    import pprint
    pprint.pprint(load(), width=120)


def load_exploratory():
    """Post-hoc exploratory results (registered plan docs/EXPLORATORY_ANALYSIS_PLAN.md), read from reports/exploratory/.

    Kept separate from load(): these are NOT frozen primary outputs and are always labelled post-hoc in the portfolio.
    """
    E = ROOT / "reports" / "exploratory"
    x1 = {(r["annotation"], r["pathway"]): r for r in csv.DictReader(open(E / "x1_molecular_anchoring" / "anchoring_tests.csv"))}
    s1 = json.loads((E / "x1_molecular_anchoring" / "x1_summary.json").read_text())
    s2 = json.loads((E / "x2_competing_risks" / "x2_summary.json").read_text())
    cif = {(r["group"], r["cause"], int(float(r["horizon"]))): r for r in csv.DictReader(open(E / "x2_competing_risks" / "cumulative_incidence.csv"))}
    s3 = json.loads((E / "x3_batch" / "x3_summary.json").read_text())
    s4 = json.loads((E / "x4_mmr_provenance" / "x4_summary.json").read_text())
    s6 = json.loads((E / "x6_scores_vs_composition" / "x6_summary.json").read_text())
    cor = {(r["score"], r["composition"]): float(r["rho"]) for r in csv.DictReader(open(E / "x6_scores_vs_composition" / "score_composition_spearman.csv"))}
    s7 = json.loads((E / "x7_composition_adjusted_cox" / "x7_summary.json").read_text())
    t8 = list(csv.DictReader(open(E / "x8_external_audit" / "eligibility_table.csv")))
    s8 = json.loads((E / "x8_external_audit" / "x8_summary.json").read_text())
    g = s3["frailty"]["gamma"]
    mmr_e2f = next(r for r in s7["mmr_adjusted_sensitivity"]["rows"] if r["model"] == "R0 + dMMR" and r["term"] == "E2F")
    return {
        "tp53_p53_diff": float(x1[("tp53.mutation", "p53")]["median_diff_sd"]), "tp53_p53_q": float(x1[("tp53.mutation", "p53")]["q_BH"]),
        "tp53_pmmr_diff": s1["secondary_tp53_within_pMMR"]["mean_diff_sd_mutant_minus_wt"], "tp53_pmmr_p": s1["secondary_tp53_within_pMMR"]["batch_stratified_permutation_p"],
        "braf_wnt_diff": float(x1[("braf.mutation", "Wnt")]["median_diff_sd"]), "cimp_wnt_diff": float(x1[("cimp.status", "Wnt")]["median_diff_sd"]),
        "cif_rec60": float(cif[("all", "recurrence", 60)]["cif"]), "cif_death60": float(cif[("all", "death without recurrence", 60)]["cif"]),
        "naive60": s2["naive_one_minus_KM_recurrence"]["60"],
        "theta": g["theta_hat"], "theta_ci": g["theta_95_profile_interval_on_grid"],
        "frailty_hr": {k: g["coefficients"][k]["HR"] for k in ("EMT", "E2F", "Wnt", "p53")},
        "lobo_within": sum(abs(float(r["delta_in_se"])) < 1.96 for r in csv.DictReader(open(E / "x3_batch" / "lobo_diagnostics.csv"))),
        "lobo_total": sum(1 for _ in csv.DictReader(open(E / "x3_batch" / "lobo_diagnostics.csv"))),
        "p53_within_c2": next(r for r in s4["robustness"] if r["pathway"] == "p53")["within_cit"]["C2"]["eps2"],
        "emt_fib": cor[("EMT", "MCP_Fibroblasts")], "emt_stromal": cor[("EMT", "StromalScore")], "emtex_fib": cor[("EMT_ex", "MCP_Fibroblasts")],
        "emt_overlap": s6["overlap"]["EMT"]["stromal_union"]["n"], "vif_emt": s7["vif"]["R1"]["EMT"], "vif_str": s7["vif"]["R1"]["StromalScore"],
        "c0": s7["discrimination"]["C0"]["corrected_C"], "c1": s7["discrimination"]["C1"]["corrected_C"],
        "e2f_mmr": (mmr_e2f["HR"], mmr_e2f["lo"], mmr_e2f["hi"]),
        "ext_hits": s8["search"]["hits"], "ext_assessed": len(t8), "ext_eligible": s8["n_eligible_for_frozen_model_validation"],
        "ext_max_events": max(int(float(r["events_II_III"])) for r in t8 if r["events_II_III"] and r["stage_II_III_n"]),
        "ext_needed": s8["events_needed_for_relaxed_slope_target"],
    }


def load_extensions():
    """Registered research extensions E1–E3 + TCGA addendum and the pathway-audit case study (post-hoc, outcome-free)."""
    R = ROOT / "results" / "research_extensions"
    P = ROOT / "reports" / "research_extensions"
    ch = json.loads((P / "challenge_tests.json").read_text())
    tc = json.loads((P / "e3_tcga_summary.json").read_text())
    qc = {r["cohort"]: float(r["median_r"]) for r in csv.DictReader(open(R / "e3_positive_control_qc.csv"))}
    e3 = {(r["hypothesis"], r["cohort"]): r for r in csv.DictReader(open(R / "e3_replication_tests.csv"))}
    t3 = {r["hypothesis"]: r for r in csv.DictReader(open(R / "e3_tcga_replication_tests.csv"))}
    gt = {(r["annotation"], r["score"]): r for r in csv.DictReader(open(R / "e2_genotype_adjusted.csv"))}
    cs = list(csv.DictReader(open(ROOT / "reports" / "pathway_audit_case_study" / "summary.csv")))
    surv = next(float(r["survival_screen_corrected_dC"]) for r in cs if r["cohort"] == "GSE39582")
    dis = {"GSE39582": ch["C2_R6_gene_disjoint"]["GSE39582"], "GSE41258": ch["C2_R6_gene_disjoint"]["GSE41258"],
           "GSE143985": ch["C2_R6_gene_disjoint"]["GSE143985"], "TCGA-COAD": tc["R6_gene_disjoint"]}
    return {"disjoint": dis, "qc": {**qc, "TCGA-COAD": tc["qc"]["median_r"]}, "verdicts": tc["verdicts_qc_passing_plus_tcga"],
            "r1_41258": float(e3[("R1", "GSE41258")]["effect"]), "r1_tcga": float(t3["R1"]["effect"]), "r1_tcga_ci": (float(t3["R1"]["lo"]), float(t3["R1"]["hi"])),
            "r3_tcga": float(t3["R3"]["effect"]), "r4_tcga": float(t3["R4"]["effect"]), "n_tcga": tc["n_primary_tumours"],
            "braf_p53": (float(gt[("braf.mutation", "p53")]["diff_unadj"]), float(gt[("braf.mutation", "p53")]["diff_adj"])),
            "braf_wnt": (float(gt[("braf.mutation", "Wnt")]["diff_unadj"]), float(gt[("braf.mutation", "Wnt")]["diff_adj"])),
            "survival_screen": surv, "n_shared_fib": ch["C2_R6_gene_disjoint"]["n_shared_with_EMT"],
            "baseline": json.loads((P / "tcga_challenge_tests.json").read_text())["C6_random_baseline_R6"],
            "tcga_hyper": json.loads((P / "tcga_challenge_tests.json").read_text())["C5_tcga_hypermutation"]}

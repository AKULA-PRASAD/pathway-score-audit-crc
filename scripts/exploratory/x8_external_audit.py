"""X8. External-cohort eligibility audit. Registered plan §X8.

Inputs: the recorded GEO search (data/external/exploratory/geo_search/query.txt, esearch.json, esummary.json) and
per-sample characteristics text downloaded for screened candidates (chars/*.txt). No expression data are downloaded.
Expected precision replaces a fixed event threshold: widths are scaled from the Phase 5 simulation at 136 events
(C-index 95 % CI width 0.0945; calibration-slope 95 % CI width 0.63 at C = 0.65) by sqrt(136 / events), the
event-count scaling that Phase 5 showed governs precision. Riley et al. (Stat Med 2022) relaxed target: slope CI width 0.4.
"""
import json
import re
from collections import Counter

import numpy as np
import pandas as pd

import xcommon as X

GS = X.EXT / "geo_search"


def parse(g):
    txt = open(GS / "chars" / f"{g}.txt", errors="replace").read()
    out = []
    for s in txt.split("^SAMPLE")[1:]:
        d = {"gsm": s.split("\n")[0].strip(" =")}
        for m in re.finditer(r"!Sample_characteristics_ch1 = ([^\n]+)", s):
            kv = m.group(1)
            if ":" in kv:
                k, v = kv.split(":", 1)
                d[k.strip().lower()] = v.strip()
        m = re.search(r"!Sample_platform_id = (\S+)", s)
        d["_gpl"] = m.group(1) if m else ""
        m = re.search(r"!Sample_source_name_ch1 = ([^\n]+)", s)
        d["_source"] = m.group(1) if m else ""
        m = re.search(r"!Sample_data_processing = ([^\n]+)", s)
        d["_proc"] = m.group(1)[:100] if m else ""
        out.append(d)
    return pd.DataFrame(out)


def counts_14333():
    d = parse("GSE14333")
    f = d["location"].str.split(";", expand=False)
    rec = []
    for s in d["location"]:
        parts = [p.strip() for p in s.split(";")]
        kv = {"site": parts[0]}
        for p in parts[1:]:
            k, v = p.split(":")
            kv[k.strip()] = v.strip()
        rec.append(kv)
    r = pd.DataFrame(rec)
    r["stage23"] = r.DukesStage.isin(["B", "C"])
    r["event"] = (r.DFS_Cens == "0") & r.DFS_Time.notna()   # Phase 4: 0 = event
    colon = r.site.isin(["Left", "Right", "Colon"])
    return {"stage_II_III_n": int(r.stage23.sum()), "events_II_III": int((r.stage23 & r.event).sum()),
            "colon_only_II_III_n": int((r.stage23 & colon).sum()), "colon_only_events": int((r.stage23 & colon & r.event).sum()),
            "sites": dict(Counter(r.site))}


def simple(g, stage_col, stage_ok, ev_col, ev_yes, tissue_filter=None):
    d = parse(g)
    if tissue_filter:
        d = d[d.apply(tissue_filter, axis=1)]
    st = d[stage_col].astype(str).isin(stage_ok)
    ev = d[ev_col].astype(str).isin(ev_yes)
    known = ~d[ev_col].astype(str).isin(["NA", "na", "nan", ""])
    return {"stage_II_III_n": int(st.sum()), "stage_II_III_with_outcome": int((st & known).sum()), "events_II_III": int((st & ev).sum())}


def main():
    out = X.outdir("x8_external_audit")
    search = json.loads((GS / "esearch.json").read_text())
    summ = json.loads((GS / "esummary.json").read_text())
    C = {}
    C["GSE14333"] = counts_14333()
    C["GSE161158"] = simple("GSE161158", "stage", ["stage 2", "stage 3"], "dfs_event", ["1"])
    C["GSE17536"] = simple("GSE17536", "ajcc_stage", ["2", "3"], "dfs_event (disease free survival; cancer recurrence)", ["recurrence"])
    C["GSE17537"] = simple("GSE17537", "ajcc_stage", ["2", "3"], "dfs_event (disease free survival; cancer recurrence)", ["recurrence"])
    d37 = parse("GSE37892")
    C["GSE37892"] = {"stage_II_III_n": int(d37.stage.isin(["2", "3"]).sum()),
                     "events_II_III": int((d37["date at distant metastasis"] != "NA").sum())}
    d33 = parse("GSE33113")
    C["GSE33113"] = {"stage_II_III_n": int((d33.tissue == "primary tumor resection").sum()),
                     "events_II_III": int((d33["meta or recurrence within 3 years"] == "yes").sum())}
    C["GSE41258"] = simple("GSE41258", "group stage", ["II", "III"], "recurrence", ["Yes"], lambda r: r.get("tissue") == "Primary Tumor")
    C["GSE103479"] = simple("GSE103479", "dukes", ["B", "C"], "recurrence", ["Yes"])
    C["GSE106584"] = simple("GSE106584", "clinical stage", ["2", "2A", "2B", "3", "3A", "3B", "3C"], "recurrence", ["1"])
    C["GSE143985"] = simple("GSE143985", "stage", ["2", "3"], "dfs_event", ["1"])
    d24 = parse("GSE24551")
    d24 = d24[d24.tissue == "colorectal cancer biopsy"]
    # labels are swapped in the deposit: 'disease_free_survival_years' holds 0/1 and 'disease_free_survival_event' holds times
    C["GSE24551"] = {"stage_II_III_n": int(d24["tumor stage"].isin(["2", "3"]).sum()),
                     "events_II_III": int((d24["tumor stage"].isin(["2", "3"]) & (d24["disease_free_survival_years"] == "1")).sum()),
                     "note": "event and time field labels appear swapped in the deposit"}
    C["GSE38832"] = simple("GSE38832", "ajcc_stage", ["2", "3"], "dfs_event (disease free survival)", ["1 (cancer recurrence)"])
    d816 = parse("GSE81653")
    C["GSE81653"] = {"stage_II_III_n": None, "events_II_III": int((d816["recurrence status"] == "Had Recurrence").sum()),
                     "note": "stage and time not deposited; event count is all stages"}
    d58 = parse("GSE58394")
    C["GSE58394"] = {"stage_II_III_n": int(len(d58)), "events_II_III": int((d58.recurrency == "recurrent").sum()), "note": "stage II only; no time"}
    # ---- eligibility table (assessments from the deposited fields; see docs/EXTERNAL_COHORT_ELIGIBILITY_AUDIT.md)
    A = [
        ("GSE14333", "Peter MacCallum (Melbourne) + Moffitt", "GPL570 (MAS5)", "DFS (DFS_Cens 0 = event)", "yes", "age, sex, Dukes stage, site, AdjCTX: all five", "Independent of GSE39582. Overlaps GSE17536 (Moffitt) — Phase 4", "includes rectal tumours; MAS5 not RMA; DFS time unit undocumented"),
        ("GSE161158", "Moffitt + Vanderbilt", "GPL570 (RMA+ComBat)", "DFS", "yes", "age, stage only", "Independent of GSE39582; re-deposit overlapping GSE17536/GSE17537 institutions", "sex, location, chemotherapy missing"),
        ("GSE17536", "Moffitt", "GPL570 (RMA)", "DFS (recurrence)", "yes", "age, sex, stage", "Overlaps GSE14333 (Phase 4) and GSE161158", "location, chemotherapy missing"),
        ("GSE17537", "Vanderbilt", "GPL570 (RMA)", "DFS (recurrence)", "yes", "age, sex, stage", "Overlaps GSE161158 institution", "location, chemotherapy missing"),
        ("GSE37892", "Six French centres", "GPL570", "distant metastasis (dates)", "yes (dates)", "age, sex, stage, location", "French centres; overlap with the CIT programme cannot be excluded from the deposit", "endpoint is distant metastasis only; chemotherapy missing"),
        ("GSE33113", "AMC Amsterdam", "GPL570 (MAS5)", "metastasis or recurrence within 3 years", "yes", "age, sex (stage II only)", "Independent", "composite endpoint truncated at 3 years; stage II only; location, chemotherapy missing"),
        ("GSE41258", "MSKCC", "GPL96 (U133A, PLIER)", "recurrence with interval", "yes", "age, sex, stage, location", "Independent", "U133A lacks part of the 607 genes; chemotherapy missing"),
        ("GSE103479", "Aberdeen, Barcelona, Florence, Dublin", "GPL23985 (Almac Xcel, RMA+ComBat)", "recurrence; progression-free survival time", "yes", "age, sex, stage, tumour site, adjuvant treatment: all five", "Independent", "different array platform: frozen scores cannot be reproduced exactly"),
        ("GSE106584", "US (Mayo-type cohort)", "GPL17586 (HTA 2.0)", "recurrence; DFS days", "yes", "stage, adjuvant therapy", "Independent", "age, sex, location missing; different platform; includes stage I and IV"),
        ("GSE143985", "Japan", "GPL570 (MAS5)", "DFS", "yes", "stage, chemotherapy", "Independent", "age, sex, location missing"),
        ("GSE24551", "Oslo", "GPL5175/GPL11028 (HuEx exon)", "DFS (fields swapped)", "yes", "stage, MSI", "Independent", "age, sex, location, chemotherapy missing; exon platform"),
        ("GSE38832", "Vanderbilt and others", "GPL570 (MAS5)", "DFS", "yes", "stage", "Overlaps Vanderbilt deposits", "few events; covariates missing"),
        ("GSE81653", "Unknown", "GPL16686 (HuGene 2.0)", "recurrence status only", "no", "none", "Unknown", "no time-to-event, no stage, no covariates"),
        ("GSE58394", "Japan", "GPL18709 (Agilent)", "recurrence status", "no", "age, sex (stage II only)", "Independent", "no time; Agilent platform"),
    ]
    INDEP = {g: True for g, *_ in A}
    INDEP["GSE37892"] = None          # French centres: overlap with CIT cannot be excluded from the deposit
    INDEP["GSE81653"] = None          # institution not stated
    ENDPOINT_OK = {"GSE14333": "unverified (DFS; death handling not stated)", "GSE161158": "unverified (DFS; death handling not stated)",
                   "GSE17536": "partly (recurrence-coded DFS)", "GSE17537": "partly (recurrence-coded DFS)",
                   "GSE37892": "no (distant metastasis only)", "GSE33113": "no (composite, truncated at 3 years)",
                   "GSE41258": "partly (recurrence with interval)", "GSE103479": "unverified (PFS time with recurrence flag)",
                   "GSE106584": "partly (recurrence; DFS days)", "GSE143985": "unverified (DFS)", "GSE24551": "unverified (DFS; fields swapped)",
                   "GSE38832": "unverified (DFS)", "GSE81653": "no (no time)", "GSE58394": "no (no time)"}
    rows = []
    for g, inst, plat, endp, time_, preds, indep, issues in A:
        c = C.get(g, {})
        ev = c.get("events_II_III")
        w_c = 0.0945 * np.sqrt(136 / ev) if ev else np.nan
        w_s = 0.63 * np.sqrt(136 / ev) if ev else np.nan
        crit = {"c1_independent_of_GSE39582": INDEP[g],
                "c2_compatible_recurrence_endpoint": ENDPOINT_OK[g],
                "c3_all_five_M1_predictors": "all five" in preds,
                "c4_exact_score_reproduction_from_deposit": plat.startswith("GPL570 (RMA"),
                "c5_relaxed_precision_target_slope_width_le_0_4": bool(ev and w_s <= 0.4)}
        eligible = (crit["c1_independent_of_GSE39582"] is True and crit["c2_compatible_recurrence_endpoint"].startswith("yes")
                    and crit["c3_all_five_M1_predictors"] and crit["c4_exact_score_reproduction_from_deposit"])
        rows.append(dict(cohort=g, institutions=inst, platform=plat, endpoint=endp, time_to_event=time_, m1_predictors=preds,
                         independence=indep, main_issues=issues, stage_II_III_n=c.get("stage_II_III_n"), events_II_III=ev,
                         approx_C_CI_width=round(w_c, 3) if ev else None, approx_slope_CI_width=round(w_s, 2) if ev else None,
                         **crit, frozen_model_validation_eligible=bool(eligible), count_notes=c.get("note", "")))
    tab = pd.DataFrame(rows)
    tab.to_csv(out / "eligibility_table.csv", index=False)
    X.write_json(out / "x8_summary.json", {"plan_section": "X8", "search": {"database": "NCBI GEO DataSets (gds) via E-utilities",
                                          "query": (GS / "query.txt").read_text().strip(), "date": "2026-11-08",
                                          "hits": int(search["esearchresult"]["count"]), "screened_n_ge_80": 85},
                                          "counts": C, "events_needed_for_relaxed_slope_target": round(136 * (0.63 / 0.4) ** 2),
                                          "n_eligible_for_frozen_model_validation": int(tab.frozen_model_validation_eligible.sum())})
    print(tab[["cohort", "events_II_III", "approx_slope_CI_width", "c1_independent_of_GSE39582", "c2_compatible_recurrence_endpoint", "c3_all_five_M1_predictors", "c4_exact_score_reproduction_from_deposit", "frozen_model_validation_eligible"]].to_string())
    print(C["GSE14333"])
    print("events for slope width 0.4 ~", round(136 * (0.63 / 0.4) ** 2))


if __name__ == "__main__":
    main()

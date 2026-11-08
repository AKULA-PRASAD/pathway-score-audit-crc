"""Phase 1: cohort and endpoint INVESTIGATION (read-only).

Regenerates every quantitative claim in
docs/PHASE1_COHORT_ENDPOINT_INVESTIGATION.md.

This script investigates and reports. It makes NO decision: nothing is
excluded, recoded, imputed or censored, and no cohort or endpoint is locked.
Scripts 01-14, raw data, the cached matrix and the clinical metadata are
untouched. Outputs are written only to reports/phase1/.

Supplementary Table S1 (per-patient annotation from the source publication)
is fetched on demand to data/external/ and cross-checked against the GEO
metadata. It is never written into metadata/ and never overrides GEO.

Usage:
    .venv/bin/python scripts/15_investigate_cohort_endpoint.py
    .venv/bin/python scripts/15_investigate_cohort_endpoint.py --no-fetch
"""

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths  # noqa: E402

OUT = paths.ROOT / "reports" / "phase1"

TUMOUR_PATTERN = "primary colorectal Adenocarcinoma"

# Marisa et al. 2013, PLoS Med 10(5):e1001453, Table S1 "Patient and Tumor
# Characteristics" -- per-patient annotation, 566 CIT identifiers.
TABLE_S1_URL = (
    "https://journals.plos.org/plosmedicine/article/file"
    "?id=10.1371/journal.pmed.1001453.s009&type=supplementary"
)
TABLE_S1_PATH = Path("data/external/marisa2013_tableS1.xls")


def fetch_table_s1(root: Path, allow_fetch: bool = True) -> Path | None:
    """Fetch Table S1 if absent. Read-only external evidence; never merged."""
    dest = root / TABLE_S1_PATH
    if dest.exists():
        return dest
    if not allow_fetch:
        return None
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        req = urllib.request.Request(
            TABLE_S1_URL, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            blob = r.read()
        if not blob.startswith(b"\xd0\xcf\x11\xe0"):   # OLE2 magic
            print("  WARNING: Table S1 download is not an .xls file; skipping.")
            return None
        dest.write_bytes(blob)
        print(f"  fetched Table S1 -> {TABLE_S1_PATH} "
              f"(sha256 {hashlib.sha256(blob).hexdigest()[:16]}...)")
        return dest
    except Exception as exc:                     # noqa: BLE001
        print(f"  WARNING: could not fetch Table S1 ({type(exc).__name__}). "
              "Cross-check skipped.")
        return None


def load() -> pd.DataFrame:
    return pd.read_csv(paths.SAMPLE_METADATA)


def tumours(meta: pd.DataFrame) -> pd.DataFrame:
    return meta[meta["source_name"].str.contains(
        TUMOUR_PATTERN, case=False, na=False)].copy()


def section(title: str) -> None:
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--no-fetch", action="store_true",
                    help="Skip the Table S1 download.")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    meta = load()
    tum = tumours(meta)
    s23 = tum[tum["tnm.stage"].isin([2, 3])].copy()
    complete = s23[s23[["rfs.event", "rfs.delay"]].notna().all(axis=1)].copy()

    findings: dict = {}

    # ---------------------------------------------------------------
    section("1. COHORT FLOW (counts only; nothing excluded)")
    # ---------------------------------------------------------------
    flow = [
        ("deposited samples", len(meta)),
        ("primary tumours", len(tum)),
        ("stage II-III primary tumours", len(s23)),
        ("with documented RFS (event and time)", len(complete)),
        ("  of those, rfs.delay > 0", int((complete["rfs.delay"] > 0).sum())),
    ]
    prev = None
    for label, n in flow:
        delta = "" if prev is None else f"  (-{prev - n})"
        print(f"  {label:42s} {n:4d}{delta}")
        prev = n
    findings["cohort_flow"] = dict(flow)

    # ---------------------------------------------------------------
    section("2. PUBLICATION CROSS-CHECK")
    # ---------------------------------------------------------------
    disc = complete[complete["dataset"] == "discovery"]
    print(f"  discovery stage II-III with documented RFS : {len(disc)}")
    print(f"  Marisa et al. 2013 reports                 : 359")
    print(f"  match                                      : {len(disc) == 359}")
    findings["publication_crosscheck"] = {
        "our_discovery_count": len(disc), "published": 359,
        "match": bool(len(disc) == 359),
    }

    # ---------------------------------------------------------------
    section("3. TIME UNIT EVIDENCE")
    # ---------------------------------------------------------------
    both = tum[tum[["rfs.delay", "os.delay_(months)"]].notna().all(axis=1)]
    viol = int((both["rfs.delay"] > both["os.delay_(months)"]).sum())
    eq = int((both["rfs.delay"] == both["os.delay_(months)"]).sum())
    print(f"  records with both times          : {len(both)}")
    print(f"  rfs.delay range                  : {both['rfs.delay'].min():.0f}"
          f" - {both['rfs.delay'].max():.0f}")
    print(f"  os.delay range                   : "
          f"{both['os.delay_(months)'].min():.0f}"
          f" - {both['os.delay_(months)'].max():.0f}")
    print(f"  rfs.delay > os.delay (violations): {viol}")
    print(f"  rfs.delay == os.delay            : {eq} "
          f"({eq / len(both) * 100:.1f}%)")
    print(f"  both integer-valued              : "
          f"{bool((both['rfs.delay'] % 1 == 0).all())}")
    print("  NOTE: supplementary methods say 'five years/60 months',")
    print("        equating the analysis horizon to months.")
    findings["time_unit"] = {
        "n": len(both), "violations_rfs_gt_os": viol, "equal": eq,
        "rfs_max": float(both["rfs.delay"].max()),
        "os_max": float(both["os.delay_(months)"].max()),
    }

    # ---------------------------------------------------------------
    section("4. DEATH WITHOUT RECURRENCE: EMPIRICAL ENCODING TEST")
    # ---------------------------------------------------------------
    full = s23[s23[["rfs.event", "rfs.delay", "os.event",
                    "os.delay_(months)"]].notna().all(axis=1)]
    rows = []
    for (r, o), g in full.groupby([full["rfs.event"].astype(int),
                                   full["os.event"].astype(int)]):
        e = int((g["rfs.delay"] == g["os.delay_(months)"]).sum())
        rows.append({"rfs.event": r, "os.event": o, "n": len(g),
                     "rfs==os": e, "pct": round(e / len(g) * 100, 1)})
    table = pd.DataFrame(rows)
    print(table.to_string(index=False))
    dead = full[(full["rfs.event"] == 0) & (full["os.event"] == 1)]
    eq_dead = int((dead["rfs.delay"] == dead["os.delay_(months)"]).sum())
    print(f"\n  died without recorded recurrence : {len(dead)}")
    print(f"  rfs.delay == os.delay            : {eq_dead}/{len(dead)}")
    print("  => RFS time censored at death; record not counted as an event.")
    print("  => death without recurrence = CENSORING (inferred from data;")
    print("     NOT stated in the publication or supplementary methods).")
    findings["death_without_recurrence"] = {
        "n_died_no_recurrence": len(dead), "n_time_equals_os": eq_dead,
        "inference": "censoring", "source": "data pattern, not publication",
    }
    table.to_csv(OUT / "event_os_crosstab.csv", index=False)

    # ---------------------------------------------------------------
    section("5. ZERO-TIME RECORDS (reported verbatim, NOT excluded)")
    # ---------------------------------------------------------------
    cols = ["sample_id", "title", "dataset", "tnm.stage", "tnm.t", "tnm.n",
            "tnm.m", "rfs.event", "rfs.delay", "os.event",
            "os.delay_(months)", "chemotherapy.adjuvant",
            "age.at.diagnosis_(year)", "sex", "mmr.status",
            "cit.molecularsubtype"]
    zero = s23[s23["rfs.delay"] == 0]
    print(zero[cols].to_string(index=False))
    zero[cols].to_csv(OUT / "zero_time_records.csv", index=False)

    grp_a = zero[zero["rfs.event"] == 1]
    grp_b = zero[zero["rfs.event"] == 0]
    print(f"\n  Pattern A -- rfs.event=1, rfs.delay=0 : {len(grp_a)}"
          f"  (os.delay>0 in {int((grp_a['os.delay_(months)'] > 0).sum())})")
    print(f"  Pattern B -- rfs.event=0, rfs.delay=0 : {len(grp_b)}"
          f"  (os.delay==0 in {int((grp_b['os.delay_(months)'] == 0).sum())},"
          f" os.event==1 in {int((grp_b['os.event'] == 1).sum())})")

    print("\n  -- EVIDENCE: is 'rfs.delay=0 + rfs.event=1' a coding convention? --")
    s4 = tum[tum["tnm.stage"] == 4]
    s4z = s4[s4["rfs.delay"] == 0]
    print(f"     stage IV tumours                     : {len(s4)}")
    print(f"     stage IV with rfs.delay == 0         : {len(s4z)}"
          f" ({len(s4z) / len(s4) * 100:.0f}%)")
    print(f"     of those, rfs.event == 1             : "
          f"{int((s4z['rfs.event'] == 1).sum())}/{len(s4z)}")
    print("     => in stage IV (metastatic at baseline) this pattern is the")
    print("        dominant coding. It is consistent with 'disease present at")
    print("        baseline', NOT with a recurrence observed in month 0.")

    print("\n  -- EVIDENCE: is os.delay==0 a true value or a missing sentinel? --")
    vc = tum["os.delay_(months)"].value_counts().sort_index()
    print(f"     os.delay == 0: {int(vc.get(0.0, 0))} | == 1: {int(vc.get(1.0, 0))}"
          f" | == 2: {int(vc.get(2.0, 0))} | == 3: {int(vc.get(3.0, 0))}")
    print("     => no spike at 0 relative to 1 and 2, so the data do NOT")
    print("        support a sentinel interpretation; but they do not rule it")
    print("        out either. UNRESOLVED.")

    findings["zero_time"] = {
        "n": len(zero),
        "pattern_a_event_at_zero": len(grp_a),
        "pattern_b_zero_followup": len(grp_b),
        "sample_ids": zero["sample_id"].tolist(),
        "stage4_total": len(s4),
        "stage4_rfs_delay_zero": len(s4z),
        "stage4_zero_with_event": int((s4z["rfs.event"] == 1).sum()),
        "os_delay_value_counts_low": {
            str(k): int(vc.get(k, 0)) for k in (0.0, 1.0, 2.0, 3.0)},
    }

    # ---------------------------------------------------------------
    section("6. TNM CONSISTENCY")
    # ---------------------------------------------------------------
    m1 = tum[(tum["tnm.m"] == "M1") & (tum["tnm.stage"] != 4)]
    print(f"  M1 but stage != IV : {len(m1)}")
    if len(m1):
        print(m1[["sample_id", "title", "dataset", "tnm.stage", "tnm.t",
                  "tnm.n", "tnm.m", "rfs.event", "rfs.delay"]]
              .to_string(index=False))
    print(f"  N1/N2 but stage not III/IV : "
          f"{len(tum[tum['tnm.n'].isin(['N1', 'N2']) & ~tum['tnm.stage'].isin([3, 4])])}")
    print(f"  N0 but stage III           : "
          f"{len(tum[(tum['tnm.n'] == 'N0') & (tum['tnm.stage'] == 3)])}")
    print(f"  missing tnm.m in stage II-III : {int(s23['tnm.m'].isna().sum())}")
    findings["tnm_consistency"] = {
        "m1_not_stage4": m1["sample_id"].tolist(),
        "missing_m_stage23": int(s23["tnm.m"].isna().sum()),
    }

    # ---------------------------------------------------------------
    section("7. PAIRED NON-TUMOURAL SAMPLES")
    # ---------------------------------------------------------------
    nt = meta[~meta["source_name"].str.contains(
        TUMOUR_PATTERN, case=False, na=False)]
    dep = meta[meta["dependancy_sample"].notna()]
    linked = set(dep["dependancy_sample"].dropna())
    tum_titles = set(tum["title"])
    print(f"  non-tumoural samples        : {len(nt)}")
    print(f"  with a dependancy link      : {int(nt['dependancy_sample'].notna().sum())}")
    print(f"  unlinked non-tumoural       : {int(nt['dependancy_sample'].isna().sum())}")
    print(f"  links resolving to a tumour : {len(linked & tum_titles)}/{len(linked)}")
    findings["paired_normals"] = {
        "n_non_tumoural": len(nt), "n_linked": len(linked),
        "links_resolve": len(linked & tum_titles),
    }

    # ---------------------------------------------------------------
    section("8. ENDPOINT COMPLETENESS BY SUBSET")
    # ---------------------------------------------------------------
    rows = []
    for ds in ["discovery", "validation"]:
        g = s23[s23["dataset"] == ds]
        cg = g[g[["rfs.event", "rfs.delay"]].notna().all(axis=1)]
        rows.append({
            "subset": ds, "stage II-III": len(g), "documented RFS": len(cg),
            "events": int((cg["rfs.event"] == 1).sum()),
            "censored": int((cg["rfs.event"] == 0).sum()),
            "event rate %": round((cg["rfs.event"] == 1).mean() * 100, 1),
            "zero-time": int((cg["rfs.delay"] == 0).sum()),
            "median FU (mo)": float(cg["rfs.delay"].median()),
            "events <=60mo": int(((cg["rfs.event"] == 1)
                                  & (cg["rfs.delay"] <= 60)).sum()),
        })
    comp = pd.DataFrame(rows)
    print(comp.to_string(index=False))
    comp.to_csv(OUT / "endpoint_completeness.csv", index=False)
    findings["completeness"] = comp.to_dict(orient="records")

    # ---------------------------------------------------------------
    section("9. COVARIATE AUDIT (stage II-III, documented RFS)")
    # ---------------------------------------------------------------
    cov = ["age.at.diagnosis_(year)", "sex", "tnm.stage", "tnm.t", "tnm.n",
           "tnm.m", "tumor.location", "chemotherapy.adjuvant",
           "chemotherapy.adjuvant.type", "mmr.status", "cimp.status",
           "cin.status", "kras.mutation", "braf.mutation", "tp53.mutation",
           "cit.molecularsubtype", "rna.extraction.method"]
    rows = []
    for c in cov:
        nn = int(complete[c].notna().sum())
        rows.append({
            "covariate": c, "complete": nn, "missing": len(complete) - nn,
            "missing %": round((len(complete) - nn) / len(complete) * 100, 1),
            "levels": int(complete[c].nunique(dropna=True)),
            "events with data": int(
                complete.loc[complete[c].notna(), "rfs.event"].sum()),
        })
    covt = pd.DataFrame(rows)
    print(covt.to_string(index=False))
    covt.to_csv(OUT / "covariate_audit.csv", index=False)
    n_events = int((complete["rfs.event"] == 1).sum())
    print(f"\n  total events = {n_events} -> ~{n_events // 10} parameters at 10 EPV")
    findings["covariates"] = {"n_events": n_events,
                              "epv10_budget": n_events // 10}

    # ---------------------------------------------------------------
    section("10. ADJUVANT CHEMOTHERAPY BY STAGE (confounding by indication)")
    # ---------------------------------------------------------------
    ev = complete.groupby(
        [complete["chemotherapy.adjuvant"].fillna("MISSING"),
         complete["tnm.stage"]])["rfs.event"].agg(["size", "sum"])
    ev.columns = ["n", "events"]
    ev["event %"] = (ev["events"] / ev["n"] * 100).round(1)
    print(ev.to_string())
    ev.to_csv(OUT / "chemo_by_stage.csv")
    findings["chemo_by_stage"] = ev.reset_index().to_dict(orient="records")

    # ---------------------------------------------------------------
    section("11. CROSS-CHECK vs PUBLICATION TABLE S1 (per-patient)")
    # ---------------------------------------------------------------
    s1_path = fetch_table_s1(paths.ROOT, allow_fetch=not args.no_fetch)
    if s1_path is None:
        print("  Table S1 unavailable; cross-check skipped.")
        findings["table_s1"] = {"status": "unavailable"}
    else:
        try:
            s1 = pd.read_excel(s1_path, sheet_name=0, header=1)
            s1 = s1[s1["id"].astype(str).str.startswith("CIT")].copy()
            mg = meta.merge(s1, left_on="title", right_on="id", how="inner")
            print(f"  Table S1 per-patient rows      : {len(s1)}")
            print(f"  matched to GEO by CIT id       : {len(mg)}")
            print(f"  Table S1 columns               : {list(s1.columns)}")
            print("  NOTE: Table S1 carries TNM.Stage but NO TNM.T/N/M and no OS,")
            print("        so it cannot arbitrate a stage-vs-M conflict.")

            mg["stage_s1"] = pd.to_numeric(mg["TNM.Stage"], errors="coerce")
            sd = mg[(mg["tnm.stage"] != mg["stage_s1"])
                    & ~(mg["tnm.stage"].isna() & mg["stage_s1"].isna())]
            ed = mg[(mg["rfs.event"] != mg["RFS.event"])
                    & ~(mg["rfs.event"].isna() & mg["RFS.event"].isna())]
            dd = mg[(mg["rfs.delay"] != mg["RFS.delay"])
                    & ~(mg["rfs.delay"].isna() & mg["RFS.delay"].isna())]
            print(f"\n  stage disagreements GEO vs S1  : {len(sd)}")
            if len(sd):
                print(sd[["title", "sample_id", "tnm.stage", "stage_s1",
                          "tnm.m"]].to_string(index=False))
            print(f"  rfs.event disagreements        : {len(ed)}")
            print(f"  rfs.delay disagreements        : {len(dd)}")

            zt = mg[mg["title"].isin(
                zero["title"].tolist())][
                ["title", "sample_id", "tnm.stage", "tnm.m", "TNM.Stage",
                 "rfs.event", "RFS.event", "rfs.delay", "RFS.delay"]]
            print("\n  the 8 zero-time records as published:")
            print(zt.to_string(index=False))
            zt.to_csv(OUT / "zero_time_vs_tableS1.csv", index=False)

            findings["table_s1"] = {
                "status": "ok", "rows": len(s1), "matched": len(mg),
                "has_tnm_mnm": False, "has_os": False,
                "stage_disagreements": sd["title"].tolist(),
                "rfs_event_disagreements": len(ed),
                "rfs_delay_disagreements": len(dd),
            }
        except Exception as exc:                       # noqa: BLE001
            print(f"  Table S1 parse failed ({type(exc).__name__}: {exc}).")
            print("  (.xls needs the optional 'xlrd' reader.)")
            findings["table_s1"] = {"status": f"parse-failed: {type(exc).__name__}"}

    # ---------------------------------------------------------------
    section("12. FIVE-YEAR ADMINISTRATIVE CENSORING -- EXACT MECHANICS")
    # ---------------------------------------------------------------
    print("  Rule (proposed, NOT applied here): for every record,")
    print("     if rfs.delay > 60:  time := 60  AND  event := 0")
    print("     else:               time and event unchanged")
    print("  An event at month 75 becomes a CENSORED observation at month 60.")
    print("  Original columns are never overwritten; derived columns only.\n")
    for ds in ["discovery", "validation"]:
        g = complete[complete["dataset"] == ds]
        ev_all = int((g["rfs.event"] == 1).sum())
        ev_60 = int(((g["rfs.event"] == 1) & (g["rfs.delay"] <= 60)).sum())
        reclass = int(((g["rfs.event"] == 1) & (g["rfs.delay"] > 60)).sum())
        trunc = int((g["rfs.delay"] > 60).sum())
        print(f"  {ds:11s} n={len(g):4d} | events {ev_all:3d} -> {ev_60:3d} "
              f"| events reclassified as censored: {reclass:2d} "
              f"| times truncated to 60: {trunc:3d}")
    findings["five_year_censoring"] = {
        ds: {
            "events_uncensored": int(
                (complete[complete["dataset"] == ds]["rfs.event"] == 1).sum()),
            "events_within_60": int(
                ((complete[complete["dataset"] == ds]["rfs.event"] == 1)
                 & (complete[complete["dataset"] == ds]["rfs.delay"] <= 60)).sum()),
            "times_truncated": int(
                (complete[complete["dataset"] == ds]["rfs.delay"] > 60).sum()),
        } for ds in ["discovery", "validation"]
    }

    # ---------------------------------------------------------------
    section("13. FIVE-YEAR CUTOFF -- RECORD-LEVEL VERIFICATION")
    # ---------------------------------------------------------------
    print("  PRECONDITION: a '60' cutoff presumes rfs.delay is in MONTHS.")
    print("  That unit is PROVISIONAL (see section 15). The rule below cannot")
    print("  be confirmed as the primary analysis until the unit is settled.\n")
    vcols = ["sample_id", "title", "tnm.stage", "rfs.event", "rfs.delay",
             "os.event", "os.delay_(months)", "chemotherapy.adjuvant"]
    recl = complete[(complete["dataset"] == "discovery")
                    & (complete["rfs.event"] == 1)
                    & (complete["rfs.delay"] > 60)]
    print("  Discovery records whose EVENT STATUS would be reclassified:")
    print(recl[vcols].sort_values("rfs.delay").to_string(index=False))
    print(f"\n    count                       : {len(recl)}")
    print(f"    all rfs.event == 1          : {bool((recl['rfs.event'] == 1).all())}")
    print(f"    all rfs.delay > 60          : {bool((recl['rfs.delay'] > 60).all())}")
    print(f"    all rfs.delay <= os.delay   : "
          f"{bool((recl['rfs.delay'] <= recl['os.delay_(months)']).all())}")
    recl[vcols].to_csv(OUT / "five_year_reclassified.csv", index=False)

    print("\n  Why validation reclassifies ZERO despite records beyond 60:")
    for ds in ["discovery", "validation"]:
        g = complete[complete["dataset"] == ds]
        over = g[g["rfs.delay"] > 60]
        print(f"    {ds:11s} >60: {len(over):3d}"
              f" | events among them: {int((over['rfs.event'] == 1).sum()):2d}"
              f" | censored among them: {int((over['rfs.event'] == 0).sum()):3d}")
        print(f"    {'':11s} latest EVENT time: "
              f"{g.loc[g['rfs.event'] == 1, 'rfs.delay'].max():.0f}"
              f" | latest CENSORED time: "
              f"{g.loc[g['rfs.event'] == 0, 'rfs.delay'].max():.0f}")
    print("    => truncation changes TIME for every record >60, but changes")
    print("       EVENT STATUS only for records that were events. Validation's")
    print("       latest recurrence is well under 60, so all of its >60 records")
    print("       are already censored and truncation is a no-op for the count.")
    findings["five_year_reclassified"] = recl["sample_id"].tolist()

    # ---------------------------------------------------------------
    section("14. ZERO-TIME SCENARIOS (counts only; nothing applied)")
    # ---------------------------------------------------------------
    gA = complete[(complete["rfs.delay"] == 0) & (complete["rfs.event"] == 1)]
    gB = complete[(complete["rfs.delay"] == 0) & (complete["rfs.event"] == 0)]
    print(f"  Group A (rfs.event=1, rfs.delay=0): {len(gA)}  "
          f"{gA['sample_id'].tolist()}")
    print(f"  Group B (rfs.event=0, rfs.delay=0): {len(gB)}  "
          f"{gB['sample_id'].tolist()}\n")
    scen = []
    for name, drop in [("S0 retain all 8", pd.Index([])),
                       ("S1 exclude Group A only", gA.index),
                       ("S2 exclude Group B only", gB.index),
                       ("S3 exclude both (A+B)", gA.index.union(gB.index))]:
        k = complete.drop(index=drop)
        scen.append({
            "scenario": name, "n": len(k),
            "events": int((k["rfs.event"] == 1).sum()),
            "censored": int((k["rfs.event"] == 0).sum()),
            "disc n": int((k["dataset"] == "discovery").sum()),
            "disc events": int(((k["dataset"] == "discovery")
                                & (k["rfs.event"] == 1)).sum()),
            "val n": int((k["dataset"] == "validation").sum()),
            "val events": int(((k["dataset"] == "validation")
                               & (k["rfs.event"] == 1)).sum()),
        })
    sdf = pd.DataFrame(scen)
    print(sdf.to_string(index=False))
    sdf.to_csv(OUT / "zero_time_scenarios.csv", index=False)

    cit446 = complete[complete["title"] == "CIT446"]
    print(f"\n  CIT446 sensitivity (independent of the above):")
    print(f"    with CIT446    : n={len(complete)}, "
          f"events={int((complete['rfs.event'] == 1).sum())}")
    k = complete.drop(index=cit446.index)
    print(f"    without CIT446 : n={len(k)}, "
          f"events={int((k['rfs.event'] == 1).sum())}")
    findings["zero_time_scenarios"] = scen

    # ---------------------------------------------------------------
    section("15. TIME-UNIT DOCUMENTATION SEARCH -- SOURCES EXHAUSTED")
    # ---------------------------------------------------------------
    searched = [
        ("GEO SOFT !Sample_characteristics_ch1",
         "rfs.delay has NO unit label; os.delay IS labelled '(months)'"),
        ("Publication main text (PLoS Med)", "no unit statement for the column"),
        ("Text S1 supplementary methods",
         "'five years/60 months' -- the ANALYSIS HORIZON, not the column"),
        ("Table S1 per-patient annotation",
         "single sheet, no legend, no units, RFS.delay unlabelled"),
        ("ArrayExpress E-GEOD-39582 SDRF",
         "Characteristics [rfs.delay] present; NO Unit[] column anywhere; "
         "no overall-survival column at all"),
        ("GEO SuperSeries GSE40967", "no RFS metadata whatsoever"),
    ]
    for src, result in searched:
        print(f"  {src:40s} -> {result}")
    print("\n  OUTCOME: no source labels the unit. Status remains PROVISIONAL.")
    print("  Remaining untried avenue: direct correspondence with the CIT")
    print("  programme / corresponding author.")
    findings["time_unit_search"] = [
        {"source": s, "result": r} for s, r in searched]

    (OUT / "phase1_findings.json").write_text(
        json.dumps(findings, indent=2, default=str) + "\n")

    section("DONE")
    print(f"  evidence written to {OUT.relative_to(paths.ROOT)}/")
    print("  No cohort or endpoint decision was made. Nothing was excluded,")
    print("  recoded, imputed or censored. Original files unchanged.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

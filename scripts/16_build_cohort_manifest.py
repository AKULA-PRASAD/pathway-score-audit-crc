"""Phase 1 (final): build the cohort membership manifest.

Produces a per-sample manifest recording, for every deposited sample, which
eligibility criteria it meets and which analysis scenario it belongs to.

WHAT THIS IS NOT
----------------
This does not filter, modify, or create an analysis dataset. It writes a
LABEL table. Original source values are copied verbatim and never altered;
`metadata/`, `data/raw/`, the cached matrix and scripts 01-14 are untouched.

The two exclusion criteria are kept SEPARATE throughout, because they rest on
different justifications and an external cohort may contain one without the
other:

  Criterion 1 -- "recurrence-free at the time origin"
      Removes records with rfs.delay == 0 AND rfs.event == 1 (Group A).
      Conceptual: a patient recorded as recurring at surgery was never
      recurrence-free, so was never at risk.

  Criterion 2 -- "> 0 units of recorded follow-up"
      Removes records with rfs.delay == 0 AND rfs.event == 0 (Group B).
      Standard follow-up eligibility; these contribute zero person-time.

Scenarios: S0 retain all | S1 apply C1 only | S2 apply C2 only | S3 apply both.
S3 is the APPROVED PROPOSED PRIMARY cohort; S0/S1/S2 are sensitivity analyses.

Derived five-year columns are emitted as SEPARATE columns and are
CONDITIONAL on the provisional assumption that rfs.delay is in months.

Usage:
    .venv/bin/python scripts/16_build_cohort_manifest.py
"""

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths  # noqa: E402

OUT = paths.ROOT / "reports" / "phase1"
TUMOUR_PATTERN = "primary colorectal Adenocarcinoma"
HORIZON = 60  # CONDITIONAL on the provisional months assumption


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(paths.SAMPLE_METADATA)

    m = meta.copy()

    # --- eligibility flags (labels only; no row is dropped here) ---
    m["is_primary_tumour"] = m["source_name"].str.contains(
        TUMOUR_PATTERN, case=False, na=False)
    m["is_stage_2_3"] = m["tnm.stage"].isin([2, 3])
    m["has_documented_rfs"] = m[["rfs.event", "rfs.delay"]].notna().all(axis=1)

    zero = m["has_documented_rfs"] & (m["rfs.delay"] == 0)
    m["fails_criterion_1"] = zero & (m["rfs.event"] == 1)   # Group A
    m["fails_criterion_2"] = zero & (m["rfs.event"] == 0)   # Group B
    m["zero_time_group"] = ""
    m.loc[m["fails_criterion_1"], "zero_time_group"] = "A"
    m.loc[m["fails_criterion_2"], "zero_time_group"] = "B"

    base = m["is_primary_tumour"] & m["is_stage_2_3"] & m["has_documented_rfs"]
    m["in_S0"] = base
    m["in_S1"] = base & ~m["fails_criterion_1"]
    m["in_S2"] = base & ~m["fails_criterion_2"]
    m["in_S3"] = base & ~m["fails_criterion_1"] & ~m["fails_criterion_2"]
    m["is_primary_cohort"] = m["in_S3"]

    # --- derived endpoint columns: NEW columns, originals untouched ---
    m["rfs_time_full"] = m["rfs.delay"]
    m["rfs_event_full"] = m["rfs.event"]
    over = m["has_documented_rfs"] & (m["rfs.delay"] > HORIZON)
    m["rfs_time_5y"] = m["rfs.delay"].where(~over, HORIZON)
    m["rfs_event_5y"] = m["rfs.event"].where(~over, 0.0)
    m["reclassified_by_horizon"] = over & (m["rfs.event"] == 1)

    # verify originals really are untouched
    src = pd.read_csv(paths.SAMPLE_METADATA)
    for col in ["rfs.delay", "rfs.event", "os.event", "os.delay_(months)",
                "tnm.stage", "tnm.m"]:
        assert m[col].equals(src[col]), f"source column {col} was modified"
    print("  source columns verified byte-identical to metadata/ ✓")

    cols = ["sample_id", "title", "dataset", "source_name", "tnm.stage",
            "tnm.m", "rfs.event", "rfs.delay", "os.event",
            "os.delay_(months)", "is_primary_tumour", "is_stage_2_3",
            "has_documented_rfs", "zero_time_group", "fails_criterion_1",
            "fails_criterion_2", "in_S0", "in_S1", "in_S2", "in_S3",
            "is_primary_cohort", "rfs_time_full", "rfs_event_full",
            "rfs_time_5y", "rfs_event_5y", "reclassified_by_horizon"]
    m[cols].to_csv(OUT / "cohort_manifest.csv", index=False)

    # ---------------- cohort flow ----------------
    print("\n" + "=" * 72)
    print("COHORT FLOW (approved primary = S3)")
    print("=" * 72)
    flow = [
        ("deposited samples", len(m)),
        ("primary tumours", int(m["is_primary_tumour"].sum())),
        ("stage II-III primary tumours",
         int((m["is_primary_tumour"] & m["is_stage_2_3"]).sum())),
        ("with documented RFS (event and time)", int(base.sum())),
        ("minus Criterion 1 failures (Group A, n=3)",
         int((base & ~m["fails_criterion_1"]).sum())),
        ("minus Criterion 2 failures (Group B, n=5)", int(m["in_S3"].sum())),
    ]
    prev = None
    for label, n in flow:
        d = "" if prev is None else f"   (-{prev - n})"
        print(f"  {label:44s} {n:4d}{d}")
        prev = n

    # ---------------- counts per scenario ----------------
    print("\n" + "=" * 72)
    print("SAMPLE AND EVENT COUNTS BY SCENARIO")
    print("=" * 72)
    rows = []
    for s, desc in [("S3", "PRIMARY (C1+C2 applied)"),
                    ("S0", "sensitivity: retain all 8"),
                    ("S1", "sensitivity: C1 only"),
                    ("S2", "sensitivity: C2 only")]:
        g = m[m[f"in_{s}"]]
        r = {"scenario": s, "description": desc, "n": len(g),
             "events_full": int((g["rfs_event_full"] == 1).sum()),
             "censored_full": int((g["rfs_event_full"] == 0).sum()),
             "events_5y": int((g["rfs_event_5y"] == 1).sum()),
             "median_fu_full": float(g["rfs_time_full"].median())}
        for ds in ["discovery", "validation"]:
            h = g[g["dataset"] == ds]
            r[f"{ds[:4]}_n"] = len(h)
            r[f"{ds[:4]}_ev_full"] = int((h["rfs_event_full"] == 1).sum())
            r[f"{ds[:4]}_ev_5y"] = int((h["rfs_event_5y"] == 1).sum())
        rows.append(r)
    counts = pd.DataFrame(rows)
    print(counts.to_string(index=False))
    counts.to_csv(OUT / "cohort_counts.csv", index=False)

    print("\n  NOTE: *_5y columns are CONDITIONAL on the provisional assumption")
    print("        that rfs.delay is in months. Not a verified fact.")

    # ---------------- excluded records, named ----------------
    print("\n" + "=" * 72)
    print("EXCLUDED RECORDS (named, values verbatim)")
    print("=" * 72)
    ex = m[base & (m["fails_criterion_1"] | m["fails_criterion_2"])]
    print(ex[["sample_id", "title", "dataset", "tnm.stage", "tnm.m",
              "rfs.event", "rfs.delay", "os.event", "os.delay_(months)",
              "zero_time_group"]].to_string(index=False))
    ex.to_csv(OUT / "excluded_records.csv", index=False)

    print("\n  CIT446 (GSM972402) is excluded by Criterion 1 only.")
    print("  Its stage/M conflict is NOT recoded and NOT used as a reason.")

    summary = {
        "approved_primary_scenario": "S3",
        "criteria": {
            "criterion_1": "recurrence-free at time origin (excludes Group A)",
            "criterion_2": "> 0 units of recorded follow-up (excludes Group B)",
        },
        "horizon_months": HORIZON,
        "horizon_status": "CONDITIONAL on provisional months assumption",
        "flow": dict(flow),
        "counts": counts.to_dict(orient="records"),
        "excluded_sample_ids": ex["sample_id"].tolist(),
        "source_values_modified": False,
    }
    (OUT / "cohort_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")

    print("\n" + "=" * 72)
    print("  manifest -> reports/phase1/cohort_manifest.csv")
    print("  counts   -> reports/phase1/cohort_counts.csv")
    print("  excluded -> reports/phase1/excluded_records.csv")
    print("  Labels only. No source value altered; no dataset filtered.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

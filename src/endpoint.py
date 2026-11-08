"""Endpoint derivation.

Hierarchy (approved):
  PRIMARY     rfs_event_full / rfs_time_full      -- unconditional
  CONDITIONAL rfs_event_60units / rfs_time_60units -- depends on the
              PROVISIONAL assumption that rfs.delay is in months
  SEPARATE    os_event_source / os_time_source     -- preserved verbatim,
              NOT combined with RFS here

Source columns are copied verbatim and never altered. Derived columns are
additional; nothing is overwritten.

The conditional columns are named `60units`, not `5y`. Naming them `5y` would
assert a time unit that no source documents.
"""

from __future__ import annotations

import pandas as pd

from src.provenance import (CENSORING_INTERPRETATION, ENDPOINT_HIERARCHY,
                            TIME_UNIT_ASSUMPTION)

HORIZON_UNITS = 60

SRC_RFS_EVENT = "rfs.event"
SRC_RFS_TIME = "rfs.delay"
SRC_OS_EVENT = "os.event"
SRC_OS_TIME = "os.delay_(months)"


def build_endpoint(meta: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    """Endpoint table for `sample_ids`, aligned by label, in the given order."""
    idx = meta.set_index("sample_id")
    missing = [s for s in sample_ids if s not in idx.index]
    if missing:
        raise KeyError(f"{len(missing)} samples absent from metadata: "
                       f"{missing[:5]}")
    m = idx.loc[sample_ids]

    out = pd.DataFrame(index=pd.Index(sample_ids, name="sample_id"))

    # ---- verbatim source copies -------------------------------------
    out["rfs_event_source"] = m[SRC_RFS_EVENT].to_numpy()
    out["rfs_time_source"] = m[SRC_RFS_TIME].to_numpy()
    out["os_event_source"] = m[SRC_OS_EVENT].to_numpy()
    out["os_time_source"] = m[SRC_OS_TIME].to_numpy()

    # ---- PRIMARY: full follow-up ------------------------------------
    out["rfs_event_full"] = out["rfs_event_source"]
    out["rfs_time_full"] = out["rfs_time_source"]

    # ---- CONDITIONAL: 60-unit administrative horizon ----------------
    over = out["rfs_time_source"] > HORIZON_UNITS
    out["rfs_event_60units"] = out["rfs_event_source"].where(~over, 0.0)
    out["rfs_time_60units"] = out["rfs_time_source"].where(
        ~over, float(HORIZON_UNITS))
    out["reclassified_by_horizon"] = (
        over & (out["rfs_event_source"] == 1)).astype(bool)

    return out.reset_index()


def endpoint_metadata() -> dict:
    return {
        "hierarchy": ENDPOINT_HIERARCHY,
        "primary_endpoint": "rfs_event_full / rfs_time_full",
        "conditional_endpoint": "rfs_event_60units / rfs_time_60units",
        "conditional_on": TIME_UNIT_ASSUMPTION,
        "horizon_units": HORIZON_UNITS,
        "horizon_rule": ("if rfs_time_source > 60: time := 60 and event := 0; "
                         "else unchanged"),
        "time_unit_assumption": TIME_UNIT_ASSUMPTION,
        "censoring_interpretation": CENSORING_INTERPRETATION,
        "os_fields": "preserved separately and verbatim; NOT combined with RFS",
        "unresolved": [
            "time unit of rfs.delay (PROVISIONAL: months)",
            "death-without-recurrence censoring (PROVISIONAL)",
            "zero-time records (UNRESOLVED; see cohort scenarios)",
        ],
    }

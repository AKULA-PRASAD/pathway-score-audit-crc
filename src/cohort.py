"""Cohort materialisation from the approved Phase 1 manifest.

This module READS the frozen Phase 1 cohort manifest. It never re-implements
eligibility logic and never recomputes membership: the manifest is the single
source of truth, approved in Phase 1.

Outcome-blind in the sense that matters here: no endpoint value is inspected
to decide membership. Membership is looked up, not derived.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from src import paths

MANIFEST = paths.ROOT / "reports" / "phase1" / "cohort_manifest.csv"
COHORTS = ("S0", "S1", "S2", "S3")
PRIMARY_COHORT = "S3"

CRITERION_LABELS = {
    "fails_criterion_1": "C1: not recurrence-free at the time origin "
                         "(rfs.delay == 0 AND rfs.event == 1)",
    "fails_criterion_2": "C2: zero units of recorded follow-up "
                         "(rfs.delay == 0 AND rfs.event == 0)",
}


@dataclass
class Cohort:
    name: str
    sample_ids: list[str]
    manifest_rows: pd.DataFrame

    def __len__(self) -> int:
        return len(self.sample_ids)


def load_manifest() -> pd.DataFrame:
    if not MANIFEST.exists():
        raise FileNotFoundError(
            f"Phase 1 cohort manifest not found at {MANIFEST}. "
            "Run scripts/16_build_cohort_manifest.py first.")
    return pd.read_csv(MANIFEST)


def canonical_order(manifest: pd.DataFrame,
                    reference_columns: list[str]) -> list[str]:
    """Canonical sample order = the probe cache's column order, restricted.

    Fixing the order against the cache (not against the manifest's row order)
    keeps every downstream artifact consistent with the expression source.
    """
    wanted = set(manifest["sample_id"])
    return [s for s in reference_columns if s in wanted]


def materialise(name: str, manifest: pd.DataFrame,
                reference_columns: list[str]) -> Cohort:
    """Return the cohort's sample ids in canonical order, by lookup only."""
    if name not in COHORTS:
        raise ValueError(f"unknown cohort {name!r}; expected {COHORTS}")
    col = f"in_{name}"
    if col not in manifest.columns:
        raise KeyError(f"manifest lacks column {col!r}")
    rows = manifest[manifest[col].astype(bool)].copy()
    ids = canonical_order(rows, reference_columns)
    if len(ids) != len(rows):
        missing = set(rows["sample_id"]) - set(reference_columns)
        raise ValueError(
            f"{len(missing)} {name} samples absent from the expression "
            f"matrix: {sorted(missing)[:5]}")
    rows = rows.set_index("sample_id").loc[ids].reset_index()
    return Cohort(name=name, sample_ids=ids, manifest_rows=rows)


def exclusion_table(manifest: pd.DataFrame) -> pd.DataFrame:
    """Every excluded sample with the named criterion that removed it."""
    # Restrict to the cohort BASE (stage II-III primary tumours with a
    # documented RFS). Zero-time records outside that base -- e.g. stage IV,
    # where rfs.delay == 0 is the coding convention for disease at baseline --
    # are not cohort exclusions and must not be listed as such.
    base = manifest[
        manifest["is_primary_tumour"].astype(bool)
        & manifest["is_stage_2_3"].astype(bool)
        & manifest["has_documented_rfs"].astype(bool)
    ]
    out = []
    for _, r in base.iterrows():
        for flag, label in CRITERION_LABELS.items():
            if bool(r.get(flag, False)):
                out.append({"sample_id": r["sample_id"], "title": r["title"],
                            "dataset": r["dataset"], "criterion": flag,
                            "criterion_description": label,
                            "rfs_event": r["rfs.event"],
                            "rfs_delay": r["rfs.delay"]})
    return pd.DataFrame(out).sort_values(["criterion", "sample_id"])

"""External sample annotation recovered from ArrayExpress E-GEOD-39582.

The GEO deposit does NOT expose the ComBat batch assignments used by the
depositors. The ArrayExpress SDRF for the same study does, as
`Characteristics [normalizationcombatbatch]` and `[normalizationcombatmod]`.

This module fetches and reads that file. It is DESCRIPTIVE ONLY: recovering
the batch labels lets us ASK whether a probe-selection criterion is tracking
batch structure. It performs no batch correction and no normalization, and
nothing here modifies the expression data.
"""

from __future__ import annotations

import hashlib
import urllib.request
from pathlib import Path

import pandas as pd

from src import paths

SDRF_URL = "https://www.ebi.ac.uk/biostudies/files/E-GEOD-39582/E-GEOD-39582.sdrf.txt"
SDRF_PATH = Path("data/external/E-GEOD-39582.sdrf.txt")

BATCH_COL = "Characteristics [normalizationcombatbatch]"
MOD_COL = "Characteristics [normalizationcombatmod]"


def fetch_sdrf(allow_fetch: bool = True) -> Path | None:
    dest = paths.ROOT / SDRF_PATH
    if dest.exists():
        return dest
    if not allow_fetch:
        return None
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        req = urllib.request.Request(SDRF_URL,
                                     headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            blob = r.read()
        if b"normalizationcombatbatch" not in blob:
            print("  WARNING: SDRF lacks the ComBat batch column; not saved.")
            return None
        dest.write_bytes(blob)
        print(f"  fetched SDRF -> {SDRF_PATH} "
              f"(sha256 {hashlib.sha256(blob).hexdigest()[:16]}...)")
        return dest
    except Exception as exc:                        # noqa: BLE001
        print(f"  WARNING: SDRF fetch failed ({type(exc).__name__}).")
        return None


def load_batch_assignments(allow_fetch: bool = True) -> pd.DataFrame | None:
    """Return sample_id -> combat_batch, combat_mod. None if unavailable."""
    path = fetch_sdrf(allow_fetch)
    if path is None:
        return None
    sdrf = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)

    src_col = next((c for c in sdrf.columns
                    if c.strip().lower() in ("source name", "sample name")), None)
    if src_col is None or BATCH_COL not in sdrf.columns:
        return None

    out = pd.DataFrame({
        "sdrf_sample": sdrf[src_col].astype(str).str.strip(),
        "combat_batch": sdrf[BATCH_COL].astype(str).str.strip(),
        "combat_mod": sdrf[MOD_COL].astype(str).str.strip()
        if MOD_COL in sdrf.columns else None,
    }).drop_duplicates("sdrf_sample")

    # SDRF source names carry the GSM accession; extract it.
    gsm = out["sdrf_sample"].str.extract(r"(GSM\d+)", expand=False)
    out["sample_id"] = gsm
    return out.dropna(subset=["sample_id"])[
        ["sample_id", "combat_batch", "combat_mod"]].reset_index(drop=True)

"""Resolve a recorded M2 specification hash against the amendment chain.

Earlier phases embedded the M2 panel hashes in their own records:

  * ``config/analysis/A19_FREEZE_RECORD.json``  -> ``upstream_frozen_hashes``
  * ``config/analysis/QC_DECISIONS_RECORD.json`` -> ``frozen_artifacts_unmodified``

Those are **historical attestations**: a statement that at Phase 6D / Phase 7C,
these were the bytes that phase actually read. They are deliberately NOT
rewritten when a later amendment changes the file, because rewriting them would
destroy the record of what was verified at the time.

A checker that compares such an attestation to *today's* disk therefore conflates
two different questions:

  1. "is this record intact, and did the file change in a way nobody approved?"
     -- the question worth asking, and
  2. "is the file byte-identical to what it was months ago?"
     -- false after ANY approved amendment, however well documented.

``resolve_hash`` answers (1). A recorded hash is accepted when it is either the
current file, or a version the M2 freeze record documents as superseded by an
approved amendment. Anything else is undocumented drift and still fails.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from src import paths

FREEZE_RECORD = paths.ROOT / "config" / "pathways" / "M2_FREEZE_RECORD.json"

#: files this module can resolve, by basename. The freeze record is included
#: because an amendment necessarily rewrites it too -- it is the document that
#: records the amendment -- so earlier phases' attestations about IT also need
#: resolving against the chain.
M2_SPEC_FILES = ("m2_panel_4set.json", "m2_panel_5set.json",
                 "M2_FREEZE_RECORD.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _record() -> dict:
    return json.loads(FREEZE_RECORD.read_text())


def _superseded_entries(filename: str) -> list[dict]:
    """Every superseded version of `filename`, from both places they are kept.

    Per-spec lists carry the panel history. ``original_freeze.sha256`` carries
    the pre-amendment bytes of all three files, including the freeze record
    itself -- which has no entry under ``specifications`` because it is the
    record, not a specification.
    """
    rec = _record()
    out: list[dict] = []
    for spec in rec["specifications"].values():
        if Path(spec["file"]).name != filename:
            continue
        for entry in spec.get("superseded_sha256", []):
            out.append(entry if isinstance(entry, dict)
                       else {"sha256": entry, "superseded_by": None})
    orig = rec.get("original_freeze", {}).get("sha256", {})
    if filename in orig:
        seen = {e["sha256"] for e in out}
        if orig[filename] not in seen:
            out.append({"sha256": orig[filename],
                        "superseded_by": rec["amendment_history"][0]["amendment_id"]
                        if rec.get("amendment_history") else None})
    # the freeze record itself: every amendment rewrites it, so it keeps its own
    # superseded list rather than appearing under "specifications".
    if filename == FREEZE_RECORD.name:
        seen = {e["sha256"] for e in out}
        for entry in rec.get("record_superseded_sha256", []):
            if entry["sha256"] not in seen:
                out.append(entry)
    return out


def superseded_hashes(filename: str) -> list[str]:
    """Hashes of `filename` that an approved amendment has superseded."""
    return [e["sha256"] for e in _superseded_entries(filename)]


def amendment_for(filename: str, digest: str) -> str | None:
    """Id of the amendment that superseded `digest`, or None."""
    for e in _superseded_entries(filename):
        if e["sha256"] == digest:
            return e.get("superseded_by")
    return None


def resolve_hash(path: Path, recorded: str) -> tuple[bool, str]:
    """Is `recorded` the current file, or a documented superseded version?

    Returns ``(ok, explanation)``. ``ok`` is False only for undocumented drift.
    """
    path = Path(path)
    name = path.name
    current = sha256(path)
    if recorded == current:
        return True, "current"
    if name in M2_SPEC_FILES and recorded in superseded_hashes(name):
        return True, f"superseded by {amendment_for(name, recorded)}"
    return False, (
        f"UNDOCUMENTED DRIFT in {name}: recorded {recorded[:16]} is neither the "
        f"current file ({current[:16]}) nor any hash the M2 freeze record "
        f"documents as superseded by an approved amendment")

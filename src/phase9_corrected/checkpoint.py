"""Append-only, self-verifying JSONL checkpoints for the corrected rerun.

File: <output_dir>/runs/run<r>.jsonl
  line 1     {"kind":"header", schema, run_id, config_sha256, manifest_sha256}
  line 2..   {"kind":"record", run_id, k, record_sha256, record}
record_sha256 = sha256(canonical_json(record)). Each line is written with one os.write
on an O_APPEND descriptor followed by fsync, so a crash can leave at most one partial
FINAL line. Loading is strict: any malformed line, hash mismatch, identity mismatch or
duplicate k raises. A partial final line is never repaired silently; `recover_partial_tail`
must be called explicitly and preserves the original file in quarantine first.

A per-run lock file (O_CREAT|O_EXCL) prevents two writers on one run.
"""
from __future__ import annotations

import json
import os
import shutil
import socket
from datetime import datetime, timezone
from pathlib import Path

from . import SCHEMA_VERSION
from .config import canonical_json, assert_output_dir_allowed

import hashlib


class CheckpointError(RuntimeError):
    pass


class CheckpointCorruption(CheckpointError):
    pass


class DuplicateReplicateError(CheckpointError):
    pass


class LockHeldError(CheckpointError):
    pass


def _sha(obj) -> str:
    return hashlib.sha256(canonical_json(obj)).hexdigest()


class CheckpointStore:
    def __init__(self, output_dir: Path, run_id: int, config_sha256: str, manifest_sha256: str,
                 root: Path):
        assert_output_dir_allowed(output_dir, root)
        self.dir = Path(output_dir) / "runs"
        self.run_id = int(run_id)
        self.config_sha256, self.manifest_sha256 = config_sha256, manifest_sha256
        self.path = self.dir / f"run{self.run_id}.jsonl"
        self.lock_path = self.dir / f"run{self.run_id}.lock"
        self._lock_fd = None
        self.done = {}

    # ------------------------------------------------------------------ locking
    def acquire(self):
        self.dir.mkdir(parents=True, exist_ok=True)
        try:
            fd = os.open(self.lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            raise LockHeldError(f"{self.lock_path} exists: another writer is active or a previous "
                                "run ended abnormally; inspect and remove it manually") from None
        os.write(fd, json.dumps({"pid": os.getpid(), "host": socket.gethostname(),
                                 "utc": datetime.now(timezone.utc).isoformat()}).encode())
        os.fsync(fd)
        self._lock_fd = fd

    def release(self):
        if self._lock_fd is not None:
            os.close(self._lock_fd)
            self._lock_fd = None
            self.lock_path.unlink(missing_ok=True)

    def __enter__(self):
        self.acquire()
        try:
            self.open()
        except BaseException:
            self.release()
            raise
        return self

    def __exit__(self, *exc):
        self.release()

    # ------------------------------------------------------------------ io
    def _header(self) -> dict:
        return {"kind": "header", "schema": SCHEMA_VERSION, "run_id": self.run_id,
                "config_sha256": self.config_sha256, "manifest_sha256": self.manifest_sha256}

    def _append_line(self, obj):
        data = canonical_json(obj) + b"\n"
        fd = os.open(self.path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
        try:
            n = os.write(fd, data)
            if n != len(data):
                raise CheckpointError("short write")
            os.fsync(fd)
        finally:
            os.close(fd)

    def open(self):
        """Create (with header) or strictly load an existing checkpoint."""
        if not self.path.exists():
            self.dir.mkdir(parents=True, exist_ok=True)
            self._append_line(self._header())
            self.done = {}
            return self
        self.done = self.load_records(self.path, expect_header=self._header())
        return self

    @staticmethod
    def load_records(path: Path, expect_header: dict | None = None) -> dict:
        raw = Path(path).read_bytes()
        if raw and not raw.endswith(b"\n"):
            raise CheckpointCorruption(f"{path}: partial final line (interrupted write)")
        lines = raw.split(b"\n")[:-1] if raw else []
        if not lines:
            raise CheckpointCorruption(f"{path}: empty checkpoint (no header)")
        try:
            objs = [json.loads(l) for l in lines]
        except json.JSONDecodeError as ex:
            raise CheckpointCorruption(f"{path}: malformed line: {ex}") from None
        head = objs[0]
        if head.get("kind") != "header":
            raise CheckpointCorruption(f"{path}: first line is not a header")
        if expect_header is not None and head != expect_header:
            raise CheckpointCorruption(f"{path}: header mismatch (different config/manifest/run)")
        done = {}
        for i, o in enumerate(objs[1:], 2):
            if o.get("kind") != "record":
                raise CheckpointCorruption(f"{path}:{i}: unexpected line kind")
            r = o["record"]
            if _sha(r) != o["record_sha256"]:
                raise CheckpointCorruption(f"{path}:{i}: record hash mismatch")
            if r["run_id"] != head["run_id"] or o["run_id"] != r["run_id"] or o["k"] != r["k"]:
                raise CheckpointCorruption(f"{path}:{i}: identity mismatch")
            if r["k"] in done:
                raise DuplicateReplicateError(f"{path}:{i}: duplicate k={r['k']}")
            done[r["k"]] = r
        return done

    def append(self, rec: dict):
        if self._lock_fd is None:
            raise CheckpointError("append without holding the run lock")
        if rec["run_id"] != self.run_id:
            raise CheckpointError("record run_id does not match this checkpoint")
        if rec["k"] in self.done:
            raise DuplicateReplicateError(f"k={rec['k']} already recorded for run {self.run_id}")
        # JSON round-trip first, so the in-memory copy equals what is on disk
        rec = json.loads(canonical_json(rec))
        self._append_line({"kind": "record", "run_id": self.run_id, "k": rec["k"],
                           "record_sha256": _sha(rec), "record": rec})
        self.done[rec["k"]] = rec

    def ordered(self, B: int | None = None) -> list:
        """Records in canonical k order; with B, requires exactly k = 0..B-1."""
        ks = sorted(self.done)
        if B is not None and ks != list(range(B)):
            raise CheckpointError(f"run {self.run_id}: k set is not 0..{B - 1}")
        return [self.done[k] for k in ks]

    def file_sha256(self) -> str:
        return hashlib.sha256(self.path.read_bytes()).hexdigest()

    # ------------------------------------------------------------------ recovery
    def recover_partial_tail(self, quarantine_dir: Path) -> dict:
        """Explicit recovery from an interrupted final write. Preserves the original file
        byte-for-byte in quarantine, then atomically rewrites the valid prefix. Refuses
        if any defect other than a single partial final line is present."""
        raw = self.path.read_bytes()
        if raw.endswith(b"\n"):
            raise CheckpointError("no partial final line: nothing to recover")
        cut = raw.rfind(b"\n") + 1
        prefix = raw[:cut]
        tmp = self.path.with_suffix(".jsonl.tmp")
        tmp.write_bytes(prefix)
        self.load_records(tmp, expect_header=self._header())   # prefix must be fully valid
        q = Path(quarantine_dir)
        q.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        saved = q / f"{self.path.name}.{stamp}.orig"
        shutil.copy2(self.path, saved)
        fd = os.open(tmp, os.O_RDONLY)
        os.fsync(fd); os.close(fd)
        os.replace(tmp, self.path)
        return {"quarantined_original": str(saved), "bytes_dropped": len(raw) - cut,
                "original_sha256": hashlib.sha256(raw).hexdigest()}

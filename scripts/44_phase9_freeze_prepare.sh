#!/bin/bash
# Phase 9 corrected pipeline -- Decision 13 P2: NON-SCIENTIFIC freeze-candidate preparation.
#
# Usage:  scripts/44_phase9_freeze_prepare.sh <candidate-output.json>
#
# Resolves the initial-freeze bootstrap: the scientific launcher (scripts/42) refuses to
# start Python until an owner-approved freeze record and anchor exist, so the FIRST record
# is produced here instead. This script:
#   * never invokes the scientific driver, the launcher, the verification gate, the
#     post-run runner or any scientific operation; it starts only the capture-only entry
#     point scripts/44_phase9_freeze_capture.py;
#   * uses the same interpreter and flags as the scientific launcher:
#       .venv/bin/python -I -B -X pycache_prefix=<fresh, empty, external directory>
#   * writes the candidate OUTSIDE the repository (it must not already exist) and never
#     touches config/analysis/PHASE9_CORRECTED_FREEZE_RECORD.json or the owner config;
#   * independently re-checks every inventory hash in the candidate with /usr/bin/shasum
#     and prints the candidate's SHA-256.
# It does NOT establish a trust anchor. The owner must independently re-hash and review the
# candidate before copying it into the repository and approving Decision 13
# (docs/PHASE9_DECISION13_FREEZE_PROCEDURE.md).
#
# Exit codes: 2 usage; 10 tools/location; 11 output exists; 16 inventory re-check failed;
# 19 output or bytecode directory refused; 20 bytecode directory not empty; 21 interpreter
# missing; 22 bytecode written; otherwise the capture entry point's own status.
set -euo pipefail
umask 077
export LC_ALL=C

die() {
  local code="$1"; shift
  printf 'p9c-freeze-prepare: REFUSED (code %s): %s\n' "$code" "$*" >&2
  exit "$code"
}

WORK=""
PYC=""
cleanup() {
  if [ -n "$WORK" ]; then /bin/rm -rf "$WORK"; fi
  if [ -n "$PYC" ]; then /bin/rm -rf "$PYC"; fi
  return 0
}
trap cleanup EXIT

JQ=/usr/bin/jq
SHASUM=/usr/bin/shasum
[ -x "$JQ" ] || die 10 "/usr/bin/jq is required"
[ -x "$SHASUM" ] || die 10 "/usr/bin/shasum is required"

SELF_REL="scripts/44_phase9_freeze_prepare.sh"
CAPTURE_REL="scripts/44_phase9_freeze_capture.py"

[ "$#" -eq 1 ] || die 2 "usage: $0 <candidate-output.json>"
OUT_ARG="$1"
OUT_DIR="$(cd "$(dirname "$OUT_ARG")" 2>/dev/null && pwd -P)" || die 19 "output directory does not exist"
OUT="$OUT_DIR/$(basename "$OUT_ARG")"
[ ! -e "$OUT" ] || die 11 "candidate output already exists: $OUT"

SELF_DIR="$(cd "$(dirname "$0")" && pwd -P)" || die 10 "cannot resolve the script directory"
ROOT="$(cd "$SELF_DIR/.." && pwd -P)" || die 10 "cannot resolve the repository root"
[ "$SELF_DIR/$(basename "$0")" = "$ROOT/$SELF_REL" ] || die 10 "script is not at $SELF_REL"
case "$OUT_DIR/" in
  "$ROOT"/*) die 19 "the candidate must be written outside the repository" ;;
esac
cd "$ROOT"

WORK="$(/usr/bin/mktemp -d "${TMPDIR:-/tmp}/p9c_prep.XXXXXX")" || die 19 "mktemp failed"
PYC="$(/usr/bin/mktemp -d "${TMPDIR:-/tmp}/p9c_pyc.XXXXXX")" || die 19 "mktemp failed"
PYC="$(cd "$PYC" && pwd -P)" || die 19 "cannot resolve the bytecode directory"
case "$PYC/" in
  "$ROOT"/*) die 19 "bytecode directory must be outside the repository" ;;
esac
[ -z "$(/bin/ls -A "$PYC")" ] || die 20 "bytecode directory is not empty"
PY="$ROOT/.venv/bin/python"
[ -x "$PY" ] || die 21 "virtual-environment interpreter not found: $PY"

set +e
"$PY" -I -B -X pycache_prefix="$PYC" "$ROOT/$CAPTURE_REL" "$OUT"
status=$?
set -e
if [ -n "$(/bin/ls -A "$PYC")" ]; then
  printf 'p9c-freeze-prepare: bytecode was written despite -B (policy violation)\n' >&2
  [ "$status" -ne 0 ] || status=22
fi
[ "$status" -eq 0 ] || exit "$status"

# ---- independent re-check of the candidate's inventory with /usr/bin/shasum
"$JQ" -er '.code | to_entries[] | "\(.value)  \(.key)"' "$OUT" > "$WORK/inventory.sha256" \
  || die 16 "cannot extract the candidate inventory"
if ! "$SHASUM" -a 256 --strict -c "$WORK/inventory.sha256" > "$WORK/check.txt" 2>&1; then
  /usr/bin/grep -v ': OK$' "$WORK/check.txt" >&2 || true
  die 16 "candidate inventory does not match the files on disk"
fi
CAND_SHA="$("$SHASUM" -a 256 "$OUT")"
printf 'p9c-freeze-prepare: CANDIDATE written (NOT approved, NOT a trust anchor)\n'
printf '  path:   %s\n  sha256: %s\n  files:  %s (all re-checked with shasum)\n' \
  "$OUT" "${CAND_SHA%% *}" "$(/usr/bin/wc -l < "$WORK/inventory.sha256" | /usr/bin/tr -d ' ')"
printf '  next:   owner re-hashes and reviews it (docs/PHASE9_DECISION13_FREEZE_PROCEDURE.md)\n'
exit 0

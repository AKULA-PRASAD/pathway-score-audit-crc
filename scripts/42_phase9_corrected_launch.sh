#!/bin/bash
# Phase 9 corrected-pipeline launcher -- owner decision DEC-E1 (2026-11-06).
#
# Usage:  scripts/42_phase9_corrected_launch.sh <validate|run> <owner-config.json> [driver args...]
#         scripts/42_phase9_corrected_launch.sh verify <owner-config.json>
#
# `validate` and `run` start the driver (scripts/41); `verify` starts the post-run
# verification runner (scripts/43, Decision 13 P1) under exactly the same checks and flags.
#
# Before Python starts, this launcher verifies -- using only sealed macOS system tools
# (/usr/bin/jq, /usr/bin/shasum) and never the scientific Python package -- that:
#   * the owner-approved Decision 13 entry is approved and both hashes are well formed;
#   * the freeze record's bytes match the owner-approved `freeze_record_sha256`;
#   * every file in the freeze record's code inventory (including THIS launcher and the
#     driver) matches its frozen SHA-256, and no inventory file is missing;
#   * the driver's inventory entry equals the owner-approved driver `sha256`;
#   * every Decision 13 required file (driver, launcher, runner, freeze-preparation scripts,
#     historical driver, src helpers, frozen pipeline modules) is in the inventory;
#   * the package directory contains exactly the frozen set of modules.
# It then runs the virtual-environment interpreter as
#   python -I -B -X pycache_prefix=<fresh, empty, external directory> <driver> ...
# so that no repository-local bytecode cache can be read and none is written. Python's
# exit status is preserved; the temporary directories are always removed.
#
# TRUST LIMITATION -- MANUAL ANCHOR STEP. This launcher cannot authenticate itself: a
# modified launcher could skip these checks. Before each execution campaign the owner must
# independently confirm (1) shasum -a 256 config/analysis/PHASE9_CORRECTED_FREEZE_RECORD.json
# equals the approved `freeze_record_sha256`, and (2) shasum -a 256 of this file equals its
# entry in that freeze record. The launcher re-verifies its own entry on every invocation.
# These are source-file and launch-flag checks, not an attestation of in-memory code; a
# short time-of-check/time-of-use window and the unhashed interpreter binary remain.
set -euo pipefail
umask 077
export LC_ALL=C

die() {
  local code="$1"; shift
  printf 'p9c-launch: REFUSED (code %s): %s\n' "$code" "$*" >&2
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

SELF_REL="scripts/42_phase9_corrected_launch.sh"
DRIVER_REL="scripts/41_phase9_corrected_driver.py"
VERIFY_REL="scripts/43_phase9_corrected_verify.py"
# Decision 13: must equal orchestration.REQUIRED_FILES (checked by the test suite)
REQUIRED_REL=("$DRIVER_REL" "$SELF_REL" "$VERIFY_REL" "scripts/44_phase9_freeze_prepare.sh"
              "scripts/44_phase9_freeze_capture.py" "scripts/40_phase9_dc_interval.py"
              "src/__init__.py" "src/paths.py" "src/reproducibility.py" "src/phase9_execution.py"
              "src/pathway_scoring.py" "src/score_standardisation.py")
FREEZE_REL="config/analysis/PHASE9_CORRECTED_FREEZE_RECORD.json"
PKG_DIR="src/phase9_corrected"
D13='.owner_decisions_outside_schema.D13_driver_freeze'

# ---- arguments (the config path is resolved against the CALLER's working directory)
[ "$#" -ge 2 ] || die 2 "usage: $0 <validate|run|verify> <owner-config.json> [driver args...]"
CMD="$1"
CONFIG_ARG="$2"
shift 2
case "$CMD" in
  validate|run) ;;
  verify) [ "$#" -eq 0 ] || die 2 "verify takes no further arguments" ;;
  *) die 2 "command must be 'validate', 'run' or 'verify'" ;;
esac
[ -f "$CONFIG_ARG" ] || die 11 "owner configuration not found: $CONFIG_ARG"
CONFIG="$(cd "$(dirname "$CONFIG_ARG")" && pwd -P)/$(basename "$CONFIG_ARG")" \
  || die 11 "cannot resolve the owner configuration path"

# ---- repository root, independent of the caller's working directory
SELF_DIR="$(cd "$(dirname "$0")" && pwd -P)" || die 10 "cannot resolve the launcher directory"
ROOT="$(cd "$SELF_DIR/.." && pwd -P)" || die 10 "cannot resolve the repository root"
[ "$SELF_DIR/$(basename "$0")" = "$ROOT/$SELF_REL" ] || die 10 "launcher is not at $SELF_REL"
cd "$ROOT"

WORK="$(/usr/bin/mktemp -d "${TMPDIR:-/tmp}/p9c_launch.XXXXXX")" || die 19 "mktemp failed"

# ---- 1. owner-approved anchor (read with /usr/bin/jq; no project Python)
[ -f "$FREEZE_REL" ] || die 11 "freeze record absent: $FREEZE_REL"
"$JQ" -e "$D13.approved == true" "$CONFIG" >/dev/null 2>&1 || die 12 "Decision 13 is not approved"
FRH="$("$JQ" -er "$D13.freeze_record_sha256" "$CONFIG" 2>/dev/null)" \
  || die 12 "owner freeze_record_sha256 is missing"
DRH="$("$JQ" -er "$D13.sha256" "$CONFIG" 2>/dev/null)" || die 12 "owner driver sha256 is missing"
DRP="$("$JQ" -er "$D13.driver_path" "$CONFIG" 2>/dev/null)" || die 12 "owner driver_path is missing"
[[ "$FRH" =~ ^[0-9a-f]{64}$ ]] || die 13 "freeze_record_sha256 is not 64 lowercase hex characters"
[[ "$DRH" =~ ^[0-9a-f]{64}$ ]] || die 13 "driver sha256 is not 64 lowercase hex characters"
[ "$DRP" = "$DRIVER_REL" ] || die 13 "driver_path must be $DRIVER_REL"

# ---- 2. freeze-record bytes against the anchor, BEFORE its inventory is used
ACT="$("$SHASUM" -a 256 "$FREEZE_REL")" || die 14 "cannot hash the freeze record"
ACT="${ACT%% *}"
[ "$ACT" = "$FRH" ] || die 14 "freeze record does not match the owner-approved freeze_record_sha256"

# ---- 3. inventory shape: relative, safe paths; 64-hex values; no configuration files
"$JQ" -e '.code | type == "object" and length > 0' "$FREEZE_REL" >/dev/null 2>&1 \
  || die 15 "freeze record has no code inventory"
"$JQ" -e '[.code | to_entries[]
           | select((.key | test("^[A-Za-z0-9_./-]+$") | not) or (.key | startswith("/"))
                    or (.key | contains("..")) or (.key | startswith("config/"))
                    or ((.value | type) != "string") or (.value | test("^[0-9a-f]{64}$") | not))]
          | length == 0' "$FREEZE_REL" >/dev/null 2>&1 \
  || die 15 "freeze inventory contains an unsafe path, a configuration file or a malformed hash"
"$JQ" -er --arg p "$SELF_REL" '.code[$p]' "$FREEZE_REL" >/dev/null 2>&1 \
  || die 15 "the launcher is not in the freeze inventory"
INV_DRV="$("$JQ" -er --arg p "$DRIVER_REL" '.code[$p]' "$FREEZE_REL" 2>/dev/null)" \
  || die 15 "the driver is not in the freeze inventory"
[ "$INV_DRV" = "$DRH" ] || die 17 "owner driver sha256 differs from the freeze inventory's driver entry"
for rel in "${REQUIRED_REL[@]}"; do
  "$JQ" -er --arg p "$rel" '.code[$p]' "$FREEZE_REL" >/dev/null 2>&1 \
    || die 15 "required file is not in the freeze inventory: $rel"
done

# ---- 4. every frozen file (including this launcher): changed or MISSING -> refuse
"$JQ" -er '.code | to_entries[] | "\(.value)  \(.key)"' "$FREEZE_REL" > "$WORK/inventory.sha256" \
  || die 15 "cannot extract the freeze inventory"
if ! "$SHASUM" -a 256 --strict -c "$WORK/inventory.sha256" > "$WORK/check.txt" 2>&1; then
  /usr/bin/grep -v ': OK$' "$WORK/check.txt" >&2 || true
  die 16 "frozen inventory mismatch (a file changed or is missing)"
fi

# ---- 5. package file set: no unexpected and no missing module (explicit comparison)
/usr/bin/find "$PKG_DIR" -maxdepth 1 -type f -name '*.py' > "$WORK/pkg_disk.unsorted" \
  || die 18 "cannot list $PKG_DIR"
/usr/bin/sort "$WORK/pkg_disk.unsorted" > "$WORK/pkg_disk.txt"
"$JQ" -r --arg d "$PKG_DIR/" '.code | keys[] | select(startswith($d))' "$FREEZE_REL" \
  > "$WORK/pkg_frozen.unsorted" || die 18 "cannot list the frozen package"
/usr/bin/sort "$WORK/pkg_frozen.unsorted" > "$WORK/pkg_frozen.txt"
set +e
/usr/bin/cmp -s "$WORK/pkg_disk.txt" "$WORK/pkg_frozen.txt"
cmp_status=$?
set -e
if [ "$cmp_status" -ne 0 ]; then
  /usr/bin/diff "$WORK/pkg_disk.txt" "$WORK/pkg_frozen.txt" >&2 || true
  die 18 "package file set differs from the freeze (unexpected or missing module)"
fi

# ---- 6. fresh, empty, external bytecode directory; interpreter
PYC="$(/usr/bin/mktemp -d "${TMPDIR:-/tmp}/p9c_pyc.XXXXXX")" || die 19 "mktemp failed"
PYC="$(cd "$PYC" && pwd -P)" || die 19 "cannot resolve the bytecode directory"
case "$PYC/" in
  "$ROOT"/*) die 19 "bytecode directory must be outside the repository" ;;
esac
[ -z "$(/bin/ls -A "$PYC")" ] || die 20 "bytecode directory is not empty"
PY="$ROOT/.venv/bin/python"
[ -x "$PY" ] || die 21 "virtual-environment interpreter not found: $PY"

# ---- 7. launch; preserve Python's exit status; cleanup never masks it
set +e
if [ "$CMD" = "verify" ]; then
  "$PY" -I -B -X pycache_prefix="$PYC" "$ROOT/$VERIFY_REL" "$CONFIG"
else
  "$PY" -I -B -X pycache_prefix="$PYC" "$ROOT/$DRIVER_REL" "$CMD" "$CONFIG" ${1+"$@"}
fi
status=$?
set -e
if [ -n "$(/bin/ls -A "$PYC")" ]; then
  printf 'p9c-launch: bytecode was written despite -B (policy violation)\n' >&2
  if [ "$status" -eq 0 ]; then status=22; fi
fi
exit "$status"

"""Canonical project paths.

Single source of truth, anchored to the repository root via __file__ so that
every path resolves identically regardless of the current working directory.

Scripts 01-12 are unchanged and keep their own path handling; this module is
for new code only.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# --- raw inputs (read-only; never written by this project) ---
RAW_DIR = ROOT / "data" / "raw"
SOFT_SERIES = RAW_DIR / "gse39582" / "GSE39582_family.soft.gz"
GPL570_ANNOTATION = RAW_DIR / "gpl570" / "GPL570_platform.soft"
RAW_CHECKSUMS = RAW_DIR / "CHECKSUMS.sha256"

# --- curated inputs produced by scripts 01-02 ---
METADATA_DIR = ROOT / "metadata"
SAMPLE_METADATA = METADATA_DIR / "GSE39582_sample_metadata.csv"
ACQUISITION_INFO = METADATA_DIR / "GSE39582_acquisition_info.json"

# --- generated outputs (the only directory new code writes to) ---
PROCESSED_DIR = ROOT / "data" / "processed"
PROBE_MATRIX = PROCESSED_DIR / "probe_matrix.parquet"
PROBE_MATRIX_MANIFEST = PROCESSED_DIR / "probe_matrix_manifest.json"

# --- existing audit outputs (read-only for new code) ---
REPORTS_DIR = ROOT / "reports"
RESULTS_DIR = ROOT / "results"

# --- project scaffolding ---
SCRIPTS_DIR = ROOT / "scripts"
CONFIG_FILE = ROOT / "config" / "analysis_config.yaml"
TESTS_DIR = ROOT / "tests"
BASELINE_CHECKSUMS = TESTS_DIR / "baseline_checksums.json"

# Platform constant for GPL570 (Affymetrix HG-U133 Plus 2.0).
# This is a property of the array design, not an analysis choice.
GPL570_PROBE_COUNT = 54675

# Sample count of the full GSE39582 series as deposited. The full series is
# cached; restricting to a subset is a cohort-eligibility decision reserved
# for Phase 1.
SERIES_SAMPLE_COUNT = 585


def ensure_processed_dir() -> Path:
    """Create the generated-output directory if absent. Returns the path."""
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    return PROCESSED_DIR

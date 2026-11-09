"""Skip tests that need the local data bundle when it is absent, for example in a fresh clone.

The raw GEO files (data/raw, 356 MB) and the derived caches (data/processed) are not in Git.
The tests listed in data_dependent_tests.txt read them, or launch scripts that expect the
project virtual environment. The list was built by running the whole suite in a clean clone.
When data/processed exists, nothing is skipped and the full suite runs unchanged.
"""
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA_PRESENT = (ROOT / "data" / "processed").is_dir()
_LISTED = {
    line.strip()
    for line in Path(__file__).with_name("data_dependent_tests.txt").read_text().splitlines()
    if line.strip() and not line.startswith("#")
}


def pytest_collection_modifyitems(config, items):
    if DATA_PRESENT:
        return
    skip = pytest.mark.skip(reason="needs the local data bundle (data/raw, data/processed); see README, 'Reproduce the analysis'")
    for item in items:
        if item.nodeid in _LISTED:
            item.add_marker(skip)

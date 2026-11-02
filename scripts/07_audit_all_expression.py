import gzip
import numpy as np
import pandas as pd
from pathlib import Path

soft_file = Path("data/raw/gse39582/GSE39582_family.soft.gz")
metadata_file = Path("metadata/GSE39582_sample_metadata.csv")

# Load metadata and identify primary tumor samples
metadata = pd.read_csv(metadata_file)

tumor_ids = set(
    metadata.loc[
        metadata["source_name"].str.contains(
            "primary colorectal Adenocarcinoma",
            case=False,
            na=False
        ),
        "sample_id"
    ]
)

expected_samples = len(tumor_ids)

print("Expected primary tumor samples:", expected_samples)

# Track expression structure
sample_count = 0
probe_reference = None
probe_count_by_sample = []
missing_values_by_sample = []
non_numeric_by_sample = []
sample_ranges = []

inside_sample = False
inside_table = False
current_sample = None
current_probe_count = 0
current_missing = 0
current_non_numeric = 0
current_min = float("inf")
current_max = float("-inf")
current_probe_ids = []

def finish_sample():
    global sample_count

    if current_sample not in tumor_ids:
        return

    sample_count += 1
    probe_count_by_sample.append((current_sample, current_probe_count))
    missing_values_by_sample.append((current_sample, current_missing))
    non_numeric_by_sample.append((current_sample, current_non_numeric))

    sample_ranges.append(
        (current_sample, current_min, current_max)
    )

with gzip.open(
    soft_file, "rt", encoding="utf-8", errors="replace"
) as f:

    for line in f:

        line = line.rstrip("\n")

        if line.startswith("^SAMPLE = "):
            if inside_sample and inside_table:
                finish_sample()

            current_sample = line.split("=", 1)[1].strip()
            inside_sample = current_sample in tumor_ids
            inside_table = False

            current_probe_count = 0
            current_missing = 0
            current_non_numeric = 0
            current_min = float("inf")
            current_max = float("-inf")
            current_probe_ids = []

        elif inside_sample and line == "!sample_table_begin":
            inside_table = True
            continue

        elif inside_sample and line == "!sample_table_end":
            inside_table = False
            finish_sample()
            inside_sample = False
            current_sample = None

        elif inside_sample and inside_table:

            if current_probe_count == 0 and line.startswith("ID_REF"):
                continue

            parts = line.split("\t")

            if len(parts) < 2:
                current_non_numeric += 1
                continue

            probe_id = parts[0]
            current_probe_ids.append(probe_id)

            try:
                value = float(parts[1])

                if not np.isfinite(value):
                    current_missing += 1
                else:
                    current_min = min(current_min, value)
                    current_max = max(current_max, value)

            except ValueError:
                current_non_numeric += 1

            current_probe_count += 1

        if sample_count % 100 == 0 and sample_count > 0:
            # Progress is printed only at sample boundaries below.
            pass

    if inside_sample and inside_table:
        finish_sample()

print("\n========== SAMPLE COMPLETENESS ==========")
print("Primary tumor samples expected:", expected_samples)
print("Primary tumor expression tables found:", sample_count)

probe_counts = pd.Series(
    dict(probe_count_by_sample),
    name="probe_count"
)

print("\nProbe counts across samples:")
print(probe_counts.describe())

print("\nSamples with unexpected probe counts:")
print(probe_counts[probe_counts != 54675])

missing_counts = pd.Series(dict(missing_values_by_sample))
non_numeric_counts = pd.Series(dict(non_numeric_by_sample))

print("\n========== VALUE QUALITY ==========")
print("Total non-finite values:", int(missing_counts.sum()))
print("Total non-numeric values:", int(non_numeric_counts.sum()))

print("\nSamples with non-finite values:")
print(missing_counts[missing_counts > 0])

print("\nSamples with non-numeric values:")
print(non_numeric_counts[non_numeric_counts > 0])

ranges = pd.DataFrame(
    sample_ranges,
    columns=["sample_id", "minimum", "maximum"]
)

print("\n========== EXPRESSION VALUE RANGE ==========")
print(ranges[["minimum", "maximum"]].describe())

print("\n========== AUDIT STATUS ==========")

if sample_count == expected_samples:
    print("All primary tumor samples have expression tables.")
else:
    print("WARNING: Some primary tumor expression tables are missing.")

if (probe_counts == 54675).all():
    print("All samples have 54,675 probe rows.")
else:
    print("WARNING: Probe counts differ between samples.")

if missing_counts.sum() == 0 and non_numeric_counts.sum() == 0:
    print("All expression values are finite numeric values.")
else:
    print("WARNING: Expression values require further investigation.")

print("\nAudit complete. Original data were not modified.")

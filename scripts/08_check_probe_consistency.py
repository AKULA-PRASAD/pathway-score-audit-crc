import gzip
import pandas as pd
from pathlib import Path

soft_file = Path("data/raw/gse39582/GSE39582_family.soft.gz")
metadata_file = Path("metadata/GSE39582_sample_metadata.csv")

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

reference_probes = None
sample_count = 0
mismatch_samples = []
duplicate_probe_samples = []
missing_probe_samples = []
current_probes = []
current_sample = None
inside_sample = False
inside_table = False

with gzip.open(
    soft_file, "rt", encoding="utf-8", errors="replace"
) as f:

    for line in f:

        line = line.rstrip("\n")

        if line.startswith("^SAMPLE = "):
            current_sample = line.split("=", 1)[1].strip()
            inside_sample = current_sample in tumor_ids
            inside_table = False
            current_probes = []

        elif inside_sample and line == "!sample_table_begin":
            inside_table = True
            continue

        elif inside_sample and line == "!sample_table_end":
            inside_table = False
            sample_count += 1

            if len(current_probes) != 54675:
                missing_probe_samples.append(current_sample)

            if len(current_probes) != len(set(current_probes)):
                duplicate_probe_samples.append(current_sample)

            if reference_probes is None:
                reference_probes = current_probes.copy()
            elif current_probes != reference_probes:
                mismatch_samples.append(current_sample)

            inside_sample = False

        elif inside_sample and inside_table:

            if line.startswith("ID_REF"):
                continue

            parts = line.split("\t")

            if parts:
                current_probes.append(parts[0])

print("\n========== PROBE CONSISTENCY AUDIT ==========")
print("Samples checked:", sample_count)

if reference_probes is not None:
    print("Reference probe count:", len(reference_probes))
    print("Unique reference probes:", len(set(reference_probes)))

print("\nSamples with probe-order or probe-ID mismatches:",
      len(mismatch_samples))

print(mismatch_samples[:20])

print("\nSamples with duplicate probe IDs:",
      len(duplicate_probe_samples))

print(duplicate_probe_samples[:20])

print("\nSamples with unexpected probe counts:",
      len(missing_probe_samples))

print(missing_probe_samples[:20])

print("\n========== FINAL STATUS ==========")

if (
    sample_count == len(tumor_ids)
    and len(mismatch_samples) == 0
    and len(duplicate_probe_samples) == 0
    and len(missing_probe_samples) == 0
):
    print("PASS: All primary tumor samples have identical probe IDs in identical order.")
else:
    print("WARNING: Probe consistency issues require investigation.")

print("\nAudit complete. Original data were not modified.")


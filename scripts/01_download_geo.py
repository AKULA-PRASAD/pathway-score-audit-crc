from pathlib import Path
import GEOparse
import pandas as pd
import json
from datetime import datetime, timezone

# Project directories
PROJECT = Path(__file__).resolve().parents[1]

RAW_DIR = PROJECT / "data" / "raw" / "gse39582"
METADATA_DIR = PROJECT / "metadata"

RAW_DIR.mkdir(parents=True, exist_ok=True)
METADATA_DIR.mkdir(parents=True, exist_ok=True)

# Download GEO series
print("Retrieving GSE39582 metadata...")

gse = GEOparse.get_GEO(
    geo="GSE39582",
    destdir=str(RAW_DIR),
    silent=False
)

print("\nGEO accession:", gse.name)
print("Number of samples:", len(gse.gsms))
print("Platform IDs:", list(gse.gpls.keys()))

# Extract sample-level metadata
records = []

for sample_id, sample in gse.gsms.items():

    record = {
        "sample_id": sample_id,
        "title": sample.metadata.get("title", [None])[0],
        "platform_id": sample.metadata.get("platform_id", [None])[0],
        "source_name": sample.metadata.get("source_name_ch1", [None])[0],
    }

    # Preserve all sample characteristics
    for characteristic in sample.metadata.get(
        "characteristics_ch1", []
    ):
        if ":" in characteristic:
            key, value = characteristic.split(":", 1)
            key = key.strip().lower().replace(" ", "_")

            record[key] = value.strip()

    records.append(record)

sample_metadata = pd.DataFrame(records)

# Save metadata
output_file = METADATA_DIR / "GSE39582_sample_metadata.csv"

sample_metadata.to_csv(output_file, index=False)

# Save acquisition information
acquisition_info = {
    "accession": gse.name,
    "download_date_utc": datetime.now(timezone.utc).isoformat(),
    "number_of_samples": len(gse.gsms),
    "platform_ids": list(gse.gpls.keys()),
    "metadata_columns": sample_metadata.columns.tolist(),
}

with open(
    METADATA_DIR / "GSE39582_acquisition_info.json",
    "w"
) as file:
    json.dump(acquisition_info, file, indent=4)

print("\nMetadata dimensions:", sample_metadata.shape)
print("\nMetadata columns:")
print(sample_metadata.columns.tolist())

print("\nFirst five records:")
print(sample_metadata.head().to_string(index=False))

print("\nSaved metadata to:", output_file)

import gzip
from pathlib import Path

soft_file = Path(
    "data/raw/gse39582/GSE39582_family.soft.gz"
)

target_sample = "GSM971957"

if not soft_file.exists():
    raise FileNotFoundError(f"File not found: {soft_file}")

sample_found = False
inside_target = False
inside_table = False

metadata_lines = []
table_header = None
first_rows = []
row_count = 0

with gzip.open(soft_file, "rt", encoding="utf-8", errors="replace") as f:

    for line in f:

        line = line.rstrip("\n")

        if line.startswith("^SAMPLE = "):
            inside_target = line.strip() == f"^SAMPLE = {target_sample}"

            if inside_target:
                sample_found = True

        elif inside_target:

            if line == "!sample_table_begin":
                inside_table = True
                continue

            if line == "!sample_table_end":
                inside_table = False
                break

            if inside_table:
                if table_header is None:
                    table_header = line
                else:
                    row_count += 1
                    if len(first_rows) < 10:
                        first_rows.append(line)

            elif line.startswith("!"):
                metadata_lines.append(line)

if not sample_found:
    raise ValueError(f"Sample {target_sample} was not found.")

print("\n========== SAMPLE METADATA ==========")

for line in metadata_lines:
    if any(term in line.lower() for term in [
        "title",
        "source_name",
        "characteristics",
        "description",
        "data_processing",
        "platform_id",
        "supplementary_file"
    ]):
        print(line)

print("\n========== EXPRESSION TABLE ==========")
print("Sample:", target_sample)
print("Column header:", table_header)
print("Number of probe rows:", row_count)

print("\nFirst 10 expression records:")
for row in first_rows:
    print(row)

print("\nInspection complete. No original data were modified.")

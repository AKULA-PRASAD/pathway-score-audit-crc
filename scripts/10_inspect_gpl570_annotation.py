import pandas as pd
from pathlib import Path

file_path = Path("data/raw/gpl570/GPL570_platform.soft")

if not file_path.exists():
    raise FileNotFoundError(file_path)

# Locate the beginning of the platform annotation table
header_index = None

with open(file_path, "r", encoding="utf-8", errors="replace") as f:
    for i, line in enumerate(f):
        if line.startswith("ID\t"):
            header_index = i
            break

if header_index is None:
    raise ValueError("Could not locate the GPL570 annotation table.")

print("Annotation table begins at line:", header_index)

# Read annotation table
annotation = pd.read_csv(
    file_path,
    sep="\t",
    skiprows=header_index,
    dtype=str,
    low_memory=False,
    on_bad_lines="warn"
)

# Remove any trailing non-table records
if "ID" not in annotation.columns:
    raise ValueError("Annotation table does not contain an ID column.")

annotation = annotation[
    annotation["ID"].notna()
    & annotation["ID"].str.strip().ne("")
].copy()

print("\n========== ANNOTATION STRUCTURE ==========")
print("Rows:", len(annotation))
print("Columns:", len(annotation.columns))
print("Unique probe IDs:", annotation["ID"].nunique())
print("Duplicate probe IDs:", annotation["ID"].duplicated().sum())

print("\n========== AVAILABLE ANNOTATION COLUMNS ==========")

for column in annotation.columns:
    print("-", column)

# Identify likely gene annotation columns
gene_columns = [
    column for column in annotation.columns
    if any(term in column.lower() for term in [
        "gene symbol",
        "gene title",
        "gene assignment",
        "entrez",
        "refseq"
    ])
]

print("\n========== GENE-RELATED COLUMNS ==========")

for column in gene_columns:
    print(f"\nColumn: {column}")
    print(annotation[column].head(10).to_string(index=False))

print("\n========== EXAMPLE PROBE RECORDS ==========")

example_probes = [
    "1007_s_at",
    "1053_at",
    "117_at",
    "121_at",
    "1255_g_at"
]

print(
    annotation[annotation["ID"].isin(example_probes)]
    .to_string(index=False)
)

print("\nInspection complete. Original annotation unchanged.")

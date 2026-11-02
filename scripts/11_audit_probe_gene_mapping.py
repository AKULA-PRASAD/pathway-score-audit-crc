import gzip
import pandas as pd
from pathlib import Path
from collections import Counter

soft_file = Path("data/raw/gse39582/GSE39582_family.soft.gz")
annotation_file = Path("data/raw/gpl570/GPL570_platform.soft")

# --------------------------------------------------
# 1. Read GPL570 annotation table
# --------------------------------------------------

with open(annotation_file, "r", encoding="utf-8", errors="replace") as f:
    lines = f.readlines()

header_index = next(
    (i for i, line in enumerate(lines) if line.startswith("ID\t")),
    None
)

if header_index is None:
    raise ValueError("Annotation table header not found.")

annotation = pd.read_csv(
    annotation_file,
    sep="\t",
    skiprows=header_index,
    dtype=str,
    low_memory=False
)

annotation = annotation.dropna(subset=["ID"]).copy()

# --------------------------------------------------
# 2. Extract expression probe IDs from one sample
# --------------------------------------------------

target_sample = "GSM971957"
expression_probes = []
inside_target = False
inside_table = False

with gzip.open(soft_file, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        line = line.rstrip("\n")

        if line.startswith("^SAMPLE = "):
            inside_target = (
                line.split("=", 1)[1].strip() == target_sample
            )
            inside_table = False

        elif inside_target and line == "!sample_table_begin":
            inside_table = True
            continue

        elif inside_target and line == "!sample_table_end":
            break

        elif inside_target and inside_table:
            if line.startswith("ID_REF"):
                continue

            parts = line.split("\t")
            if parts:
                expression_probes.append(parts[0])

# --------------------------------------------------
# 3. Compare expression and annotation probe IDs
# --------------------------------------------------

expression_set = set(expression_probes)
annotation_set = set(annotation["ID"])

print("\n========== PROBE ID COMPARISON ==========")
print("Expression probes:", len(expression_set))
print("Annotation probes:", len(annotation_set))

missing_annotation = expression_set - annotation_set
extra_annotation = annotation_set - expression_set

print("Expression probes missing from annotation:", len(missing_annotation))
print(sorted(missing_annotation)[:30])

print("Annotation probes absent from expression:", len(extra_annotation))
print(sorted(extra_annotation)[:30])

# --------------------------------------------------
# 4. Classify probe-to-gene mappings
# --------------------------------------------------

def parse_symbols(value):
    if pd.isna(value):
        return []

    symbols = [
        x.strip()
        for x in str(value).split("///")
    ]

    # Remove empty and placeholder annotations
    symbols = [
        x for x in symbols
        if x and x not in {"---", "NA", "N/A"}
    ]

    return list(dict.fromkeys(symbols))


annotation["parsed_symbols"] = annotation["Gene Symbol"].apply(
    parse_symbols
)

unmapped = annotation[
    annotation["parsed_symbols"].apply(len) == 0
]

single_mapped = annotation[
    annotation["parsed_symbols"].apply(len) == 1
]

multi_mapped = annotation[
    annotation["parsed_symbols"].apply(len) > 1
]

print("\n========== PROBE MAPPING CATEGORIES ==========")
print("Unmapped probes:", len(unmapped))
print("Single-gene probes:", len(single_mapped))
print("Multi-gene probes:", len(multi_mapped))

print("\nExamples of multi-gene probes:")
print(
    multi_mapped[["ID", "Gene Symbol"]]
    .head(15)
    .to_string(index=False)
)

# --------------------------------------------------
# 5. Count genes represented by multiple probes
# --------------------------------------------------

gene_probe_counts = Counter()

for symbols in single_mapped["parsed_symbols"]:
    gene_probe_counts[symbols[0]] += 1

genes_with_multiple_probes = {
    gene: count
    for gene, count in gene_probe_counts.items()
    if count > 1
}

print("\n========== GENE-LEVEL MAPPING ==========")
print("Unique unambiguous genes:", len(gene_probe_counts))
print("Genes represented by multiple probes:",
      len(genes_with_multiple_probes))

print("\nGenes with the most probes:")
for gene, count in sorted(
    genes_with_multiple_probes.items(),
    key=lambda x: x[1],
    reverse=True
)[:20]:
    print(gene, count)

print("\n========== SUMMARY ==========")
print("Expression probes:", len(expression_set))
print("Annotation probes:", len(annotation_set))
print("Unambiguous gene-mapped probes:", len(single_mapped))
print("Ambiguous probes:", len(multi_mapped))
print("Unmapped probes:", len(unmapped))
print("Unique unambiguous genes:", len(gene_probe_counts))

print("\nAudit complete. No original data were modified.")

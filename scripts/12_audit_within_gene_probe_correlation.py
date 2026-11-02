
from pathlib import Path
from itertools import combinations
import gzip
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# --------------------------------------------------
# 1. PROJECT PATHS
# --------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]

SOFT_FILE = ROOT / "data/raw/gse39582/GSE39582_family.soft.gz"
ANNOTATION_FILE = ROOT / "data/raw/gpl570/GPL570_platform.soft"

METADATA_FILE = ROOT / "metadata/GSE39582_sample_metadata.csv"

OUTPUT_DIR = ROOT / "results/probe_correlation_audit"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------
# 2. LOAD SAMPLE METADATA
# --------------------------------------------------

metadata = pd.read_csv(METADATA_FILE)

primary = metadata[
    metadata["source_name"].str.contains(
        "primary colorectal adenocarcinoma",
        case=False,
        na=False
    )
].copy()

sample_ids = primary["sample_id"].tolist()

print("Primary tumor samples:", len(sample_ids))

# --------------------------------------------------
# 3. LOAD GPL570 ANNOTATION
# --------------------------------------------------

with open(ANNOTATION_FILE, "r", encoding="utf-8", errors="replace") as f:
    lines = f.readlines()

header_index = next(
    i for i, line in enumerate(lines)
    if line.startswith("ID\tGB_ACC")
)

annotation = pd.read_csv(
    ANNOTATION_FILE,
    sep="\t",
    skiprows=header_index,
    dtype=str,
    comment=None
)

annotation = annotation[
    annotation["ID"].notna()
    & ~annotation["ID"].str.startswith("!")
].copy()

# --------------------------------------------------
# 4. IDENTIFY UNAMBIGUOUS PROBE-GENE MAPPINGS
# --------------------------------------------------

def parse_gene_symbol(value):

    if pd.isna(value):
        return None

    symbols = [
        x.strip()
        for x in str(value).split("///")
        if x.strip() not in {"", "---", "NA", "N/A"}
    ]

    symbols = list(dict.fromkeys(symbols))

    if len(symbols) == 1:
        return symbols[0]

    return None


annotation["gene_symbol_clean"] = (
    annotation["Gene Symbol"].apply(parse_gene_symbol)
)

mapping = annotation[
    annotation["gene_symbol_clean"].notna()
][["ID", "gene_symbol_clean"]].copy()

print("Unambiguous probes:", len(mapping))
print(
    "Unique genes:",
    mapping["gene_symbol_clean"].nunique()
)

# --------------------------------------------------
# 5. EXTRACT EXPRESSION VALUES FROM SOFT FILE
# --------------------------------------------------

probe_to_gene = dict(
    zip(mapping["ID"], mapping["gene_symbol_clean"])
)

sample_set = set(sample_ids)

expression_by_sample = {}

current_sample = None
current_rows = []
inside_table = False

with gzip.open(
    SOFT_FILE,
    "rt",
    encoding="utf-8",
    errors="replace"
) as f:

    for line in f:

        if line.startswith("^SAMPLE = "):

            if current_sample in sample_set and current_rows:
                expression_by_sample[current_sample] = dict(current_rows)

            current_sample = line.strip().split("=", 1)[1].strip()

            current_rows = []
            inside_table = False

        elif (
            line.startswith("ID_REF\t")
            and current_sample in sample_set
        ):

            inside_table = True

        elif inside_table:

            if line.startswith("!sample_table_end"):
                inside_table = False
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 2:
                continue

            probe_id = parts[0]

            if probe_id not in probe_to_gene:
                continue

            try:
                value = float(parts[1])
            except ValueError:
                continue

            if np.isfinite(value):
                current_rows.append((probe_id, value))

    if current_sample in sample_set and current_rows:
        expression_by_sample[current_sample] = dict(current_rows)

print("Samples extracted:", len(expression_by_sample))

if len(expression_by_sample) != len(sample_ids):
    raise ValueError(
        "Not all primary tumor samples were extracted. "
        "Check the SOFT parsing logic before proceeding."
    )

# --------------------------------------------------
# 6. BUILD PROBE-LEVEL EXPRESSION MATRIX
# --------------------------------------------------

expression = pd.DataFrame(expression_by_sample)

expression.index.name = "probe_id"

print("Probe-level matrix:", expression.shape)

# --------------------------------------------------
# 7. GROUP PROBES BY GENE
# --------------------------------------------------

gene_groups = mapping.groupby("gene_symbol_clean")["ID"]

multi_probe_genes = gene_groups.filter(
    lambda x: len(x) > 1
)

genes_to_test = {
    gene: list(probes)
    for gene, probes in gene_groups
    if len(probes) > 1
}

print("Genes with multiple probes:", len(genes_to_test))

# --------------------------------------------------
# 8. CALCULATE WITHIN-GENE SPEARMAN CORRELATIONS
# --------------------------------------------------

correlation_results = []

for gene, probes in genes_to_test.items():

    available_probes = [
        probe for probe in probes
        if probe in expression.index
    ]

    if len(available_probes) < 2:
        continue

    for probe_a, probe_b in combinations(available_probes, 2):

        x = expression.loc[probe_a].values
        y = expression.loc[probe_b].values

        result = spearmanr(x, y, nan_policy="omit")

        correlation_results.append({
            "gene_symbol": gene,
            "probe_a": probe_a,
            "probe_b": probe_b,
            "spearman_rho": result.statistic,
            "p_value": result.pvalue
        })

correlations = pd.DataFrame(correlation_results)

if correlations.empty:
    raise ValueError("No within-gene probe correlations were calculated.")

# --------------------------------------------------
# 9. SUMMARIZE CORRELATIONS BY GENE
# --------------------------------------------------

gene_summary = (
    correlations.groupby("gene_symbol")
    .agg(
        number_of_probe_pairs=("spearman_rho", "count"),
        median_spearman_rho=("spearman_rho", "median"),
        minimum_spearman_rho=("spearman_rho", "min"),
        maximum_spearman_rho=("spearman_rho", "max")
    )
    .reset_index()
)

gene_summary = gene_summary.sort_values(
    "median_spearman_rho"
)

# --------------------------------------------------
# 10. SAVE RESULTS
# --------------------------------------------------

correlations.to_csv(
    OUTPUT_DIR / "within_gene_probe_correlations.csv",
    index=False
)

gene_summary.to_csv(
    OUTPUT_DIR / "within_gene_correlation_summary.csv",
    index=False
)

# --------------------------------------------------
# 11. PRINT AUDIT SUMMARY
# --------------------------------------------------

print("\n========== WITHIN-GENE CORRELATION AUDIT ==========")

print("Genes evaluated:", len(gene_summary))
print("Probe pairs evaluated:", len(correlations))

print("\nCorrelation distribution:")

print(
    correlations["spearman_rho"].describe(
        percentiles=[0.10, 0.25, 0.50, 0.75, 0.90]
    )
)

print("\nGene-level median correlation distribution:")

print(
    gene_summary["median_spearman_rho"].describe(
        percentiles=[0.10, 0.25, 0.50, 0.75, 0.90]
    )
)

print("\nLowest-agreement genes:")

print(
    gene_summary.head(20).to_string(index=False)
)

print("\nHighest-agreement genes:")

print(
    gene_summary.sort_values(
        "median_spearman_rho",
        ascending=False
    ).head(20).to_string(index=False)
)

print("\nResults saved to:", OUTPUT_DIR)

print("\nAudit complete. Original data were not modified.")


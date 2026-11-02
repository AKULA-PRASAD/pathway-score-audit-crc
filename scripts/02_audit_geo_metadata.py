from pathlib import Path
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]

INPUT = PROJECT / "metadata" / "GSE39582_sample_metadata.csv"
OUTPUT = PROJECT / "reports"

OUTPUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)

print("\n========== BASIC DATASET INFORMATION ==========")

print("Total samples:", len(df))
print("Total columns:", len(df.columns))
print("Duplicate sample IDs:", df["sample_id"].duplicated().sum())

print("\n========== STUDY SUBSETS ==========")

print(df["dataset"].value_counts(dropna=False))

print("\n========== TUMOR STAGE ==========")

print(df["tnm.stage"].value_counts(dropna=False))

print("\n========== RECURRENCE INFORMATION ==========")

print("RFS event distribution:")
print(df["rfs.event"].value_counts(dropna=False))

print("\nRFS delay distribution:")
print(df["rfs.delay"].describe())

print("\n========== OVERALL SURVIVAL ==========")

print("OS event distribution:")
print(df["os.event"].value_counts(dropna=False))

print("\nOS delay distribution:")
print(df["os.delay_(months)"].describe())

print("\n========== MISSING DATA ==========")

missing = pd.DataFrame({
    "missing_count": df.isna().sum(),
    "missing_percentage": df.isna().mean() * 100
}).sort_values("missing_percentage", ascending=False)

print(missing.to_string())

print("\n========== SAMPLE SOURCE ==========")

print(df["source_name"].value_counts(dropna=False))

print("\n========== CLINICAL VARIABLES ==========")

clinical_columns = [
    "sex",
    "age.at.diagnosis_(year)",
    "tnm.stage",
    "tnm.t",
    "tnm.n",
    "tnm.m",
    "tumor.location",
    "chemotherapy.adjuvant",
    "mmr.status"
]

for col in clinical_columns:
    if col in df.columns:
        print(f"\n{col}")
        print(df[col].value_counts(dropna=False).head(15))

# Save missingness report
missing.to_csv(OUTPUT / "GSE39582_missingness.csv")

# Save stage and endpoint summaries
summary = {
    "total_samples": len(df),
    "unique_sample_ids": df["sample_id"].nunique(),
    "study_subsets": df["dataset"].value_counts(dropna=False).to_dict(),
    "tumor_stages": df["tnm.stage"].value_counts(dropna=False).to_dict(),
    "rfs_events": df["rfs.event"].value_counts(dropna=False).to_dict(),
    "os_events": df["os.event"].value_counts(dropna=False).to_dict()
}

pd.Series(summary).to_json(
    OUTPUT / "GSE39582_clinical_summary.json",
    indent=4
)

print("\nAudit reports saved in:", OUTPUT)

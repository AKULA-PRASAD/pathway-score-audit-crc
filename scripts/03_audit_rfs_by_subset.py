import pandas as pd

# Load the original metadata without modifying it
df = pd.read_csv("metadata/GSE39582_sample_metadata.csv")

# Keep primary colorectal adenocarcinoma samples only
tumor = df[
    df["source_name"].str.contains(
        "primary colorectal Adenocarcinoma",
        case=False,
        na=False
    )
].copy()

# Define stage I-III eligibility
stage_1_3 = tumor[tumor["tnm.stage"].isin([1, 2, 3])].copy()

print("\n========== PRIMARY TUMOR COHORT ==========")
print("Total primary tumor samples:", len(tumor))

print("\nSamples by study subset:")
print(tumor["dataset"].value_counts(dropna=False))

print("\n========== STAGE I-III COHORT ==========")
print("Total stage I-III samples:", len(stage_1_3))

print("\nStage distribution by subset:")
print(pd.crosstab(
    stage_1_3["dataset"],
    stage_1_3["tnm.stage"],
    margins=True
))

print("\n========== RFS AUDIT BY SUBSET ==========")

for subset in ["discovery", "validation"]:
    group = tumor[tumor["dataset"] == subset]

    print(f"\n--- {subset.upper()} ---")
    print("Total primary tumor samples:", len(group))

    print("RFS event counts:")
    print(group["rfs.event"].value_counts(dropna=False).sort_index())

    print("RFS complete cases:", group["rfs.event"].notna().sum()
          , "event and time both present:",
          group[["rfs.event", "rfs.delay"]].notna().all(axis=1).sum())

    print("RFS zero-time records:",
          (group["rfs.delay"] == 0).sum())

    print("\nStage I-III RFS event counts:")
    eligible = group[group["tnm.stage"].isin([1, 2, 3])]

    print("Eligible samples:", len(eligible))
    print(eligible["rfs.event"].value_counts(dropna=False).sort_index())

    print("Eligible complete RFS records:",
          eligible[["rfs.event", "rfs.delay"]].notna().all(axis=1).sum())

    print("Eligible zero-time records:",
          (eligible["rfs.delay"] == 0).sum())

print("\n========== RFS FOLLOW-UP SUMMARY ==========")

print("\nAll primary tumors:")
print(tumor.groupby("dataset")["rfs.delay"].describe())

print("\nStage I-III primary tumors:")
print(stage_1_3.groupby("dataset")["rfs.delay"].describe())

print("\n========== SAMPLE TITLE UNIQUENESS ==========")
print("Unique titles:", tumor["title"].nunique())
print("Total primary tumor samples:", len(tumor))

print("\nAudit complete. No original data were modified.")

import pandas as pd

df = pd.read_csv("metadata/GSE39582_sample_metadata.csv")

# Keep primary colorectal tumor samples, stage I-III
tumor = df[
    df["source_name"].str.contains(
        "primary colorectal Adenocarcinoma",
        case=False,
        na=False
    )
    & df["tnm.stage"].isin([1, 2, 3])
].copy()

# Identify records with zero RFS time
zero = tumor[tumor["rfs.delay"] == 0].copy()

# Identify records with missing RFS information
missing = tumor[
    tumor[["rfs.event", "rfs.delay"]].isna().any(axis=1)
].copy()

columns = [
    "sample_id",
    "title",
    "dataset",
    "tnm.stage",
    "tnm.t",
    "tnm.n",
    "tnm.m",
    "rfs.event",
    "rfs.delay",
    "os.event",
    "os.delay_(months)",
    "chemotherapy.adjuvant",
    "age.at.diagnosis_(year)"
]

print("\n========== ZERO RFS-TIME RECORDS ==========")
print("Number of records:", len(zero))
print(zero[columns].to_string(index=False))

print("\n========== MISSING RFS RECORDS ==========")
print("Number of records:", len(missing))
print(missing[columns].to_string(index=False))

print("\n========== ZERO-TIME RECORDS BY SUBSET ==========")
print(pd.crosstab(
    zero["dataset"],
    zero["rfs.event"],
    margins=True,
    dropna=False
))

print("\n========== ZERO RFS TIME VS OS ==========")
print(zero[[
    "sample_id",
    "dataset",
    "rfs.event",
    "rfs.delay",
    "os.event",
    "os.delay_(months)"
]].to_string(index=False))

print("\nAudit complete. No original data were modified.")

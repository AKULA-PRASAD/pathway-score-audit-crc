import pandas as pd

df = pd.read_csv("metadata/GSE39582_sample_metadata.csv")

# Primary colorectal adenocarcinoma samples, stage II-III
cohort = df[
    df["source_name"].str.contains(
        "primary colorectal Adenocarcinoma",
        case=False,
        na=False
    )
    & df["tnm.stage"].isin([2, 3])
].copy()

print("\n========== STAGE II-III COHORT ==========")
print("Total samples:", len(cohort))

for subset in ["discovery", "validation"]:
    group = cohort[cohort["dataset"] == subset].copy()

    print(f"\n========== {subset.upper()} ==========")
    print("Total samples:", len(group))

    print("\nRFS event distribution:")
    print(group["rfs.event"].value_counts(dropna=False))

    print("\nComplete RFS records:")
    print(
        group[["rfs.event", "rfs.delay"]]
        .notna()
        .all(axis=1)
        .sum()
    )

    print("\nRFS zero-time records:")
    print((group["rfs.delay"] == 0).sum())

    print("\nRFS event/time combinations:")
    print(pd.crosstab(
        group["rfs.event"],
        group["rfs.delay"].eq(0),
        dropna=False
    ))

    print("\nFollow-up summary:")
    print(group["rfs.delay"].describe())

print("\n========== STAGE II-III OVERALL ==========")
print("RFS events:", (cohort["rfs.event"] == 1).sum())
print("RFS censored:", (cohort["rfs.event"] == 0).sum())
print("Missing RFS event:", cohort["rfs.event"].isna().sum())
print("Missing RFS time:", cohort["rfs.delay"].isna().sum())

print("\nAudit complete. Original metadata unchanged.")

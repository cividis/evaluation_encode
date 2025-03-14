import pandas as pd


def prepare_metadata(input_file="resources/metadata.tsv"):
    # Read metadata and filter required experiments
    df = (
        pd.read_csv(input_file, sep="\t")
        .query("`File format` == 'bed narrowPeak'")
        .query(
            "(`Assay` in ['ATAC-seq', 'TF ChIP-seq']) or (`Assay` == 'Histone ChIP-seq' and `Experiment target`.str.startswith('H3K27ac'))"
        )
    )

    # Process metadata
    df.set_index("File accession", inplace=True)
    df.index.name = None
    df.sort_values("Assay", inplace=True)
    return df


def all_human_transcription_factors(filepath="resources/mmc2.xlsx"):
    # The Human Transcription Factors - Lambert et al. 2018
    # https://www.cell.com/cell/fulltext/S0092-8674(18)30106-5
    df = pd.read_excel(filepath, sheet_name="Table S1. Related to Figure 1B")
    human_tfs = df[df["Is TF?"].str.strip().str.lower() == "yes"]["Unnamed: 1"]
    return human_tfs


# Load metadata
df = prepare_metadata()

# Step 1: Find valid TFs (human TFs with >= 5 cell types)
valid_tfs = (
    df[df["Assay"] == "TF ChIP-seq"]
    .groupby("Experiment target")["Biosample term name"]
    .nunique()
    .reset_index(name="unique_cell_types")
    .query("unique_cell_types >= 5")["Experiment target"]
)

human_tfs = all_human_transcription_factors()
valid_tfs = valid_tfs[valid_tfs.str.split("-").str[0].isin(human_tfs)]

# Step 2: Filter the DataFrame for valid TFs and relevant assays
filtered_df = df[
    ((df["Assay"] == "TF ChIP-seq") & (df["Experiment target"].isin(valid_tfs)))
    | (df["Assay"].isin(["Histone ChIP-seq", "ATAC-seq"]))
]

# Step 3: Identify valid cell types
valid_cell_types = (
    filtered_df.groupby("Biosample term name")["Assay"]
    .apply(set)
    .reset_index(name="assays")
    .loc[
        lambda df: df["assays"].apply(
            lambda x: "TF ChIP-seq" in x
            and ("Histone ChIP-seq" in x or "ATAC-seq" in x)
        ),
        "Biosample term name",
    ]
)

# Step 4: Filter `filtered_df` for only valid cell types
filtered_df = filtered_df[filtered_df["Biosample term name"].isin(valid_cell_types)]

# Step 7: Separate `query_df` and `target_df`
query_df = filtered_df[filtered_df["Assay"].isin(["Histone ChIP-seq", "ATAC-seq"])]
target_df = filtered_df[filtered_df["Assay"] == "TF ChIP-seq"]

# Get relevant columns
query_df = query_df[["Biosample term name", "Assay"]].reset_index()
target_df = target_df[
    ["Biosample term name", "Experiment target", "Assay"]
].reset_index()

# Step 8: Merge on `Biosample term name`
merged_df = query_df.merge(
    target_df, on="Biosample term name", suffixes=("_query", "_target")
)

# Step 9: Save the result
merged_df.to_csv("resources/benchmark_set.csv", index=False)

# Print summary
print(f"Query samples: {query_df.shape[0]}")
print(f"Target samples: {target_df.shape[0]}")
print(f"Merged samples: {merged_df.shape[0]}")

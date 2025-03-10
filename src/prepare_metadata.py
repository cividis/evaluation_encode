import os
import pandas as pd

df = (
    pd.read_csv("metadata/metadata.tsv", sep="\t")
    .query("`File format` == 'bed narrowPeak'")
    .query(
        "(`Assay` in ['ATAC-seq', 'TF ChIP-seq']) or (`Assay` == 'Histone ChIP-seq' and `Experiment target`.str.startswith('H3K27ac'))",
    )
)
df.set_index("File accession", inplace=True)
df.sort_values("Assay", inplace=True)

os.makedirs("data", exist_ok=True)
df.to_parquet("data/metadata.parquet")

import pandas as pd

df = pd.read_csv("metadata/metadata.tsv", sep="\t")
df.query("`File format` == 'bed narrowPeak'", inplace=True)
df.query(
    "(`Assay` in ['ATAC-seq', 'TF ChIP-seq']) or (`Assay` == 'Histone ChIP-seq' and `Experiment target`.str.startswith('H3K27ac'))",
    inplace=True,
)
experiments = df["File accession"].tolist()
print(experiments)

print(df["Assay"].value_counts())

import pandas as pd
from tfsage.download import download_encode
from tfsage.features import extract_features_parallel, load_region_set

df = (
    pd.read_csv("metadata/metadata.tsv", sep="\t")
    .query("`File format` == 'bed narrowPeak'")
    .query(
        "(`Assay` in ['ATAC-seq', 'TF ChIP-seq']) or (`Assay` == 'Histone ChIP-seq' and `Experiment target`.str.startswith('H3K27ac'))",
    )
)
experiments = df["File accession"].tolist()


rule all:
    input:
        "data/rp_matrix.parquet",


rule download_experiment:
    output:
        "work/experiments/{experiment}.bed",
    run:
        download_encode(wildcards.experiment, output[0])


rule extract_features:
    input:
        expand("work/experiments/{experiment}.bed", experiment=experiments),
    output:
        "data/rp_matrix.parquet",
    run:
        from pathlib import Path

        gene_loc_set = load_region_set("hg38")

        df = extract_features_parallel(input, gene_loc_set)
        df.columns = [f"{Path(f).stem}" for f in input]
        df.to_parquet(output[0])

import pandas as pd
from tfsage.download import download_encode
from tfsage.features import extract_features_parallel, load_region_set

df = (
    pd.read_csv("metadata/metadata.tsv", sep="\t")
    .query("`File format` == 'bed narrowPeak'")
    .query(
        "(`Assay` in ['ATAC-seq', 'TF ChIP-seq']) or (`Assay` == 'Histone ChIP-seq' and `Experiment target`.str.startswith('H3K27ac'))",
    )
    .query("`Biosample treatments`.isnull()")
    .query("`Biosample genetic modifications methods`.isnull()")
)
experiments = df["File accession"].tolist()


rule all:
    input:
        # "data/rp_matrix.parquet",
        expand("tmp/{experiment}.done", experiment=experiments),


rule download_experiment:
    output:
        "downloads/{experiment}.bed",
    run:
        download_encode(wildcards.experiment, output[0])


rule extract_features_test:
    input:
        "downloads/{experiment}.bed",
    output:
        "tmp/{experiment}.done",
    run:
        from pathlib import Path
        from tfsage.features import extract_features

        gene_loc_set = load_region_set("hg38")
        df = extract_features(input[0], gene_loc_set)
        Path(output[0]).touch()


rule extract_features:
    input:
        expand("downloads/{experiment}.bed", experiment=experiments),
    output:
        "data/rp_matrix.parquet",
    threads: 16
    run:
        from pathlib import Path

        gene_loc_set = load_region_set("hg38")

        df = extract_features_parallel(input, gene_loc_set, max_workers=threads)
        df.columns = [f"{Path(f).stem}" for f in input]
        df.to_parquet(output[0])

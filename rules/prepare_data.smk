import pandas as pd
from tfsage.download import download_encode
from tfsage.features import extract_features_parallel, load_region_set


rule download_experiment:
    output:
        "downloads/{experiment}.bed",
    run:
        download_encode(wildcards.experiment, output[0])


checkpoint prepare_metadata:
    params:
        input_file="metadata/metadata.tsv",
    output:
        "data/metadata.parquet",
    run:
        df = (
            pd.read_csv(params.input_file, sep="\t")
            .query("`File format` == 'bed narrowPeak'")
            .query(
                "(`Assay` in ['ATAC-seq', 'TF ChIP-seq']) or (`Assay` == 'Histone ChIP-seq' and `Experiment target`.str.startswith('H3K27ac'))",
            )
        )
        df.set_index("File accession", inplace=True)
        df.sort_values("Assay", inplace=True)
        df.to_parquet(output[0])


def aggregate_input(wildcards):
    df = pd.read_parquet(checkpoints.prepare_metadata.get().output[0])
    experiments = df.index.tolist()
    return expand("downloads/{experiment}.bed", experiment=experiments)


rule extract_features:
    input:
        aggregate_input,
    output:
        "data/rp_matrix.parquet",
    threads: 48
    resources:
        mem_mb=96000,
    run:
        from pathlib import Path

        gene_loc_set = load_region_set("hg38")

        df = extract_features_parallel(input, gene_loc_set, max_workers=threads)
        df.columns = [f"{Path(f).stem}" for f in input]
        df.to_parquet(output[0])


rule aggregate_genes:
    input:
        "data/rp_matrix.parquet",
    output:
        "data/rp_matrix_agg.parquet",
    run:
        df = pd.read_parquet(input[0])
        df["gene"] = df.index.str.split(":").str[1]
        df = df.groupby("gene").mean()
        df.index.name = None
        df.to_parquet(output[0])

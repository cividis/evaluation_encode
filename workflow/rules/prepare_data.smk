import os
import pandas as pd
from pathlib import Path
from tfsage.download import download_encode
from tfsage.features import extract_features_parallel, load_region_set
from tfsage.embedding import run_seurat_integration, compute_distances


rule download_experiment:
    output:
        "downloads/{experiment}.bed",
    retries: 5
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
        df.index.name = None
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
        "data/rp_matrix_tss.parquet",
    threads: 48
    resources:
        mem_mb=96000,
    run:
        gene_loc_set = load_region_set("hg38")
        df = extract_features_parallel(input, gene_loc_set, max_workers=threads)
        df.columns = [f"{Path(f).stem}" for f in input]
        df.to_parquet(output[0])


rule aggregate_genes:
    input:
        "data/rp_matrix_tss.parquet",
    output:
        "data/rp_matrix_gene.parquet",
    run:
        df = pd.read_parquet(input[0])
        df["gene"] = df.index.str.split(":").str[1]
        df = df.groupby("gene").mean()
        df.index.name = None
        df.to_parquet(output[0])


rule embed_and_integrate:
    input:
        "data/rp_matrix_{suffix}.parquet",
        "data/metadata.parquet",
    output:
        directory("data/embeddings_{suffix}"),
    params:
        methods=[
            "CCAIntegration",
            "HarmonyIntegration",
            "JointPCAIntegration",
            "RPCAIntegration",
            "FastMNNIntegration",
            "none",
        ],
    run:
        run_seurat_integration(input[0], input[1], output[0], "Assay", params.methods)


rule compute_pairwse_distances:
    input:
        "data/embeddings_{suffix}",
    output:
        directory("data/distances_{suffix}/{method}"),
    params:
        distance_metrics=["euclidean", "cosine", "correlation"],
    run:
        os.makedirs(output[0], exist_ok=True)
        input_file = f"{input[0]}/{wildcards.method}.parquet"
        df = pd.read_parquet(input_file)
        df.set_index("__index_level_0__", inplace=True)
        for metric in params.distance_metrics:
            output_file = f"{output[0]}/{metric}.parquet"
            distances = compute_distances(df, metric)
            distances.to_parquet(output_file)

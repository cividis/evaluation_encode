import pandas as pd
from tfsage.download import download_encode
from tfsage.embedding import run_seurat_integration


rule download_experiment:
    output:
        config["results_dir"] + "downloads/{experiment}.bed",
    retries: 5
    run:
        download_encode(wildcards.experiment, output[0])


checkpoint prepare_metadata:
    params:
        input_file=config["encode_metadata"],
    output:
        config["results_dir"] + "data/metadata.parquet",
    script:
        "../scripts/prepare_metadata.py"


def aggregate_input(wildcards):
    df = pd.read_parquet(checkpoints.prepare_metadata.get().output[0])
    experiments = df.index.tolist()
    return expand(
        config["results_dir"] + "downloads/{experiment}.bed", experiment=experiments
    )


rule extract_features:
    input:
        aggregate_input,
    output:
        config["results_dir"] + "data/rp_matrix_tss.parquet",
    threads: 48
    resources:
        mem_mb=96000,
    script:
        "../scripts/extract_features.py"


rule aggregate_genes:
    input:
        config["results_dir"] + "data/rp_matrix_tss.parquet",
    output:
        config["results_dir"] + "data/rp_matrix_gene.parquet",
    script:
        "../scripts/aggregate_genes.py"


rule embed_and_integrate:
    input:
        config["results_dir"] + "data/rp_matrix_{suffix}.parquet",
        config["results_dir"] + "data/metadata.parquet",
    output:
        directory(config["results_dir"] + "data/embeddings_{suffix}"),
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


rule pairwise_distances:
    input:
        config["results_dir"] + "data/embeddings_{suffix}",
    output:
        directory(config["results_dir"] + "data/distances_{suffix}/{method}"),
    params:
        distance_metrics=["euclidean", "cosine", "correlation"],
    script:
        "../scripts/pairwise_distances.py"

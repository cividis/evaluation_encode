import pandas as pd
from tfsage.download import download_encode
from tfsage.embedding import run_seurat_integration


rule download:
    output:
        config["results_dir"] + "downloads/{experiment}.bed",
    retries: 5
    run:
        download_encode(wildcards.experiment, output[0])


checkpoint metadata:
    input:
        config["encode_metadata"],
    output:
        config["results_dir"] + "data/metadata.parquet",
    script:
        "../scripts/prepare_data/metadata.py"


def all_experiments(wildcards):
    experiments = pd.read_parquet(checkpoints.metadata.get().output[0]).index.tolist()
    template = config["results_dir"] + "downloads/{experiment}.bed"
    return expand(template, experiment=experiments)


rule features_tss:
    input:
        all_experiments,
    output:
        config["results_dir"] + "data/rp_matrix/tss.parquet",
    threads: 48
    resources:
        mem_mb=96000,
    script:
        "../scripts/prepare_data/features_tss.py"


rule features_gene:
    input:
        config["results_dir"] + "data/rp_matrix/tss.parquet",
    output:
        config["results_dir"] + "data/rp_matrix/gene.parquet",
    script:
        "../scripts/prepare_data/features_gene.py"


rule embeddings:
    input:
        rp_matrix=config["results_dir"] + "data/rp_matrix/{features}.parquet",
        metadata=config["results_dir"] + "data/metadata.parquet",
    output:
        directory(config["results_dir"] + "data/embeddings/{features}"),
    params:
        methods=config["embedding_methods"],
    run:
        run_seurat_integration(
            rp_matrix=input.rp_matrix,
            metadata=input.metadata,
            output_dir=output[0],
            align_key="Assay",
            methods=params.methods,
        )


rule distances:
    input:
        config["results_dir"] + "data/embeddings/{features}",
    output:
        directory(config["results_dir"] + "data/distances/{features}/{method}"),
    params:
        metrics=config["distance_metrics"],
    script:
        "../scripts/prepare_data/distances.py"

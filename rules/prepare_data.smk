import pandas as pd
from tfsage.download import download_encode
from tfsage.features import extract_features_parallel, load_region_set

df = pd.read_parquet("data/metadata.parquet")
experiments = df.index.tolist()


rule download_experiment:
    output:
        "downloads/{experiment}.bed",
    run:
        download_encode(wildcards.experiment, output[0])


rule extract_features:
    input:
        expand("downloads/{experiment}.bed", experiment=experiments),
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

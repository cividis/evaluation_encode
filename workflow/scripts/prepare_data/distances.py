from snakemake.script import snakemake
import os
import pandas as pd
from tfsage import search


def compute_distances(input_dir, output_dir, method, metrics):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Read input embedding data
    input_file = f"{input_dir}/{method}.parquet"
    df = pd.read_parquet(input_file)
    df.set_index("__index_level_0__", inplace=True)

    # Compute and save distances for each metric
    for metric in metrics:
        output_file = f"{output_dir}/{metric}.parquet"
        distances = search.compute_distances(df, metric)
        distances.to_parquet(output_file)


compute_distances(
    snakemake.input[0],
    snakemake.output[0],
    snakemake.wildcards.method,
    snakemake.params.metrics,
)

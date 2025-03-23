import pandas as pd
import numpy as np
from tfsage.search import generate_ranked_list
from tfsage import generate
from snakemake.script import snakemake


def synthesize_experiments(
    distances_file, metadata_file, query_id, factor, top_n, data_dir, output_file
):
    # Read distances
    distances_df = pd.read_parquet(distances_file)

    # Generate scoring function
    p = distances_df.to_numpy().var()
    scoring_function = lambda x: np.exp(-(x**2) / p)

    # Read metadata
    metadata = pd.read_parquet(metadata_file).filter(
        items=["Biosample term name", "Assay", "Experiment target"],
        axis=1,
    )

    # Get cell type of the query
    cell_type = metadata.loc[query_id, "Biosample term name"]

    # Generate ranked list for the query, removing the cell type
    ranked_list = (
        generate_ranked_list(
            distances_df,
            query_id,
            metadata=metadata,
            scoring_function=scoring_function,
        )
        .query("`Experiment target` == @factor")
        .query("`Biosample term name` != @cell_type")
        .filter(["score", "Experiment target", "Biosample term name"], axis=1)
        .sort_values("score", ascending=False)
        .head(top_n)
    )

    bed_files = [data_dir + x + ".bed" for x in ranked_list.index]
    weights = ranked_list["score"].tolist()
    result = generate.synthesize_experiments(bed_files, weights)
    result.to_parquet(output_file)


synthesize_experiments(
    snakemake.params.distances,
    snakemake.input.metadata,
    snakemake.wildcards.query_id,
    snakemake.wildcards.factor,
    int(snakemake.wildcards.n),
    snakemake.params.data_dir,
    snakemake.output[0],
)

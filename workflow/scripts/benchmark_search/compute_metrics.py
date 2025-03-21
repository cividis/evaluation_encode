import pandas as pd
import numpy as np
from tqdm import tqdm
from sklearn import metrics
from tfsage.search import generate_ranked_list
from typing import Callable
from snakemake.script import snakemake


def ranked_list_metrics(ranked_list: pd.DataFrame) -> tuple:
    ranked_list = ranked_list.assign(
        rank=lambda x: x["score"].rank(ascending=False),
        precision=lambda x: x["hit"].cumsum() / x["rank"],
    )
    reciprocal_rank = 1 / ranked_list.query("hit")["rank"].values[0]
    average_precision = ranked_list.query("hit")["precision"].mean()
    ndcg_score = metrics.ndcg_score([ranked_list["hit"]], [ranked_list["score"]])
    return reciprocal_rank, average_precision, ndcg_score


def generate_results(
    distances_df: pd.DataFrame,
    metadata: pd.DataFrame,
    benchmark_df: pd.DataFrame,
    scoring_function: Callable,
) -> pd.DataFrame:
    results = []
    for query_id, cell_type, factor in tqdm(
        benchmark_df.itertuples(index=False), total=len(benchmark_df)
    ):
        ranked_list = (
            generate_ranked_list(
                distances_df,
                query_id,
                metadata=metadata,
                scoring_function=scoring_function,
            )
            .query("`Experiment target` == @factor")
            .filter(["score", "Experiment target", "Biosample term name"], axis=1)
            .assign(hit=lambda x: x["Biosample term name"] == cell_type)
        )
        metrics = ranked_list_metrics(ranked_list)
        results.append([query_id, cell_type, factor, *metrics])

    columns = [
        "query_id",
        "cell_type",
        "factor",
        "reciprocal_rank",
        "average_precision",
        "ndcg_score",
    ]
    df = pd.DataFrame(results, columns=columns)
    return df


def compute_metrics(distances_file, metadata_file, benchmark_file, output_file):
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

    # Read benchmark
    benchmark_df = (
        pd.read_csv(benchmark_file)
        .filter(
            items=["index_query", "Biosample term name", "Experiment target"],
            axis=1,
        )
        .drop_duplicates()
    )

    df = generate_results(distances_df, metadata, benchmark_df, scoring_function)
    df.to_csv(output_file, index=False)


compute_metrics(
    snakemake.params.distances,
    snakemake.input.metadata,
    snakemake.input.benchmark,
    snakemake.output[0],
)

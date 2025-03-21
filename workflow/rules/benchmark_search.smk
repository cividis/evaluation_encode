rule compute_metrics:
    input:
        distances_dir=config["results_dir"] + "data/distances/{features}/{method}",
        metadata=config["results_dir"] + "data/metadata.parquet",
        benchmark=config["benchmark_set"],
    output:
        config["results_dir"] + "benchmark/search/{features}/{method}/{metric}.csv",
    params:
        distances=lambda w, input: input.distances_dir + f"/{w.metric}.parquet",
    script:
        "../scripts/benchmark_search/compute_metrics.py"


rule aggregate_metrics:
    input:
        expand(
            config["results_dir"]
            + "benchmark/search/{features}/{method}/{metric}.csv",
            features=["tss", "gene"],
            method=config["embedding_methods"],
            metric=config["distance_metrics"],
        ),
    output:
        config["results_dir"] + "benchmark/search.csv",
    script:
        "../scripts/benchmark_search/aggregate_metrics.py"

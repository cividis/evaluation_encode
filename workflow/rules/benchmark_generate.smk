import pandas as pd

benchmark_df = pd.read_csv(config["benchmark_set"])


rule synthesize_experiments:
    input:
        distances_dir=config["results_dir"] + "data/distances/gene/FastMNNIntegration",
        metadata=config["results_dir"] + "data/metadata.parquet",
    output:
        config["results_dir"]
        + "benchmark/generate/tfsage_experiments/head_{n}/{query_id}_{factor}.parquet",
    params:
        distances=lambda w, input: input.distances_dir + "/cosine.parquet",
        data_dir=config["results_dir"] + "downloads/",
    script:
        "../scripts/benchmark_generate/synthesize_experiments.py"


rule generate_test_set:
    input:
        query_file=config["results_dir"] + "downloads/{query_id}.bed",
        target_file=config["results_dir"] + "downloads/{target_id}.bed",
    output:
        config["results_dir"]
        + "benchmark/generate/test_sets/{test_set_name}/{query_id}_{target_id}.parquet",
    params:
        config_test_set=lambda w: config["test_sets"][w.test_set_name],
    script:
        "../scripts/benchmark_generate/generate_test_set.py"


rule compute_metrics_tfsage:
    input:
        predictions=config["results_dir"]
        + "benchmark/generate/tfsage_experiments/head_{n}/{query_id}_{factor}.parquet",
        test_set=config["results_dir"]
        + "benchmark/generate/test_sets/{test_set_name}/{query_id}_{target_id}.parquet",
    output:
        config["results_dir"]
        + "benchmark/generate/metrics/{test_set_name}/tfsage/head_{n}/{query_id}_{target_id}_{factor}.json",
    params:
        method_class="tfsage",
        method_name=lambda w: f"tfsage-{w.n}",
    script:
        "../scripts/benchmark_generate/compute_classification_metrics.py"


rule compute_metrics_motif_scan:
    input:
        predictions=config["results_dir"] + "motif_scan/{motif_db}/{query_id}.bed",
        test_set=config["results_dir"]
        + "benchmark/generate/test_sets/{test_set_name}/{query_id}_{target_id}.parquet",
        m2f_path="resources/motif_databases_filtered/{motif_db}_filtered.motif2factors.txt",
    output:
        config["results_dir"]
        + "benchmark/generate/metrics/{test_set_name}/motif_scan/{motif_db}/{query_id}_{target_id}_{factor}.json",
    params:
        method_class="motif scan",
        method_name="{motif_db}",
    script:
        "../scripts/benchmark_generate/compute_classification_metrics.py"


rule aggregate_metrics_generate:
    input:
        all_tfsage=expand(
            expand(
                config["results_dir"]
                + "benchmark/generate/metrics/{test_set_name}/tfsage/head_{n}/{{query_id}}_{{target_id}}_{{factor}}.json",
                n=config["tfsage_n"],
                test_set_name=list(config["test_sets"].keys()),
            ),
            zip,
            query_id=benchmark_df["index_query"],
            target_id=benchmark_df["index_target"],
            factor=benchmark_df["Experiment target"],
        ),
        all_motif_scan=expand(
            expand(
                config["results_dir"]
                + "benchmark/generate/metrics/{test_set_name}/motif_scan/{motif_db}/{{query_id}}_{{target_id}}_{{factor}}.json",
                motif_db=config["motif_dbs"],
                test_set_name=list(config["test_sets"].keys()),
            ),
            zip,
            query_id=benchmark_df["index_query"],
            target_id=benchmark_df["index_target"],
            factor=benchmark_df["Experiment target"],
        ),
    output:
        config["results_dir"] + "benchmark/generate.csv",
    script:
        "../scripts/benchmark_generate/aggregate_metrics.py"

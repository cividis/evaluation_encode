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
        + "benchmark/generate/test_sets/{query_id}_{target_id}.parquet",
    params:
        peak_width=200,
        blacklist_chroms=["chrM", "chrY"],
        n_samples=1000,
        p_positive=0.5,
        random_state=42,
    script:
        "../scripts/benchmark_generate/generate_test_set.py"


rule compute_metrics_tfsage:
    input:
        predictions=config["results_dir"]
        + "benchmark/generate/tfsage_experiments/head_{n}/{query_id}_{factor}.parquet",
        test_set=config["results_dir"]
        + "benchmark/generate/test_sets/{query_id}_{target_id}.parquet",
    output:
        config["results_dir"]
        + "benchmark/generate/metrics/tfsage/head_{n}/{query_id}_{target_id}_{factor}.json",
    params:
        method_class="tfsage",
        method_name=lambda w: f"tfsage-{w.n}",
    script:
        "../scripts/benchmark_generate/compute_classification_metrics.py"


rule compute_metrics_motif_scan:
    input:
        predictions=config["results_dir"] + "motif_scan/{motif_db}/{query_id}.bed",
        test_set=config["results_dir"]
        + "benchmark/generate/test_sets/{query_id}_{target_id}.parquet",
        m2f_path="resources/motif_databases_filtered/{motif_db}_filtered.motif2factors.txt",
    output:
        config["results_dir"]
        + "benchmark/generate/metrics/motif_scan/{motif_db}/{query_id}_{target_id}_{factor}.json",
    params:
        method_class="motif scan",
        method_name="{motif_db}",
    script:
        "../scripts/benchmark_generate/compute_classification_metrics.py"

import pandas as pd


rule filter_motifs:
    input:
        motif2factors="resources/motif_databases/{motif_db}.motif2factors.txt",
        pfm="resources/motif_databases/{motif_db}.pfm",
    output:
        motif2factors="resources/motif_databases_filtered/{motif_db}_filtered.motif2factors.txt",
        pfm="resources/motif_databases_filtered/{motif_db}_filtered.pfm",
    params:
        tf_list=(
            pd.read_csv(config["benchmark_set"])["Experiment target"]
            .str.split("-")
            .str[0]
            .unique()
            .tolist()
        ),
    script:
        "../scripts/filter_motifs.py"


rule motif_scan:
    input:
        config["results_dir"] + "downloads/{experiment}.bed",
        pfm="resources/motif_databases_filtered/{motif_db}_filtered.pfm",
    output:
        config["results_dir"] + "motif_scan/{motif_db}/{experiment}.bed",
    params:
        opts=lambda w, input: f"-g GRCh38 -p {input.pfm}",
    threads: 12
    resources:
        mem_mb=24000,
    conda:
        "gimme_env"
    shell:
        """
        gimme scan {input[0]} {params.opts} -N {threads} -b > {output}
        """

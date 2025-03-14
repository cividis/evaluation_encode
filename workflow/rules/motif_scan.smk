rule motif_scan:
    input:
        config["results_dir"] + "downloads/{experiment}.bed",
    output:
        config["results_dir"] + "motif_scan/{experiment}.bed",
    params:
        opts="-g GRCh38 -p HOCOMOCOv11_HUMAN",
    threads: 12
    resources:
        mem_mb=24000,
    conda:
        "gimme_env"
    shell:
        """
        gimme scan {input} {params.opts} -N {threads} -b > {output}
        """

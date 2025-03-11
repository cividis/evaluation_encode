include: "rules/prepare_data.smk"


rule all:
    input:
        "data/rp_matrix.parquet",
        "data/rp_matrix_agg.parquet",

import pandas as pd
from snakemake.script import snakemake


def prepare_metadata(input_file, output_file):
    # Read metadata and filter required experiments
    df = (
        pd.read_csv(input_file, sep="\t")
        .query("`File format` == 'bed narrowPeak'")
        .query(
            "(`Assay` in ['ATAC-seq', 'TF ChIP-seq']) or (`Assay` == 'Histone ChIP-seq' and `Experiment target`.str.startswith('H3K27ac'))"
        )
    )

    # Process metadata and save as parquet
    df.set_index("File accession", inplace=True)
    df.index.name = None
    df.sort_values("Assay", inplace=True)
    df.to_parquet(output_file)


prepare_metadata(snakemake.input[0], snakemake.output[0])

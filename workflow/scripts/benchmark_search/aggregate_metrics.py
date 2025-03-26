import pandas as pd
from snakemake.script import snakemake


def aggregate_metrics(input_files, output_file):
    df_list = [pd.read_csv(f) for f in input_files]
    df = pd.concat(df_list)
    df.to_csv(output_file, index=False)


aggregate_metrics(snakemake.input, snakemake.output[0])

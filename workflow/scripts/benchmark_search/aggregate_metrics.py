import pandas as pd
from snakemake.script import snakemake


def return_dataframe(filepath):
    *_, features, method, metric = filepath.rsplit(".", 1)[0].split("/")
    df = pd.read_csv(filepath).assign(features=features, method=method, metric=metric)
    return df


def aggregate_metrics(input_files, output_file):
    df_list = [return_dataframe(f) for f in input_files]
    df = pd.concat(df_list)
    df.to_csv(output_file, index=False)


aggregate_metrics(snakemake.input, snakemake.output[0])

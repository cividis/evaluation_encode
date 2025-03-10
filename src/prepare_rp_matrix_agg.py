import pandas as pd

df = pd.read_parquet("data/rp_matrix.parquet")
df["gene"] = df.index.str.split(":").str[1]
df = df.groupby("gene").mean()
df.index.name = None
df.to_parquet("data/rp_matrix_agg.parquet")

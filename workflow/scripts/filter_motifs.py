import pandas as pd
from tqdm import tqdm
from snakemake.script import snakemake


def filter_motifs(
    motif2factors_path,
    pfm_path,
    tf_list,
    motif2factors_path_out,
    pfm_path_out,
):
    # Filter motif2factors database
    motif_db = pd.read_csv(motif2factors_path, sep="\t").query("Factor in @tf_list")
    motif_list = motif_db["Motif"].unique()

    # Filter PFM file
    pfm_filtered = filter_pfm(pfm_path, motif_list)

    # Save filtered motif2factors and PFM
    motif_db.to_csv(motif2factors_path_out, sep="\t", index=False)
    with open(pfm_path_out, "w") as f:
        f.writelines(pfm_filtered)


def filter_pfm(pfm_path, motif_list):
    with open(pfm_path, "r") as f:
        pfm_lines = f.readlines()

    include = False
    lines = []
    for line in tqdm(pfm_lines):
        if line.startswith(">"):
            this_motif = line.replace(">", "").strip()
            include = any(this_motif == motif for motif in motif_list)
        if include or line.startswith("#"):
            lines.append(line)

    return lines


filter_motifs(
    snakemake.input.motif2factors,
    snakemake.input.pfm,
    snakemake.params.tf_list,
    snakemake.output.motif2factors,
    snakemake.output.pfm,
)

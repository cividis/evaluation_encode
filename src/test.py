from tfsage.features import extract_features, load_region_set, prepare_region_set

gene_loc_set = load_region_set("hg38")
bed_file = "downloads/ENCFF624NUZ.bed"
region_set = prepare_region_set(bed_file, gene_loc_set.genome)


features = extract_features(bed_file, gene_loc_set)
print(features)

import pandas as pd
df = (
    pd.read_csv("metadata/metadata.tsv", sep="\t")
    .query("`File format` == 'bed narrowPeak'")
    .query(
        "(`Assay` in ['ATAC-seq', 'TF ChIP-seq']) or (`Assay` == 'Histone ChIP-seq' and `Experiment target`.str.startswith('H3K27ac'))",
    )
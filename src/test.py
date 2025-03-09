from tfsage.features import extract_features, load_region_set, prepare_region_set

gene_loc_set = load_region_set("hg38")
bed_file = "downloads/ENCFF624NUZ.bed"
region_set = prepare_region_set(bed_file, gene_loc_set.genome)


features = extract_features(bed_file, gene_loc_set)
print(features)

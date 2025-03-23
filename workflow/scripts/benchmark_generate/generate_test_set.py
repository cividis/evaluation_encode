from tfsage import utils
from snakemake.script import snakemake


def generate_test_set(
    query_file,
    target_file,
    peak_width,
    blacklist_chroms,
    n_samples,
    p_positive,
    random_state,
    output_file,
):
    df = utils.generate_test_set(query_file, target_file, peak_width)
    df = df.query("chrom not in @blacklist_chroms")
    df = utils.stratified_sample(df, n_samples, p_positive, random_state)
    df.to_parquet(output_file)


generate_test_set(
    snakemake.input.query_file,
    snakemake.input.target_file,
    snakemake.params.peak_width,
    snakemake.params.blacklist_chroms,
    snakemake.params.n_samples,
    snakemake.params.p_positive,
    snakemake.params.random_state,
    snakemake.output[0],
)

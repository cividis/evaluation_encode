from tfsage.embedding import run_seurat_integration

methods = [
    "CCAIntegration",
    "HarmonyIntegration",
    "JointPCAIntegration",
    "RPCAIntegration",
    "FastMNNIntegration",
    "none",
]

run_seurat_integration(
    rp_matrix="data/rp_matrix.parquet",
    metadata="data/metadata.parquet",
    output_dir="data/embeddings",
    methods=methods,
)

run_seurat_integration(
    rp_matrix="data/rp_matrix_agg.parquet",
    metadata="data/metadata.parquet",
    output_dir="data/embeddings_agg",
    methods=methods,
)

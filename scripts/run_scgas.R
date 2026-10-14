log_con <- file(snakemake@log[[1]], "w")
sink(log_con, type = "message")
on.exit({ sink(type = "message"); close(log_con) })
Sys.setenv(OPENBLAS_NUM_THREADS = snakemake@threads)

fragment_path <- snakemake@input[["fragment"]]
barcodes_path <- snakemake@input[["barcodes"]]
rna_counts_path <- snakemake@input[["rna_counts"]]
ref_dir       <- snakemake@params[["ref_dir"]]
genome        <- snakemake@params[["genome"]]
cfg        <- snakemake@params[["cfg"]]
n_cores    <- snakemake@params[["n_cores"]]
chunk_size <- snakemake@params[["chunk_size"]]
pkg_dir    <- snakemake@params[["pkg_dir"]]
out_path      <- snakemake@output[["activity"]]
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

# Load scGAS (from local checkout if not installed in the env)
if (!requireNamespace("scGAS", quietly = TRUE)) {
  pkg_dir <- normalizePath(pkg_dir, mustWork = FALSE)
  if (!file.exists(file.path(pkg_dir, "DESCRIPTION"))) {
    stop("[scgas] scGAS is not installed and no package DESCRIPTION was found at: ",
         pkg_dir, "\n  Either install scGAS into the activated r_benchmark env:\n",
         "    Rscript -e 'devtools::install(\"/path/to/scGAS-R\")'\n",
         "  or set `scgas.pkg_dir` in config.yaml to the scGAS-R checkout.")
  }
  message("[scgas] scGAS not installed; loading from local checkout: ", pkg_dir)
  devtools::load_all(pkg_dir, quiet = TRUE)
} else {
  library(scGAS)
}

barcodes   <- readLines(barcodes_path)
rna_counts <- readRDS(rna_counts_path)
shared_bc  <- intersect(barcodes, colnames(rna_counts))
message("[scgas] Cells: ", length(shared_bc))

lsi_dims <- seq(cfg$lsi_dims_start, cfg$lsi_dims_end)

# Intermediate results (Seurat objects + scGAS checkpoints) are kept in
# <outdir>/methods/scgas/scgas_run/. scGAS resumes from these checkpoints if
# they exist WITHOUT checking whether inputs or parameters changed: delete that
# folder after changing QC thresholds, input data or scgas parameters.
run_name <- "scgas_run"
message("[scgas] Intermediate results: ", file.path(dirname(out_path), run_name))

# ── Step 1: Preprocess ────────────────────────────────────────────────────────
message("[scgas] Step 1: Preprocess")
obj <- scgas_preprocess(
  fragment_path  = fragment_path,
  data_dir       = ref_dir,
  genome_version = genome,
  rna_matrix     = rna_counts[, shared_bc, drop = FALSE],
  cell_barcodes  = shared_bc,
  lsi_dims       = lsi_dims,
  n_cores        = n_cores,
  verbose        = TRUE,
  out_dir        = dirname(out_path),   # inherited by all later steps via obj@misc
  run_name       = run_name,
  save_obj       = TRUE                 # -> seurat_preprocessed.rds
)

# ── Step 2: Metacells ─────────────────────────────────────────────────────────
message("[scgas] Step 2: Metacells")
obj <- scgas_metacell(
  obj,
  gamma    = cfg$gamma,
  min_size = cfg$min_size,
  max_size = cfg$max_size,
  n_cores  = n_cores,
  save_obj = TRUE                       # -> seurat_metacell.rds
)

# ── Step 3: Train models ──────────────────────────────────────────────────────
message("[scgas] Step 3: Train Lasso models")
obj <- scgas_train_models(obj, n_cores = n_cores,
                          save_obj = TRUE)   # -> seurat_models_trained.rds

# ── Step 4: Compute scGAS ─────────────────────────────────────────────────────
message("[scgas] Step 4: Compute single-cell GAS")
obj <- scgas_compute(
  obj,
  lsi_dims          = lsi_dims,
  knn_k             = cfg$knn_k,
  seed_fraction     = cfg$seed_fraction,
  chunk_size        = chunk_size,
  run_dim_reduction = TRUE,             # scgaspca / scgasumap in the saved object
  save_obj          = TRUE              # -> seurat_scgas_computed.rds
)

# ── Extract and save ──────────────────────────────────────────────────────────
scgas_mat <- Seurat::GetAssayData(obj, assay = "scGAS", layer = "data")
scgas_mat <- as(scgas_mat, "sparseMatrix")
saveRDS(scgas_mat, out_path)
message("[scgas] Saved: ", out_path, "  (", nrow(scgas_mat), " genes x ", ncol(scgas_mat), " cells)")

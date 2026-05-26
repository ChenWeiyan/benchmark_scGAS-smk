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

# Load scGAS (from local package if not installed)
if (!requireNamespace("scGAS", quietly = TRUE)) {
  message("[scgas] Loading from local package: ", pkg_dir)
  devtools::load_all(pkg_dir, quiet = TRUE)
} else {
  library(scGAS)
}

barcodes   <- readLines(barcodes_path)
rna_counts <- readRDS(rna_counts_path)
shared_bc  <- intersect(barcodes, colnames(rna_counts))
message("[scgas] Cells: ", length(shared_bc))

lsi_dims <- seq(cfg$lsi_dims_start, cfg$lsi_dims_end)

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
  verbose        = TRUE
)

# ── Step 2: Metacells ─────────────────────────────────────────────────────────
message("[scgas] Step 2: Metacells")
obj <- scgas_metacell(
  obj,
  gamma    = cfg$gamma,
  min_size = cfg$min_size,
  max_size = cfg$max_size,
  n_cores  = n_cores
)

# ── Step 3: Train models ──────────────────────────────────────────────────────
message("[scgas] Step 3: Train Lasso models")
obj <- scgas_train_models(obj, n_cores = n_cores)

# ── Step 4: Compute scGAS ─────────────────────────────────────────────────────
message("[scgas] Step 4: Compute single-cell GAS")
obj <- scgas_compute(
  obj,
  lsi_dims          = lsi_dims,
  knn_k             = cfg$knn_k,
  seed_fraction     = cfg$seed_fraction,
  chunk_size        = chunk_size,
  run_dim_reduction = FALSE
)

# ── Extract and save ──────────────────────────────────────────────────────────
scgas_mat <- Seurat::GetAssayData(obj, assay = "scGAS", layer = "data")
scgas_mat <- as(scgas_mat, "sparseMatrix")
saveRDS(scgas_mat, out_path)
message("[scgas] Saved: ", out_path, "  (", nrow(scgas_mat), " genes x ", ncol(scgas_mat), " cells)")

log_con <- file(snakemake@log[[1]], "w")
sink(log_con, type = "message")
on.exit({ sink(type = "message"); close(log_con) })
Sys.setenv(OPENBLAS_NUM_THREADS = 1L)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

methods     <- snakemake@params[["methods"]]
methods_dir <- snakemake@params[["methods_dir"]]

for (m in methods) {
  raw_path <- file.path(methods_dir, m, "activity_raw.rds")
  out_path <- file.path(methods_dir, m, "activity_lognorm.rds")

  message("[normalize] Loading raw matrix for: ", m)
  mat <- readRDS(raw_path)
  if (!inherits(mat, "sparseMatrix")) mat <- as(mat, "CsparseMatrix")
  message("[normalize] ", m, ": ", nrow(mat), " genes x ", ncol(mat), " cells")

  so <- CreateSeuratObject(counts = mat, assay = "GA")
  so <- NormalizeData(so, assay = "GA", verbose = FALSE)
  norm_mat <- GetAssayData(so, assay = "GA", layer = "data")

  saveRDS(norm_mat, out_path)
  message("[normalize] Saved: ", out_path)
}

message("[normalize] Done.")

log_con <- file(snakemake@log[[1]], "w")
sink(log_con, type = "message")
on.exit({ sink(type = "message"); close(log_con) })
Sys.setenv(OPENBLAS_NUM_THREADS = snakemake@threads)

# TODO: Install ArchR: devtools::install_github("GreenleafLab/ArchR", ref="master", repos=BiocManager::repositories())
# ArchR >= 1.0.2 required

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
})

fragment_path <- snakemake@input[["fragment"]]
barcodes_path <- snakemake@input[["barcodes"]]
genome        <- snakemake@params[["genome"]]
cfg           <- snakemake@params[["cfg"]]
n_cores       <- snakemake@params[["n_cores"]]
archr_dir     <- snakemake@params[["archr_dir"]]
out_path      <- snakemake@output[["activity"]]
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
dir.create(archr_dir, showWarnings = FALSE, recursive = TRUE)

barcodes <- readLines(barcodes_path)
addArchRGenome(genome)
addArchRThreads(threads = n_cores)

# ── Create Arrow file ─────────────────────────────────────────────────────────
arrow_files <- createArrowFiles(
  inputFiles     = fragment_path,
  sampleNames    = "sample",
  validBarcodes  = barcodes,
  minTSS         = 0,       # QC already applied upstream
  minFrags       = 0,
  addTileMat     = TRUE,
  addGeneScoreMat= TRUE,
  force          = TRUE,
  outputDirectory= archr_dir
)

# ── Create ArchR project ──────────────────────────────────────────────────────
proj <- ArchRProject(ArrowFiles = arrow_files, outputDirectory = archr_dir, copyArrows = FALSE)
proj <- proj[barcodes[barcodes %in% getCellNames(proj)], ]

# ── Extract GeneScoreMatrix ───────────────────────────────────────────────────
message("[archr] Extracting GeneScoreMatrix ...")
gs_se <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
gs_mat <- assay(gs_se, "GeneScoreMatrix")
rownames(gs_mat) <- rowData(gs_se)$name

saveRDS(gs_mat, out_path)
message("[archr] Saved raw gene score matrix: ", out_path)

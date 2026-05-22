log_con <- file(snakemake@log[[1]], "w")
sink(log_con, type = "message")
on.exit({ sink(type = "message"); close(log_con) })

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(Matrix)
})

fragment_path <- snakemake@input[["fragment"]]
rna_path      <- snakemake@input[["rna"]]
genome        <- snakemake@params[["genome"]]
qc            <- snakemake@params[["qc"]]
n_cores       <- snakemake@params[["n_cores"]]

# ── Load RNA ──────────────────────────────────────────────────────────────────
message("[preprocess] Loading RNA: ", rna_path)
if (grepl("\\.h5$", rna_path, ignore.case = TRUE)) {
  rna_counts <- Read10X_h5(rna_path)
  if (is.list(rna_counts)) rna_counts <- rna_counts[["Gene Expression"]]
} else if (dir.exists(rna_path)) {
  rna_counts <- Read10X(data.dir = rna_path)
  if (is.list(rna_counts)) rna_counts <- rna_counts[["Gene Expression"]]
} else {
  stop("[preprocess] rna must be a .h5 file or MTX directory")
}
message("[preprocess] RNA: ", nrow(rna_counts), " genes x ", ncol(rna_counts), " cells")

# ── RNA QC ────────────────────────────────────────────────────────────────────
mt_pattern   <- "^MT-"
mt_genes     <- grep(mt_pattern, rownames(rna_counts), value = TRUE, ignore.case = TRUE)
lib_sz_rna   <- Matrix::colSums(rna_counts)
pct_mt       <- if (length(mt_genes) > 0)
  Matrix::colSums(rna_counts[mt_genes, ]) / lib_sz_rna * 100
else
  setNames(rep(0, ncol(rna_counts)), colnames(rna_counts))
n_feat_rna   <- Matrix::colSums(rna_counts > 0)

# ── Load gene annotation ──────────────────────────────────────────────────────
message("[preprocess] Loading annotation for genome: ", genome)
ann <- .load_annotation(genome)

# ── ATAC QC via coarse 10kb tiling ────────────────────────────────────────────
message("[preprocess] Building coarse tiling matrix for ATAC QC ...")
sl <- seqlengths(ann)
sl <- sl[!is.na(sl) & sl > 0]
tiles <- tileGenome(sl, tilewidth = 10000, cut.last.tile.in.chrom = TRUE)
tiles <- keepStandardChromosomes(tiles, pruning.mode = "coarse")

frag_obj <- CreateFragmentObject(fragment_path,
                                 cells = intersect(colnames(rna_counts),
                                                   .read_fragment_barcodes(fragment_path)))
cm_qc <- FeatureMatrix(fragments = frag_obj, features = tiles,
                        cells = colnames(rna_counts))
# drop all-zero cells
cm_qc <- cm_qc[, Matrix::colSums(cm_qc) > 0, drop = FALSE]

chrom_assay <- CreateChromatinAssay(counts = cm_qc, fragments = frag_obj,
                                    annotation = ann, min.cells = 0, min.features = 0)
obj_qc <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")

message("[preprocess] Computing NucleosomeSignal ...")
obj_qc <- NucleosomeSignal(obj_qc)

message("[preprocess] Computing TSSEnrichment ...")
obj_qc <- TSSEnrichment(obj_qc, fast = FALSE)

# ── Merge QC metrics ──────────────────────────────────────────────────────────
atac_cells <- colnames(obj_qc)
rna_cells  <- colnames(rna_counts)
shared     <- intersect(atac_cells, rna_cells)
message("[preprocess] Cells with both ATAC and RNA: ", length(shared))

qc_df <- data.frame(
  cell              = shared,
  nFeature_ATAC     = as.numeric(obj_qc$nFeature_ATAC[shared]),
  TSS_enrichment    = as.numeric(obj_qc$TSS.enrichment[shared]),
  nucleosome_signal = as.numeric(obj_qc$nucleosome_signal[shared]),
  nFeature_RNA      = as.numeric(n_feat_rna[shared]),
  pct_mt            = as.numeric(pct_mt[shared]),
  stringsAsFactors  = FALSE
)
saveRDS(qc_df, snakemake@output[["qc_metrics"]])

# ── Filter cells ──────────────────────────────────────────────────────────────
pass <- with(qc_df,
  nFeature_ATAC     >= qc$min_features_atac     &
  nFeature_ATAC     <= qc$max_features_atac     &
  TSS_enrichment    >= qc$min_tss_enrichment    &
  nucleosome_signal <= qc$max_nucleosome_signal &
  nFeature_RNA      >= qc$min_features_rna      &
  nFeature_RNA      <= qc$max_features_rna      &
  pct_mt            <= qc$max_pct_mt
)
filtered_cells <- qc_df$cell[pass]
message("[preprocess] Cells passing QC: ", length(filtered_cells), " / ", nrow(qc_df))
if (length(filtered_cells) < 50)
  warning("[preprocess] Very few cells passed QC — check thresholds")

writeLines(filtered_cells, snakemake@output[["barcodes"]])

# ── Save filtered RNA ─────────────────────────────────────────────────────────
rna_filt   <- rna_counts[, filtered_cells, drop = FALSE]
gene_n     <- Matrix::rowSums(rna_filt > 0)
keep_genes <- gene_n >= qc$min_cells_per_gene
rna_filt   <- rna_filt[keep_genes, , drop = FALSE]
message("[preprocess] Genes retained: ", nrow(rna_filt))
saveRDS(rna_filt, snakemake@output[["rna_counts"]])

lib_sz      <- Matrix::colSums(rna_filt)
rna_lognorm <- log1p(Matrix::t(Matrix::t(rna_filt) / lib_sz) * 1e4)
saveRDS(rna_lognorm, snakemake@output[["rna_lognorm"]])

message("[preprocess] Done.")


# ── Helpers ───────────────────────────────────────────────────────────────────

.load_annotation <- function(genome) {
  pkg_map <- c(
    hg19 = "EnsDb.Hsapiens.v75",
    hg38 = "EnsDb.Hsapiens.v86"
  )
  pkg <- pkg_map[[genome]]
  if (is.null(pkg)) stop("[preprocess] No EnsDb mapped for genome: ", genome)
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("[preprocess] Install: BiocManager::install('", pkg, "')")
  ann <- Signac::GetGRangesFromEnsDb(
    ensdb   = get(pkg, envir = asNamespace(pkg)),
    verbose = FALSE
  )
  GenomeInfoDb::seqlevelsStyle(ann) <- "UCSC"
  GenomeInfoDb::genome(ann)         <- genome
  ann
}

.read_fragment_barcodes <- function(frag_path, n_lines = 5e5) {
  con <- if (grepl("\\.gz$", frag_path)) gzcon(file(frag_path, "rb")) else file(frag_path, "r")
  on.exit(close(con))
  bcs <- character(0)
  repeat {
    lines <- readLines(con, n = 10000L, warn = FALSE)
    if (length(lines) == 0) break
    data_lines <- lines[!startsWith(lines, "#")]
    if (length(data_lines)) {
      fields <- do.call(rbind, strsplit(data_lines, "\t"))
      bcs    <- union(bcs, fields[, 4])
    }
    if (length(bcs) >= n_lines) break
  }
  bcs
}

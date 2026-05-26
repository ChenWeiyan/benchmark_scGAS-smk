log_con <- file(snakemake@log[[1]], "w")
sink(log_con, type = "message")
on.exit({ sink(type = "message"); close(log_con) })
Sys.setenv(OPENBLAS_NUM_THREADS = snakemake@threads)

# TODO: Install cicero: BiocManager::install("cicero")
# Requires monocle3: devtools::install_github("cole-trapnell-lab/monocle3")
# cicero >= 1.3 required

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(cicero)
  library(monocle3)
  library(GenomicRanges)
  library(Matrix)
})

fragment_path <- snakemake@input[["fragment"]]
barcodes_path <- snakemake@input[["barcodes"]]
genome        <- snakemake@params[["genome"]]
cfg           <- snakemake@params[["cfg"]]
n_cores       <- snakemake@params[["n_cores"]]
outdir        <- snakemake@params[["outdir"]]
out_path      <- snakemake@output[["activity"]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

barcodes <- readLines(barcodes_path)

# ── Load annotation & build peak matrix ──────────────────────────────────────
pkg_map <- c(hg19 = "EnsDb.Hsapiens.v75", hg38 = "EnsDb.Hsapiens.v86")
ann <- {
  pkg <- pkg_map[[genome]]
  a   <- Signac::GetGRangesFromEnsDb(ensdb = get(pkg, envir = asNamespace(pkg)), verbose = FALSE)
  GenomeInfoDb::seqlevelsStyle(a) <- "UCSC"
  a
}

frag_obj <- CreateFragmentObject(fragment_path, cells = barcodes)

# Use called peaks or tiles as features for cicero
# TODO: replace tile_gr with called peaks (from MACS2/ArchR) for better results
sl    <- seqlengths(ann); sl <- sl[!is.na(sl) & sl > 0]
tiles <- tileGenome(sl, tilewidth = 500, cut.last.tile.in.chrom = TRUE)
tiles <- keepStandardChromosomes(tiles, pruning.mode = "coarse")

message("[cicero] Building peak matrix ...")
peak_mat <- FeatureMatrix(fragments = frag_obj, features = tiles, cells = barcodes)

# ── Convert to CellDataSet for cicero ────────────────────────────────────────
message("[cicero] Building CellDataSet ...")
peak_mat <- peak_mat[Matrix::rowSums(peak_mat) > 0, , drop = FALSE]
cds <- new_cell_data_set(
  expression_data    = peak_mat,
  cell_metadata      = data.frame(cell = barcodes, row.names = barcodes),
  gene_metadata      = data.frame(
    gene_short_name  = rownames(peak_mat),
    chr              = as.character(GenomicRanges::seqnames(
                         GRanges(sub("-[0-9]+-[0-9]+$", "", rownames(peak_mat)),
                                 IRanges(1, 1)))),
    bp1              = as.numeric(gsub(".*-([0-9]+)-[0-9]+$", "\\1", rownames(peak_mat))),
    bp2              = as.numeric(gsub(".*-[0-9]+-([0-9]+)$", "\\1", rownames(peak_mat))),
    row.names        = rownames(peak_mat)
  )
)
cds <- detect_genes(cds)
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, method = "LSI")
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "LSI")

# ── Run cicero ────────────────────────────────────────────────────────────────
message("[cicero] Running cicero (window = ", cfg$window, ") ...")
umap_coords <- reducedDims(cds)[["UMAP"]]
cicero_out  <- run_cicero(cds, genomic_coords = umap_coords,
                          window = cfg$window, silent = TRUE)

# ── Gene activity: annotate connections, aggregate ───────────────────────────
message("[cicero] Annotating connections and building gene activity matrix ...")
gene_annot <- as.data.frame(ann)[, c("seqnames", "start", "end", "gene_name")]
gene_annot$chr <- as.character(gene_annot$seqnames)
unnorm_ga <- build_gene_activity_matrix(cds, cicero_out, gene_annotation_sub = gene_annot,
                                        window = cfg$distance_constraint)
unnorm_ga <- unnorm_ga[!grepl("^NA", rownames(unnorm_ga)), , drop = FALSE]

saveRDS(unnorm_ga, out_path)
message("[cicero] Saved raw gene activity: ", out_path)

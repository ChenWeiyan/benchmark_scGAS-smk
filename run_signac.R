log_con <- file(snakemake@log[[1]], "w")
sink(log_con, type = "message")
on.exit({ sink(type = "message"); close(log_con) })
Sys.setenv(OPENBLAS_NUM_THREADS = snakemake@threads)

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(Matrix)
})

fragment_path <- snakemake@input[["fragment"]]
barcodes_path <- snakemake@input[["barcodes"]]
genome        <- snakemake@params[["genome"]]
upstream      <- snakemake@params[["upstream"]]
downstream    <- snakemake@params[["downstream"]]
out_path      <- snakemake@output[["activity"]]
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

barcodes <- readLines(barcodes_path)
message("[signac] Cells: ", length(barcodes))

# ── Gene annotation ───────────────────────────────────────────────────────────
pkg_map <- c(hg19 = "EnsDb.Hsapiens.v75", hg38 = "EnsDb.Hsapiens.v86")
pkg <- pkg_map[[genome]]
if (is.null(pkg)) stop("[signac] No EnsDb for genome: ", genome)
if (!requireNamespace(pkg, quietly = TRUE))
  stop("[signac] Install: BiocManager::install('", pkg, "')")

genes_gr <- genes(get(pkg, envir = asNamespace(pkg)),
                  filter  = ~ gene_biotype == "protein_coding",
                  columns = c("gene_id", "gene_name"))
seqlevelsStyle(genes_gr) <- "UCSC"
genome(genes_gr)          <- genome
genes_gr <- keepStandardChromosomes(genes_gr, pruning.mode = "coarse")

message("[signac] Genes: ", length(genes_gr))

# ── Extend gene bodies ────────────────────────────────────────────────────────
genes_ext <- Extend(genes_gr, upstream = upstream, downstream = downstream)

# ── Fragment object with filtered cells ──────────────────────────────────────
frag_obj <- CreateFragmentObject(fragment_path, cells = barcodes)

# ── Gene activity matrix ──────────────────────────────────────────────────────
message("[signac] Computing FeatureMatrix (", length(genes_ext), " gene bodies) ...")
ga_raw <- FeatureMatrix(
  fragments = frag_obj,
  features  = genes_ext,
  cells     = barcodes
)

# Map region strings to gene names
region_str    <- GRangesToString(genes_ext)
gene_name_lut <- setNames(genes_ext$gene_name, region_str)
rownames(ga_raw) <- gene_name_lut[rownames(ga_raw)]
ga_raw <- ga_raw[!is.na(rownames(ga_raw)), , drop = FALSE]
ga_raw <- ga_raw[!duplicated(rownames(ga_raw)), , drop = FALSE]

ga_raw <- ga_raw[, barcodes, drop = FALSE]   # ensure column order matches barcodes
message("[signac] Gene activity matrix: ", nrow(ga_raw), " genes x ", ncol(ga_raw), " cells")

saveRDS(ga_raw, out_path)
message("[signac] Saved raw gene activity: ", out_path)

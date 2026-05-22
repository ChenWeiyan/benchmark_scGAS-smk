log_con <- file(snakemake@log[[1]], "w")
sink(log_con, type = "message")
on.exit({ sink(type = "message"); close(log_con) })

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(Rsamtools)
})

frag_path  <- snakemake@input[["fragment"]]
chain_file <- snakemake@input[["chain"]]
out_gz     <- snakemake@output[["lifted"]]
n_cores    <- snakemake@params[["n_cores"]]
outdir     <- snakemake@params[["outdir"]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("[liftover] Reading fragment file: ", frag_path)
frags <- fread(frag_path, sep = "\t", header = FALSE,
               col.names = c("chr", "start", "end", "barcode", "score"),
               nThread   = n_cores)
message("[liftover] Fragments loaded: ", nrow(frags))

message("[liftover] Converting to GRanges (BED 0-based → 1-based) ...")
gr <- GRanges(
  seqnames = frags$chr,
  ranges   = IRanges(start = frags$start + 1L, end = frags$end),
  barcode  = frags$barcode,
  score    = frags$score
)
rm(frags); gc()

message("[liftover] Loading chain: ", chain_file)
chain <- import.chain(chain_file)

message("[liftover] Running liftOver ...")
gr_lifted <- unlist(liftOver(gr, chain))
n_in   <- length(gr)
n_out  <- length(gr_lifted)
rm(gr, chain); gc()
message(sprintf("[liftover] Retained %d / %d fragments (%.1f%%)",
                n_out, n_in, 100 * n_out / n_in))

out_df <- data.table(
  chr     = as.character(seqnames(gr_lifted)),
  start   = start(gr_lifted) - 1L,   # back to 0-based BED
  end     = end(gr_lifted),
  barcode = gr_lifted$barcode,
  score   = gr_lifted$score
)
setorder(out_df, chr, start)
rm(gr_lifted); gc()

out_tsv <- sub("\\.gz$", "", out_gz)
message("[liftover] Writing sorted TSV: ", out_tsv)
fwrite(out_df, out_tsv, sep = "\t", col.names = FALSE, nThread = n_cores)
rm(out_df); gc()

message("[liftover] bgzip + tabix indexing ...")
bgzip(out_tsv, dest = out_gz, overwrite = TRUE)
indexTabix(out_gz, format = "bed")
file.remove(out_tsv)

message("[liftover] Done: ", out_gz)

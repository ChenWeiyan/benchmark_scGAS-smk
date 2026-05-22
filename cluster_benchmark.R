log_con <- file(snakemake@log[[1]], "w")
sink(log_con, type = "message")
on.exit({ sink(type = "message"); close(log_con) })
Sys.setenv(OPENBLAS_NUM_THREADS = 1L)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(Matrix)
  library(igraph)
  library(aricode)
})

methods       <- snakemake@params[["methods"]]
methods_dir   <- snakemake@params[["methods_dir"]]
cell_type_csv <- snakemake@params[["cell_type_csv"]]
cell_type_col <- snakemake@params[["cell_type_col"]]
resolutions   <- snakemake@params[["resolutions"]]
n_pcs         <- snakemake@params[["n_pcs"]]
out_csv       <- snakemake@output[["metrics_csv"]]
out_pdf       <- snakemake@output[["plot_pdf"]]
out_png       <- snakemake@output[["plot_png"]]
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)

# ── Cell type labels ───────────────────────────────────────────────────────────
message("[cluster] Loading cell type labels from: ", cell_type_csv)
ct_dt  <- fread(cell_type_csv, select = c("barcode", cell_type_col))
ct_map <- setNames(ct_dt[[cell_type_col]], ct_dt$barcode)

# ── Per-method clustering sweep ───────────────────────────────────────────────
all_metrics <- list()

for (m in methods) {
  path <- file.path(methods_dir, m, "activity_lognorm.rds")
  if (!file.exists(path)) { warning("[cluster] Missing: ", path); next }
  message("[cluster] Processing method: ", m)
  mat <- readRDS(path)
  if (!inherits(mat, "sparseMatrix")) mat <- as(mat, "CsparseMatrix")

  cells       <- colnames(mat)
  true_labels <- ct_map[cells]
  keep        <- !is.na(true_labels) & nchar(true_labels) > 0
  message(sprintf("[cluster] %s: %d labelled cells / %d total", m, sum(keep), length(cells)))
  if (sum(keep) < 50) { warning("[cluster] Too few labelled cells, skipping: ", m); next }

  mat_k  <- mat[, keep, drop = FALSE]
  true_k <- as.character(true_labels[keep])

  so <- CreateSeuratObject(counts = mat_k, assay = "GA")
  n_var  <- min(3000, nrow(mat_k))
  n_pc   <- min(n_pcs, n_var - 1L, sum(keep) - 1L)
  so <- FindVariableFeatures(so, nfeatures = n_var, verbose = FALSE)
  so <- ScaleData(so, features = VariableFeatures(so), verbose = FALSE)
  so <- RunPCA(so, features = VariableFeatures(so), npcs = n_pc, verbose = FALSE)
  so <- FindNeighbors(so, reduction = "pca", dims = seq_len(n_pc), verbose = FALSE)

  for (res in resolutions) {
    so_r  <- FindClusters(so, resolution = res, verbose = FALSE)
    pred  <- as.character(Idents(so_r))
    true_int <- as.integer(as.factor(true_k))
    pred_int <- as.integer(as.factor(pred))
    n_cl  <- length(unique(pred))
    ari   <- igraph::compare(true_int, pred_int, method = "adjusted.rand")
    nmi   <- igraph::compare(true_int, pred_int, method = "nmi")
    ami   <- aricode::AMI(true_k, pred)
    all_metrics[[length(all_metrics) + 1L]] <- data.frame(
      method = m, resolution = res, n_clusters = n_cl,
      ARI = ari, NMI = nmi, AMI = ami, stringsAsFactors = FALSE
    )
    message(sprintf("[cluster] %-12s  res=%.2f  k=%d  ARI=%.3f  NMI=%.3f  AMI=%.3f",
                    m, res, n_cl, ari, nmi, ami))
  }
}

if (length(all_metrics) == 0) stop("[cluster] No results produced.")
metrics_df <- do.call(rbind, all_metrics)
fwrite(metrics_df, out_csv)

# ── Colour palette ─────────────────────────────────────────────────────────────
base_cols <- c(scgas="#E64B35", signac="#4DBBD5", archr="#00A087",
               cicero="#3C5488", linger="#F39B7F", cisformer="#8491B4", scarlink="#91D1C2")
method_lvls    <- unique(metrics_df$method)
method_colours <- base_cols[method_lvls]
method_colours[is.na(method_colours)] <- "#999999"
names(method_colours) <- method_lvls
metrics_df$method <- factor(metrics_df$method, levels = method_lvls)

# ── Long format for ARI / NMI / AMI line plot ─────────────────────────────────
metrics_long <- reshape(
  metrics_df[, c("method", "resolution", "ARI", "NMI", "AMI")],
  varying   = c("ARI", "NMI", "AMI"),
  v.names   = "value",
  timevar   = "metric",
  times     = c("ARI", "NMI", "AMI"),
  direction = "long",
  idvar     = c("method", "resolution")
)
metrics_long$metric <- factor(metrics_long$metric, levels = c("ARI", "NMI", "AMI"))
metrics_long$method <- factor(metrics_long$method, levels = method_lvls)

p_line <- ggplot(metrics_long, aes(x = resolution, y = value,
                                    colour = method, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = method_colours) +
  labs(title = "Clustering benchmark: activity matrix → Leiden resolution sweep",
       x = "Resolution", y = "Score", colour = "Method") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

n_cl_df <- unique(metrics_df[, c("method", "resolution", "n_clusters")])
p_ncl <- ggplot(n_cl_df, aes(x = resolution, y = n_clusters,
                               colour = method, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_manual(values = method_colours) +
  labs(title = "Number of clusters by resolution", x = "Resolution", y = "N clusters") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

combined <- p_line / p_ncl +
  plot_layout(heights = c(3, 1)) +
  plot_annotation(
    title    = "Clustering benchmark vs cell-type ground truth",
    subtitle = sprintf("Methods: %s  |  %d resolutions  |  cell type column: '%s'",
                       paste(method_lvls, collapse = ", "),
                       length(resolutions), cell_type_col),
    theme = theme(plot.title    = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 11))
  )

ggsave(out_pdf, combined, width = 14, height = 10)
ggsave(out_png, combined, width = 14, height = 10, dpi = 200)
message("[cluster] Done. Metrics saved: ", out_csv)

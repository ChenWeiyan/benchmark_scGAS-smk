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
})

rna_lognorm_path <- snakemake@input[["rna_lognorm"]]
methods          <- snakemake@params[["methods"]]
methods_dir      <- snakemake@params[["methods_dir"]]
min_expr_frac    <- snakemake@params[["min_expr_frac"]]
n_hvg            <- snakemake@params[["n_hvg"]]
n_top_deg        <- snakemake@params[["n_top_deg"]]
cell_type_csv    <- snakemake@params[["cell_type_csv"]]
cell_type_col    <- snakemake@params[["cell_type_col"]]
out_csv          <- snakemake@output[["pcc_csv"]]
out_pdf          <- snakemake@output[["plot_pdf"]]
out_png          <- snakemake@output[["plot_png"]]
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)

# ── Load RNA log-norm ──────────────────────────────────────────────────────────
message("[evaluate] Loading RNA log-norm ...")
rna_mat   <- readRDS(rna_lognorm_path)
if (!inherits(rna_mat, "sparseMatrix")) rna_mat <- as(rna_mat, "CsparseMatrix")
expr_frac <- Matrix::rowMeans(rna_mat > 0)
rna_genes <- rownames(rna_mat)[expr_frac >= min_expr_frac]
message("[evaluate] RNA genes retained (expr_frac >= ", min_expr_frac, "): ", length(rna_genes))

# ── Load activity matrices ─────────────────────────────────────────────────────
activity_list <- list()
for (m in methods) {
  path <- file.path(methods_dir, m, "activity_lognorm.rds")
  if (!file.exists(path)) { warning("[evaluate] Missing: ", path); next }
  mat <- readRDS(path)
  if (!inherits(mat, "sparseMatrix")) mat <- as(mat, "CsparseMatrix")
  activity_list[[m]] <- mat
  message("[evaluate] Loaded '", m, "': ", nrow(mat), " genes x ", ncol(mat), " cells")
}
if (length(activity_list) == 0) stop("[evaluate] No activity matrices loaded.")

# ── Common genes / cells ───────────────────────────────────────────────────────
common_genes <- Reduce(intersect, c(list(rna_genes), lapply(activity_list, rownames)))
common_cells <- Reduce(intersect, c(list(colnames(rna_mat)), lapply(activity_list, colnames)))
message("[evaluate] Common genes: ", length(common_genes),
        "  |  Common cells: ", length(common_cells))
rna_sub <- as.matrix(rna_mat[common_genes, common_cells, drop = FALSE])

# ── HVG gene set (top n_hvg from RNA) ─────────────────────────────────────────
message("[evaluate] Selecting top ", n_hvg, " HVGs from RNA ...")
so_hvg <- CreateSeuratObject(counts = rna_mat[rna_genes, common_cells, drop = FALSE])
so_hvg <- SetAssayData(so_hvg, assay = "RNA", layer = "data",
                       new.data = rna_mat[rna_genes, common_cells, drop = FALSE])
so_hvg <- FindVariableFeatures(so_hvg, nfeatures = n_hvg, verbose = FALSE)
hvg_genes <- intersect(common_genes, VariableFeatures(so_hvg))
rm(so_hvg)
message("[evaluate] HVG gene set: ", length(hvg_genes), " genes")

# ── DEG gene set (optional – requires cell_type_csv) ──────────────────────────
deg_genes <- NULL
if (nchar(cell_type_csv) > 0 && file.exists(cell_type_csv)) {
  message("[evaluate] Computing DEGs from: ", cell_type_csv)
  ct_dt  <- fread(cell_type_csv, select = c("barcode", cell_type_col))
  ct_map <- setNames(ct_dt[[cell_type_col]], ct_dt$barcode)
  ct_vec <- ct_map[common_cells]
  ct_vec[is.na(ct_vec) | ct_vec == ""] <- "Unknown"
  message("[evaluate] Cell types: ", paste(sort(unique(ct_vec)), collapse = ", "))

  so_deg <- CreateSeuratObject(counts = rna_mat[rna_genes, common_cells, drop = FALSE])
  so_deg <- SetAssayData(so_deg, assay = "RNA", layer = "data",
                         new.data = rna_mat[rna_genes, common_cells, drop = FALSE])
  Idents(so_deg) <- ct_vec
  markers <- FindAllMarkers(so_deg, assay = "RNA", only.pos = TRUE,
                            min.pct = 0.10, logfc.threshold = 0.25,
                            test.use = "wilcox", verbose = FALSE)
  fwrite(markers, file.path(dirname(out_csv), "deg_markers.csv"))
  markers_sig  <- markers[markers$p_val_adj < 0.05, ]
  top_per_type <- do.call(rbind, lapply(split(markers_sig, markers_sig$cluster), function(x)
    head(x[order(-x$avg_log2FC), ], n_top_deg)))
  deg_genes <- unique(intersect(top_per_type$gene, common_genes))
  rm(so_deg)
  message("[evaluate] DEG gene set: ", length(deg_genes), " genes")
} else {
  message("[evaluate] No cell_type_csv — skipping DEG gene set")
}

# ── Fast row-wise PCC ──────────────────────────────────────────────────────────
fast_pcc <- function(A, B) {
  A <- as.matrix(A); B <- as.matrix(B)
  A_c <- A - rowMeans(A); B_c <- B - rowMeans(B)
  num <- rowSums(A_c * B_c)
  den <- sqrt(rowSums(A_c^2)) * sqrt(rowSums(B_c^2))
  ifelse(den == 0, NA_real_, num / den)
}

compute_pcc_set <- function(genes, label) {
  message("[evaluate] PCC [", label, "] — ", length(genes), " genes ...")
  rna_g <- rna_sub[genes, , drop = FALSE]
  do.call(rbind, lapply(names(activity_list), function(m) {
    act_g <- as.matrix(activity_list[[m]][genes, common_cells, drop = FALSE])
    pcc   <- fast_pcc(act_g, rna_g)
    data.frame(filter = label, method = m, gene = genes, pcc = pcc,
               stringsAsFactors = FALSE)
  }))
}

gene_sets <- list(all = common_genes, hvg = hvg_genes)
if (!is.null(deg_genes)) gene_sets[["deg"]] <- deg_genes

pcc_raw <- do.call(rbind, lapply(names(gene_sets), function(f)
  compute_pcc_set(gene_sets[[f]], f)))
pcc_raw <- pcc_raw[!is.na(pcc_raw$pcc), ]
fwrite(pcc_raw, out_csv)

for (f in names(gene_sets)) {
  for (m in names(activity_list)) {
    v <- pcc_raw$pcc[pcc_raw$filter == f & pcc_raw$method == m]
    message(sprintf("[evaluate] [%s] %-12s  median=%.3f  n=%d", f, m, median(v), length(v)))
  }
}

# ── Colour palette ─────────────────────────────────────────────────────────────
base_cols <- c(scgas="#E64B35", signac="#4DBBD5", archr="#00A087",
               cicero="#3C5488", linger="#F39B7F", cisformer="#8491B4", scarlink="#91D1C2")
method_lvls    <- unique(pcc_raw$method)
method_colours <- base_cols[method_lvls]
method_colours[is.na(method_colours)] <- "#999999"
names(method_colours) <- method_lvls

# ── Build factor labels with gene counts ──────────────────────────────────────
gene_n <- vapply(names(gene_sets), function(f) length(gene_sets[[f]]), integer(1))
label_map <- c(
  all = sprintf("All\n(%d genes)", gene_n["all"]),
  hvg = sprintf("HVG\n(%d genes)", gene_n["hvg"])
)
if (!is.null(deg_genes)) label_map["deg"] <- sprintf("DEG\n(%d genes)", gene_n["deg"])

pcc_plot          <- pcc_raw
pcc_plot$filter   <- factor(pcc_plot$filter, levels = names(label_map),
                             labels = unname(label_map))
pcc_plot$method   <- factor(pcc_plot$method, levels = method_lvls)

med_df <- aggregate(pcc ~ filter + method, data = pcc_plot, FUN = median)

# ── Panel A: density faceted by filter ────────────────────────────────────────
p_density <- ggplot(pcc_plot, aes(x = pcc, colour = method, fill = method)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  geom_vline(data = med_df, aes(xintercept = pcc, colour = method),
             linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ filter, nrow = 1) +
  scale_colour_manual(values = method_colours) +
  scale_fill_manual(values = method_colours) +
  labs(title = "PCC distribution", x = "Pearson corr. (activity vs RNA)", y = "Density") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

# ── Panel B: violin + box faceted by filter ───────────────────────────────────
p_box <- ggplot(pcc_plot, aes(x = method, y = pcc, fill = method)) +
  geom_violin(alpha = 0.4, colour = NA) +
  geom_boxplot(width = 0.18, outlier.size = 0.3, outlier.alpha = 0.2, colour = "grey20") +
  stat_summary(fun = median, geom = "point", size = 2, colour = "white") +
  facet_wrap(~ filter, nrow = 1) +
  scale_fill_manual(values = method_colours) +
  labs(x = NULL, y = "Pearson corr. (activity vs RNA)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1))

# ── Panel C: scatter or median-bar for "all" gene set ─────────────────────────
ref_method    <- if ("scgas" %in% method_lvls) "scgas" else method_lvls[1]
other_methods <- setdiff(method_lvls, ref_method)
all_label     <- label_map["all"]

ref_pcc <- setNames(
  pcc_plot$pcc[pcc_plot$filter == all_label & pcc_plot$method == ref_method],
  pcc_plot$gene[pcc_plot$filter == all_label & pcc_plot$method == ref_method]
)

if (length(other_methods) == 1) {
  m2      <- other_methods[1]
  cmp_pcc <- setNames(
    pcc_plot$pcc[pcc_plot$filter == all_label & pcc_plot$method == m2],
    pcc_plot$gene[pcc_plot$filter == all_label & pcc_plot$method == m2]
  )
  shared  <- intersect(names(ref_pcc), names(cmp_pcc))
  sdf     <- na.omit(data.frame(x = cmp_pcc[shared], y = ref_pcc[shared]))
  lim     <- range(c(sdf$x, sdf$y))
  wt      <- suppressWarnings(wilcox.test(sdf$y, sdf$x, paired = TRUE))
  p_lbl   <- if (wt$p.value < 2.2e-16) "p < 2.2e-16" else sprintf("p = %.2e", wt$p.value)
  p_scatter <- ggplot(sdf, aes(x = x, y = y)) +
    geom_point(alpha = 0.15, size = 0.6, colour = "grey40") +
    geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
    geom_smooth(method = "lm", se = TRUE, colour = method_colours[ref_method], linewidth = 0.8) +
    coord_equal(xlim = lim, ylim = lim) +
    annotate("text", x = lim[1], y = lim[2], label = p_lbl, hjust = 0, vjust = 1, size = 3.5) +
    labs(title = sprintf("Per-gene PCC: %s vs %s (all genes)", ref_method, m2),
         x = paste(m2, "PCC"), y = paste(ref_method, "PCC")) +
    theme_classic(base_size = 12)
} else {
  rank_df <- aggregate(pcc ~ method,
                       data = pcc_plot[pcc_plot$filter == all_label, ], FUN = median)
  rank_df <- rank_df[order(rank_df$pcc, decreasing = TRUE), ]
  rank_df$method <- factor(rank_df$method, levels = rank_df$method)
  p_scatter <- ggplot(rank_df, aes(x = method, y = pcc, fill = method)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = method_colours) +
    labs(title = "Median PCC ranking (all genes)", x = NULL, y = "Median Pearson corr.") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
}

# ── Panel D: ECDF faceted by filter ───────────────────────────────────────────
p_ecdf <- ggplot(pcc_plot, aes(x = pcc, colour = method)) +
  stat_ecdf(linewidth = 0.8) +
  geom_vline(data = med_df, aes(xintercept = pcc, colour = method),
             linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ filter, nrow = 1) +
  scale_colour_manual(values = method_colours) +
  labs(x = "Pearson corr. (activity vs RNA)", y = "Fraction of genes") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

# ── Assemble & save ────────────────────────────────────────────────────────────
combined <- (p_density / p_box) | (p_scatter / p_ecdf)
combined <- combined +
  plot_annotation(
    title    = "Gene activity benchmark: chromatin accessibility vs RNA expression",
    subtitle = sprintf("Gene sets: %s  |  %d common cells",
                       paste(names(label_map), collapse = " / "), length(common_cells)),
    theme = theme(plot.title    = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 11))
  )

ggsave(out_pdf, combined, width = 18, height = 12)
ggsave(out_png, combined, width = 18, height = 12, dpi = 200)
message("[evaluate] Plots saved. Done.")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("PROJECT_ROOT")
if (BASE == "") stop("PROJECT_ROOT is not set.")
RUN <- file.path(BASE, "bampe_peak_validation_20260622")

OUTDIR <- file.path(RUN, "plots_validation_extra")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

log2cpm_file <- file.path(RUN, "counts", "bampe_consensus_log2cpm.clean.tsv")
metadata_file <- file.path(RUN, "differential_atac", "bampe_sample_metadata.tsv")
published_file <- file.path(BASE, "data", "raw", "Supplementary_data_3_ATAC_stimulation_DA_peaks.txt.gz")
pairs_file <- file.path(RUN, "correlation", "bampe_to_geo_best_overlap_pairs.tsv")
da_dir <- file.path(RUN, "differential_atac")

message("Output directory: ", OUTDIR)

# -------------------------------
# 1. SPEARMAN SAMPLE HEATMAP
# -------------------------------

message("Reading BAMPE log2CPM matrix...")
mat_dt <- fread(log2cpm_file)

# First column should be peak_id
peak_col <- names(mat_dt)[1]
peak_ids <- mat_dt[[peak_col]]

mat <- as.matrix(mat_dt[, -1, with = FALSE])
rownames(mat) <- peak_ids

# Remove peaks with zero variance
message("Filtering zero-variance peaks...")
vars <- apply(mat, 1, var, na.rm = TRUE)
mat <- mat[vars > 0, , drop = FALSE]
vars <- vars[vars > 0]

# To keep this fast and interpretable, use top variable peaks
top_n <- min(50000, nrow(mat))
message("Using top variable peaks for sample correlation: ", top_n)

top_idx <- order(vars, decreasing = TRUE)[seq_len(top_n)]
mat_top <- mat[top_idx, , drop = FALSE]

message("Computing sample-wise Spearman correlation matrix...")
cor_mat <- cor(mat_top, method = "spearman", use = "pairwise.complete.obs")

cor_dt <- as.data.table(as.table(cor_mat))
setnames(cor_dt, c("sample_1", "sample_2", "spearman"))

# Read metadata if available
if (file.exists(metadata_file)) {
  meta <- fread(metadata_file)
  
  # Try to detect sample column
  sample_candidates <- c("sample", "sample_id", "sample_name")
  sample_col <- sample_candidates[sample_candidates %in% names(meta)][1]
  
  if (!is.na(sample_col)) {
    meta_small <- unique(meta[, .(
      sample = get(sample_col),
      cell_type = if ("cell_type" %in% names(meta)) cell_type else NA_character_,
      condition = if ("condition" %in% names(meta)) condition else NA_character_,
      donor = if ("donor" %in% names(meta)) donor else NA_character_
    )])
    
    # Keep matrix sample order
    sample_order <- colnames(cor_mat)
    cor_dt[, sample_1 := factor(sample_1, levels = sample_order)]
    cor_dt[, sample_2 := factor(sample_2, levels = rev(sample_order))]
  }
} else {
  sample_order <- colnames(cor_mat)
  cor_dt[, sample_1 := factor(sample_1, levels = sample_order)]
  cor_dt[, sample_2 := factor(sample_2, levels = rev(sample_order))]
}

p_heat <- ggplot(cor_dt, aes(x = sample_1, y = sample_2, fill = spearman)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0.5,
    limits = c(0, 1),
    name = "Spearman"
  ) +
  labs(
    title = "BAMPE log2CPM sample-wise Spearman correlation",
    subtitle = paste0("Top ", format(top_n, big.mark = ","), " most variable peaks"),
    x = "Sample",
    y = "Sample"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(OUTDIR, "08_BAMPE_log2CPM_sample_spearman_heatmap.png"),
       p_heat, width = 8, height = 7, dpi = 300)

ggsave(file.path(OUTDIR, "08_BAMPE_log2CPM_sample_spearman_heatmap.pdf"),
       p_heat, width = 8, height = 7)

message("Saved Spearman heatmap.")

# -------------------------------
# 2. RECONSTRUCT DA LOGFC MATCHED TABLE
# -------------------------------

message("Reading BAMPE-to-GEO best-overlap pairs...")
pairs <- fread(pairs_file)

# Make column names robust
# Expected columns should include BAMPE peak and GEO peak IDs.
print(names(pairs))

# Detect BAMPE/GEO peak ID columns.
# If the file has no header, fread names columns V1, V2, V3.
if (all(c("V1", "V2") %in% names(pairs))) {
  message("Pairs file has no header. Assuming V1 = BAMPE peak ID, V2 = GEO peak ID.")
  bampe_col <- "V1"
  geo_col <- "V2"
} else {
  bampe_col_candidates <- c("bampe_peak_id", "peak_id", "query_peak_id", "bampe_peak")
  geo_col_candidates <- c("geo_peak_id", "reference_peak_id", "subject_peak_id", "geo_peak")

  bampe_col <- bampe_col_candidates[bampe_col_candidates %in% names(pairs)][1]
  geo_col <- geo_col_candidates[geo_col_candidates %in% names(pairs)][1]

  if (is.na(bampe_col) || is.na(geo_col)) {
    stop("Could not detect BAMPE/GEO peak ID columns in pairs file. Check names(pairs).")
  }
}

pairs <- pairs[, .(
  bampe_peak_id = get(bampe_col),
  geo_peak_id = get(geo_col)
)]
pairs <- unique(pairs)

message("Reading Calderon published stimulation DA peaks...")
pub <- fread(published_file)

# Published file expected columns: logFC, AveExpr, t, P.Value, adj.P.Val, B, peak_id, contrast
required_pub <- c("logFC", "peak_id", "contrast")
missing_pub <- setdiff(required_pub, names(pub))
if (length(missing_pub) > 0) {
  stop("Published DA file missing required columns: ", paste(missing_pub, collapse = ", "))
}

pub <- pub[, .(
  geo_peak_id = peak_id,
  contrast,
  calderon_logFC = as.numeric(logFC),
  calderon_adjP = if ("adj.P.Val" %in% names(pub)) as.numeric(adj.P.Val) else NA_real_
)]

# Convert contrast like Bulk_B_S-Bulk_B_U into cell_type
pub[, cell_type := sub("_S-.*$", "", contrast)]

message("Reading reprocessed BAMPE DA files...")
da_files <- list.files(
  da_dir,
  pattern = "^DA_BAMPE_.*_S_vs_U\\.tsv$",
  full.names = TRUE
)

if (length(da_files) == 0) {
  stop("No BAMPE DA files found in: ", da_dir)
}

read_one_da <- function(f) {
  x <- fread(f)
  if (!all(c("peak_id", "logFC") %in% names(x))) {
    stop("DA file missing peak_id/logFC: ", f)
  }
  if (!"cell_type" %in% names(x)) {
    cell_type_from_file <- basename(f)
    cell_type_from_file <- sub("^DA_BAMPE_", "", cell_type_from_file)
    cell_type_from_file <- sub("_S_vs_U\\.tsv$", "", cell_type_from_file)
    x[, cell_type := cell_type_from_file]
  }
  x[, .(
    bampe_peak_id = peak_id,
    cell_type,
    bampe_logFC = as.numeric(logFC),
    bampe_adjP = if ("adj.P.Val" %in% names(x)) as.numeric(adj.P.Val) else NA_real_,
    bampe_significant = if ("significant" %in% names(x)) as.character(significant) else NA_character_
  )]
}

our_da <- rbindlist(lapply(da_files, read_one_da), use.names = TRUE, fill = TRUE)

message("Merging BAMPE DA with GEO peak matches...")
our_da_geo <- merge(our_da, pairs, by = "bampe_peak_id", allow.cartesian = TRUE)

message("Merging with Calderon published DA...")
matched <- merge(
  our_da_geo,
  pub,
  by = c("geo_peak_id", "cell_type"),
  allow.cartesian = TRUE
)

matched <- matched[!is.na(bampe_logFC) & !is.na(calderon_logFC)]

message("Matched DA rows: ", nrow(matched))
fwrite(matched, file.path(OUTDIR, "BAMPE_vs_Calderon_DA_logFC_matched_points.tsv"), sep = "\t")

# -------------------------------
# 3. DA LOGFC SCATTERPLOT - ALL POINTS
# -------------------------------

# To prevent huge overplotting, sample points for plotting only.
# Correlation is calculated on full matched table.
set.seed(1)
plot_n <- min(100000, nrow(matched))
matched_plot <- matched[sample(.N, plot_n)]

rho_all <- cor(matched$bampe_logFC, matched$calderon_logFC, method = "spearman", use = "complete.obs")
pearson_all <- cor(matched$bampe_logFC, matched$calderon_logFC, method = "pearson", use = "complete.obs")

p_scatter_all <- ggplot(matched_plot, aes(x = calderon_logFC, y = bampe_logFC)) +
  geom_point(alpha = 0.15, size = 0.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  labs(
    title = "BAMPE DA logFC concordance with Calderon",
    subtitle = paste0(
      "Matched published stimulation DA peaks; Spearman = ",
      round(rho_all, 3),
      ", Pearson = ",
      round(pearson_all, 3),
      "; plotted points = ",
      format(plot_n, big.mark = ",")
    ),
    x = "Calderon published DA logFC",
    y = "BAMPE limma-voom DA logFC"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(OUTDIR, "09_BAMPE_vs_Calderon_DA_logFC_scatter_all.png"),
       p_scatter_all, width = 7, height = 6, dpi = 300)

ggsave(file.path(OUTDIR, "09_BAMPE_vs_Calderon_DA_logFC_scatter_all.pdf"),
       p_scatter_all, width = 7, height = 6)

message("Saved all-point DA scatterplot.")

# -------------------------------
# 4. DA LOGFC SCATTERPLOT - FACETED BY CELL TYPE
# -------------------------------

# Downsample within each cell type for display
matched_plot_facet <- matched[, {
  n_keep <- min(.N, 4000)
  .SD[sample(.N, n_keep)]
}, by = cell_type]

cell_cor <- matched[, .(
  n = .N,
  spearman = cor(bampe_logFC, calderon_logFC, method = "spearman", use = "complete.obs")
), by = cell_type]

matched_plot_facet <- merge(matched_plot_facet, cell_cor, by = "cell_type")
matched_plot_facet[, facet_label := paste0(cell_type, "\nρ=", round(spearman, 3), ", n=", n)]

p_scatter_facet <- ggplot(matched_plot_facet, aes(x = calderon_logFC, y = bampe_logFC)) +
  geom_point(alpha = 0.18, size = 0.25) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.35) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2) +
  facet_wrap(~ facet_label, scales = "free", ncol = 4) +
  labs(
    title = "DA logFC concordance by cell type",
    subtitle = "BAMPE limma-voom vs Calderon published stimulation DA peaks",
    x = "Calderon published DA logFC",
    y = "BAMPE DA logFC"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(size = 7)
  )

ggsave(file.path(OUTDIR, "10_BAMPE_vs_Calderon_DA_logFC_scatter_by_cell_type.png"),
       p_scatter_facet, width = 12, height = 10, dpi = 300)

ggsave(file.path(OUTDIR, "10_BAMPE_vs_Calderon_DA_logFC_scatter_by_cell_type.pdf"),
       p_scatter_facet, width = 12, height = 10)

message("Saved faceted DA scatterplot.")

# -------------------------------
# 5. CELL TYPE x CELL TYPE SPEARMAN HEATMAP OF DA LOGFC
# -------------------------------
# This heatmap checks whether DA effect profiles are similar between cell types
# using BAMPE logFC values on shared BAMPE peaks.

message("Creating DA logFC cell-type Spearman heatmap...")

wide_da <- dcast(
  our_da[, .(bampe_peak_id, cell_type, bampe_logFC)],
  bampe_peak_id ~ cell_type,
  value.var = "bampe_logFC"
)

da_mat <- as.matrix(wide_da[, -1, with = FALSE])
rownames(da_mat) <- wide_da$bampe_peak_id

# Remove rows with too much missingness
keep <- rowSums(!is.na(da_mat)) >= 3
da_mat <- da_mat[keep, , drop = FALSE]

da_cor <- cor(da_mat, method = "spearman", use = "pairwise.complete.obs")

da_cor_dt <- as.data.table(as.table(da_cor))
setnames(da_cor_dt, c("cell_type_1", "cell_type_2", "spearman"))

ct_order <- colnames(da_cor)
da_cor_dt[, cell_type_1 := factor(cell_type_1, levels = ct_order)]
da_cor_dt[, cell_type_2 := factor(cell_type_2, levels = rev(ct_order))]

p_da_heat <- ggplot(da_cor_dt, aes(x = cell_type_1, y = cell_type_2, fill = spearman)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman"
  ) +
  labs(
    title = "Cell-type Spearman heatmap of BAMPE DA logFC profiles",
    subtitle = "Correlation of stimulation logFC profiles across cell types",
    x = "Cell type",
    y = "Cell type"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(OUTDIR, "11_BAMPE_DA_logFC_celltype_spearman_heatmap.png"),
       p_da_heat, width = 8, height = 7, dpi = 300)

ggsave(file.path(OUTDIR, "11_BAMPE_DA_logFC_celltype_spearman_heatmap.pdf"),
       p_da_heat, width = 8, height = 7)

message("Saved DA logFC cell-type Spearman heatmap.")

message("Done. Files saved in: ", OUTDIR)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

MOTIF_OUT <- args[1]
BAMPE_DIR <- args[2]
BAMPE_Q05_OUT <- args[3]

PLOT_DIR <- file.path(BAMPE_Q05_OUT, "plots_interpretation_qc")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

message("writing plots to: ", PLOT_DIR)

pick_col <- function(dt, candidates, required = TRUE) {
  hit <- intersect(candidates, names(dt))
  if (length(hit) == 0) {
    if (required) stop("Could not find any of these columns: ", paste(candidates, collapse = ", "))
    return(NA_character_)
  }
  hit[1]
}

count_bed_rows <- function(f) {
  as.integer(system(sprintf("wc -l < %s", shQuote(f)), intern = TRUE))
}

save_plot <- function(p, file, width = 9, height = 6) {
  ggsave(file, p, width = width, height = height, dpi = 300)
}

# ============================================================
# Files
# ============================================================

da_summary_file <- file.path(BAMPE_DIR, "differential_atac", "BAMPE_limma_voom_DA_summary.tsv")
general_summary_file <- file.path(MOTIF_OUT, "tables", "all_motif_general_annotation_summary.tsv")

peak_manifest_hg19 <- file.path(BAMPE_Q05_OUT, "tables", "q05_celltype_accessible_peak_files.hg19.tsv")
counts_by_motif_celltype_file <- file.path(BAMPE_Q05_OUT, "tables", "q05_accessible_counts_by_motif_celltype.tsv")
orientation_file <- file.path(BAMPE_Q05_OUT, "tables", "q05_gene_body_accessible_counts_by_orientation.tsv")
region_file <- file.path(BAMPE_Q05_OUT, "tables", "q05_genic_region_breakdown_all_rows.tsv")
tss_file <- file.path(BAMPE_Q05_OUT, "tables", "q05_tss_distance_bins_by_orientation.tsv")

# ============================================================
# 1. q05 significant peak counts per cell type
# ============================================================

if (file.exists(peak_manifest_hg19)) {
  message("making q05 significant peak count plots.")

  peak_manifest <- fread(peak_manifest_hg19)
  stopifnot(all(c("cell_type", "bed_file") %in% names(peak_manifest)))

  peak_manifest[, n_q05_significant_peaks := vapply(bed_file, count_bed_rows, integer(1))]

  fwrite(peak_manifest, file.path(PLOT_DIR, "table_q05_peak_counts_by_celltype.tsv"), sep = "\t")

  p <- ggplot(peak_manifest[order(n_q05_significant_peaks)], aes(
    x = reorder(cell_type, n_q05_significant_peaks),
    y = n_q05_significant_peaks
  )) +
    geom_col(fill = "grey35") +
    coord_flip() +
    theme_bw(base_size = 14) +
    labs(
      title = "q<0.05 significant ATAC peaks per cell type",
      x = "Cell type",
      y = "Number of significant peaks"
    )

  save_plot(p, file.path(PLOT_DIR, "01_q05_significant_peak_counts_per_celltype.png"), 8, 8)
}

# ============================================================
# 2. Differential ATAC summary plots
# ============================================================

if (file.exists(da_summary_file)) {
  message("making DA summary plots.")

  da <- fread(da_summary_file)

  da[, frac_significant := significant_peaks / tested_peaks]
  da[, frac_up_in_stim := up_in_stim / significant_peaks]
  da[, frac_down_in_stim := down_in_stim / significant_peaks]

  fwrite(da, file.path(PLOT_DIR, "table_bampe_da_summary_with_fractions.tsv"), sep = "\t")

  p1 <- ggplot(da[order(significant_peaks)], aes(
    x = reorder(cell_type, significant_peaks),
    y = significant_peaks
  )) +
    geom_col(fill = "grey35") +
    coord_flip() +
    theme_bw(base_size = 14) +
    labs(
      title = "Differentially accessible peaks per cell type",
      x = "Cell type",
      y = "Significant DA peaks"
    )

  save_plot(p1, file.path(PLOT_DIR, "02_da_significant_peak_counts_per_celltype.png"), 8, 8)

  p2 <- ggplot(da[order(frac_significant)], aes(
    x = reorder(cell_type, frac_significant),
    y = frac_significant
  )) +
    geom_col(fill = "grey35") +
    coord_flip() +
    theme_bw(base_size = 14) +
    labs(
      title = "Fraction of tested peaks that are DA",
      x = "Cell type",
      y = "DA peaks / tested peaks"
    )

  save_plot(p2, file.path(PLOT_DIR, "03_da_fraction_of_tested_peaks_per_celltype.png"), 8, 8)

  da_long <- melt(
    da[, .(cell_type, up_in_stim, down_in_stim)],
    id.vars = "cell_type",
    variable.name = "direction",
    value.name = "n_peaks"
  )

  p3 <- ggplot(da_long, aes(
    x = reorder(cell_type, n_peaks, FUN = sum),
    y = n_peaks,
    fill = direction
  )) +
    geom_col(position = "stack") +
    coord_flip() +
    theme_bw(base_size = 14) +
    labs(
      title = "Direction of differential accessibility",
      x = "Cell type",
      y = "Number of significant DA peaks"
    )

  save_plot(p3, file.path(PLOT_DIR, "04_da_up_vs_down_in_stim_per_celltype.png"), 8, 8)
}

# ============================================================
# 3. Read general motif summary
# ============================================================

if (!file.exists(general_summary_file)) {
  stop("Missing general motif summary file: ", general_summary_file)
}

general_summary <- fread(general_summary_file)

motif_col_general <- pick_col(general_summary, c("archetype", "motif_archetype", "motif"))
gene_body_total_col <- pick_col(general_summary, c(
  "protein_coding_or_lncRNA_rows",
  "rows_in_protein_coding_or_lncRNA",
  "gene_body_rows"
))

setnames(general_summary, motif_col_general, "motif")
setnames(general_summary, gene_body_total_col, "gene_body_total_rows")

general_summary <- general_summary[, .(motif, gene_body_total_rows)]

# ============================================================
# 4. Read motif × cell-type accessibility counts
# ============================================================

if (!file.exists(counts_by_motif_celltype_file)) {
  stop("Missing counts-by-motif-celltype file: ", counts_by_motif_celltype_file)
}

acc <- fread(counts_by_motif_celltype_file)

motif_col <- pick_col(acc, c("archetype", "motif", "motif_archetype"))
celltype_col <- pick_col(acc, c("cell_type"))
accessible_col <- pick_col(acc, c(
  "accessible_tfbs_count",
  "n_accessible",
  "accessible_rows",
  "n_accessible_rows"
), required = FALSE)

gene_body_accessible_col <- pick_col(acc, c(
  "gene_body_accessible_tfbs_count",
  "n_gene_body_accessible",
  "gene_body_accessible_rows",
  "n_gene_body_accessible_rows"
), required = FALSE)

if (is.na(gene_body_accessible_col) && !is.na(accessible_col)) {
  gene_body_accessible_col <- accessible_col
}

setnames(acc, motif_col, "motif")
setnames(acc, celltype_col, "cell_type")
setnames(acc, gene_body_accessible_col, "gene_body_accessible_rows")

acc <- merge(acc, general_summary, by = "motif", all.x = TRUE)

acc[, gene_body_accessible_fraction := gene_body_accessible_rows / gene_body_total_rows]

# ============================================================
# 5. Cell-type level plots
# ============================================================

celltype_summary <- acc[, .(
  mean_gene_body_accessible_rows = mean(gene_body_accessible_rows, na.rm = TRUE),
  median_gene_body_accessible_rows = median(gene_body_accessible_rows, na.rm = TRUE)
), by = cell_type]

fwrite(celltype_summary, file.path(PLOT_DIR, "table_celltype_accessibility_summary.tsv"), sep = "\t")

p <- ggplot(celltype_summary[order(mean_gene_body_accessible_rows)], aes(
  x = reorder(cell_type, mean_gene_body_accessible_rows),
  y = mean_gene_body_accessible_rows
)) +
  geom_col(fill = "grey35") +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "Mean gene-body accessible TFBS count per motif archetype",
    subtitle = "Computed separately within each cell type",
    x = "Cell type",
    y = "Mean gene-body accessible TFBS rows per motif"
  )

save_plot(p, file.path(PLOT_DIR, "05_mean_gene_body_accessible_tfbs_per_celltype.png"), 8, 8)

# ============================================================
# 6. Distribution excluding giant outlier motifs
# ============================================================

acc_no_big <- acc[!motif %in% c("GC-tract", "KLF_SP_2")]

p <- ggplot(acc_no_big, aes(x = gene_body_accessible_rows)) +
  geom_histogram(bins = 50, fill = "grey35", color = "white", linewidth = 0.2) +
  scale_x_log10() +
  theme_bw(base_size = 14) +
  labs(
    title = "Distribution of gene-body accessible TFBS counts",
    subtitle = "Excluding GC-tract and KLF_SP_2",
    x = "Gene-body accessible TFBS rows per motif-cell type pair, log10 scale",
    y = "Number of motif-cell type pairs"
  )

save_plot(p, file.path(PLOT_DIR, "06_distribution_gene_body_accessible_tfbs_excluding_GC_KLF.png"), 9, 6)

# ============================================================
# 7. Motif-level summaries
# ============================================================

motif_summary <- acc[, .(
  mean_gene_body_accessible_rows = mean(gene_body_accessible_rows, na.rm = TRUE),
  median_gene_body_accessible_rows = median(gene_body_accessible_rows, na.rm = TRUE),
  mean_gene_body_accessible_fraction = mean(gene_body_accessible_fraction, na.rm = TRUE)
), by = .(motif, gene_body_total_rows)]

motif_summary <- motif_summary[order(-mean_gene_body_accessible_rows)]

fwrite(motif_summary, file.path(PLOT_DIR, "table_motif_accessibility_summary.tsv"), sep = "\t")

top25 <- motif_summary[1:25]
top25[, motif := factor(motif, levels = rev(motif))]

p <- ggplot(top25, aes(x = motif, y = mean_gene_body_accessible_rows)) +
  geom_col(fill = "grey35") +
  coord_flip() +
  scale_y_log10() +
  theme_bw(base_size = 14) +
  labs(
    title = "Top 25 motifs by mean gene-body accessible TFBS count",
    x = "Motif archetype",
    y = "Mean accessible TFBS rows, log10 scale"
  )

save_plot(p, file.path(PLOT_DIR, "07_top25_motifs_by_mean_gene_body_accessible_count.png"), 9, 8)

p <- ggplot(motif_summary, aes(
  x = gene_body_total_rows,
  y = mean_gene_body_accessible_rows
)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw(base_size = 14) +
  labs(
    title = "Motif abundance vs accessible gene-body TFBS count",
    subtitle = "Checks whether large accessible counts are mostly driven by motif abundance",
    x = "Total gene-body TFBS rows per motif, log10 scale",
    y = "Mean accessible gene-body TFBS rows, log10 scale"
  )

save_plot(p, file.path(PLOT_DIR, "08_motif_abundance_vs_accessible_count_scatter.png"), 8, 6)

motif_summary_frac <- motif_summary[is.finite(mean_gene_body_accessible_fraction)]
p <- ggplot(motif_summary_frac, aes(
  x = gene_body_total_rows,
  y = mean_gene_body_accessible_fraction
)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_x_log10() +
  theme_bw(base_size = 14) +
  labs(
    title = "Motif abundance vs accessible fraction",
    subtitle = "Normalizes for motif abundance",
    x = "Total gene-body TFBS rows per motif, log10 scale",
    y = "Mean accessible fraction"
  )

save_plot(p, file.path(PLOT_DIR, "09_motif_abundance_vs_accessible_fraction_scatter.png"), 8, 6)

# ============================================================
# 8. Heatmap of motif profiles across cell types
# ============================================================

heat_dt <- dcast(
  acc[, .(motif, cell_type, gene_body_accessible_rows)],
  motif ~ cell_type,
  value.var = "gene_body_accessible_rows",
  fill = 0
)

heat_mat <- as.matrix(heat_dt[, -1])
rownames(heat_mat) <- heat_dt$motif

rv <- apply(heat_mat, 1, var)
top_n <- min(100, nrow(heat_mat))
heat_mat_top <- heat_mat[order(rv, decreasing = TRUE)[1:top_n], , drop = FALSE]

heat_long <- as.data.table(log10(heat_mat_top + 1), keep.rownames = "motif")
heat_long <- melt(
  heat_long,
  id.vars = "motif",
  variable.name = "cell_type",
  value.name = "log10_accessible_rows"
)

p <- ggplot(heat_long, aes(
  x = cell_type,
  y = motif,
  fill = log10_accessible_rows
)) +
  geom_tile() +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 5)
  ) +
  labs(
    title = "Top variable accessible motif profiles across cell types",
    x = "Cell type",
    y = "Motif archetype",
    fill = "log10(count + 1)"
  )

save_plot(
  p,
  file.path(PLOT_DIR, "10_heatmap_top_variable_motif_profiles_by_celltype.png"),
  width = 12,
  height = 14
)

# ============================================================
# 9. Genic region breakdown
# ============================================================

if (file.exists(region_file)) {
  region_dt <- fread(region_file)

  region_col <- pick_col(region_dt, c("genic_region_priority", "genic_region", "region"))
  n_col <- pick_col(region_dt, c("N", "count", "n"))

  setnames(region_dt, region_col, "region")
  setnames(region_dt, n_col, "N")

  region_summary <- region_dt[, .(N = sum(N, na.rm = TRUE)), by = region]
  region_summary[, fraction := N / sum(N)]

  fwrite(region_summary, file.path(PLOT_DIR, "table_genic_region_breakdown.tsv"), sep = "\t")

  p <- ggplot(region_summary[order(-N)], aes(
    x = reorder(region, N),
    y = N
  )) +
    geom_col(fill = "grey35") +
    coord_flip() +
    theme_bw(base_size = 14) +
    labs(
      title = "Genic region breakdown of annotated TFBS rows",
      x = "Genic region",
      y = "Number of rows"
    )

  save_plot(p, file.path(PLOT_DIR, "11_genic_region_breakdown.png"), 8, 6)
}

# ============================================================
# 10. Orientation balance
# ============================================================

if (file.exists(orientation_file)) {
  orient_dt <- fread(orientation_file)

  orient_col <- pick_col(orient_dt, c("motif_orientation", "orientation"))
  n_col <- pick_col(orient_dt, c("N", "count", "n"))
  motif_col2 <- pick_col(orient_dt, c("archetype", "motif", "motif_archetype"), required = FALSE)

  setnames(orient_dt, orient_col, "motif_orientation")
  setnames(orient_dt, n_col, "N")

  orient_summary <- orient_dt[, .(N = sum(N, na.rm = TRUE)), by = motif_orientation]

  fwrite(orient_summary, file.path(PLOT_DIR, "table_orientation_balance.tsv"), sep = "\t")

  p <- ggplot(orient_summary, aes(x = motif_orientation, y = N)) +
    geom_col(fill = "grey35") +
    theme_bw(base_size = 14) +
    labs(
      title = "Orientation balance among accessible gene-body TFBSs",
      x = "Motif orientation",
      y = "Number of rows"
    )

  save_plot(p, file.path(PLOT_DIR, "12_orientation_balance.png"), 7, 5)
}

# ============================================================
# 11. TSS distance by orientation
# ============================================================

if (file.exists(tss_file)) {
  tss_dt <- fread(tss_file)

  orient_col <- pick_col(tss_dt, c("motif_orientation", "orientation"))
  bin_col <- pick_col(tss_dt, c("tss_distance_bin", "distance_bin", "tss_bin"))
  n_col <- pick_col(tss_dt, c("N", "count", "n"))

  setnames(tss_dt, orient_col, "motif_orientation")
  setnames(tss_dt, bin_col, "tss_distance_bin")
  setnames(tss_dt, n_col, "N")

  tss_summary <- tss_dt[, .(N = sum(N, na.rm = TRUE)), by = .(motif_orientation, tss_distance_bin)]

  fwrite(tss_summary, file.path(PLOT_DIR, "table_tss_distance_by_orientation.tsv"), sep = "\t")

  p <- ggplot(tss_summary, aes(
    x = tss_distance_bin,
    y = N,
    fill = motif_orientation
  )) +
    geom_col(position = "dodge") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "TSS distance distribution by motif orientation",
      x = "Distance to nearest upstream TSS bin",
      y = "Number of rows"
    )

  save_plot(p, file.path(PLOT_DIR, "13_tss_distance_distribution_by_orientation.png"), 10, 6)
}

message("Done.")

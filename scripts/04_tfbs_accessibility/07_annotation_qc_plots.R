#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
BAMPE_Q05_OUT <- args[1]

PLOT_DIR <- file.path(BAMPE_Q05_OUT, "plots_interpretation_qc")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

region_file <- file.path(BAMPE_Q05_OUT, "tables", "q05_genic_region_breakdown_all_rows.tsv")
orientation_file <- file.path(BAMPE_Q05_OUT, "tables", "q05_gene_body_accessible_counts_by_orientation.tsv")
tss_file <- file.path(BAMPE_Q05_OUT, "tables", "q05_tss_distance_bins_by_orientation.tsv")

save_plot <- function(p, file, width = 8, height = 6) {
  ggsave(file, p, width = width, height = height, dpi = 300)
}

# ------------------------------------------------------------
# 11. Genic-region breakdown
# ------------------------------------------------------------

if (file.exists(region_file)) {
  region_dt <- fread(region_file)

  region_col <- intersect(c("genic_region_priority", "genic_region", "region", "region_class"), names(region_dt))[1]
  count_col <- intersect(c("n_rows", "N", "count", "n"), names(region_dt))[1]

  if (is.na(region_col)) stop("Could not find region column in region file.")

  if (!is.na(count_col)) {
    region_summary <- region_dt[, .(N = sum(get(count_col), na.rm = TRUE)), by = region_col]
  } else {
    region_summary <- region_dt[, .N, by = region_col]
  }

  setnames(region_summary, region_col, "genic_region")
  region_summary <- region_summary[!is.na(genic_region) & genic_region != ""]
  region_summary <- region_summary[order(-N)]
  region_summary[, fraction := N / sum(N)]

  fwrite(region_summary, file.path(PLOT_DIR, "table_genic_region_breakdown_for_plot.tsv"), sep = "\t")

  p <- ggplot(region_summary, aes(x = reorder(genic_region, N), y = N)) +
    geom_col(fill = "grey35") +
    coord_flip() +
    scale_y_log10() +
    theme_bw(base_size = 14) +
    labs(
      title = "Genic-region breakdown of annotated TFBS rows",
      subtitle = "Sanity check for gene annotation; intronic signal should usually dominate gene-body sequence",
      x = "Genic region",
      y = "TFBS rows, log10 scale"
    )

  save_plot(p, file.path(PLOT_DIR, "11_genic_region_breakdown.png"), 8, 6)
}

# ------------------------------------------------------------
# 12. Orientation balance
# ------------------------------------------------------------

if (file.exists(orientation_file)) {
  orient_dt <- fread(orientation_file)

  orient_col <- intersect(c("motif_orientation", "orientation"), names(orient_dt))[1]
  count_col <- intersect(c("gene_body_accessible_rows_orientation", "n_rows", "N", "count", "n"), names(orient_dt))[1]

  if (is.na(orient_col)) stop("Could not find orientation column in orientation file.")
  if (is.na(count_col)) stop("Could not find count column in orientation file.")

  orient_summary <- orient_dt[, .(N = sum(get(count_col), na.rm = TRUE)), by = orient_col]
  setnames(orient_summary, orient_col, "motif_orientation")

  fwrite(orient_summary, file.path(PLOT_DIR, "table_orientation_balance.tsv"), sep = "\t")

  p <- ggplot(orient_summary, aes(x = motif_orientation, y = N)) +
    geom_col(fill = "grey35") +
    theme_bw(base_size = 14) +
    labs(
      title = "Default vs reverse-complement orientation balance",
      subtitle = "Checks whether duplicated motif orientations are represented sensibly",
      x = "Motif orientation",
      y = "Accessible gene-body TFBS rows"
    )

  save_plot(p, file.path(PLOT_DIR, "12_orientation_balance.png"), 7, 5)
}

# ------------------------------------------------------------
# 13. TSS distance by orientation
# ------------------------------------------------------------

if (file.exists(tss_file)) {
  tss_dt <- fread(tss_file)

  orient_col <- intersect(c("motif_orientation", "orientation"), names(tss_dt))[1]
  bin_col <- intersect(c("tss_distance_bin", "distance_bin", "tss_bin"), names(tss_dt))[1]
  count_col <- intersect(c("n_rows", "N", "count", "n"), names(tss_dt))[1]

  if (is.na(orient_col)) stop("Could not find orientation column in TSS file.")
  if (is.na(bin_col)) stop("Could not find TSS distance bin column.")
  if (is.na(count_col)) stop("Could not find count column in TSS file.")

  tss_summary <- tss_dt[, .(N = sum(get(count_col), na.rm = TRUE)), by = c(orient_col, bin_col)]
  setnames(tss_summary, c(orient_col, bin_col), c("motif_orientation", "tss_distance_bin"))

  tss_summary[, tss_distance_bin := factor(
    tss_distance_bin,
    levels = c("0", "1-1kb", "1kb-10kb", "10kb-100kb", "100kb-1Mb", ">1Mb", "negative_or_other")
  )]

  fwrite(tss_summary, file.path(PLOT_DIR, "table_tss_distance_by_orientation.tsv"), sep = "\t")

  p <- ggplot(tss_summary, aes(x = tss_distance_bin, y = N, fill = motif_orientation)) +
    geom_col(position = "dodge") +
    scale_y_log10() +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "TSS-distance distribution by motif orientation",
      subtitle = "Strand-aware sanity check for upstream TSS annotation",
      x = "Distance to upstream TSS bin",
      y = "TFBS rows, log10 scale",
      fill = "Orientation"
    )

  save_plot(p, file.path(PLOT_DIR, "13_tss_distance_distribution_by_orientation.png"), 10, 6)
}

message("Done remaining plots.")

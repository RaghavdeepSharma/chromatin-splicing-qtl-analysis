suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("PROJECT_ROOT")
if (BASE == "") stop("PROJECT_ROOT is not set.")
RUN  <- file.path(BASE, "bampe_peak_validation_20260622")
OUT  <- file.path(RUN, "published_DA_correlation_v1")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

published_file <- file.path(BASE, "data/raw/Supplementary_data_3_ATAC_stimulation_DA_peaks.txt.gz")
pairs_file <- file.path(RUN, "correlation/bampe_to_geo_best_overlap_pairs.tsv")
our_da_dir <- file.path(RUN, "differential_atac")

if (!file.exists(published_file)) stop("Missing published DA file: ", published_file)
if (!file.exists(pairs_file)) stop("Missing overlap-pair file: ", pairs_file)

message("Reading best-overlap pairs...")
pairs <- fread(pairs_file, header = FALSE)
setnames(pairs, c("bam_peak_id", "geo_peak_id", "overlap_bp"))

message("Reading Calderon published stimulation ATAC DA...")
pub <- fread(cmd = paste("zcat", published_file))

needed <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "peak_id", "contrast")
missing <- setdiff(needed, names(pub))
if (length(missing) > 0) stop("Missing columns in published table: ", paste(missing, collapse = ", "))

pub <- pub[grepl("_S-.*_U$", contrast)]
pub[, published_cell_type := sub("_S-.*$", "", contrast)]

cell_types <- c(
  "Bulk_B",
  "CD8pos_T",
  "Central_memory_CD8pos_T",
  "Effector_CD4pos_T",
  "Effector_memory_CD8pos_T",
  "Follicular_T_Helper",
  "Gamma_delta_T",
  "Mature_NK",
  "Mem_B",
  "Memory_Teffs",
  "Memory_Tregs",
  "Monocytes",
  "Naive_B",
  "Naive_CD8_T",
  "Naive_Teffs",
  "Regulatory_T",
  "Th1_precursors",
  "Th17_precursors",
  "Th2_precursors"
)

results <- list()

for (ct in cell_types) {
  safe_ct <- gsub("[^A-Za-z0-9_]+", "_", ct)
  our_file <- file.path(our_da_dir, paste0("DA_BAMPE_", safe_ct, "_S_vs_U.tsv"))

  if (!file.exists(our_file)) {
    message("Skipping missing reprocessed DA file: ", our_file)
    next
  }

  message("Processing ", ct)

  our <- fread(our_file)
  our <- our[, .(
    bam_peak_id = peak_id,
    our_cell_type = cell_type,
    our_logFC = logFC,
    our_P.Value = P.Value,
    our_adj.P.Val = adj.P.Val,
    our_significant = significant
  )]

  our_mapped <- pairs[our, on = "bam_peak_id"]

  pub_ct <- pub[published_cell_type == ct, .(
    geo_peak_id = peak_id,
    published_contrast = contrast,
    published_logFC = logFC,
    published_P.Value = P.Value,
    published_adj.P.Val = adj.P.Val
  )]

  merged <- pub_ct[our_mapped, on = "geo_peak_id", nomatch = 0]
  merged <- merged[!is.na(our_logFC) & !is.na(published_logFC)]

  rho_all <- suppressWarnings(cor(
    merged$our_logFC,
    merged$published_logFC,
    method = "spearman",
    use = "complete.obs"
  ))

  pearson_all <- suppressWarnings(cor(
    merged$our_logFC,
    merged$published_logFC,
    method = "pearson",
    use = "complete.obs"
  ))

  both_sig <- merged[
    our_significant == TRUE &
      published_adj.P.Val < 0.01 &
      abs(published_logFC) > 1
  ]

  rho_both_sig <- if (nrow(both_sig) >= 100) {
    suppressWarnings(cor(
      both_sig$our_logFC,
      both_sig$published_logFC,
      method = "spearman",
      use = "complete.obs"
    ))
  } else {
    NA_real_
  }

  fwrite(
    merged,
    file.path(OUT, paste0("matched_DA_", safe_ct, "_BAMPE_vs_Calderon.tsv")),
    sep = "\t"
  )

  results[[ct]] <- data.table(
    cell_type = ct,
    n_matched_peaks = nrow(merged),
    spearman_logFC_all = rho_all,
    pearson_logFC_all = pearson_all,
    n_our_significant = sum(merged$our_significant == TRUE, na.rm = TRUE),
    n_published_significant = sum(
      merged$published_adj.P.Val < 0.01 &
        abs(merged$published_logFC) > 1,
      na.rm = TRUE
    ),
    n_both_significant = nrow(both_sig),
    spearman_logFC_both_significant = rho_both_sig
  )
}

res <- rbindlist(results, fill = TRUE)

fwrite(
  res,
  file.path(OUT, "BAMPE_vs_Calderon_published_DA_logFC_correlation_summary.tsv"),
  sep = "\t"
)

overall <- data.table(
  metric = c(
    "n_cell_types",
    "median_spearman_logFC_all",
    "mean_spearman_logFC_all",
    "min_spearman_logFC_all",
    "max_spearman_logFC_all"
  ),
  value = c(
    nrow(res),
    median(res$spearman_logFC_all, na.rm = TRUE),
    mean(res$spearman_logFC_all, na.rm = TRUE),
    min(res$spearman_logFC_all, na.rm = TRUE),
    max(res$spearman_logFC_all, na.rm = TRUE)
  )
)

fwrite(
  overall,
  file.path(OUT, "BAMPE_vs_Calderon_published_DA_overall_summary.tsv"),
  sep = "\t"
)

message("Done.")
print(res)
print(overall)

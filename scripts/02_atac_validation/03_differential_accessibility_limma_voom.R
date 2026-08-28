suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

BASE <- Sys.getenv("PROJECT_ROOT")
if (BASE == "") stop("PROJECT_ROOT is not set.")
RUN  <- file.path(BASE, "bampe_peak_validation_20260622")
OUT  <- file.path(RUN, "differential_atac")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

count_file <- file.path(RUN, "counts/bampe_consensus_counts.clean.tsv")
tss_file <- file.path(BASE, "qc/biosample_tss_enrichment.real.tsv")

message("Reading count matrix...")
counts_dt <- fread(count_file)
setnames(counts_dt, 1, "peak_id")

peak_ids <- counts_dt$peak_id
count_mat <- as.matrix(counts_dt[, -1, with = FALSE])
storage.mode(count_mat) <- "integer"
rownames(count_mat) <- peak_ids

samples <- colnames(count_mat)

parse_sample <- function(x) {
  donor <- sub("-.*$", "", x)
  condition <- ifelse(grepl("-no_trea", x), "U", "S")
  cell_type <- sub("^[^-]+-", "", x)
  cell_type <- sub("-no_treament$", "", cell_type)
  cell_type <- sub("-no_treatment$", "", cell_type)
  cell_type <- sub("-treatment[0-9]+$", "", cell_type)
  data.table(sample = x, donor = donor, cell_type = cell_type, condition = condition)
}

meta <- rbindlist(lapply(samples, parse_sample))

if (file.exists(tss_file)) {
  message("Reading TSS enrichment...")
  tss <- fread(tss_file)
  setnames(tss, 1, "sample")
  tss_col <- grep("tss", names(tss), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(tss_col)) {
    tss <- tss[, .(sample, TSS = get(tss_col))]
    meta <- merge(meta, tss, by = "sample", all.x = TRUE)
  } else {
    meta[, TSS := NA_real_]
  }
} else {
  meta[, TSS := NA_real_]
}

fwrite(meta, file.path(OUT, "bampe_sample_metadata.tsv"), sep = "\t")

cell_types <- sort(unique(meta$cell_type))
summary_rows <- list()

for (ct in cell_types) {
  m <- meta[cell_type == ct]
  
  donors_with_both <- m[, .N, by = .(donor, condition)][, .N, by = donor][N == 2, donor]
  m <- m[donor %in% donors_with_both]
  
  if (length(unique(m$donor)) < 3) {
    message("Skipping ", ct, ": fewer than 3 paired donors")
    next
  }
  
  idx <- match(m$sample, colnames(count_mat))
  mat <- count_mat[, idx, drop = FALSE]
  colnames(mat) <- m$sample
  
  # Calderon-style expressed/accessibility filter: >=1 CPM in at least two samples
  y <- DGEList(counts = mat)
  cpm_mat <- cpm(y)
  keep <- rowSums(cpm_mat >= 1) >= 2
  
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  
  m[, condition := factor(condition, levels = c("U", "S"))]
  m[, donor := factor(donor)]
  
  # Include TSS if available and non-missing for this subset
  use_tss <- "TSS" %in% names(m) && all(!is.na(m$TSS)) && length(unique(m$TSS)) > 1
  
  if (use_tss) {
    design <- model.matrix(~ donor + TSS + condition, data = m)
  } else {
    design <- model.matrix(~ donor + condition, data = m)
  }
  
  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  coef_name <- grep("conditionS", colnames(design), value = TRUE)
  if (length(coef_name) != 1) {
    stop("Could not find conditionS coefficient for ", ct)
  }
  
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  tt <- data.table(
    peak_id = rownames(tt),
    cell_type = ct,
    n_samples = nrow(m),
    n_donors = length(unique(m$donor)),
    logFC = tt$logFC,
    AveExpr = tt$AveExpr,
    t = tt$t,
    P.Value = tt$P.Value,
    adj.P.Val = tt$adj.P.Val,
    B = tt$B
  )
  
  tt[, significant := adj.P.Val < 0.01 & abs(logFC) > 1]
  
  safe_ct <- gsub("[^A-Za-z0-9_]+", "_", ct)
  out_file <- file.path(OUT, paste0("DA_BAMPE_", safe_ct, "_S_vs_U.tsv"))
  fwrite(tt, out_file, sep = "\t")
  
  summary_rows[[ct]] <- data.table(
    cell_type = ct,
    n_samples = nrow(m),
    n_donors = length(unique(m$donor)),
    tested_peaks = nrow(tt),
    significant_peaks = sum(tt$significant),
    up_in_stim = sum(tt$significant & tt$logFC > 0),
    down_in_stim = sum(tt$significant & tt$logFC < 0),
    used_TSS_covariate = use_tss
  )
  
  message("Finished ", ct, ": ", sum(tt$significant), " significant peaks")
}

summary <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary, file.path(OUT, "BAMPE_limma_voom_DA_summary.tsv"), sep = "\t")

message("Done.")
print(summary)

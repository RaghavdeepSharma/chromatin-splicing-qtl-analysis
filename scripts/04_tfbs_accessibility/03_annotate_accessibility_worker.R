library(data.table)
library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)

motif_file <- args[1]
peak_manifest <- args[2]
outdir <- args[3]

dir.create(file.path(outdir, "annotated_with_accessibility"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "tables"), showWarnings = FALSE, recursive = TRUE)

message("Reading motif annotation: ", motif_file)
dt <- fread(cmd = paste("gzip -dc", shQuote(motif_file)))

required_cols <- c("chr", "start", "end")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required motif columns: ", paste(missing_cols, collapse = ", "))
}

motif_name <- basename(motif_file)
motif_name <- sub("\\.general_annotated\\.tsv\\.gz$", "", motif_name)

message("Annotating accessibility for motif: ", motif_name)

# The motif table keeps BED-style 0-based start and 1-based end.
# Convert BED starts to 1-based closed coordinates for GenomicRanges.
motif_gr <- GRanges(
  seqnames = dt$chr,
  ranges = IRanges(start = dt$start + 1, end = dt$end)
)

manifest <- fread(peak_manifest)

if (!all(c("cell_type", "bed_file") %in% names(manifest))) {
  stop("Peak manifest must have columns: cell_type and bed_file")
}

safe_colname <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+$", "", x)
  x <- gsub("^_+", "", x)
  paste0("accessible_", x)
}

summary_rows <- list()

for (i in seq_len(nrow(manifest))) {
  cell_type <- manifest$cell_type[i]
  bed_file <- manifest$bed_file[i]
  colname <- safe_colname(cell_type)

  message("Reading accessible peaks for cell type: ", cell_type)

  if (!file.exists(bed_file)) {
    stop("Cannot find peak BED file: ", bed_file)
  }

  peaks <- fread(bed_file, header = FALSE)
  if (ncol(peaks) < 3) {
    stop("Peak BED has fewer than 3 columns: ", bed_file)
  }

  setnames(peaks, 1:3, c("chr", "start", "end"))

  peak_gr <- GRanges(
    seqnames = peaks$chr,
    ranges = IRanges(start = peaks$start + 1, end = peaks$end)
  )

  hits <- findOverlaps(motif_gr, peak_gr, ignore.strand = TRUE)

  acc <- integer(nrow(dt))
  if (length(hits) > 0) {
    acc[unique(queryHits(hits))] <- 1L
  }

  dt[, (colname) := acc]

  summary_rows[[i]] <- data.table(
    motif = motif_name,
    cell_type = cell_type,
    peak_file = bed_file,
    n_accessible_peaks = length(peak_gr),
    n_tfbs_rows = nrow(dt),
    n_accessible_tfbs_rows = sum(acc),
    fraction_accessible_tfbs_rows = sum(acc) / nrow(dt),
    accessibility_column = colname
  )
}

out_file <- file.path(outdir, "annotated_with_accessibility", paste0(motif_name, ".general_annotated.celltype_accessibility.tsv.gz"))
summary_file <- file.path(outdir, "tables", paste0(motif_name, ".celltype_accessibility_summary.tsv"))

message("Writing annotated output: ", out_file)
fwrite(dt, out_file, sep = "\t")

summary_dt <- rbindlist(summary_rows)
fwrite(summary_dt, summary_file, sep = "\t")

message("Summary written: ", summary_file)
message("Completed motif: ", motif_name)

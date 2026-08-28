#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

base_dir <- Sys.getenv("KERIMOV_BASE")

if (base_dir == "") {
  stop("KERIMOV_BASE is not set.")
}

input_file <- file.path(
  base_dir,
  "05_processed",
  "mapped_sqtl_bed_manifest.tsv"
)

output_file <- file.path(
  base_dir,
  "02_manifests",
  "sqtl_atac_overlap_manifest.tsv"
)

manifest <- read.delim(
  input_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

required <- c(
  "dataset_id",
  "calderon_celltype",
  "match_quality",
  "include_strict_primary",
  "include_expanded_primary",
  "include_sensitivity",
  "high_pip_bed",
  "low_pip_bed",
  "calderon_peak_file"
)

missing_columns <- setdiff(required, colnames(manifest))

if (length(missing_columns) > 0) {
  stop(
    "Missing manifest columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

overlap_manifest <- manifest[, required]

files_to_check <- c(
  overlap_manifest$high_pip_bed,
  overlap_manifest$low_pip_bed,
  overlap_manifest$calderon_peak_file
)

missing_files <- unique(
  files_to_check[!file.exists(files_to_check)]
)

if (length(missing_files) > 0) {
  stop(
    "Missing input files:\n",
    paste(missing_files, collapse = "\n")
  )
}

if (anyDuplicated(overlap_manifest$dataset_id)) {
  stop("Duplicate dataset IDs detected.")
}

write.table(
  overlap_manifest,
  output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

cat("ATAC-overlap manifest created.\n")
cat("Datasets: ", nrow(overlap_manifest), "\n", sep = "")
cat("Output: ", output_file, "\n", sep = "")

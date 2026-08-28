#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

base_dir <- Sys.getenv("KERIMOV_BASE")

if (base_dir == "") {
  stop("KERIMOV_BASE is not set. Source config.sh first.")
}

inventory_file <- file.path(
  base_dir,
  "05_processed",
  "immune_leafcutter_dataset_inventory.tsv"
)

manifest_dir <- file.path(base_dir, "02_manifests")
raw_dir      <- file.path(base_dir, "03_raw_susie")
results_dir  <- file.path(base_dir, "06_results")

dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

inventory <- read.delim(
  inventory_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  check.names = FALSE
)

required <- c(
  "study_id",
  "dataset_id",
  "study_label",
  "sample_group",
  "tissue_label",
  "condition_label",
  "sample_size",
  "quant_method",
  "preliminary_cell_class",
  "calderon_priority"
)

missing <- setdiff(required, names(inventory))

if (length(missing) > 0) {
  stop("Missing columns: ", paste(missing, collapse = ", "))
}

# Download direct matches plus broad lineage matches.
manifest <- inventory[
  inventory$calderon_priority %in%
    c("priority", "broad_lineage_match"),
]

manifest$download_tier <- ifelse(
  manifest$calderon_priority == "priority",
  "primary",
  "secondary_lineage"
)

# Use sample_group to infer Treg state because tissue_label appears
# reversed for some Treg metadata records.
manifest$state_interpretation <- manifest$condition_label

is_treg <- grepl(
  "^Treg_(naive|memory)$",
  manifest$sample_group,
  ignore.case = TRUE
)

manifest$state_interpretation[is_treg] <- sub(
  "^Treg_",
  "",
  manifest$sample_group[is_treg],
  ignore.case = TRUE
)

manifest$manual_note <- ""

manifest$manual_note[is_treg] <-
  "Treg state inferred from sample_group; tissue_label appears inconsistent"

manifest$manual_note[
  manifest$calderon_priority == "broad_lineage_match"
] <- paste0(
  manifest$manual_note[
    manifest$calderon_priority == "broad_lineage_match"
  ],
  ifelse(
    manifest$manual_note[
      manifest$calderon_priority == "broad_lineage_match"
    ] == "",
    "",
    "; "
  ),
  "Broad lineage match only; do not treat as an exact Calderon cell type"
)

ftp_root <- "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie"

manifest$url <- paste0(
  ftp_root,
  "/",
  manifest$study_id,
  "/",
  manifest$dataset_id,
  "/",
  manifest$dataset_id,
  ".credible_sets.tsv.gz"
)

manifest$local_file <- file.path(
  raw_dir,
  manifest$study_id,
  manifest$dataset_id,
  paste0(manifest$dataset_id, ".credible_sets.tsv.gz")
)

manifest <- manifest[
  order(
    manifest$download_tier,
    manifest$study_label,
    manifest$sample_group
  ),
]

manifest_columns <- c(
  "study_id",
  "dataset_id",
  "study_label",
  "sample_group",
  "tissue_label",
  "condition_label",
  "state_interpretation",
  "sample_size",
  "preliminary_cell_class",
  "calderon_priority",
  "download_tier",
  "manual_note",
  "url",
  "local_file"
)

manifest <- manifest[, manifest_columns]

manifest_file <- file.path(
  manifest_dir,
  "susie_credible_sets_manifest.tsv"
)

write.table(
  manifest,
  manifest_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

summary_table <- as.data.frame(
  table(
    manifest$download_tier,
    manifest$preliminary_cell_class
  )
)

summary_table <- summary_table[summary_table$Freq > 0, ]

names(summary_table) <- c(
  "download_tier",
  "cell_class",
  "n_datasets"
)

write.table(
  summary_table,
  file.path(results_dir, "susie_download_manifest_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nSuSiE download manifest created.\n")
cat("Total datasets:     ", nrow(manifest), "\n", sep = "")
cat(
  "Primary datasets:   ",
  sum(manifest$download_tier == "primary"),
  "\n",
  sep = ""
)
cat(
  "Secondary datasets: ",
  sum(manifest$download_tier == "secondary_lineage"),
  "\n",
  sep = ""
)

cat("\nManifest:\n")
cat(manifest_file, "\n")

cat("\nCounts by download tier and cell class:\n")
print(summary_table, row.names = FALSE)

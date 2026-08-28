#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

base_dir <- Sys.getenv("KERIMOV_BASE")

if (base_dir == "") {
  stop("KERIMOV_BASE is not defined. Source config.sh first.")
}

metadata_file <- file.path(
  base_dir,
  "01_metadata",
  "dataset_metadata.tsv"
)

processed_dir <- file.path(base_dir, "05_processed")
results_dir   <- file.path(base_dir, "06_results")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

metadata <- read.delim(
  metadata_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  check.names = FALSE
)

required_columns <- c(
  "study_id",
  "dataset_id",
  "study_label",
  "sample_group",
  "tissue_id",
  "tissue_label",
  "condition_label",
  "sample_size",
  "quant_method",
  "pmid",
  "study_type"
)

missing_columns <- setdiff(required_columns, colnames(metadata))

if (length(missing_columns) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Keep one row per bulk LeafCutter dataset.
bulk_leafcutter <- metadata[
  tolower(metadata$quant_method) == "leafcutter" &
    tolower(metadata$study_type) == "bulk",
]

# Search both the tissue name and sample-group description.
search_text <- tolower(
  paste(
    bulk_leafcutter$study_label,
    bulk_leafcutter$sample_group,
    bulk_leafcutter$tissue_label,
    bulk_leafcutter$condition_label
  )
)

immune_pattern <- paste(
  c(
    "immune",
    "blood",
    "lymph",
    "monocyte",
    "macrophage",
    "neutrophil",
    "eosinophil",
    "basophil",
    "dendritic",
    "microglia",
    "mast.cell",
    "b.cell",
    "b_cell",
    "bcell",
    "plasmablast",
    "plasma.cell",
    "t.cell",
    "t_cell",
    "tcell",
    "cd4",
    "cd8",
    "treg",
    "regulatory.t",
    "th1",
    "th2",
    "th17",
    "tfh",
    "follicular.helper",
    "gamma.delta",
    "natural.killer",
    "nk.cell",
    "nk_cell",
    "lcl",
    "lymphoblast"
  ),
  collapse = "|"
)

bulk_leafcutter$is_immune_candidate <- grepl(
  immune_pattern,
  search_text,
  perl = TRUE
)

immune_candidates <- bulk_leafcutter[
  bulk_leafcutter$is_immune_candidate,
]

# Preliminary classification used only to help review the inventory.
classify_cell <- function(text) {
  text <- tolower(text)

  if (grepl("macrophage", text)) return("macrophage")
  if (grepl("monocyte", text)) return("monocyte")
  if (grepl("neutrophil", text)) return("neutrophil")
  if (grepl("dendritic", text)) return("dendritic_cell")
  if (grepl("microglia", text)) return("microglia")
  if (grepl("plasmablast|plasma.cell", text)) return("plasmablast")
  if (grepl("treg|regulatory.t", text)) return("regulatory_T_cell")
  if (grepl("tfh|follicular.helper", text)) return("follicular_helper_T_cell")
  if (grepl("th17", text)) return("TH17_cell")
  if (grepl("th2", text)) return("TH2_cell")
  if (grepl("th1", text)) return("TH1_cell")
  if (grepl("gamma.delta", text)) return("gamma_delta_T_cell")
  if (grepl("cd8", text)) return("CD8_T_cell")
  if (grepl("cd4", text)) return("CD4_T_cell")
  if (grepl("natural.killer|nk.cell|nk_cell", text)) return("NK_cell")
  if (grepl("b.cell|b_cell|bcell", text)) return("B_cell")
  if (grepl("lcl|lymphoblast", text)) return("LCL")
  if (grepl("whole.blood|blood", text)) return("blood")
  if (grepl("t.cell|t_cell|tcell", text)) return("T_cell")

  return("other_immune")
}

immune_text <- paste(
  immune_candidates$sample_group,
  immune_candidates$tissue_label,
  immune_candidates$condition_label
)

immune_candidates$preliminary_cell_class <- vapply(
  immune_text,
  classify_cell,
  character(1)
)

# Mark likely relevance to the Calderon immune-cell ATAC resource.
immune_candidates$calderon_priority <- ifelse(
  immune_candidates$preliminary_cell_class %in%
    c(
      "monocyte",
      "neutrophil",
      "dendritic_cell",
      "plasmablast",
      "regulatory_T_cell",
      "follicular_helper_T_cell",
      "TH17_cell",
      "TH2_cell",
      "TH1_cell",
      "gamma_delta_T_cell",
      "CD8_T_cell",
      "CD4_T_cell",
      "NK_cell",
      "B_cell",
      "T_cell"
    ),
  "priority",
  ifelse(
    immune_candidates$preliminary_cell_class == "macrophage",
    "broad_lineage_match",
    "no_direct_match"
  )
)

# Sort inventory for easier review.
immune_candidates <- immune_candidates[
  order(
    immune_candidates$study_label,
    immune_candidates$preliminary_cell_class,
    immune_candidates$condition_label
  ),
]

inventory_file <- file.path(
  processed_dir,
  "immune_leafcutter_dataset_inventory.tsv"
)

all_leafcutter_file <- file.path(
  processed_dir,
  "all_bulk_leafcutter_datasets.tsv"
)

write.table(
  bulk_leafcutter,
  all_leafcutter_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

write.table(
  immune_candidates,
  inventory_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

# Summary by study.
study_summary <- aggregate(
  dataset_id ~ study_id + study_label,
  data = immune_candidates,
  FUN = length
)

colnames(study_summary)[
  colnames(study_summary) == "dataset_id"
] <- "n_leafcutter_datasets"

study_summary <- study_summary[
  order(-study_summary$n_leafcutter_datasets),
]

write.table(
  study_summary,
  file.path(results_dir, "immune_leafcutter_summary_by_study.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Distinct combinations needed for Calderon mapping.
mapping_review <- unique(
  immune_candidates[
    ,
    c(
      "study_id",
      "study_label",
      "sample_group",
      "tissue_label",
      "condition_label",
      "preliminary_cell_class",
      "calderon_priority"
    )
  ]
)

mapping_review <- mapping_review[
  order(
    mapping_review$preliminary_cell_class,
    mapping_review$study_label,
    mapping_review$condition_label
  ),
]

write.table(
  mapping_review,
  file.path(results_dir, "immune_cell_condition_review.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nCompleted immune LeafCutter inventory.\n\n")
cat("All bulk LeafCutter datasets: ", nrow(bulk_leafcutter), "\n", sep = "")
cat("Immune candidates:            ", nrow(immune_candidates), "\n", sep = "")
cat("Unique studies:               ", length(unique(immune_candidates$study_id)), "\n", sep = "")
cat("Unique cell classes:          ", length(unique(immune_candidates$preliminary_cell_class)), "\n", sep = "")

cat("\nDataset counts by preliminary cell class:\n")
print(
  sort(
    table(immune_candidates$preliminary_cell_class),
    decreasing = TRUE
  )
)

cat("\nDataset counts by study:\n")
print(study_summary, row.names = FALSE)

cat("\nFiles written:\n")
cat(all_leafcutter_file, "\n")
cat(inventory_file, "\n")
cat(file.path(results_dir, "immune_leafcutter_summary_by_study.tsv"), "\n")
cat(file.path(results_dir, "immune_cell_condition_review.tsv"), "\n")

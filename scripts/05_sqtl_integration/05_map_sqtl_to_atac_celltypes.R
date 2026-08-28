#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

base_dir <- Sys.getenv("KERIMOV_BASE")

if (base_dir == "") {
  stop("KERIMOV_BASE is not set. Source config.sh first.")
}

manifest_file <- file.path(
  base_dir,
  "02_manifests",
  "susie_credible_sets_manifest.tsv"
)

inventory_file <- file.path(
  base_dir,
  "01_metadata",
  "approved_bampe_q05_hg38_peak_inventory.tsv"
)

processed_dir <- file.path(base_dir, "05_processed")
results_dir   <- file.path(base_dir, "06_results")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.delim(
  manifest_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

inventory <- read.delim(
  inventory_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

required_manifest <- c(
  "study_id",
  "dataset_id",
  "study_label",
  "sample_group",
  "tissue_label",
  "condition_label",
  "state_interpretation",
  "sample_size",
  "preliminary_cell_class"
)

required_inventory <- c(
  "celltype_label",
  "peak_file",
  "n_rows",
  "coordinate_check",
  "use_for_analysis"
)

missing_manifest <- setdiff(
  required_manifest,
  colnames(manifest)
)

missing_inventory <- setdiff(
  required_inventory,
  colnames(inventory)
)

if (length(missing_manifest) > 0) {
  stop(
    "Manifest missing columns: ",
    paste(missing_manifest, collapse = ", ")
  )
}

if (length(missing_inventory) > 0) {
  stop(
    "Inventory missing columns: ",
    paste(missing_inventory, collapse = ", ")
  )
}

usable_inventory <- inventory[
  inventory$coordinate_check == "PASS" &
    inventory$use_for_analysis == "YES",
]

if (anyDuplicated(usable_inventory$celltype_label)) {
  stop("Duplicate Calderon cell-type labels found in inventory.")
}

# ------------------------------------------------------------------
# Curated mapping specification
#
# Match levels:
# exact_cell_type = direct equivalent
# close_subtype   = biologically close but subtype is not identical
# broad_lineage   = lineage-level approximation only
# no_match        = no defensible Calderon BED
#
# The approved final BED files are cell-type level and do not retain
# resting/stimulated condition as separate BED files.
# ------------------------------------------------------------------

mapping_spec <- data.frame(
  dataset_id = c(
    "QTD000025",
    "QTD000030",
    "QTD000035",
    "QTD000040",
    "QTD000109",
    "QTD000433",
    "QTD000418",
    "QTD000413",
    "QTD000423",
    "QTD000428",
    "QTD000478",
    "QTD000488",
    "QTD000483",
    "QTD000498",
    "QTD000493",
    "QTD000503",
    "QTD000508",
    "QTD000513",
    "QTD000443",
    "QTD000453",
    "QTD000463",
    "QTD000448",
    "QTD000458",
    "QTD000468",
    "QTD000473",
    "QTD000010",
    "QTD000020",
    "QTD000005",
    "QTD000015",
    "QTD000383",
    "QTD000388",
    "QTD000393"
  ),

  calderon_celltype = c(
    "Monocytes",
    NA,
    "Naive_Teffs",
    "Naive_Tregs",
    "Naive_Teffs",
    "Monocytes",
    "Monocytes",
    "Monocytes",
    "Monocytes",
    "Monocytes",
    "Naive_B",
    "Effector_CD4pos_T",
    "Naive_Teffs",
    "CD8pos_T",
    "Naive_CD8_T",
    "Monocytes",
    "Monocytes",
    "Mature_NK",
    "Follicular_T_Helper",
    "Th1_precursors",
    "Th17_precursors",
    "Th17_precursors",
    "Th2_precursors",
    "Memory_Tregs",
    "Naive_Tregs",
    "Monocytes",
    "Monocytes",
    "Monocytes",
    "Monocytes",
    "Monocytes",
    "Monocytes",
    "Monocytes"
  ),

  match_quality = c(
    "exact_cell_type",
    "no_match",
    "close_subtype",
    "exact_cell_type",
    "broad_lineage",
    "exact_cell_type",
    "exact_cell_type",
    "exact_cell_type",
    "exact_cell_type",
    "exact_cell_type",
    "exact_cell_type",
    "close_subtype",
    "close_subtype",
    "close_subtype",
    "exact_cell_type",
    "close_subtype",
    "exact_cell_type",
    "close_subtype",
    "exact_cell_type",
    "broad_lineage",
    "broad_lineage",
    "broad_lineage",
    "broad_lineage",
    "exact_cell_type",
    "exact_cell_type",
    "broad_lineage",
    "broad_lineage",
    "broad_lineage",
    "broad_lineage",
    "broad_lineage",
    "broad_lineage",
    "broad_lineage"
  ),

  mapping_reason = c(
    "Direct monocyte match",
    "No neutrophil BED in approved Calderon pipeline",
    "Generic CD4 T-cell mapped to closest naive CD4 effector population",
    "Naive Treg mapped using sample_group rather than inconsistent tissue_label",
    "Generic T-cell population lacks CD4/CD8 subtype; exploratory CD4 effector mapping",
    "Direct monocyte cell-type match; IAV condition is not separated in final BED",
    "Direct monocyte cell-type match; LPS condition is not separated in final BED",
    "Direct resting monocyte match",
    "Direct monocyte cell-type match; Pam3CSK4 condition is not separated in final BED",
    "Direct monocyte cell-type match; R848 condition is not separated in final BED",
    "Direct naive B-cell match",
    "Activated generic CD4 T cells mapped to closest effector CD4 population",
    "Naive generic CD4 T cells mapped to naive CD4 effector population",
    "Activated generic CD8 T cells mapped to bulk CD8-positive T population",
    "Direct naive CD8 T-cell match",
    "CD16-positive monocytes mapped to broader Calderon monocyte population",
    "Direct resting monocyte match",
    "Peripheral NK population mapped to mature NK population",
    "Direct follicular helper T-cell match; memory state not separated in final BED",
    "Memory Th1 cells mapped to available Th1 precursor population",
    "Memory Th1/17 cells mapped to available Th17 precursor population",
    "Memory Th17 cells mapped to available Th17 precursor population",
    "Memory Th2 cells mapped to available Th2 precursor population",
    "Direct memory Treg match using sample_group",
    "Direct naive Treg match using sample_group",
    "Macrophage mapped to monocyte lineage for sensitivity analysis only",
    "Macrophage mapped to monocyte lineage for sensitivity analysis only",
    "Macrophage mapped to monocyte lineage for sensitivity analysis only",
    "Macrophage mapped to monocyte lineage for sensitivity analysis only",
    "Macrophage mapped to monocyte lineage for sensitivity analysis only",
    "Macrophage mapped to monocyte lineage for sensitivity analysis only",
    "Macrophage mapped to monocyte lineage for sensitivity analysis only"
  ),

  stringsAsFactors = FALSE
)

# Validate one mapping row per manifest dataset.
missing_mapping <- setdiff(
  manifest$dataset_id,
  mapping_spec$dataset_id
)

extra_mapping <- setdiff(
  mapping_spec$dataset_id,
  manifest$dataset_id
)

if (length(missing_mapping) > 0) {
  stop(
    "Manifest datasets missing from mapping: ",
    paste(missing_mapping, collapse = ", ")
  )
}

if (length(extra_mapping) > 0) {
  stop(
    "Mapping contains unexpected datasets: ",
    paste(extra_mapping, collapse = ", ")
  )
}

if (anyDuplicated(mapping_spec$dataset_id)) {
  stop("Duplicate dataset IDs in mapping specification.")
}

# Preserve manifest order.
mapping_index <- match(
  manifest$dataset_id,
  mapping_spec$dataset_id
)

mapping <- cbind(
  manifest[
    ,
    c(
      "study_id",
      "dataset_id",
      "study_label",
      "sample_group",
      "tissue_label",
      "condition_label",
      "state_interpretation",
      "sample_size",
      "preliminary_cell_class"
    )
  ],
  mapping_spec[
    mapping_index,
    c(
      "calderon_celltype",
      "match_quality",
      "mapping_reason"
    )
  ]
)

mapping$calderon_condition_scope <- ifelse(
  is.na(mapping$calderon_celltype),
  "not_applicable",
  "cell_type_level_BED_condition_collapsed"
)

# Three increasingly permissive analysis definitions.
mapping$include_strict_primary <- ifelse(
  mapping$match_quality == "exact_cell_type",
  "YES",
  "NO"
)

mapping$include_expanded_primary <- ifelse(
  mapping$match_quality %in%
    c("exact_cell_type", "close_subtype"),
  "YES",
  "NO"
)

mapping$include_sensitivity <- ifelse(
  mapping$match_quality != "no_match",
  "YES",
  "NO"
)

# Attach approved peak-file information.
inventory_index <- match(
  mapping$calderon_celltype,
  usable_inventory$celltype_label
)

mapping$calderon_peak_file <- usable_inventory$peak_file[
  inventory_index
]

mapping$calderon_peak_count <- usable_inventory$n_rows[
  inventory_index
]

# Validate all non-null mappings resolve to approved files.
mapped_rows <- !is.na(mapping$calderon_celltype)

if (any(is.na(mapping$calderon_peak_file[mapped_rows]))) {
  bad_labels <- unique(
    mapping$calderon_celltype[
      mapped_rows &
        is.na(mapping$calderon_peak_file)
    ]
  )

  stop(
    "Mapped Calderon labels missing approved BED files: ",
    paste(bad_labels, collapse = ", ")
  )
}

if (any(!file.exists(mapping$calderon_peak_file[mapped_rows]))) {
  bad_files <- mapping$calderon_peak_file[
    mapped_rows &
      !file.exists(mapping$calderon_peak_file)
  ]

  stop(
    "Mapped peak files do not exist: ",
    paste(bad_files, collapse = ", ")
  )
}

mapping$mapping_status <- ifelse(
  mapping$match_quality == "no_match",
  "EXCLUDED_NO_APPROVED_CELLTYPE",
  "MAPPED_TO_APPROVED_BAMPE_HG38"
)

mapping_output <- file.path(
  processed_dir,
  "calderon_sqtl_dataset_mapping.tsv"
)

write.table(
  mapping,
  mapping_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

# Summaries.
quality_summary <- as.data.frame(
  table(mapping$match_quality),
  stringsAsFactors = FALSE
)

colnames(quality_summary) <- c(
  "match_quality",
  "n_datasets"
)

quality_summary <- quality_summary[
  order(-quality_summary$n_datasets),
]

tier_summary <- data.frame(
  analysis_set = c(
    "strict_primary",
    "expanded_primary",
    "sensitivity"
  ),
  n_datasets = c(
    sum(mapping$include_strict_primary == "YES"),
    sum(mapping$include_expanded_primary == "YES"),
    sum(mapping$include_sensitivity == "YES")
  )
)

calderon_summary <- aggregate(
  dataset_id ~ calderon_celltype,
  data = mapping[!is.na(mapping$calderon_celltype), ],
  FUN = length
)

colnames(calderon_summary)[
  colnames(calderon_summary) == "dataset_id"
] <- "n_sqtl_datasets"

calderon_summary <- calderon_summary[
  order(
    -calderon_summary$n_sqtl_datasets,
    calderon_summary$calderon_celltype
  ),
]

write.table(
  quality_summary,
  file.path(
    results_dir,
    "calderon_mapping_summary_by_quality.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  tier_summary,
  file.path(
    results_dir,
    "calderon_mapping_summary_by_analysis_set.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  calderon_summary,
  file.path(
    results_dir,
    "calderon_mapping_summary_by_celltype.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nCalderon–sQTL mapping completed.\n\n")

cat("Datasets in manifest:       ", nrow(manifest), "\n", sep = "")
cat(
  "Datasets mapped to a BED:  ",
  sum(!is.na(mapping$calderon_peak_file)),
  "\n",
  sep = ""
)
cat(
  "Datasets without a match:  ",
  sum(mapping$match_quality == "no_match"),
  "\n",
  sep = ""
)

cat("\nMapping quality:\n")
print(quality_summary, row.names = FALSE)

cat("\nAnalysis sets:\n")
print(tier_summary, row.names = FALSE)

cat("\nDatasets per Calderon cell type:\n")
print(calderon_summary, row.names = FALSE)

cat("\nOutput:\n")
cat(mapping_output, "\n")

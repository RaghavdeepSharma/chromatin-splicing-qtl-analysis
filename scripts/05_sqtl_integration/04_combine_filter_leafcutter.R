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

processed_dir <- file.path(base_dir, "05_processed")
results_dir   <- file.path(base_dir, "06_results")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.delim(
  manifest_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  fill = TRUE
)

required_manifest_columns <- c(
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
  "local_file"
)

missing_manifest_columns <- setdiff(
  required_manifest_columns,
  colnames(manifest)
)

if (length(missing_manifest_columns) > 0) {
  stop(
    "Manifest is missing columns: ",
    paste(missing_manifest_columns, collapse = ", ")
  )
}

missing_files <- manifest$local_file[
  !file.exists(manifest$local_file)
]

invalid_files <- manifest$local_file[
  file.exists(manifest$local_file) &
    file.info(manifest$local_file)$size == 0
]

if (length(missing_files) > 0) {
  writeLines(
    missing_files,
    file.path(results_dir, "missing_credible_set_files.txt")
  )

  stop(
    length(missing_files),
    " credible-set files are missing. See missing_credible_set_files.txt"
  )
}

if (length(invalid_files) > 0) {
  stop(
    length(invalid_files),
    " credible-set files are empty."
  )
}

required_qtl_columns <- c(
  "molecular_trait_id",
  "gene_id",
  "cs_id",
  "variant",
  "rsid",
  "cs_size",
  "pip",
  "pvalue",
  "beta",
  "se",
  "z",
  "cs_min_r2",
  "region"
)

read_credible_set <- function(i) {

  file_path <- manifest$local_file[i]

  message(
    "[", i, "/", nrow(manifest), "] Reading ",
    manifest$dataset_id[i],
    " — ",
    manifest$study_label[i],
    " / ",
    manifest$sample_group[i]
  )

  connection <- gzfile(file_path, open = "rt")

  qtl <- tryCatch(
    read.delim(
      connection,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    finally = close(connection)
  )

  missing_qtl_columns <- setdiff(
    required_qtl_columns,
    colnames(qtl)
  )

  if (length(missing_qtl_columns) > 0) {
    stop(
      "File ",
      file_path,
      " is missing columns: ",
      paste(missing_qtl_columns, collapse = ", ")
    )
  }

  qtl$pip <- suppressWarnings(as.numeric(qtl$pip))

  if (anyNA(qtl$pip)) {
    stop("Non-numeric PIP values found in: ", file_path)
  }

  if (any(qtl$pip < 0 | qtl$pip > 1)) {
    stop("PIP outside the range 0–1 in: ", file_path)
  }

  # Add dataset metadata.
  qtl$study_id              <- manifest$study_id[i]
  qtl$dataset_id            <- manifest$dataset_id[i]
  qtl$study_label           <- manifest$study_label[i]
  qtl$sample_group          <- manifest$sample_group[i]
  qtl$tissue_label          <- manifest$tissue_label[i]
  qtl$condition_label       <- manifest$condition_label[i]
  qtl$state_interpretation  <- manifest$state_interpretation[i]
  qtl$sample_size           <- manifest$sample_size[i]
  qtl$cell_class            <- manifest$preliminary_cell_class[i]
  qtl$calderon_priority     <- manifest$calderon_priority[i]
  qtl$download_tier         <- manifest$download_tier[i]
  qtl$quant_method          <- "leafcutter"
  qtl$source_file           <- file_path

  # Parse variant: chr1_111108395_A_G
  variant_parts <- strsplit(
    qtl$variant,
    split = "_",
    fixed = TRUE
  )

  qtl$chromosome <- vapply(
    variant_parts,
    function(x) if (length(x) >= 1) x[1] else NA_character_,
    character(1)
  )

  qtl$position <- suppressWarnings(
    as.integer(
      vapply(
        variant_parts,
        function(x) if (length(x) >= 2) x[2] else NA_character_,
        character(1)
      )
    )
  )

  qtl$ref <- vapply(
    variant_parts,
    function(x) if (length(x) >= 3) x[3] else NA_character_,
    character(1)
  )

  qtl$alt <- vapply(
    variant_parts,
    function(x) if (length(x) >= 4) x[4] else NA_character_,
    character(1)
  )

  if (anyNA(qtl$position)) {
    stop("Unable to parse one or more variant positions in: ", file_path)
  }

  # Parse LeafCutter trait:
  # 1:111139666:111140038:clu_35622_+
  trait_parts <- strsplit(
    qtl$molecular_trait_id,
    split = ":",
    fixed = TRUE
  )

  qtl$junction_chromosome <- vapply(
    trait_parts,
    function(x) if (length(x) >= 1) x[1] else NA_character_,
    character(1)
  )

  qtl$junction_start <- suppressWarnings(
    as.integer(
      vapply(
        trait_parts,
        function(x) if (length(x) >= 2) x[2] else NA_character_,
        character(1)
      )
    )
  )

  qtl$junction_end <- suppressWarnings(
    as.integer(
      vapply(
        trait_parts,
        function(x) if (length(x) >= 3) x[3] else NA_character_,
        character(1)
      )
    )
  )

  qtl$leafcutter_cluster <- vapply(
    trait_parts,
    function(x) if (length(x) >= 4) x[4] else NA_character_,
    character(1)
  )

  qtl$significant_qtl <- qtl$pip > 0.90

  qtl$pip_class <- ifelse(
    qtl$significant_qtl,
    "high_PIP",
    "low_PIP"
  )

  metadata_columns <- c(
    "study_id",
    "dataset_id",
    "study_label",
    "sample_group",
    "tissue_label",
    "condition_label",
    "state_interpretation",
    "sample_size",
    "cell_class",
    "calderon_priority",
    "download_tier",
    "quant_method"
  )

  parsed_columns <- c(
    "chromosome",
    "position",
    "ref",
    "alt",
    "junction_chromosome",
    "junction_start",
    "junction_end",
    "leafcutter_cluster",
    "significant_qtl",
    "pip_class"
  )

  original_columns <- required_qtl_columns

  qtl[
    ,
    c(
      metadata_columns,
      original_columns,
      parsed_columns,
      "source_file"
    )
  ]
}

association_list <- lapply(
  seq_len(nrow(manifest)),
  read_credible_set
)

association_data <- do.call(
  rbind,
  association_list
)

row.names(association_data) <- NULL

# Remove only exact duplicate association records.
association_key <- paste(
  association_data$dataset_id,
  association_data$gene_id,
  association_data$molecular_trait_id,
  association_data$cs_id,
  association_data$variant,
  sep = "\r"
)

duplicate_associations <- duplicated(association_key)
n_duplicate_associations <- sum(duplicate_associations)

if (n_duplicate_associations > 0) {
  association_data <- association_data[
    !duplicate_associations,
  ]
}

# Variant-level table:
# one record per dataset and variant, retaining the row with maximum PIP.
variant_order <- order(
  association_data$dataset_id,
  association_data$variant,
  -association_data$pip,
  association_data$molecular_trait_id
)

ordered_data <- association_data[
  variant_order,
]

dataset_variant_key <- paste(
  ordered_data$dataset_id,
  ordered_data$variant,
  sep = "\r"
)

variant_data <- ordered_data[
  !duplicated(dataset_variant_key),
]

row.names(variant_data) <- NULL

variant_data$max_pip <- variant_data$pip

variant_data$significant_qtl <- (
  variant_data$max_pip > 0.90
)

variant_data$pip_class <- ifelse(
  variant_data$significant_qtl,
  "high_PIP",
  "low_PIP"
)

high_pip_associations <- association_data[
  association_data$pip > 0.90,
]

high_pip_variants <- variant_data[
  variant_data$max_pip > 0.90,
]

low_pip_variants <- variant_data[
  variant_data$max_pip <= 0.90,
]

write_gzip_table <- function(data, output_file) {

  connection <- gzfile(
    output_file,
    open = "wt"
  )

  tryCatch(
    write.table(
      data,
      connection,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE,
      na = ""
    ),
    finally = close(connection)
  )
}

association_output <- file.path(
  processed_dir,
  "all_leafcutter_credible_set_associations.tsv.gz"
)

high_association_output <- file.path(
  processed_dir,
  "high_pip_leafcutter_associations.tsv.gz"
)

variant_output <- file.path(
  processed_dir,
  "leafcutter_dataset_variant_maxpip.tsv.gz"
)

high_variant_output <- file.path(
  processed_dir,
  "high_pip_leafcutter_dataset_variants.tsv.gz"
)

low_variant_output <- file.path(
  processed_dir,
  "low_pip_leafcutter_dataset_variants.tsv.gz"
)

write_gzip_table(
  association_data,
  association_output
)

write_gzip_table(
  high_pip_associations,
  high_association_output
)

write_gzip_table(
  variant_data,
  variant_output
)

write_gzip_table(
  high_pip_variants,
  high_variant_output
)

write_gzip_table(
  low_pip_variants,
  low_variant_output
)

# Dataset-level summary.
dataset_ids <- unique(association_data$dataset_id)

summary_list <- lapply(
  dataset_ids,
  function(current_dataset) {

    associations <- association_data[
      association_data$dataset_id == current_dataset,
    ]

    variants <- variant_data[
      variant_data$dataset_id == current_dataset,
    ]

    data.frame(
      study_id = associations$study_id[1],
      dataset_id = current_dataset,
      study_label = associations$study_label[1],
      sample_group = associations$sample_group[1],
      tissue_label = associations$tissue_label[1],
      condition_label = associations$condition_label[1],
      cell_class = associations$cell_class[1],
      sample_size = associations$sample_size[1],
      n_association_rows = nrow(associations),
      n_unique_splice_traits = length(
        unique(associations$molecular_trait_id)
      ),
      n_unique_genes = length(
        unique(associations$gene_id)
      ),
      n_unique_variants = nrow(variants),
      n_high_pip_associations = sum(
        associations$pip > 0.90
      ),
      n_high_pip_variants = sum(
        variants$max_pip > 0.90
      ),
      maximum_pip = max(
        associations$pip,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }
)

dataset_summary <- do.call(
  rbind,
  summary_list
)

dataset_summary <- dataset_summary[
  order(
    dataset_summary$study_label,
    dataset_summary$sample_group,
    dataset_summary$condition_label
  ),
]

summary_output <- file.path(
  results_dir,
  "leafcutter_sqtl_summary_by_dataset.tsv"
)

write.table(
  dataset_summary,
  summary_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

overall_summary <- data.frame(
  metric = c(
    "datasets_processed",
    "association_rows",
    "exact_duplicate_associations_removed",
    "unique_dataset_variants",
    "high_pip_association_rows",
    "high_pip_dataset_variants",
    "low_pip_dataset_variants",
    "unique_genes",
    "unique_splice_traits"
  ),
  value = c(
    length(unique(association_data$dataset_id)),
    nrow(association_data),
    n_duplicate_associations,
    nrow(variant_data),
    nrow(high_pip_associations),
    nrow(high_pip_variants),
    nrow(low_pip_variants),
    length(unique(association_data$gene_id)),
    length(unique(association_data$molecular_trait_id))
  )
)

write.table(
  overall_summary,
  file.path(
    results_dir,
    "leafcutter_sqtl_overall_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nLeafCutter processing complete.\n\n")

print(
  overall_summary,
  row.names = FALSE
)

cat("\nOutputs:\n")
cat(association_output, "\n")
cat(high_association_output, "\n")
cat(variant_output, "\n")
cat(high_variant_output, "\n")
cat(low_variant_output, "\n")
cat(summary_output, "\n")

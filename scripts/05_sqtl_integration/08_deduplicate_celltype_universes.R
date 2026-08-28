#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

base_dir <- Sys.getenv("KERIMOV_BASE")

if (base_dir == "") {
  stop("KERIMOV_BASE is not set. Source config.sh first.")
}

variant_file <- file.path(
  base_dir,
  "05_processed",
  "leafcutter_dataset_variants_protein_coding_lncRNA.tsv.gz"
)

mapping_file <- file.path(
  base_dir,
  "05_processed",
  "gene_filtered_sqtl_bed_manifest.tsv"
)

expected_counts_file <- file.path(
  base_dir,
  "06_results",
  "gene_biotype_filter_by_analysis_set.tsv"
)

output_bed_dir <- file.path(
  base_dir,
  "05_processed",
  "gene_filtered_sqtl_beds",
  "deduplicated_by_analysis_set"
)

manifest_output <- file.path(
  base_dir,
  "02_manifests",
  "motif_sqtl_deduplicated_celltype_manifest.tsv"
)

summary_output <- file.path(
  base_dir,
  "06_results",
  "motif_sqtl_deduplication_by_analysis_set.tsv"
)

dir.create(
  output_bed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

required_files <- c(
  variant_file,
  mapping_file,
  expected_counts_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing required files:\n",
    paste(missing_files, collapse = "\n")
  )
}

message("Reading gene-filtered sQTL variants...")

variant_connection <- gzfile(
  variant_file,
  open = "rt"
)

variants <- tryCatch(
  read.delim(
    variant_connection,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
  ),
  finally = close(variant_connection)
)

message("Reading dataset mapping manifest...")

mapping <- read.delim(
  mapping_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

message("Reading expected analysis-set counts...")

expected_counts <- read.delim(
  expected_counts_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

required_variant_columns <- c(
  "dataset_id",
  "chromosome",
  "bed_start",
  "bed_end",
  "variant",
  "gene_id_unversioned",
  "max_pip",
  "calderon_celltype"
)

required_mapping_columns <- c(
  "dataset_id",
  "calderon_celltype",
  "match_quality",
  "include_strict_primary",
  "include_expanded_primary",
  "include_sensitivity"
)

required_expected_columns <- c(
  "analysis_set",
  "retained_high_pip_rows",
  "retained_low_pip_rows"
)

missing_variant_columns <- setdiff(
  required_variant_columns,
  colnames(variants)
)

missing_mapping_columns <- setdiff(
  required_mapping_columns,
  colnames(mapping)
)

missing_expected_columns <- setdiff(
  required_expected_columns,
  colnames(expected_counts)
)

if (length(missing_variant_columns) > 0) {
  stop(
    "Variant table missing columns: ",
    paste(missing_variant_columns, collapse = ", ")
  )
}

if (length(missing_mapping_columns) > 0) {
  stop(
    "Mapping table missing columns: ",
    paste(missing_mapping_columns, collapse = ", ")
  )
}

if (length(missing_expected_columns) > 0) {
  stop(
    "Expected-count table missing columns: ",
    paste(missing_expected_columns, collapse = ", ")
  )
}

if (anyDuplicated(mapping$dataset_id)) {
  stop("Duplicate dataset IDs detected in the mapping manifest.")
}

variants$max_pip <- suppressWarnings(
  as.numeric(variants$max_pip)
)

variants$bed_start <- suppressWarnings(
  as.integer(variants$bed_start)
)

variants$bed_end <- suppressWarnings(
  as.integer(variants$bed_end)
)

if (anyNA(variants$max_pip)) {
  stop("Missing or non-numeric max_pip values detected.")
}

if (anyNA(variants$bed_start) || anyNA(variants$bed_end)) {
  stop("Missing or non-numeric BED coordinates detected.")
}

mapping_index <- match(
  variants$dataset_id,
  mapping$dataset_id
)

if (anyNA(mapping_index)) {
  stop(
    "Some filtered variant datasets are absent from the mapping manifest."
  )
}

variants$match_quality <- mapping$match_quality[
  mapping_index
]

analysis_sets <- list(
  strict_primary = "include_strict_primary",
  expanded_primary = "include_expanded_primary",
  sensitivity = "include_sensitivity"
)

chromosome_order <- c(
  paste0("chr", 1:22),
  "chrX",
  "chrY",
  "chrM"
)

manifest_rows <- list()
summary_rows <- list()

manifest_counter <- 0L
summary_counter <- 0L

for (analysis_name in names(analysis_sets)) {

  inclusion_column <- analysis_sets[[analysis_name]]

  selected_mapping <- mapping[
    mapping[[inclusion_column]] == "YES",
  ]

  selected_dataset_ids <- selected_mapping$dataset_id

  selected_variants <- variants[
    variants$dataset_id %in% selected_dataset_ids,
  ]

  expected_row <- expected_counts[
    expected_counts$analysis_set == analysis_name,
  ]

  if (nrow(expected_row) != 1) {
    stop(
      "Expected exactly one count-summary row for ",
      analysis_name
    )
  }

  expected_raw_rows <- (
    expected_row$retained_high_pip_rows +
      expected_row$retained_low_pip_rows
  )

  if (nrow(selected_variants) != expected_raw_rows) {
    stop(
      analysis_name,
      ": expected ",
      expected_raw_rows,
      " raw rows, found ",
      nrow(selected_variants)
    )
  }

  analysis_output_dir <- file.path(
    output_bed_dir,
    analysis_name
  )

  dir.create(
    analysis_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  celltypes <- sort(
    unique(selected_variants$calderon_celltype)
  )

  analysis_unique_total <- 0L
  analysis_high_total <- 0L
  analysis_low_total <- 0L
  analysis_duplicates_total <- 0L

  message("")
  message("Analysis set: ", analysis_name)

  for (celltype in celltypes) {

    cell_data <- selected_variants[
      selected_variants$calderon_celltype == celltype,
    ]

    raw_rows <- nrow(cell_data)

    cell_data$variant_gene_key <- paste(
      cell_data$variant,
      cell_data$gene_id_unversioned,
      sep = "|"
    )

    order_index <- order(
      cell_data$variant_gene_key,
      -cell_data$max_pip,
      cell_data$dataset_id
    )

    cell_data <- cell_data[
      order_index,
    ]

    # Ensure a variant–gene pair never has conflicting coordinates.
    first_index <- match(
      cell_data$variant_gene_key,
      cell_data$variant_gene_key
    )

    coordinate_string <- paste(
      cell_data$chromosome,
      cell_data$bed_start,
      cell_data$bed_end,
      sep = ":"
    )

    coordinate_conflict <- (
      coordinate_string !=
        coordinate_string[first_index]
    )

    if (any(coordinate_conflict)) {

      conflicting_keys <- unique(
        cell_data$variant_gene_key[
          coordinate_conflict
        ]
      )

      stop(
        analysis_name,
        " / ",
        celltype,
        ": conflicting coordinates for: ",
        paste(head(conflicting_keys, 20), collapse = ", ")
      )
    }

    # Rows are ordered by descending PIP, so the first row for each
    # variant–gene pair carries the maximum PIP.
    keep_row <- !duplicated(
      cell_data$variant_gene_key
    )

    deduplicated <- cell_data[
      keep_row,
    ]

    deduplicated$pip_class <- ifelse(
      deduplicated$max_pip > 0.90,
      "high_PIP",
      "low_PIP"
    )

    deduplicated$bed_score <- as.integer(
      round(deduplicated$max_pip * 1000)
    )

    deduplicated$bed_score <- pmax(
      0L,
      pmin(1000L, deduplicated$bed_score)
    )

    deduplicated$bed_name <- paste(
      analysis_name,
      celltype,
      deduplicated$variant,
      deduplicated$gene_id_unversioned,
      sep = "|"
    )

    chromosome_rank <- match(
      deduplicated$chromosome,
      chromosome_order
    )

    chromosome_rank[
      is.na(chromosome_rank)
    ] <- length(chromosome_order) + 1L

    deduplicated <- deduplicated[
      order(
        chromosome_rank,
        deduplicated$chromosome,
        deduplicated$bed_start,
        deduplicated$bed_end,
        deduplicated$variant,
        deduplicated$gene_id_unversioned
      ),
    ]

    output_bed <- file.path(
      analysis_output_dir,
      paste0(
        celltype,
        ".",
        analysis_name,
        ".unique_variant_gene.hg38.bed.gz"
      )
    )

    bed <- data.frame(
      chromosome = deduplicated$chromosome,
      start = deduplicated$bed_start,
      end = deduplicated$bed_end,
      name = deduplicated$bed_name,
      score = deduplicated$bed_score,
      strand = ".",
      analysis_set = analysis_name,
      calderon_celltype = celltype,
      max_pip = deduplicated$max_pip,
      pip_class = deduplicated$pip_class,
      gene_id = deduplicated$gene_id_unversioned,
      variant = deduplicated$variant,
      representative_dataset_id =
        deduplicated$dataset_id,
      representative_match_quality =
        deduplicated$match_quality,
      stringsAsFactors = FALSE
    )

    output_connection <- gzfile(
      output_bed,
      open = "wt"
    )

    tryCatch(
      write.table(
        bed,
        output_connection,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        na = "."
      ),
      finally = close(output_connection)
    )

    n_unique <- nrow(deduplicated)
    n_duplicates_collapsed <- raw_rows - n_unique

    n_high <- sum(
      deduplicated$pip_class == "high_PIP"
    )

    n_low <- sum(
      deduplicated$pip_class == "low_PIP"
    )

    if ((n_high + n_low) != n_unique) {
      stop(
        analysis_name,
        " / ",
        celltype,
        ": high and low counts do not sum."
      )
    }

    selected_celltype_mapping <- selected_mapping[
      selected_mapping$calderon_celltype == celltype,
    ]

    manifest_counter <- manifest_counter + 1L

    manifest_rows[[manifest_counter]] <- data.frame(
      analysis_set = analysis_name,
      calderon_celltype = celltype,
      accessibility_column =
        paste0("accessible_", celltype),
      n_datasets =
        nrow(selected_celltype_mapping),
      dataset_ids =
        paste(
          selected_celltype_mapping$dataset_id,
          collapse = ","
        ),
      raw_dataset_variant_gene_rows = raw_rows,
      unique_variant_gene_pairs = n_unique,
      duplicate_rows_collapsed =
        n_duplicates_collapsed,
      n_high_pip_pairs = n_high,
      n_low_pip_pairs = n_low,
      deduplicated_sqtl_bed = output_bed,
      stringsAsFactors = FALSE
    )

    analysis_unique_total <- (
      analysis_unique_total + n_unique
    )

    analysis_high_total <- (
      analysis_high_total + n_high
    )

    analysis_low_total <- (
      analysis_low_total + n_low
    )

    analysis_duplicates_total <- (
      analysis_duplicates_total +
        n_duplicates_collapsed
    )

    message(
      "  ",
      celltype,
      ": ",
      raw_rows,
      " raw → ",
      n_unique,
      " unique pairs; ",
      n_duplicates_collapsed,
      " collapsed; ",
      n_high,
      " high; ",
      n_low,
      " low"
    )
  }

  summary_counter <- summary_counter + 1L

  summary_rows[[summary_counter]] <- data.frame(
    analysis_set = analysis_name,
    n_datasets = length(selected_dataset_ids),
    n_celltypes = length(celltypes),
    raw_dataset_variant_gene_rows =
      nrow(selected_variants),
    unique_variant_gene_pairs =
      analysis_unique_total,
    duplicate_rows_collapsed =
      analysis_duplicates_total,
    unique_high_pip_pairs =
      analysis_high_total,
    unique_low_pip_pairs =
      analysis_low_total,
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(
  rbind,
  manifest_rows
)

summary <- do.call(
  rbind,
  summary_rows
)

write.table(
  manifest,
  manifest_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

write.table(
  summary,
  summary_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

cat("\nCell-type sQTL deduplication completed.\n\n")

print(
  summary,
  row.names = FALSE
)

cat("\nManifest rows: ", nrow(manifest), "\n", sep = "")
cat("Manifest: ", manifest_output, "\n", sep = "")
cat("Summary:  ", summary_output, "\n", sep = "")
cat("BED root: ", output_bed_dir, "\n", sep = "")

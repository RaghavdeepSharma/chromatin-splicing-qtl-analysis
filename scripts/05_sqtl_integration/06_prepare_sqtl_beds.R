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
  "leafcutter_dataset_variant_maxpip.tsv.gz"
)

mapping_file <- file.path(
  base_dir,
  "05_processed",
  "calderon_sqtl_dataset_mapping.tsv"
)

processed_dir <- file.path(base_dir, "05_processed")
results_dir   <- file.path(base_dir, "06_results")

bed_root <- file.path(
  processed_dir,
  "mapped_sqtl_beds"
)

high_bed_dir <- file.path(bed_root, "high_PIP")
low_bed_dir  <- file.path(bed_root, "low_PIP")

dir.create(high_bed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(low_bed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(variant_file)) {
  stop("Variant file not found: ", variant_file)
}

if (!file.exists(mapping_file)) {
  stop("Mapping file not found: ", mapping_file)
}

message("Reading dataset–variant table...")

variant_connection <- gzfile(
  variant_file,
  open = "rt"
)

variants <- tryCatch(
  read.delim(
    variant_connection,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  ),
  finally = close(variant_connection)
)

message("Reading Calderon mapping...")

mapping <- read.delim(
  mapping_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_variant_columns <- c(
  "dataset_id",
  "variant",
  "rsid",
  "gene_id",
  "molecular_trait_id",
  "cs_id",
  "chromosome",
  "position",
  "ref",
  "alt",
  "max_pip"
)

required_mapping_columns <- c(
  "dataset_id",
  "calderon_celltype",
  "match_quality",
  "mapping_reason",
  "calderon_condition_scope",
  "include_strict_primary",
  "include_expanded_primary",
  "include_sensitivity",
  "calderon_peak_file",
  "calderon_peak_count",
  "mapping_status"
)

missing_variant_columns <- setdiff(
  required_variant_columns,
  colnames(variants)
)

missing_mapping_columns <- setdiff(
  required_mapping_columns,
  colnames(mapping)
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

if (anyDuplicated(mapping$dataset_id)) {
  stop("Mapping contains duplicate dataset IDs.")
}

mapping_index <- match(
  variants$dataset_id,
  mapping$dataset_id
)

if (anyNA(mapping_index)) {
  missing_ids <- unique(
    variants$dataset_id[is.na(mapping_index)]
  )

  stop(
    "Variant datasets missing from mapping: ",
    paste(missing_ids, collapse = ", ")
  )
}

mapping_columns_to_attach <- c(
  "calderon_celltype",
  "match_quality",
  "mapping_reason",
  "calderon_condition_scope",
  "include_strict_primary",
  "include_expanded_primary",
  "include_sensitivity",
  "calderon_peak_file",
  "calderon_peak_count",
  "mapping_status"
)

for (column_name in mapping_columns_to_attach) {
  variants[[column_name]] <- mapping[[column_name]][mapping_index]
}

variants$max_pip <- suppressWarnings(
  as.numeric(variants$max_pip)
)

variants$position <- suppressWarnings(
  as.integer(variants$position)
)

if (anyNA(variants$max_pip)) {
  stop("Non-numeric max_pip values detected.")
}

if (any(variants$max_pip < 0 | variants$max_pip > 1)) {
  stop("max_pip values outside 0–1 detected.")
}

if (anyNA(variants$position)) {
  stop("Non-numeric variant positions detected.")
}

# Recalculate rather than trusting prior labels.
variants$pip_class <- ifelse(
  variants$max_pip > 0.90,
  "high_PIP",
  "low_PIP"
)

variants$significant_qtl <- (
  variants$max_pip > 0.90
)

# Convert GRCh38 1-based variant coordinates to BED coordinates.
#
# SNV:
# position 100 A>G becomes BED 99–100
#
# Deletion or multi-base reference allele:
# position 100 AT>A becomes BED 99–101
#
# Insertion:
# position 100 A>AT remains anchored at BED 99–100
variants$bed_start <- variants$position - 1L

valid_reference_allele <- grepl(
  "^[ACGTN]+$",
  variants$ref,
  ignore.case = TRUE
)

reference_width <- nchar(variants$ref)

reference_width[
  is.na(reference_width) |
    reference_width < 1 |
    !valid_reference_allele
] <- 1L

variants$bed_end <- (
  variants$bed_start +
    reference_width
)

variants$variant_span_rule <- ifelse(
  valid_reference_allele,
  "reference_allele_span",
  "fallback_one_base_anchor"
)

variants$coordinate_valid <- (
  grepl("^chr", variants$chromosome) &
    variants$bed_start >= 0 &
    variants$bed_end > variants$bed_start
)

if (any(!variants$coordinate_valid)) {
  invalid_output <- file.path(
    results_dir,
    "invalid_sqtl_variant_coordinates.tsv"
  )

  write.table(
    variants[!variants$coordinate_valid, ],
    invalid_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )

  stop(
    sum(!variants$coordinate_valid),
    " variants have invalid BED coordinates. See: ",
    invalid_output
  )
}

canonical_chromosomes <- paste0(
  "chr",
  c(1:22, "X", "Y", "M")
)

variants$canonical_chromosome <- (
  variants$chromosome %in%
    canonical_chromosomes
)

variants$bed_name <- paste(
  variants$dataset_id,
  variants$variant,
  variants$gene_id,
  variants$cs_id,
  sep = "|"
)

variants$bed_score <- as.integer(
  round(variants$max_pip * 1000)
)

variants$bed_score <- pmax(
  0L,
  pmin(1000L, variants$bed_score)
)

# Save complete mapped table, including the neutrophil no-match rows.
mapped_variant_output <- file.path(
  processed_dir,
  "leafcutter_dataset_variants_with_calderon_mapping.tsv.gz"
)

output_connection <- gzfile(
  mapped_variant_output,
  open = "wt"
)

tryCatch(
  write.table(
    variants,
    output_connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  ),
  finally = close(output_connection)
)

chromosome_order <- c(
  paste0("chr", 1:22),
  "chrX",
  "chrY",
  "chrM"
)

write_bed <- function(data, output_file) {

  if (nrow(data) > 0) {

    chromosome_rank <- match(
      data$chromosome,
      chromosome_order
    )

    chromosome_rank[
      is.na(chromosome_rank)
    ] <- length(chromosome_order) + 1L

    data <- data[
      order(
        chromosome_rank,
        data$chromosome,
        data$bed_start,
        data$bed_end,
        data$variant
      ),
    ]
  }

  bed <- data.frame(
    chromosome = data$chromosome,
    start = data$bed_start,
    end = data$bed_end,
    name = data$bed_name,
    score = data$bed_score,
    strand = ".",
    dataset_id = data$dataset_id,
    calderon_celltype = data$calderon_celltype,
    max_pip = data$max_pip,
    pip_class = data$pip_class,
    gene_id = data$gene_id,
    molecular_trait_id = data$molecular_trait_id,
    variant = data$variant,
    rsid = ifelse(is.na(data$rsid) | data$rsid == "", ".", data$rsid),
    stringsAsFactors = FALSE
  )

  connection <- gzfile(
    output_file,
    open = "wt"
  )

  tryCatch(
    write.table(
      bed,
      connection,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      na = ""
    ),
    finally = close(connection)
  )
}

mapped_dataset_ids <- mapping$dataset_id[
  mapping$mapping_status ==
    "MAPPED_TO_APPROVED_BAMPE_HG38"
]

manifest_rows <- vector(
  mode = "list",
  length = length(mapped_dataset_ids)
)

for (i in seq_along(mapped_dataset_ids)) {

  current_dataset <- mapped_dataset_ids[i]

  dataset_mapping <- mapping[
    mapping$dataset_id == current_dataset,
  ]

  dataset_variants <- variants[
    variants$dataset_id == current_dataset,
  ]

  high_variants <- dataset_variants[
    dataset_variants$max_pip > 0.90,
  ]

  low_variants <- dataset_variants[
    dataset_variants$max_pip <= 0.90,
  ]

  high_bed <- file.path(
    high_bed_dir,
    paste0(current_dataset, ".high_PIP.hg38.bed.gz")
  )

  low_bed <- file.path(
    low_bed_dir,
    paste0(current_dataset, ".low_PIP.hg38.bed.gz")
  )

  message(
    "[", i, "/", length(mapped_dataset_ids), "] ",
    current_dataset,
    " — ",
    dataset_mapping$study_label,
    " / ",
    dataset_mapping$sample_group,
    " → ",
    dataset_mapping$calderon_celltype
  )

  write_bed(
    high_variants,
    high_bed
  )

  write_bed(
    low_variants,
    low_bed
  )

  manifest_rows[[i]] <- data.frame(
    study_id = dataset_mapping$study_id,
    dataset_id = current_dataset,
    study_label = dataset_mapping$study_label,
    sample_group = dataset_mapping$sample_group,
    eqtl_condition = dataset_mapping$condition_label,
    interpreted_state = dataset_mapping$state_interpretation,
    sample_size = dataset_mapping$sample_size,
    calderon_celltype = dataset_mapping$calderon_celltype,
    match_quality = dataset_mapping$match_quality,
    calderon_condition_scope =
      dataset_mapping$calderon_condition_scope,
    include_strict_primary =
      dataset_mapping$include_strict_primary,
    include_expanded_primary =
      dataset_mapping$include_expanded_primary,
    include_sensitivity =
      dataset_mapping$include_sensitivity,
    n_high_pip_variants = nrow(high_variants),
    n_low_pip_variants = nrow(low_variants),
    high_pip_bed = high_bed,
    low_pip_bed = low_bed,
    calderon_peak_file =
      dataset_mapping$calderon_peak_file,
    calderon_peak_count =
      dataset_mapping$calderon_peak_count,
    stringsAsFactors = FALSE
  )
}

bed_manifest <- do.call(
  rbind,
  manifest_rows
)

bed_manifest_output <- file.path(
  processed_dir,
  "mapped_sqtl_bed_manifest.tsv"
)

write.table(
  bed_manifest,
  bed_manifest_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

# Excluded/no-match datasets retained separately.
excluded_mapping <- mapping[
  mapping$mapping_status !=
    "MAPPED_TO_APPROVED_BAMPE_HG38",
]

write.table(
  excluded_mapping,
  file.path(
    results_dir,
    "sqtl_datasets_excluded_from_atac_mapping.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

# Summary by mapping quality.
quality_summary <- aggregate(
  cbind(
    n_high_pip_variants,
    n_low_pip_variants
  ) ~ match_quality,
  data = bed_manifest,
  FUN = sum
)

dataset_counts <- as.data.frame(
  table(bed_manifest$match_quality),
  stringsAsFactors = FALSE
)

colnames(dataset_counts) <- c(
  "match_quality",
  "n_datasets"
)

quality_summary <- merge(
  dataset_counts,
  quality_summary,
  by = "match_quality",
  all = TRUE
)

quality_summary <- quality_summary[
  order(-quality_summary$n_datasets),
]

write.table(
  quality_summary,
  file.path(
    results_dir,
    "mapped_sqtl_counts_by_match_quality.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

analysis_summary <- data.frame(
  analysis_set = c(
    "strict_primary",
    "expanded_primary",
    "sensitivity"
  ),

  n_datasets = c(
    sum(
      bed_manifest$include_strict_primary == "YES"
    ),
    sum(
      bed_manifest$include_expanded_primary == "YES"
    ),
    sum(
      bed_manifest$include_sensitivity == "YES"
    )
  ),

  n_high_pip_variants = c(
    sum(
      bed_manifest$n_high_pip_variants[
        bed_manifest$include_strict_primary == "YES"
      ]
    ),
    sum(
      bed_manifest$n_high_pip_variants[
        bed_manifest$include_expanded_primary == "YES"
      ]
    ),
    sum(
      bed_manifest$n_high_pip_variants[
        bed_manifest$include_sensitivity == "YES"
      ]
    )
  ),

  n_low_pip_variants = c(
    sum(
      bed_manifest$n_low_pip_variants[
        bed_manifest$include_strict_primary == "YES"
      ]
    ),
    sum(
      bed_manifest$n_low_pip_variants[
        bed_manifest$include_expanded_primary == "YES"
      ]
    ),
    sum(
      bed_manifest$n_low_pip_variants[
        bed_manifest$include_sensitivity == "YES"
      ]
    )
  )
)

write.table(
  analysis_summary,
  file.path(
    results_dir,
    "mapped_sqtl_counts_by_analysis_set.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

overall_summary <- data.frame(
  metric = c(
    "input_dataset_variant_rows",
    "mapped_dataset_variant_rows",
    "excluded_no_match_variant_rows",
    "mapped_datasets",
    "excluded_datasets",
    "high_pip_bed_files",
    "low_pip_bed_files",
    "mapped_high_pip_records",
    "mapped_low_pip_records",
    "noncanonical_chromosome_records"
  ),

  value = c(
    nrow(variants),
    sum(
      variants$mapping_status ==
        "MAPPED_TO_APPROVED_BAMPE_HG38"
    ),
    sum(
      variants$mapping_status !=
        "MAPPED_TO_APPROVED_BAMPE_HG38"
    ),
    nrow(bed_manifest),
    nrow(excluded_mapping),
    length(
      list.files(
        high_bed_dir,
        pattern = "\\.bed\\.gz$"
      )
    ),
    length(
      list.files(
        low_bed_dir,
        pattern = "\\.bed\\.gz$"
      )
    ),
    sum(bed_manifest$n_high_pip_variants),
    sum(bed_manifest$n_low_pip_variants),
    sum(!variants$canonical_chromosome)
  )
)

write.table(
  overall_summary,
  file.path(
    results_dir,
    "mapped_sqtl_bed_preparation_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nMapped sQTL BED preparation complete.\n\n")

print(
  overall_summary,
  row.names = FALSE
)

cat("\nCounts by analysis set:\n")
print(
  analysis_summary,
  row.names = FALSE
)

cat("\nCounts by match quality:\n")
print(
  quality_summary,
  row.names = FALSE
)

cat("\nOutputs:\n")
cat(mapped_variant_output, "\n")
cat(bed_manifest_output, "\n")
cat(high_bed_dir, "\n")
cat(low_bed_dir, "\n")

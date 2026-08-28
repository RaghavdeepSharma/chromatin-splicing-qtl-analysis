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
  "leafcutter_dataset_variants_with_calderon_mapping.tsv.gz"
)

gene_lookup_file <- file.path(
  base_dir,
  "01_metadata",
  "gencode_v47_gene_lookup.tsv"
)

old_manifest_file <- file.path(
  base_dir,
  "05_processed",
  "mapped_sqtl_bed_manifest.tsv"
)

processed_dir <- file.path(
  base_dir,
  "05_processed"
)

results_dir <- file.path(
  base_dir,
  "06_results"
)

bed_root <- file.path(
  processed_dir,
  "gene_filtered_sqtl_beds"
)

high_bed_dir <- file.path(
  bed_root,
  "high_PIP"
)

low_bed_dir <- file.path(
  bed_root,
  "low_PIP"
)

dir.create(
  high_bed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  low_bed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

required_files <- c(
  variant_file,
  gene_lookup_file,
  old_manifest_file
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

message("Reading mapped sQTL variant table...")

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

message("Reading GENCODE v47 lookup...")

gene_lookup <- read.delim(
  gene_lookup_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

message("Reading mapped sQTL BED manifest...")

old_manifest <- read.delim(
  old_manifest_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

required_variant_columns <- c(
  "dataset_id",
  "chromosome",
  "position",
  "ref",
  "alt",
  "variant",
  "rsid",
  "gene_id",
  "molecular_trait_id",
  "cs_id",
  "max_pip",
  "calderon_celltype",
  "mapping_status"
)

required_lookup_columns <- c(
  "gene_id",
  "gene_name",
  "gene_type"
)

required_manifest_columns <- c(
  "dataset_id",
  "study_id",
  "study_label",
  "sample_group",
  "eqtl_condition",
  "calderon_celltype",
  "match_quality",
  "include_strict_primary",
  "include_expanded_primary",
  "include_sensitivity",
  "calderon_peak_file"
)

missing_variant_columns <- setdiff(
  required_variant_columns,
  colnames(variants)
)

missing_lookup_columns <- setdiff(
  required_lookup_columns,
  colnames(gene_lookup)
)

missing_manifest_columns <- setdiff(
  required_manifest_columns,
  colnames(old_manifest)
)

if (length(missing_variant_columns) > 0) {
  stop(
    "Variant table missing columns: ",
    paste(missing_variant_columns, collapse = ", ")
  )
}

if (length(missing_lookup_columns) > 0) {
  stop(
    "Gene lookup missing columns: ",
    paste(missing_lookup_columns, collapse = ", ")
  )
}

if (length(missing_manifest_columns) > 0) {
  stop(
    "BED manifest missing columns: ",
    paste(missing_manifest_columns, collapse = ", ")
  )
}

if (anyDuplicated(gene_lookup$gene_id)) {
  stop("GENCODE lookup contains duplicate unversioned gene IDs.")
}

if (anyDuplicated(old_manifest$dataset_id)) {
  stop("Mapped BED manifest contains duplicate dataset IDs.")
}

variants$max_pip <- suppressWarnings(
  as.numeric(variants$max_pip)
)

variants$position <- suppressWarnings(
  as.integer(variants$position)
)

if (anyNA(variants$max_pip)) {
  stop("Non-numeric or missing max_pip values detected.")
}

if (anyNA(variants$position)) {
  stop("Non-numeric or missing variant positions detected.")
}

variants$gene_id_unversioned <- sub(
  "\\..*$",
  "",
  variants$gene_id
)

lookup_index <- match(
  variants$gene_id_unversioned,
  gene_lookup$gene_id
)

variants$gencode_gene_name <- gene_lookup$gene_name[
  lookup_index
]

variants$gencode_gene_type <- gene_lookup$gene_type[
  lookup_index
]

variants$gencode_match <- ifelse(
  is.na(lookup_index),
  "NO",
  "YES"
)

variants$gene_biotype_filter <- ifelse(
  is.na(lookup_index),
  "unmatched_gencode_v47",
  ifelse(
    variants$gencode_gene_type %in%
      c("protein_coding", "lncRNA"),
    "retain",
    "exclude_other_biotype"
  )
)

mapped_variants <- variants[
  variants$mapping_status ==
    "MAPPED_TO_APPROVED_BAMPE_HG38",
]

retained <- mapped_variants[
  mapped_variants$gene_biotype_filter == "retain",
]

retained$pip_class <- ifelse(
  retained$max_pip > 0.90,
  "high_PIP",
  "low_PIP"
)

# Recreate BED coordinates if they are not already present.
if (!all(c("bed_start", "bed_end") %in% colnames(retained))) {

  retained$bed_start <- retained$position - 1L

  valid_ref <- grepl(
    "^[ACGTN]+$",
    retained$ref,
    ignore.case = TRUE
  )

  reference_width <- nchar(retained$ref)

  reference_width[
    is.na(reference_width) |
      reference_width < 1 |
      !valid_ref
  ] <- 1L

  retained$bed_end <- (
    retained$bed_start +
      reference_width
  )
}

retained$bed_start <- as.integer(
  retained$bed_start
)

retained$bed_end <- as.integer(
  retained$bed_end
)

invalid_coordinates <- (
  is.na(retained$bed_start) |
    is.na(retained$bed_end) |
    retained$bed_start < 0 |
    retained$bed_end <= retained$bed_start |
    !grepl("^chr", retained$chromosome)
)

if (any(invalid_coordinates)) {

  invalid_output <- file.path(
    results_dir,
    "gene_filtered_invalid_coordinates.tsv"
  )

  write.table(
    retained[invalid_coordinates, ],
    invalid_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )

  stop(
    sum(invalid_coordinates),
    " invalid BED coordinates detected. See: ",
    invalid_output
  )
}

retained$rsid_clean <- ifelse(
  is.na(retained$rsid) |
    retained$rsid == "",
  ".",
  retained$rsid
)

retained$bed_score <- as.integer(
  round(retained$max_pip * 1000)
)

retained$bed_score <- pmax(
  0L,
  pmin(1000L, retained$bed_score)
)

retained$bed_name <- paste(
  retained$dataset_id,
  retained$variant,
  retained$gene_id_unversioned,
  retained$cs_id,
  sep = "|"
)

# Save all variants with gene-biotype annotation.
annotated_output <- file.path(
  processed_dir,
  "leafcutter_dataset_variants_gene_biotype_annotated.tsv.gz"
)

annotated_connection <- gzfile(
  annotated_output,
  open = "wt"
)

tryCatch(
  write.table(
    variants,
    annotated_connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  ),
  finally = close(annotated_connection)
)

# Save only mapped protein-coding/lncRNA variants.
filtered_output <- file.path(
  processed_dir,
  "leafcutter_dataset_variants_protein_coding_lncRNA.tsv.gz"
)

filtered_connection <- gzfile(
  filtered_output,
  open = "wt"
)

tryCatch(
  write.table(
    retained,
    filtered_connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  ),
  finally = close(filtered_connection)
)

chromosome_order <- c(
  paste0("chr", 1:22),
  "chrX",
  "chrY",
  "chrM"
)

write_sqtl_bed <- function(data, output_file) {

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
    gene_id = data$gene_id_unversioned,
    molecular_trait_id = data$molecular_trait_id,
    variant = data$variant,
    rsid = data$rsid_clean,
    stringsAsFactors = FALSE
  )

  output_connection <- gzfile(
    output_file,
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
}

new_manifest <- old_manifest

new_manifest$n_high_pip_before_gene_filter <-
  old_manifest$n_high_pip_variants

new_manifest$n_low_pip_before_gene_filter <-
  old_manifest$n_low_pip_variants

new_manifest$n_high_pip_variants <- 0L
new_manifest$n_low_pip_variants <- 0L

new_manifest$high_pip_bed <- ""
new_manifest$low_pip_bed <- ""

dataset_summaries <- vector(
  mode = "list",
  length = nrow(new_manifest)
)

for (i in seq_len(nrow(new_manifest))) {

  dataset_id <- new_manifest$dataset_id[i]

  dataset_mapped <- mapped_variants[
    mapped_variants$dataset_id == dataset_id,
  ]

  dataset_retained <- retained[
    retained$dataset_id == dataset_id,
  ]

  high_data <- dataset_retained[
    dataset_retained$max_pip > 0.90,
  ]

  low_data <- dataset_retained[
    dataset_retained$max_pip <= 0.90,
  ]

  high_bed <- file.path(
    high_bed_dir,
    paste0(
      dataset_id,
      ".high_PIP.protein_coding_lncRNA.hg38.bed.gz"
    )
  )

  low_bed <- file.path(
    low_bed_dir,
    paste0(
      dataset_id,
      ".low_PIP.protein_coding_lncRNA.hg38.bed.gz"
    )
  )

  message(
    "[", i, "/", nrow(new_manifest), "] ",
    dataset_id,
    " — high PIP: ",
    nrow(high_data),
    "; low PIP: ",
    nrow(low_data)
  )

  write_sqtl_bed(
    high_data,
    high_bed
  )

  write_sqtl_bed(
    low_data,
    low_bed
  )

  new_manifest$n_high_pip_variants[i] <-
    nrow(high_data)

  new_manifest$n_low_pip_variants[i] <-
    nrow(low_data)

  new_manifest$high_pip_bed[i] <-
    high_bed

  new_manifest$low_pip_bed[i] <-
    low_bed

  dataset_summaries[[i]] <- data.frame(
    dataset_id = dataset_id,
    study_label = new_manifest$study_label[i],
    sample_group = new_manifest$sample_group[i],
    calderon_celltype =
      new_manifest$calderon_celltype[i],
    match_quality =
      new_manifest$match_quality[i],
    mapped_rows_before_filter =
      nrow(dataset_mapped),
    matched_gencode_rows =
      sum(dataset_mapped$gencode_match == "YES"),
    unmatched_gencode_rows =
      sum(dataset_mapped$gencode_match == "NO"),
    excluded_other_biotype_rows =
      sum(
        dataset_mapped$gene_biotype_filter ==
          "exclude_other_biotype"
      ),
    retained_protein_coding_rows =
      sum(
        dataset_retained$gencode_gene_type ==
          "protein_coding"
      ),
    retained_lncRNA_rows =
      sum(
        dataset_retained$gencode_gene_type ==
          "lncRNA"
      ),
    retained_high_pip_rows =
      nrow(high_data),
    retained_low_pip_rows =
      nrow(low_data),
    stringsAsFactors = FALSE
  )
}

new_manifest$gene_filter <-
  "GENCODE_v47_protein_coding_or_lncRNA"

new_manifest_output <- file.path(
  processed_dir,
  "gene_filtered_sqtl_bed_manifest.tsv"
)

write.table(
  new_manifest,
  new_manifest_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

dataset_summary <- do.call(
  rbind,
  dataset_summaries
)

write.table(
  dataset_summary,
  file.path(
    results_dir,
    "gene_biotype_filter_by_dataset.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

unmatched_rows <- mapped_variants[
  mapped_variants$gencode_match == "NO",
]

unmatched_gene_summary <- data.frame(
  gene_id = sort(
    unique(unmatched_rows$gene_id_unversioned)
  ),
  stringsAsFactors = FALSE
)

if (nrow(unmatched_gene_summary) > 0) {

  unmatched_gene_summary$n_dataset_variant_rows <- vapply(
    unmatched_gene_summary$gene_id,
    function(current_gene) {
      sum(
        unmatched_rows$gene_id_unversioned ==
          current_gene
      )
    },
    integer(1)
  )

  unmatched_gene_summary$n_high_pip_rows <- vapply(
    unmatched_gene_summary$gene_id,
    function(current_gene) {
      sum(
        unmatched_rows$gene_id_unversioned ==
          current_gene &
          unmatched_rows$max_pip > 0.90
      )
    },
    integer(1)
  )

  unmatched_gene_summary <- unmatched_gene_summary[
    order(
      -unmatched_gene_summary$n_dataset_variant_rows,
      unmatched_gene_summary$gene_id
    ),
  ]
}

write.table(
  unmatched_gene_summary,
  file.path(
    results_dir,
    "unmatched_sqtl_gene_ids_gencode_v47.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

excluded_biotypes <- mapped_variants$gencode_gene_type[
  mapped_variants$gene_biotype_filter ==
    "exclude_other_biotype"
]

excluded_biotype_summary <- as.data.frame(
  table(excluded_biotypes),
  stringsAsFactors = FALSE
)

colnames(excluded_biotype_summary) <- c(
  "gene_type",
  "n_dataset_variant_rows"
)

excluded_biotype_summary <- excluded_biotype_summary[
  order(
    -excluded_biotype_summary$n_dataset_variant_rows
  ),
]

write.table(
  excluded_biotype_summary,
  file.path(
    results_dir,
    "excluded_sqtl_gene_biotypes.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

overall_summary <- data.frame(
  metric = c(
    "all_variant_rows",
    "mapped_variant_rows_before_gene_filter",
    "matched_gencode_v47_rows",
    "unmatched_gencode_v47_rows",
    "excluded_other_gene_biotype_rows",
    "retained_protein_coding_rows",
    "retained_lncRNA_rows",
    "retained_total_rows",
    "retained_high_pip_rows",
    "retained_low_pip_rows",
    "retained_unique_genes",
    "unmatched_unique_gene_ids",
    "high_pip_bed_files",
    "low_pip_bed_files"
  ),
  value = c(
    nrow(variants),
    nrow(mapped_variants),
    sum(mapped_variants$gencode_match == "YES"),
    sum(mapped_variants$gencode_match == "NO"),
    sum(
      mapped_variants$gene_biotype_filter ==
        "exclude_other_biotype"
    ),
    sum(
      retained$gencode_gene_type ==
        "protein_coding"
    ),
    sum(
      retained$gencode_gene_type ==
        "lncRNA"
    ),
    nrow(retained),
    sum(retained$max_pip > 0.90),
    sum(retained$max_pip <= 0.90),
    length(unique(retained$gene_id_unversioned)),
    length(
      unique(unmatched_rows$gene_id_unversioned)
    ),
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
    )
  )
)

write.table(
  overall_summary,
  file.path(
    results_dir,
    "gene_biotype_filter_overall_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summarize_analysis_set <- function(
  inclusion_column,
  analysis_name
) {

  selected_datasets <- new_manifest$dataset_id[
    new_manifest[[inclusion_column]] == "YES"
  ]

  selected_rows <- retained[
    retained$dataset_id %in% selected_datasets,
  ]

  data.frame(
    analysis_set = analysis_name,
    n_datasets = length(selected_datasets),
    retained_high_pip_rows =
      sum(selected_rows$max_pip > 0.90),
    retained_low_pip_rows =
      sum(selected_rows$max_pip <= 0.90),
    retained_unique_genes =
      length(
        unique(selected_rows$gene_id_unversioned)
      ),
    stringsAsFactors = FALSE
  )
}

analysis_set_summary <- rbind(
  summarize_analysis_set(
    "include_strict_primary",
    "strict_primary"
  ),
  summarize_analysis_set(
    "include_expanded_primary",
    "expanded_primary"
  ),
  summarize_analysis_set(
    "include_sensitivity",
    "sensitivity"
  )
)

write.table(
  analysis_set_summary,
  file.path(
    results_dir,
    "gene_biotype_filter_by_analysis_set.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nGene-biotype filtering completed.\n\n")

print(
  overall_summary,
  row.names = FALSE
)

cat("\nCounts by analysis set:\n")

print(
  analysis_set_summary,
  row.names = FALSE
)

cat("\nOutputs:\n")
cat(annotated_output, "\n")
cat(filtered_output, "\n")
cat(new_manifest_output, "\n")
cat(high_bed_dir, "\n")
cat(low_bed_dir, "\n")

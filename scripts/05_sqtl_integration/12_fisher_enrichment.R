#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

base_dir <- Sys.getenv("KERIMOV_BASE")

if (base_dir == "") {
  stop("KERIMOV_BASE is not set. Source config.sh first.")
}

input_dir <- file.path(
  base_dir,
  "06_results",
  "motif_component_overlap_counts"
)

results_dir <- file.path(
  base_dir,
  "06_results"
)

component_output <- file.path(
  results_dir,
  "all_motif_component_overlap_counts.tsv.gz"
)

contingency_output <- file.path(
  results_dir,
  "motif_archetype_contingency_tables.tsv"
)

fisher_output <- file.path(
  results_dir,
  "motif_archetype_fisher_results.tsv"
)

summary_output <- file.path(
  results_dir,
  "motif_archetype_analysis_summary.tsv"
)

analysis_set_summary_output <- file.path(
  results_dir,
  "motif_archetype_fisher_summary_by_analysis_set.tsv"
)

if (!dir.exists(input_dir)) {
  stop("Input directory not found: ", input_dir)
}

result_files <- list.files(
  input_dir,
  pattern = "^task_[0-9]+\\..*\\.motif_overlap_counts\\.tsv$",
  full.names = TRUE
)

result_files <- sort(result_files)

if (length(result_files) != 508) {
  stop(
    "Expected 508 motif-component result files; found ",
    length(result_files)
  )
}

message("Reading 508 motif-component result files...")

component_list <- vector(
  mode = "list",
  length = length(result_files)
)

for (i in seq_along(result_files)) {

  current_file <- result_files[i]

  current_data <- read.delim(
    current_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
  )

  if (nrow(current_data) != 29) {
    stop(
      "Expected 29 data rows in ",
      current_file,
      "; found ",
      nrow(current_data)
    )
  }

  current_data$source_result_file <- current_file

  component_list[[i]] <- current_data
}

components <- do.call(
  rbind,
  component_list
)

required_columns <- c(
  "task_id",
  "component_name",
  "candidate_archetype",
  "is_split_component",
  "chromosome_component",
  "analysis_set",
  "calderon_celltype",
  "n_accessible_intragenic_motif_gene_rows",
  "n_high_total",
  "n_high_overlap",
  "n_high_not_overlap",
  "n_low_total",
  "n_low_overlap",
  "n_low_not_overlap",
  "motif_file",
  "sqtl_bed"
)

missing_columns <- setdiff(
  required_columns,
  colnames(components)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

numeric_columns <- c(
  "task_id",
  "n_accessible_intragenic_motif_gene_rows",
  "n_high_total",
  "n_high_overlap",
  "n_high_not_overlap",
  "n_low_total",
  "n_low_overlap",
  "n_low_not_overlap"
)

for (column_name in numeric_columns) {

  components[[column_name]] <- suppressWarnings(
    as.numeric(components[[column_name]])
  )

  if (anyNA(components[[column_name]])) {
    stop(
      "Missing or non-numeric values detected in ",
      column_name
    )
  }
}

expected_component_rows <- 508 * 29

if (nrow(components) != expected_component_rows) {
  stop(
    "Expected ",
    expected_component_rows,
    " component rows; found ",
    nrow(components)
  )
}

component_key <- paste(
  components$task_id,
  components$analysis_set,
  components$calderon_celltype,
  sep = "|"
)

if (anyDuplicated(component_key)) {
  stop(
    "Duplicate task × analysis-set × cell-type rows detected."
  )
}

task_row_counts <- table(
  components$task_id
)

if (any(task_row_counts != 29)) {
  stop(
    "At least one motif component does not have 29 rows."
  )
}

if (length(unique(components$task_id)) != 508) {
  stop("Expected 508 unique task IDs.")
}

if (length(unique(components$candidate_archetype)) != 286) {
  stop(
    "Expected 286 candidate archetypes; found ",
    length(unique(components$candidate_archetype))
  )
}

if (any(
  components$n_high_overlap +
    components$n_high_not_overlap !=
    components$n_high_total
)) {
  stop(
    "High-PIP component contingency counts do not sum."
  )
}

if (any(
  components$n_low_overlap +
    components$n_low_not_overlap !=
    components$n_low_total
)) {
  stop(
    "Low-PIP component contingency counts do not sum."
  )
}

message("Saving combined component-level table...")

component_connection <- gzfile(
  component_output,
  open = "wt"
)

tryCatch(
  write.table(
    components,
    component_connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  ),
  finally = close(component_connection)
)

message(
  "Merging chromosome components into biological archetypes..."
)

group_key <- paste(
  components$candidate_archetype,
  components$analysis_set,
  components$calderon_celltype,
  sep = "\034"
)

component_groups <- split(
  components,
  group_key
)

merged_rows <- vector(
  mode = "list",
  length = length(component_groups)
)

group_counter <- 0L

for (current_group in component_groups) {

  group_counter <- group_counter + 1L

  high_totals <- unique(
    current_group$n_high_total
  )

  low_totals <- unique(
    current_group$n_low_total
  )

  sqtl_beds <- unique(
    current_group$sqtl_bed
  )

  if (length(high_totals) != 1) {
    stop(
      "Conflicting high-PIP totals for ",
      current_group$candidate_archetype[1],
      " / ",
      current_group$analysis_set[1],
      " / ",
      current_group$calderon_celltype[1]
    )
  }

  if (length(low_totals) != 1) {
    stop(
      "Conflicting low-PIP totals for ",
      current_group$candidate_archetype[1],
      " / ",
      current_group$analysis_set[1],
      " / ",
      current_group$calderon_celltype[1]
    )
  }

  if (length(sqtl_beds) != 1) {
    stop(
      "Conflicting sQTL BED files for an archetype group."
    )
  }

  high_total <- high_totals[1]
  low_total <- low_totals[1]

  high_overlap <- sum(
    current_group$n_high_overlap
  )

  low_overlap <- sum(
    current_group$n_low_overlap
  )

  if (high_overlap > high_total) {
    stop(
      "Merged high-PIP overlap exceeds the total for ",
      current_group$candidate_archetype[1],
      " / ",
      current_group$analysis_set[1],
      " / ",
      current_group$calderon_celltype[1]
    )
  }

  if (low_overlap > low_total) {
    stop(
      "Merged low-PIP overlap exceeds the total for ",
      current_group$candidate_archetype[1],
      " / ",
      current_group$analysis_set[1],
      " / ",
      current_group$calderon_celltype[1]
    )
  }

  merged_rows[[group_counter]] <- data.frame(
    candidate_archetype =
      current_group$candidate_archetype[1],

    analysis_set =
      current_group$analysis_set[1],

    calderon_celltype =
      current_group$calderon_celltype[1],

    n_component_files =
      length(unique(current_group$task_id)),

    n_chromosome_split_components =
      sum(
        unique(
          current_group[
            c(
              "task_id",
              "is_split_component"
            )
          ]
        )$is_split_component == "YES"
      ),

    n_accessible_intragenic_motif_gene_rows =
      sum(
        current_group$
          n_accessible_intragenic_motif_gene_rows
      ),

    significant_qtl_overlap_true =
      high_overlap,

    significant_qtl_overlap_false =
      high_total - high_overlap,

    nonsignificant_qtl_overlap_true =
      low_overlap,

    nonsignificant_qtl_overlap_false =
      low_total - low_overlap,

    n_high_total =
      high_total,

    n_low_total =
      low_total,

    component_task_ids =
      paste(
        sort(unique(current_group$task_id)),
        collapse = ","
      ),

    component_names =
      paste(
        sort(unique(current_group$component_name)),
        collapse = ","
      ),

    sqtl_bed =
      sqtl_beds[1],

    stringsAsFactors = FALSE
  )
}

archetypes <- do.call(
  rbind,
  merged_rows
)

expected_archetype_rows <- 286 * 29

if (nrow(archetypes) != expected_archetype_rows) {
  stop(
    "Expected ",
    expected_archetype_rows,
    " archetype rows; found ",
    nrow(archetypes)
  )
}

archetype_key <- paste(
  archetypes$candidate_archetype,
  archetypes$analysis_set,
  archetypes$calderon_celltype,
  sep = "|"
)

if (anyDuplicated(archetype_key)) {
  stop(
    "Duplicate archetype × analysis-set × cell-type rows detected."
  )
}

if (any(
  archetypes$significant_qtl_overlap_true +
    archetypes$significant_qtl_overlap_false !=
    archetypes$n_high_total
)) {
  stop(
    "Merged significant-QTL contingency counts do not sum."
  )
}

if (any(
  archetypes$nonsignificant_qtl_overlap_true +
    archetypes$nonsignificant_qtl_overlap_false !=
    archetypes$n_low_total
)) {
  stop(
    "Merged nonsignificant-QTL contingency counts do not sum."
  )
}

analysis_order <- c(
  "strict_primary",
  "expanded_primary",
  "sensitivity"
)

archetypes$analysis_order <- match(
  archetypes$analysis_set,
  analysis_order
)

archetypes <- archetypes[
  order(
    archetypes$analysis_order,
    archetypes$calderon_celltype,
    archetypes$candidate_archetype
  ),
]

archetypes$analysis_order <- NULL

message("Writing contingency tables...")

write.table(
  archetypes,
  contingency_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

message("Running Fisher exact tests...")

n_tests <- nrow(archetypes)

fisher_odds_ratio <- numeric(n_tests)
fisher_ci_lower <- numeric(n_tests)
fisher_ci_upper <- numeric(n_tests)
fisher_p_value <- numeric(n_tests)

for (i in seq_len(n_tests)) {

  a <- archetypes$
    significant_qtl_overlap_true[i]

  b <- archetypes$
    significant_qtl_overlap_false[i]

  c <- archetypes$
    nonsignificant_qtl_overlap_true[i]

  d <- archetypes$
    nonsignificant_qtl_overlap_false[i]

  contingency_table <- matrix(
    c(a, b, c, d),
    nrow = 2,
    byrow = TRUE
  )

  fisher_result <- fisher.test(
    contingency_table,
    alternative = "two.sided"
  )

  fisher_odds_ratio[i] <- unname(
    fisher_result$estimate
  )

  fisher_ci_lower[i] <- fisher_result$conf.int[1]
  fisher_ci_upper[i] <- fisher_result$conf.int[2]
  fisher_p_value[i] <- fisher_result$p.value
}

archetypes$high_pip_overlap_fraction <- (
  archetypes$significant_qtl_overlap_true /
    archetypes$n_high_total
)

archetypes$low_pip_overlap_fraction <- (
  archetypes$nonsignificant_qtl_overlap_true /
    archetypes$n_low_total
)

archetypes$overlap_fraction_difference <- (
  archetypes$high_pip_overlap_fraction -
    archetypes$low_pip_overlap_fraction
)

archetypes$haldane_anscombe_odds_ratio <- (
  (
    archetypes$significant_qtl_overlap_true + 0.5
  ) *
    (
      archetypes$
        nonsignificant_qtl_overlap_false + 0.5
    )
) / (
  (
    archetypes$significant_qtl_overlap_false + 0.5
  ) *
    (
      archetypes$
        nonsignificant_qtl_overlap_true + 0.5
    )
)

archetypes$haldane_anscombe_log2_odds_ratio <- log2(
  archetypes$haldane_anscombe_odds_ratio
)

archetypes$fisher_odds_ratio <-
  fisher_odds_ratio

archetypes$fisher_ci_lower <-
  fisher_ci_lower

archetypes$fisher_ci_upper <-
  fisher_ci_upper

archetypes$fisher_p_value <-
  fisher_p_value

# Correction across every archetype × cell-type test
# within each mapping analysis set.
archetypes$fdr_within_analysis_set <- NA_real_

for (analysis_name in unique(archetypes$analysis_set)) {

  selected <- (
    archetypes$analysis_set == analysis_name
  )

  archetypes$fdr_within_analysis_set[selected] <- p.adjust(
    archetypes$fisher_p_value[selected],
    method = "BH"
  )
}

# Separate correction across 286 archetypes within each
# cell type and analysis set.
archetypes$fdr_within_celltype <- NA_real_

celltype_groups <- interaction(
  archetypes$analysis_set,
  archetypes$calderon_celltype,
  drop = TRUE
)

for (current_group in levels(celltype_groups)) {

  selected <- (
    celltype_groups == current_group
  )

  archetypes$fdr_within_celltype[selected] <- p.adjust(
    archetypes$fisher_p_value[selected],
    method = "BH"
  )
}

archetypes$direction <- ifelse(
  archetypes$haldane_anscombe_log2_odds_ratio > 0,
  "enrichment",
  ifelse(
    archetypes$haldane_anscombe_log2_odds_ratio < 0,
    "depletion",
    "no_difference"
  )
)

archetypes$significant_global_fdr_0_05 <- ifelse(
  archetypes$fdr_within_analysis_set < 0.05,
  "YES",
  "NO"
)

archetypes$significant_celltype_fdr_0_05 <- ifelse(
  archetypes$fdr_within_celltype < 0.05,
  "YES",
  "NO"
)

archetypes <- archetypes[
  order(
    match(
      archetypes$analysis_set,
      analysis_order
    ),
    archetypes$calderon_celltype,
    archetypes$fdr_within_celltype,
    archetypes$fisher_p_value,
    archetypes$candidate_archetype
  ),
]

write.table(
  archetypes,
  fisher_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

summary_table <- data.frame(
  metric = c(
    "component_result_files",
    "component_result_rows",
    "unique_component_tasks",
    "candidate_motif_archetypes",
    "analysis_set_celltype_combinations",
    "archetype_contingency_rows",
    "expected_archetype_contingency_rows",
    "archetype_rows_with_any_high_overlap",
    "archetype_rows_with_any_low_overlap",
    "archetype_rows_with_any_overlap"
  ),
  value = c(
    length(result_files),
    nrow(components),
    length(unique(components$task_id)),
    length(unique(archetypes$candidate_archetype)),
    length(
      unique(
        paste(
          archetypes$analysis_set,
          archetypes$calderon_celltype,
          sep = "|"
        )
      )
    ),
    nrow(archetypes),
    expected_archetype_rows,
    sum(
      archetypes$significant_qtl_overlap_true > 0
    ),
    sum(
      archetypes$nonsignificant_qtl_overlap_true > 0
    ),
    sum(
      archetypes$significant_qtl_overlap_true > 0 |
        archetypes$nonsignificant_qtl_overlap_true > 0
    )
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_table,
  summary_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

analysis_set_summary <- do.call(
  rbind,
  lapply(
    analysis_order,
    function(analysis_name) {

      selected <- archetypes[
        archetypes$analysis_set == analysis_name,
      ]

      data.frame(
        analysis_set = analysis_name,
        n_celltypes =
          length(
            unique(selected$calderon_celltype)
          ),
        n_archetype_tests =
          nrow(selected),
        tests_with_high_overlap =
          sum(
            selected$
              significant_qtl_overlap_true > 0
          ),
        tests_with_low_overlap =
          sum(
            selected$
              nonsignificant_qtl_overlap_true > 0
          ),
        positive_effect_tests =
          sum(
            selected$
              haldane_anscombe_log2_odds_ratio > 0
          ),
        negative_effect_tests =
          sum(
            selected$
              haldane_anscombe_log2_odds_ratio < 0
          ),
        global_fdr_significant_tests =
          sum(
            selected$
              fdr_within_analysis_set < 0.05
          ),
        celltype_fdr_significant_tests =
          sum(
            selected$
              fdr_within_celltype < 0.05
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.table(
  analysis_set_summary,
  analysis_set_summary_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

for (analysis_name in analysis_order) {

  selected <- archetypes[
    archetypes$analysis_set == analysis_name &
      archetypes$fdr_within_analysis_set < 0.05,
  ]

  selected <- selected[
    order(
      selected$fdr_within_analysis_set,
      selected$fisher_p_value
    ),
  ]

  output_file <- file.path(
    results_dir,
    paste0(
      "motif_archetype_fisher_significant_",
      analysis_name,
      ".tsv"
    )
  )

  write.table(
    selected,
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

cat("\nMotif-archetype contingency analysis completed.\n\n")

print(
  summary_table,
  row.names = FALSE
)

cat("\nAnalysis-set summary:\n\n")

print(
  analysis_set_summary,
  row.names = FALSE
)

cat("\nMain outputs:\n")
cat(contingency_output, "\n")
cat(fisher_output, "\n")
cat(summary_output, "\n")
cat(analysis_set_summary_output, "\n")

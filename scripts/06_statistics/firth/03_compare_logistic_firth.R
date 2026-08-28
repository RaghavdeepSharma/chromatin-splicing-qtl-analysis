#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

root <- Sys.getenv("KERIMOV_BASE")
if (root == "") stop("KERIMOV_BASE is not set.")

firth_dir <- file.path(
    root,
    "07_logistic_regression/05_firth_comparison/firth_model_results"
)

logistic_file <- file.path(
    root,
    "07_logistic_regression/04_full_analysis/combined_results",
    "all_strict_primary_logistic_results.tsv.gz"
)

out_dir <- file.path(
    root,
    "07_logistic_regression/05_firth_comparison/tables"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Read all Firth result files
# ------------------------------------------------------------

firth_files <- list.files(
    firth_dir,
    pattern = "\\.strict_primary\\.firth\\.tsv$",
    full.names = TRUE
)

if (length(firth_files) != 286) {
    warning(
        "Expected 286 Firth files but found ",
        length(firth_files)
    )
}

firth <- rbindlist(
    lapply(
        firth_files,
        function(f) {
            x <- fread(f)
            x[, source_file := basename(f)]
            x
        }
    ),
    fill = TRUE,
    use.names = TRUE
)

firth[, result_key := paste(
    motif_archetype,
    analysis_set,
    calderon_celltype,
    sep = "|"
)]

# Ensure numeric columns are numeric.
numeric_columns <- c(
    "n_total",
    "n_high",
    "n_low",
    "n_overlap",
    "n_high_overlap",
    "n_high_not_overlap",
    "n_low_overlap",
    "n_low_not_overlap",
    "iterations",
    "coefficient",
    "standard_error",
    "p_value",
    "odds_ratio",
    "ci_lower_95",
    "ci_upper_95",
    "log_likelihood_penalized"
)

for (column in intersect(numeric_columns, names(firth))) {
    set(
        firth,
        j = column,
        value = suppressWarnings(
            as.numeric(firth[[column]])
        )
    )
}

firth[, valid_firth_model :=
    firth_status == "firth_fit_ok" &
    is.finite(coefficient) &
    is.finite(odds_ratio) &
    !is.na(p_value)
]

# ------------------------------------------------------------
# Firth BH adjustment
# Only valid Firth models enter the Firth multiple-testing family.
# ------------------------------------------------------------

firth[, p_adj_firth_global_bh := NA_real_]

firth[
    valid_firth_model == TRUE,
    p_adj_firth_global_bh := p.adjust(
        p_value,
        method = "BH"
    )
]

firth[, p_adj_firth_celltype_bh := NA_real_]

firth[
    valid_firth_model == TRUE,
    p_adj_firth_celltype_bh := p.adjust(
        p_value,
        method = "BH"
    ),
    by = calderon_celltype
]

firth[, significant_firth_global_fdr :=
    valid_firth_model &
    !is.na(p_adj_firth_global_bh) &
    p_adj_firth_global_bh < 0.05
]

firth[, significant_firth_celltype_fdr :=
    valid_firth_model &
    !is.na(p_adj_firth_celltype_bh) &
    p_adj_firth_celltype_bh < 0.05
]

firth[, firth_effect_direction := fifelse(
    is.na(odds_ratio),
    NA_character_,
    fifelse(
        odds_ratio > 1,
        "enrichment",
        fifelse(
            odds_ratio < 1,
            "depletion",
            "null"
        )
    )
)]

firth[, firth_event_support := fifelse(
    n_high_overlap >= 5,
    "at_least_5_high_overlap",
    fifelse(
        n_high_overlap >= 1,
        "one_to_four_high_overlap",
        "zero_high_overlap"
    )
)]

setorder(
    firth,
    calderon_celltype,
    motif_archetype
)

fwrite(
    firth,
    file.path(out_dir, "all_strict_primary_firth_results.tsv"),
    sep = "\t",
    na = "NA"
)

fwrite(
    firth[valid_firth_model == TRUE],
    file.path(out_dir, "valid_strict_primary_firth_results.tsv"),
    sep = "\t",
    na = "NA"
)

fwrite(
    firth[significant_firth_global_fdr == TRUE],
    file.path(out_dir, "significant_global_fdr_firth_results.tsv"),
    sep = "\t",
    na = "NA"
)

# ------------------------------------------------------------
# Read ordinary logistic results and rename model fields
# ------------------------------------------------------------

logistic <- fread(logistic_file)

logistic[, result_key := paste(
    motif_archetype,
    analysis_set,
    calderon_celltype,
    sep = "|"
)]

logistic_keep <- logistic[, .(
    result_key,
    motif_archetype,
    calderon_celltype,
    analysis_set,

    logistic_model_status = model_status,
    logistic_converged = converged,
    logistic_coefficient = coefficient,
    logistic_standard_error = standard_error,
    logistic_p_value = p_value,
    logistic_odds_ratio = odds_ratio,
    logistic_ci_lower_95 = ci_lower_95,
    logistic_ci_upper_95 = ci_upper_95,
    logistic_p_adj_global_bh = p_adj_global_bh,
    logistic_p_adj_celltype_bh = p_adj_celltype_bh,
    logistic_significant_global_fdr =
        significant_global_fdr,
    logistic_significant_celltype_fdr =
        significant_celltype_fdr,
    logistic_valid_primary_model =
        valid_primary_model
)]

firth_keep <- firth[, .(
    result_key,

    n_total,
    n_high,
    n_low,
    n_overlap,
    n_high_overlap,
    n_high_not_overlap,
    n_low_overlap,
    n_low_not_overlap,

    firth_status,
    firth_converged = converged,
    firth_iterations = iterations,
    firth_coefficient = coefficient,
    firth_standard_error = standard_error,
    firth_p_value = p_value,
    firth_odds_ratio = odds_ratio,
    firth_ci_lower_95 = ci_lower_95,
    firth_ci_upper_95 = ci_upper_95,
    firth_p_adj_global_bh = p_adj_firth_global_bh,
    firth_p_adj_celltype_bh = p_adj_firth_celltype_bh,
    firth_significant_global_fdr =
        significant_firth_global_fdr,
    firth_significant_celltype_fdr =
        significant_firth_celltype_fdr,
    valid_firth_model,
    firth_effect_direction,
    firth_event_support,
    firth_warning_text = warning_text
)]

comparison <- merge(
    logistic_keep,
    firth_keep,
    by = "result_key",
    all = TRUE
)

comparison[, both_models_valid :=
    logistic_valid_primary_model == TRUE &
    valid_firth_model == TRUE
]

comparison[, same_effect_direction :=
    both_models_valid &
    sign(logistic_coefficient) ==
    sign(firth_coefficient)
]

comparison[, log_or_difference :=
    firth_coefficient - logistic_coefficient
]

comparison[, comparison_category := fifelse(
    both_models_valid,
    "both_valid",
    fifelse(
        logistic_valid_primary_model != TRUE &
        valid_firth_model == TRUE,
        "firth_rescued_or_only",
        fifelse(
            logistic_valid_primary_model == TRUE &
            valid_firth_model != TRUE,
            "logistic_only",
            "neither_valid"
        )
    )
)]

setorder(
    comparison,
    calderon_celltype,
    motif_archetype
)

fwrite(
    comparison,
    file.path(
        out_dir,
        "logistic_vs_firth_side_by_side.tsv"
    ),
    sep = "\t",
    na = "NA"
)

# ------------------------------------------------------------
# Summaries
# ------------------------------------------------------------

status_summary <- comparison[, .(
    n_models = .N,
    n_logistic_valid =
        sum(logistic_valid_primary_model == TRUE, na.rm = TRUE),
    n_firth_valid =
        sum(valid_firth_model == TRUE, na.rm = TRUE),
    n_both_valid =
        sum(both_models_valid == TRUE, na.rm = TRUE),
    n_firth_rescued_or_only =
        sum(comparison_category == "firth_rescued_or_only",
            na.rm = TRUE),
    n_same_direction_among_both =
        sum(same_effect_direction == TRUE, na.rm = TRUE),
    n_firth_global_fdr =
        sum(firth_significant_global_fdr == TRUE, na.rm = TRUE),
    n_firth_global_fdr_with_5_high_overlap =
        sum(
            firth_significant_global_fdr == TRUE &
            n_high_overlap >= 5,
            na.rm = TRUE
        )
), by = calderon_celltype]

setorder(status_summary, calderon_celltype)

fwrite(
    status_summary,
    file.path(
        out_dir,
        "logistic_vs_firth_status_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
)

agreement_summary <- comparison[
    both_models_valid == TRUE,
    .(
        n_models = .N,
        pearson_log_or = cor(
            logistic_coefficient,
            firth_coefficient,
            use = "complete.obs",
            method = "pearson"
        ),
        spearman_log_or = cor(
            logistic_coefficient,
            firth_coefficient,
            use = "complete.obs",
            method = "spearman"
        ),
        median_absolute_log_or_difference =
            median(
                abs(log_or_difference),
                na.rm = TRUE
            ),
        same_direction_fraction =
            mean(
                same_effect_direction,
                na.rm = TRUE
            )
    )
]

fwrite(
    agreement_summary,
    file.path(
        out_dir,
        "logistic_vs_firth_agreement_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
)

cat("Firth files read:", length(firth_files), "\n")
cat("Total Firth rows:", nrow(firth), "\n")
cat(
    "Valid Firth models:",
    sum(firth$valid_firth_model, na.rm = TRUE),
    "\n"
)
cat(
    "Firth global-FDR significant:",
    sum(
        firth$significant_firth_global_fdr,
        na.rm = TRUE
    ),
    "\n"
)
cat(
    "Firth global-FDR significant with >=5 high overlaps:",
    sum(
        firth$significant_firth_global_fdr &
        firth$n_high_overlap >= 5,
        na.rm = TRUE
    ),
    "\n"
)
cat(
    "Models valid in both methods:",
    sum(comparison$both_models_valid, na.rm = TRUE),
    "\n"
)

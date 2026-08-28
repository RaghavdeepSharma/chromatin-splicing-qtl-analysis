#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop(
        "Usage: Rscript 13_combine_and_summarize_full_results.R ",
        "<model_results_directory> ",
        "<fisher_results.tsv> ",
        "<output_directory>"
    )
}

model_dir <- args[1]
fisher_file <- args[2]
output_dir <- args[3]

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

############################################################
# Read all logistic-regression results
############################################################

model_files <- list.files(
    model_dir,
    pattern = "^array_[0-9]+\\..*\\.strict_primary\\.logistic\\.tsv$",
    full.names = TRUE
)

model_files <- sort(model_files)

cat("Model files found:", length(model_files), "\n")

if (length(model_files) != 286) {
    stop(
        "Expected 286 model files; found ",
        length(model_files)
    )
}

result_list <- lapply(
    model_files,
    function(file) {
        dat <- fread(file)

        if (nrow(dat) != 6) {
            stop(
                "Expected 6 rows in ",
                file,
                "; found ",
                nrow(dat)
            )
        }

        dat[, source_file := basename(file)]

        dat
    }
)

results <- rbindlist(
    result_list,
    use.names = TRUE,
    fill = TRUE
)

cat("Combined model rows:", nrow(results), "\n")

if (nrow(results) != 1716) {
    stop(
        "Expected 1716 total model rows; found ",
        nrow(results)
    )
}

############################################################
# Validate unique motif × cell-type keys
############################################################

results[
    ,
    result_key :=
        paste(
            motif_archetype,
            analysis_set,
            calderon_celltype,
            sep = "|"
        )
]

if (anyDuplicated(results$result_key)) {
    stop(
        "Duplicate motif × analysis-set × cell-type results detected."
    )
}

############################################################
# Define valid ordinary logistic-regression models
############################################################

results[
    ,
    valid_primary_model :=
        model_status == "fit_ok" &
        converged == TRUE &
        is.finite(p_value) &
        is.finite(odds_ratio) &
        is.finite(ci_lower_95) &
        is.finite(ci_upper_95)
]

############################################################
# Multiple-testing corrections
############################################################

results[, p_adj_global_bh := NA_real_]
results[, p_adj_celltype_bh := NA_real_]

results[
    valid_primary_model == TRUE,
    p_adj_global_bh := p.adjust(
        p_value,
        method = "BH"
    )
]

results[
    valid_primary_model == TRUE,
    p_adj_celltype_bh := p.adjust(
        p_value,
        method = "BH"
    ),
    by = calderon_celltype
]

results[
    ,
    significant_global_fdr :=
        valid_primary_model == TRUE &
        !is.na(p_adj_global_bh) &
        p_adj_global_bh < 0.05
]

results[
    ,
    significant_celltype_fdr :=
        valid_primary_model == TRUE &
        !is.na(p_adj_celltype_bh) &
        p_adj_celltype_bh < 0.05
]

results[
    ,
    effect_direction :=
        fifelse(
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
        )
]

results[
    ,
    max_vif := pmax(
        vif_motif_overlap,
        vif_splice_distance,
        vif_tss_distance,
        vif_gc_fraction,
        na.rm = TRUE
    )
]

results[
    !is.finite(max_vif),
    max_vif := NA_real_
]

############################################################
# Save complete logistic-regression table
############################################################

setorder(
    results,
    p_adj_global_bh,
    p_value,
    motif_archetype,
    calderon_celltype,
    na.last = TRUE
)

complete_output <- file.path(
    output_dir,
    "all_strict_primary_logistic_results.tsv.gz"
)

fwrite(
    results,
    complete_output,
    sep = "\t",
    compress = "gzip"
)

############################################################
# Save valid primary results
############################################################

valid_results <- results[
    valid_primary_model == TRUE
]

fwrite(
    valid_results,
    file.path(
        output_dir,
        "valid_strict_primary_logistic_results.tsv"
    ),
    sep = "\t"
)

############################################################
# Save global-FDR significant results
############################################################

significant_global <- results[
    significant_global_fdr == TRUE
]

setorder(
    significant_global,
    p_adj_global_bh,
    p_value
)

fwrite(
    significant_global,
    file.path(
        output_dir,
        "significant_global_fdr_logistic_results.tsv"
    ),
    sep = "\t"
)

############################################################
# Save cell-type-FDR significant results
############################################################

significant_celltype <- results[
    significant_celltype_fdr == TRUE
]

setorder(
    significant_celltype,
    calderon_celltype,
    p_adj_celltype_bh,
    p_value
)

fwrite(
    significant_celltype,
    file.path(
        output_dir,
        "significant_celltype_fdr_logistic_results.tsv"
    ),
    sep = "\t"
)

############################################################
# Model-status summary
############################################################

status_summary <- results[
    ,
    .(
        n_models = .N,
        n_unique_motifs = uniqueN(motif_archetype),
        n_significant_global_fdr = sum(
            significant_global_fdr,
            na.rm = TRUE
        ),
        n_significant_celltype_fdr = sum(
            significant_celltype_fdr,
            na.rm = TRUE
        )
    ),
    by = .(
        calderon_celltype,
        model_status
    )
]

setorder(
    status_summary,
    calderon_celltype,
    model_status
)

fwrite(
    status_summary,
    file.path(
        output_dir,
        "model_status_summary.tsv"
    ),
    sep = "\t"
)

############################################################
# VIF summary
############################################################

vif_summary <- results[
    valid_primary_model == TRUE,
    .(
        n_valid_models = .N,
        median_max_vif = median(
            max_vif,
            na.rm = TRUE
        ),
        maximum_vif = max(
            max_vif,
            na.rm = TRUE
        ),
        n_vif_above_5 = sum(
            max_vif > 5,
            na.rm = TRUE
        ),
        n_vif_above_10 = sum(
            max_vif > 10,
            na.rm = TRUE
        )
    ),
    by = calderon_celltype
]

fwrite(
    vif_summary,
    file.path(
        output_dir,
        "vif_summary_by_celltype.tsv"
    ),
    sep = "\t"
)

############################################################
# Event and overlap summary
############################################################

overlap_summary <- results[
    ,
    .(
        n_models = .N,
        median_overlap = median(n_overlap),
        minimum_overlap = min(n_overlap),
        maximum_overlap = max(n_overlap),
        n_zero_overlap = sum(n_overlap == 0),
        n_zero_high_overlap = sum(n_high_overlap == 0),
        n_with_at_least_5_high_overlap =
            sum(n_high_overlap >= 5),
        n_with_at_least_10_high_overlap =
            sum(n_high_overlap >= 10)
    ),
    by = calderon_celltype
]

fwrite(
    overlap_summary,
    file.path(
        output_dir,
        "overlap_event_summary_by_celltype.tsv"
    ),
    sep = "\t"
)

############################################################
# Optional Fisher comparison
############################################################

if (file.exists(fisher_file)) {

    fisher <- fread(fisher_file)

    cat(
        "Fisher rows read:",
        nrow(fisher),
        "\n"
    )

    candidate_motif_columns <- c(
        "candidate_archetype",
        "motif_archetype"
    )

    motif_column <- candidate_motif_columns[
        candidate_motif_columns %in% names(fisher)
    ][1]

    if (is.na(motif_column)) {
        warning(
            "Could not find motif column in Fisher table. ",
            "Skipping Fisher comparison."
        )
    } else if (
        !"calderon_celltype" %in% names(fisher) ||
        !"analysis_set" %in% names(fisher)
    ) {
        warning(
            "Fisher table lacks cell type or analysis set. ",
            "Skipping Fisher comparison."
        )
    } else {

        setnames(
            fisher,
            motif_column,
            "motif_archetype"
        )

        fisher <- fisher[
            analysis_set == "strict_primary"
        ]

        fisher[
            ,
            result_key :=
                paste(
                    motif_archetype,
                    analysis_set,
                    calderon_celltype,
                    sep = "|"
                )
        ]

        fisher_keep_candidates <- c(
            "result_key",
            "motif_archetype",
            "analysis_set",
            "calderon_celltype",
            "odds_ratio",
            "p_value",
            "p_adj",
            "fdr",
            "fisher_odds_ratio",
            "fisher_p_value",
            "p_adjusted"
        )

        fisher_keep <- intersect(
            fisher_keep_candidates,
            names(fisher)
        )

        fisher_subset <- fisher[
            ,
            ..fisher_keep
        ]

        duplicate_names <- intersect(
            setdiff(names(fisher_subset), "result_key"),
            names(results)
        )

        for (column in duplicate_names) {
            setnames(
                fisher_subset,
                column,
                paste0("fisher_", column)
            )
        }

        comparison <- merge(
            results,
            fisher_subset,
            by = "result_key",
            all.x = TRUE
        )

        fwrite(
            comparison,
            file.path(
                output_dir,
                "fisher_vs_logistic_comparison.tsv.gz"
            ),
            sep = "\t",
            compress = "gzip"
        )
    }
}

############################################################
# Console summary
############################################################

cat("\nFINAL SUMMARY\n")
cat("Total models:", nrow(results), "\n")
cat(
    "Valid fit_ok models:",
    sum(results$valid_primary_model),
    "\n"
)
cat(
    "Separation-risk models:",
    sum(
        results$model_status ==
            "fit_complete_or_quasi_separation_risk"
    ),
    "\n"
)
cat(
    "No-variation models:",
    sum(
        results$model_status ==
            "skip_motif_no_variation"
    ),
    "\n"
)
cat(
    "Global-FDR significant:",
    sum(results$significant_global_fdr),
    "\n"
)
cat(
    "Cell-type-FDR significant:",
    sum(results$significant_celltype_fdr),
    "\n"
)

cat("\nModel-status table:\n")
print(
    results[
        ,
        .N,
        by = model_status
    ][
        order(-N)
    ]
)

cat("\nVIF summary:\n")
print(vif_summary)

cat("\nOutputs written to:\n")
cat(output_dir, "\n")

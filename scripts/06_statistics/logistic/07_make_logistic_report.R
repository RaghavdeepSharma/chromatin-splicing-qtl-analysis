#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

base_dir <- "04_full_analysis"

output_dir <- file.path(
    base_dir,
    "final_report"
)

plot_dir <- file.path(
    output_dir,
    "plots"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    plot_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

############################################################
# Input files
############################################################

input_files <- c(
    full = file.path(
        base_dir,
        "combined_results",
        "all_strict_primary_logistic_results.tsv.gz"
    ),
    no_tss = file.path(
        base_dir,
        "no_tss_combined_results",
        "all_strict_primary_no_tss_logistic_results.tsv.gz"
    ),
    no_gc = file.path(
        base_dir,
        "no_gc_combined_results",
        "all_strict_primary_no_gc_logistic_results.tsv.gz"
    ),
    no_splice = file.path(
        base_dir,
        "no_splice_combined_results",
        "all_strict_primary_no_splice_logistic_results.tsv.gz"
    )
)

missing_files <- input_files[
    !file.exists(input_files)
]

if (length(missing_files) > 0) {
    stop(
        "Missing input files:\n",
        paste(missing_files, collapse = "\n")
    )
}

############################################################
# Read and standardize models
############################################################

read_model <- function(file, model_name) {

    dat <- fread(
        cmd = paste(
            "gzip -dc",
            shQuote(file)
        )
    )

    dat[
        ,
        model_specification := model_name
    ]

    dat[
        ,
        comparison_key := paste(
            motif_archetype,
            calderon_celltype,
            analysis_set,
            sep = "|"
        )
    ]

    dat[
        ,
        valid_model :=
            model_status == "fit_ok" &
            converged == TRUE &
            is.finite(odds_ratio) &
            is.finite(p_value)
    ]

    dat
}

model_list <- Map(
    read_model,
    input_files,
    names(input_files)
)

all_models_long <- rbindlist(
    model_list,
    use.names = TRUE,
    fill = TRUE
)

model_labels <- c(
    full = "Full model",
    no_tss = "No TSS",
    no_gc = "No GC",
    no_splice = "No splice distance"
)

all_models_long[
    ,
    model_label := factor(
        model_labels[model_specification],
        levels = unname(model_labels)
    )
]

############################################################
# Save long comparison table
############################################################

long_keep <- c(
    "comparison_key",
    "motif_archetype",
    "calderon_celltype",
    "analysis_set",
    "model_specification",
    "model_label",
    "model_status",
    "converged",
    "n_total",
    "n_high",
    "n_low",
    "n_overlap",
    "n_high_overlap",
    "n_low_overlap",
    "coefficient",
    "standard_error",
    "odds_ratio",
    "ci_lower_95",
    "ci_upper_95",
    "p_value",
    "p_adj_global_bh",
    "p_adj_celltype_bh",
    "significant_global_fdr",
    "significant_celltype_fdr",
    "max_vif",
    "valid_model"
)

long_keep <- intersect(
    long_keep,
    names(all_models_long)
)

fwrite(
    all_models_long[
        ,
        ..long_keep
    ],
    file.path(
        output_dir,
        "all_models_sensitivity_comparison_long.tsv.gz"
    ),
    sep = "\t",
    compress = "gzip"
)

############################################################
# Wide comparison table
############################################################

wide_source <- all_models_long[
    ,
    .(
        comparison_key,
        motif_archetype,
        calderon_celltype,
        analysis_set,
        model_specification,
        model_status,
        n_high_overlap,
        n_low_overlap,
        odds_ratio,
        ci_lower_95,
        ci_upper_95,
        p_value,
        p_adj_global_bh,
        significant_global_fdr,
        valid_model
    )
]

wide <- dcast(
    wide_source,
    comparison_key +
        motif_archetype +
        calderon_celltype +
        analysis_set +
        n_high_overlap +
        n_low_overlap ~
        model_specification,
    value.var = c(
        "model_status",
        "odds_ratio",
        "ci_lower_95",
        "ci_upper_95",
        "p_value",
        "p_adj_global_bh",
        "significant_global_fdr",
        "valid_model"
    )
)

wide[
    ,
    consistent_direction :=
        fifelse(
            is.finite(odds_ratio_full) &
            is.finite(odds_ratio_no_tss) &
            is.finite(odds_ratio_no_gc) &
            is.finite(odds_ratio_no_splice),
            (
                odds_ratio_full > 1 &
                odds_ratio_no_tss > 1 &
                odds_ratio_no_gc > 1 &
                odds_ratio_no_splice > 1
            ) |
            (
                odds_ratio_full < 1 &
                odds_ratio_no_tss < 1 &
                odds_ratio_no_gc < 1 &
                odds_ratio_no_splice < 1
            ),
            NA
        )
]

wide[
    ,
    significant_all_four :=
        significant_global_fdr_full == TRUE &
        significant_global_fdr_no_tss == TRUE &
        significant_global_fdr_no_gc == TRUE &
        significant_global_fdr_no_splice == TRUE
]

wide[
    ,
    robust_event_support :=
        n_high_overlap >= 5
]

setorder(
    wide,
    -significant_all_four,
    -robust_event_support,
    p_adj_global_bh_full,
    p_value_full
)

fwrite(
    wide,
    file.path(
        output_dir,
        "all_models_sensitivity_comparison_wide.tsv"
    ),
    sep = "\t"
)

############################################################
# Primary significant table
############################################################

primary_significant <- all_models_long[
    model_specification == "full" &
    significant_global_fdr == TRUE
]

primary_significant[
    ,
    evidence_support := fifelse(
        n_high_overlap >= 5,
        "At least 5 High-PIP overlaps",
        "Fewer than 5 High-PIP overlaps"
    )
]

setorder(
    primary_significant,
    p_adj_global_bh,
    p_value
)

fwrite(
    primary_significant,
    file.path(
        output_dir,
        "primary_full_model_significant_results.tsv"
    ),
    sep = "\t"
)

############################################################
# Plot 1: model status
############################################################

status_counts <- all_models_long[
    model_specification == "full",
    .N,
    by = model_status
]

status_labels <- c(
    skip_motif_no_variation = "No motif-overlap variation",
    fit_complete_or_quasi_separation_risk = "Separation risk",
    fit_ok = "Fit successfully"
)

status_counts[
    ,
    status_label := status_labels[model_status]
]

status_counts[
    ,
    status_label := factor(
        status_label,
        levels = c(
            "No motif-overlap variation",
            "Separation risk",
            "Fit successfully"
        )
    )
]

p_status <- ggplot(
    status_counts,
    aes(
        x = status_label,
        y = N
    )
) +
    geom_col(
        width = 0.7
    ) +
    geom_text(
        aes(label = N),
        vjust = -0.4,
        size = 4
    ) +
    labs(
        title = "Outcome of 1,716 planned logistic-regression models",
        x = NULL,
        y = "Number of motif–cell-type models"
    ) +
    theme_bw(base_size = 12) +
    theme(
        axis.text.x = element_text(
            angle = 20,
            hjust = 1
        ),
        panel.grid.major.x = element_blank()
    ) +
    expand_limits(
        y = max(status_counts$N) * 1.08
    )

ggsave(
    file.path(
        plot_dir,
        "01_model_status_summary.png"
    ),
    p_status,
    width = 8,
    height = 5,
    dpi = 300
)

ggsave(
    file.path(
        plot_dir,
        "01_model_status_summary.pdf"
    ),
    p_status,
    width = 8,
    height = 5
)

############################################################
# Plot 2: forest plot
############################################################

forest <- copy(primary_significant)

forest[
    ,
    display_label := paste0(
        motif_archetype,
        " — ",
        calderon_celltype,
        " (High-PIP overlaps: ",
        n_high_overlap,
        ")"
    )
]

forest[
    ,
    display_label := factor(
        display_label,
        levels = rev(display_label)
    )
]

p_forest <- ggplot(
    forest,
    aes(
        x = odds_ratio,
        y = display_label,
        shape = evidence_support
    )
) +
    geom_vline(
        xintercept = 1,
        linetype = "dashed"
    ) +
    geom_errorbarh(
        aes(
            xmin = ci_lower_95,
            xmax = ci_upper_95
        ),
        height = 0.2
    ) +
    geom_point(
        size = 3
    ) +
    scale_x_log10() +
    labs(
        title = "Full-model motif associations passing global FDR",
        subtitle = "Adjusted for splice distance, gene-specific nearest-TSS distance and GC content",
        x = "Adjusted odds ratio, log scale",
        y = NULL,
        shape = "Evidence support"
    ) +
    theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
    )

ggsave(
    file.path(
        plot_dir,
        "02_primary_significant_forest_plot.png"
    ),
    p_forest,
    width = 10,
    height = 6,
    dpi = 300
)

ggsave(
    file.path(
        plot_dir,
        "02_primary_significant_forest_plot.pdf"
    ),
    p_forest,
    width = 10,
    height = 6
)

############################################################
# Plot 3: sensitivity comparison
############################################################

sensitivity_keys <- primary_significant$comparison_key

sensitivity <- all_models_long[
    comparison_key %in% sensitivity_keys &
    valid_model == TRUE
]

sensitivity[
    ,
    display_label := paste0(
        motif_archetype,
        " — ",
        calderon_celltype
    )
]

sensitivity[
    ,
    display_label := factor(
        display_label,
        levels = unique(
            primary_significant[
                order(p_adj_global_bh),
                paste0(
                    motif_archetype,
                    " — ",
                    calderon_celltype
                )
            ]
        )
    )
]

p_sensitivity <- ggplot(
    sensitivity,
    aes(
        x = model_label,
        y = odds_ratio,
        group = display_label
    )
) +
    geom_hline(
        yintercept = 1,
        linetype = "dashed"
    ) +
    geom_line() +
    geom_point(
        aes(
            shape = n_high_overlap >= 5
        ),
        size = 2.7
    ) +
    scale_y_log10() +
    facet_wrap(
        ~ display_label,
        scales = "free_y",
        ncol = 2
    ) +
    labs(
        title = "Sensitivity of motif odds ratios to covariate removal",
        subtitle = "Each panel shows the same motif–cell-type association under four model specifications",
        x = NULL,
        y = "Odds ratio, log scale",
        shape = "At least 5 High-PIP overlaps"
    ) +
    theme_bw(base_size = 10) +
    theme(
        axis.text.x = element_text(
            angle = 35,
            hjust = 1
        ),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
    )

ggsave(
    file.path(
        plot_dir,
        "03_covariate_sensitivity_plot.png"
    ),
    p_sensitivity,
    width = 11,
    height = 10,
    dpi = 300
)

ggsave(
    file.path(
        plot_dir,
        "03_covariate_sensitivity_plot.pdf"
    ),
    p_sensitivity,
    width = 11,
    height = 10
)

############################################################
# Summary statistics
############################################################

model_summary <- all_models_long[
    ,
    .(
        total_models = .N,
        fit_ok = sum(model_status == "fit_ok"),
        separation_risk = sum(
            model_status ==
                "fit_complete_or_quasi_separation_risk"
        ),
        no_variation = sum(
            model_status ==
                "skip_motif_no_variation"
        ),
        global_fdr_significant = sum(
            significant_global_fdr,
            na.rm = TRUE
        ),
        celltype_fdr_significant = sum(
            significant_celltype_fdr,
            na.rm = TRUE
        ),
        significant_with_5_high_overlap = sum(
            significant_global_fdr == TRUE &
            n_high_overlap >= 5,
            na.rm = TRUE
        )
    ),
    by = .(
        model_specification,
        model_label
    )
]

fwrite(
    model_summary,
    file.path(
        output_dir,
        "model_specification_summary.tsv"
    ),
    sep = "\t"
)

############################################################
# SMAD_1 focused table
############################################################

smad <- all_models_long[
    motif_archetype == "SMAD_1" &
    calderon_celltype == "Monocytes",
    .(
        model_specification,
        model_label,
        n_high_overlap,
        n_low_overlap,
        odds_ratio,
        ci_lower_95,
        ci_upper_95,
        p_value,
        p_adj_global_bh,
        significant_global_fdr,
        model_status
    )
]

fwrite(
    smad,
    file.path(
        output_dir,
        "SMAD_1_monocytes_sensitivity.tsv"
    ),
    sep = "\t"
)

############################################################
# Written report
############################################################

full_summary <- model_summary[
    model_specification == "full"
]

smad_full <- smad[
    model_specification == "full"
]

smad_all_significant <- all(
    smad$significant_global_fdr == TRUE
)

report_lines <- c(
    "FINAL LOGISTIC-REGRESSION ANALYSIS SUMMARY",
    "",
    "Primary analysis",
    paste0(
        "A total of ",
        full_summary$total_models,
        " motif–cell-type combinations were attempted ",
        "(286 motif archetypes across six strict-primary cell types)."
    ),
    paste0(
        full_summary$fit_ok,
        " models were successfully estimated using ordinary logistic regression, ",
        full_summary$separation_risk,
        " showed complete or quasi-separation risk, and ",
        full_summary$no_variation,
        " had no motif-overlap variation."
    ),
    paste0(
        full_summary$global_fdr_significant,
        " successfully estimated models passed global Benjamini–Hochberg FDR correction."
    ),
    "",
    "Primary model",
    paste0(
        "High-PIP status was modeled as a function of accessible motif overlap, ",
        "log-transformed splice-site distance, log-transformed distance to the nearest ",
        "transcript TSS of the associated gene, and local GC fraction."
    ),
    "",
    "Main finding",
    paste0(
        "SMAD_1 in Monocytes was the only full-model significant association supported ",
        "by at least five High-PIP motif overlaps. It had an adjusted odds ratio of ",
        formatC(
            smad_full$odds_ratio,
            digits = 3,
            format = "f"
        ),
        " with a 95% confidence interval of ",
        formatC(
            smad_full$ci_lower_95,
            digits = 3,
            format = "f"
        ),
        " to ",
        formatC(
            smad_full$ci_upper_95,
            digits = 3,
            format = "f"
        ),
        "."
    ),
    paste0(
        "The association was based on ",
        smad_full$n_high_overlap,
        " High-PIP overlaps and ",
        smad_full$n_low_overlap,
        " Low-PIP overlaps."
    ),
    "",
    "Sensitivity analysis",
    paste0(
        "The SMAD_1 Monocyte association remained globally FDR-significant under all ",
        "four model specifications: full, no TSS, no GC and no splice distance. ",
        "Stability across these models indicates that the result was not driven by ",
        "a single covariate."
    ),
    "",
    "Interpretation",
    paste0(
        "Accessible SMAD_1 motif overlap in Monocytes was associated with increased odds ",
        "of High-PIP sQTL status after genomic-context adjustment. This is an association ",
        "and does not establish that SMAD_1 binding directly causes the splicing effect."
    ),
    "",
    "Limitations",
    paste0(
        "Most planned models were not estimable because accessible motif overlap was sparse. ",
        "Four cell types contained only 11 to 14 total High-PIP variants, limiting motif-level ",
        "ordinary logistic regression."
    ),
    paste0(
        "Six of the seven primary global-FDR findings were supported by only one or two ",
        "High-PIP overlaps and should therefore be treated as exploratory."
    ),
    paste0(
        "The five-High-PIP-overlap criterion was used as a practical evidence-support ",
        "threshold rather than as a prespecified formal statistical cutoff."
    ),
    paste0(
        "The TSS covariate represents absolute distance to the nearest transcript TSS ",
        "of the associated gene and is not a strict strand-aware upstream-only TSS distance."
    ),
    paste0(
        "Ordinary logistic-regression estimates from complete or quasi-separated models ",
        "were excluded from primary interpretation."
    )
)

writeLines(
    report_lines,
    file.path(
        output_dir,
        "final_logistic_regression_summary.txt"
    )
)

cat("\nFINAL OUTPUTS\n")
cat("Output directory:", output_dir, "\n")
cat("Comparison rows:", nrow(wide), "\n")
cat("Primary significant models:", nrow(primary_significant), "\n")
cat("SMAD_1 significant in all four models:", smad_all_significant, "\n")
cat("\nModel summary:\n")
print(model_summary)
cat("\nSMAD_1 sensitivity:\n")
print(smad)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop(
        "Usage: Rscript 09_qc_vif_and_full_logistic.R ",
        "<input.tsv.gz> <output_prefix> <plot_directory>"
    )
}

input_file <- args[1]
output_prefix <- args[2]
plot_dir <- args[3]

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("Reading input:\n", input_file, "\n")

dat <- fread(
    cmd = paste("gzip -dc", shQuote(input_file))
)

required <- c(
    "high_pip",
    "motif_overlap",
    "splice_distance",
    "log1p_splice_distance",
    "nearest_gene_tss_distance",
    "log1p_nearest_gene_tss_distance",
    "gc_fraction"
)

missing_columns <- setdiff(required, names(dat))

if (length(missing_columns) > 0) {
    stop(
        "Missing required columns: ",
        paste(missing_columns, collapse = ", ")
    )
}

dat[, high_pip := as.integer(high_pip)]
dat[, motif_overlap := as.integer(motif_overlap)]

numeric_columns <- c(
    "splice_distance",
    "log1p_splice_distance",
    "nearest_gene_tss_distance",
    "log1p_nearest_gene_tss_distance",
    "gc_fraction"
)

for (column in numeric_columns) {
    set(
        dat,
        j = column,
        value = as.numeric(dat[[column]])
    )
}

model_columns <- c(
    "high_pip",
    "motif_overlap",
    "log1p_splice_distance",
    "log1p_nearest_gene_tss_distance",
    "gc_fraction"
)

model_dat <- dat[
    complete.cases(dat[, ..model_columns])
]

cat("Rows in input:", nrow(dat), "\n")
cat("Rows used in model:", nrow(model_dat), "\n")
cat("Rows removed:", nrow(dat) - nrow(model_dat), "\n")

cat("\nOutcome and motif table:\n")
print(
    table(
        high_pip = model_dat$high_pip,
        motif_overlap = model_dat$motif_overlap
    )
)

############################################################
# Distribution summaries
############################################################

distribution_summary <- rbindlist(
    lapply(
        numeric_columns,
        function(column) {
            x <- dat[[column]]

            data.table(
                variable = column,
                n = sum(!is.na(x)),
                missing = sum(is.na(x)),
                minimum = min(x, na.rm = TRUE),
                q1 = quantile(x, 0.25, na.rm = TRUE),
                median = median(x, na.rm = TRUE),
                mean = mean(x, na.rm = TRUE),
                q3 = quantile(x, 0.75, na.rm = TRUE),
                maximum = max(x, na.rm = TRUE),
                standard_deviation = sd(x, na.rm = TRUE)
            )
        }
    )
)

fwrite(
    distribution_summary,
    paste0(output_prefix, ".distribution_summary.tsv"),
    sep = "\t"
)

############################################################
# Histograms
############################################################

make_histogram <- function(column, title, filename, bins = 60) {
    p <- ggplot(
        dat,
        aes(x = .data[[column]])
    ) +
        geom_histogram(
            bins = bins
        ) +
        labs(
            title = title,
            x = column,
            y = "Number of variants"
        ) +
        theme_minimal(base_size = 12)

    ggsave(
        file.path(plot_dir, filename),
        p,
        width = 7,
        height = 5,
        dpi = 300
    )
}

make_histogram(
    "splice_distance",
    "Raw distance to nearest splice boundary",
    "splice_distance_raw.png"
)

make_histogram(
    "log1p_splice_distance",
    "Log-transformed distance to nearest splice boundary",
    "splice_distance_log1p.png"
)

make_histogram(
    "nearest_gene_tss_distance",
    "Raw nearest gene-specific TSS distance",
    "nearest_gene_tss_distance_raw.png"
)

make_histogram(
    "log1p_nearest_gene_tss_distance",
    "Log-transformed nearest gene-specific TSS distance",
    "nearest_gene_tss_distance_log1p.png"
)

make_histogram(
    "gc_fraction",
    "GC content in ±100 bp window",
    "gc_fraction.png"
)

############################################################
# Correlations
############################################################

predictor_columns <- c(
    "motif_overlap",
    "log1p_splice_distance",
    "log1p_nearest_gene_tss_distance",
    "gc_fraction"
)

correlation_matrix <- cor(
    model_dat[, ..predictor_columns],
    use = "pairwise.complete.obs",
    method = "pearson"
)

fwrite(
    as.data.table(
        correlation_matrix,
        keep.rownames = "variable"
    ),
    paste0(output_prefix, ".correlation_matrix.tsv"),
    sep = "\t"
)

############################################################
# Manual VIF calculation
#
# For each predictor:
#   regress predictor against all other predictors
#   VIF = 1 / (1 - R-squared)
############################################################

calculate_vif <- function(data, predictors) {
    results <- lapply(
        predictors,
        function(response_variable) {
            other_predictors <- setdiff(
                predictors,
                response_variable
            )

            formula_text <- paste(
                response_variable,
                "~",
                paste(other_predictors, collapse = " + ")
            )

            fit <- lm(
                as.formula(formula_text),
                data = data
            )

            r_squared <- summary(fit)$r.squared

            data.table(
                variable = response_variable,
                r_squared = r_squared,
                vif = 1 / (1 - r_squared)
            )
        }
    )

    rbindlist(results)
}

vif_table <- calculate_vif(
    model_dat,
    predictor_columns
)

fwrite(
    vif_table,
    paste0(output_prefix, ".vif.tsv"),
    sep = "\t"
)

############################################################
# Logistic models
############################################################

model_unadjusted <- glm(
    high_pip ~ motif_overlap,
    data = model_dat,
    family = binomial(link = "logit")
)

model_splice <- glm(
    high_pip ~
        motif_overlap +
        log1p_splice_distance,
    data = model_dat,
    family = binomial(link = "logit")
)

model_full <- glm(
    high_pip ~
        motif_overlap +
        log1p_splice_distance +
        log1p_nearest_gene_tss_distance +
        gc_fraction,
    data = model_dat,
    family = binomial(link = "logit")
)

extract_coefficients <- function(model, model_name) {
    coefficients <- summary(model)$coefficients

    data.table(
        model = model_name,
        term = rownames(coefficients),
        coefficient = coefficients[, "Estimate"],
        standard_error = coefficients[, "Std. Error"],
        z_value = coefficients[, "z value"],
        p_value = coefficients[, "Pr(>|z|)"],
        odds_ratio = exp(coefficients[, "Estimate"]),
        ci_lower_95 = exp(
            coefficients[, "Estimate"] -
                1.96 * coefficients[, "Std. Error"]
        ),
        ci_upper_95 = exp(
            coefficients[, "Estimate"] +
                1.96 * coefficients[, "Std. Error"]
        )
    )
}

coefficient_results <- rbindlist(
    list(
        extract_coefficients(
            model_unadjusted,
            "unadjusted"
        ),
        extract_coefficients(
            model_splice,
            "splice_adjusted"
        ),
        extract_coefficients(
            model_full,
            "full_adjusted"
        )
    )
)

fwrite(
    coefficient_results,
    paste0(output_prefix, ".coefficients.tsv"),
    sep = "\t"
)

diagnostics <- data.table(
    model = c(
        "unadjusted",
        "splice_adjusted",
        "full_adjusted"
    ),
    converged = c(
        model_unadjusted$converged,
        model_splice$converged,
        model_full$converged
    ),
    iterations = c(
        model_unadjusted$iter,
        model_splice$iter,
        model_full$iter
    ),
    residual_deviance = c(
        deviance(model_unadjusted),
        deviance(model_splice),
        deviance(model_full)
    ),
    null_deviance = c(
        model_unadjusted$null.deviance,
        model_splice$null.deviance,
        model_full$null.deviance
    ),
    AIC = c(
        AIC(model_unadjusted),
        AIC(model_splice),
        AIC(model_full)
    )
)

fwrite(
    diagnostics,
    paste0(output_prefix, ".diagnostics.tsv"),
    sep = "\t"
)

capture.output(
    {
        cat("DISTRIBUTION SUMMARY\n")
        print(distribution_summary)

        cat("\n\nCORRELATION MATRIX\n")
        print(correlation_matrix)

        cat("\n\nVIF TABLE\n")
        print(vif_table)

        cat("\n\nUNADJUSTED MODEL\n")
        print(summary(model_unadjusted))

        cat("\n\nSPLICE-ADJUSTED MODEL\n")
        print(summary(model_splice))

        cat("\n\nFULLY ADJUSTED MODEL\n")
        print(summary(model_full))
    },
    file = paste0(output_prefix, ".full_summary.txt")
)

cat("\nDistribution summary:\n")
print(distribution_summary)

cat("\nCorrelation matrix:\n")
print(correlation_matrix)

cat("\nVIF table:\n")
print(vif_table)

cat("\nCoefficient results:\n")
print(coefficient_results)

cat("\nDiagnostics:\n")
print(diagnostics)

cat("\nCreated outputs with prefix:\n", output_prefix, "\n")
cat("Plots written to:\n", plot_dir, "\n")

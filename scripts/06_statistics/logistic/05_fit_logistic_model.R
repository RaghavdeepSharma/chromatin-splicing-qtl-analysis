#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
    stop(
        "Usage: Rscript 05_fit_logistic_model.R ",
        "<master_covariates.tsv.gz> ",
        "<overlap_keys.tsv> ",
        "<motif_archetype> ",
        "<celltype> ",
        "<output.tsv>"
    )
}

master_file <- args[1]
overlap_file <- args[2]
motif <- args[3]
celltype <- args[4]
output_file <- args[5]

dat <- fread(
    cmd = paste("gzip -dc", shQuote(master_file))
)

required <- c(
    "variant_gene_key",
    "high_pip",
    "log1p_splice_distance",
    "log1p_nearest_gene_tss_distance",
    "gc_fraction"
)

missing_columns <- setdiff(required, names(dat))

if (length(missing_columns) > 0) {
    stop(
        "Missing columns: ",
        paste(missing_columns, collapse = ", ")
    )
}

if (file.exists(overlap_file) && file.info(overlap_file)$size > 0) {

    # Each line is one complete variant-gene key.
    # Use readLines rather than fread because the keys contain "|",
    # which fread may incorrectly detect as a delimiter.
    overlap_set <- unique(
        trimws(
            readLines(
                overlap_file,
                warn = FALSE
            )
        )
    )

    overlap_set <- overlap_set[
        !is.na(overlap_set) &
        nzchar(overlap_set)
    ]

    dat[
        ,
        motif_overlap :=
            as.integer(
                variant_gene_key %in% overlap_set
            )
    ]

} else {

    dat[, motif_overlap := 0L]
}
dat[, high_pip := as.integer(high_pip)]
dat[, motif_overlap := as.integer(motif_overlap)]

model_columns <- c(
    "high_pip",
    "motif_overlap",
    "log1p_splice_distance",
    "log1p_nearest_gene_tss_distance",
    "gc_fraction"
)

dat <- dat[
    complete.cases(dat[, ..model_columns])
]

n_total <- nrow(dat)
n_high <- sum(dat$high_pip == 1)
n_low <- sum(dat$high_pip == 0)

n_high_overlap <- sum(
    dat$high_pip == 1 &
        dat$motif_overlap == 1
)

n_high_not_overlap <- sum(
    dat$high_pip == 1 &
        dat$motif_overlap == 0
)

n_low_overlap <- sum(
    dat$high_pip == 0 &
        dat$motif_overlap == 1
)

n_low_not_overlap <- sum(
    dat$high_pip == 0 &
        dat$motif_overlap == 0
)

n_overlap <- sum(dat$motif_overlap == 1)

base_result <- data.table(
    motif_archetype = motif,
    calderon_celltype = celltype,
    analysis_set = "strict_primary",
    n_total = n_total,
    n_high = n_high,
    n_low = n_low,
    n_overlap = n_overlap,
    n_high_overlap = n_high_overlap,
    n_high_not_overlap = n_high_not_overlap,
    n_low_overlap = n_low_overlap,
    n_low_not_overlap = n_low_not_overlap
)

empty_result <- function(status, warning_text = NA_character_) {
    cbind(
        base_result,
        data.table(
            model_status = status,
            converged = FALSE,
            iterations = NA_integer_,
            coefficient = NA_real_,
            standard_error = NA_real_,
            z_value = NA_real_,
            p_value = NA_real_,
            odds_ratio = NA_real_,
            ci_lower_95 = NA_real_,
            ci_upper_95 = NA_real_,
            vif_motif_overlap = NA_real_,
            vif_splice_distance = NA_real_,
            vif_tss_distance = NA_real_,
            vif_gc_fraction = NA_real_,
            residual_deviance = NA_real_,
            null_deviance = NA_real_,
            AIC = NA_real_,
            warning_text = warning_text
        )
    )
}

if (length(unique(dat$high_pip)) < 2) {
    fwrite(
        empty_result("skip_outcome_no_variation"),
        output_file,
        sep = "\t"
    )
    quit(save = "no", status = 0)
}

if (length(unique(dat$motif_overlap)) < 2) {
    fwrite(
        empty_result("skip_motif_no_variation"),
        output_file,
        sep = "\t"
    )
    quit(save = "no", status = 0)
}

predictors <- c(
    "motif_overlap",
    "log1p_splice_distance",
    "log1p_nearest_gene_tss_distance",
    "gc_fraction"
)

calculate_vif <- function(data, variable, predictors) {
    other_predictors <- setdiff(
        predictors,
        variable
    )

    fit <- lm(
        as.formula(
            paste(
                variable,
                "~",
                paste(other_predictors, collapse = " + ")
            )
        ),
        data = data
    )

    r_squared <- summary(fit)$r.squared

    1 / (1 - r_squared)
}

vif_values <- sapply(
    predictors,
    function(variable) {
        calculate_vif(
            dat,
            variable,
            predictors
        )
    }
)

warning_messages <- character()

model <- withCallingHandlers(
    glm(
        high_pip ~
            motif_overlap +
            log1p_splice_distance +
            log1p_nearest_gene_tss_distance +
            gc_fraction,
        data = dat,
        family = binomial(link = "logit"),
        control = glm.control(maxit = 50)
    ),
    warning = function(w) {
        warning_messages <<- c(
            warning_messages,
            conditionMessage(w)
        )

        invokeRestart("muffleWarning")
    }
)

coef_table <- summary(model)$coefficients

if (!"motif_overlap" %in% rownames(coef_table)) {
    fwrite(
        empty_result(
            "fit_missing_motif_coefficient",
            paste(unique(warning_messages), collapse = "; ")
        ),
        output_file,
        sep = "\t"
    )

    quit(save = "no", status = 0)
}

estimate <- coef_table[
    "motif_overlap",
    "Estimate"
]

standard_error <- coef_table[
    "motif_overlap",
    "Std. Error"
]

z_value <- coef_table[
    "motif_overlap",
    "z value"
]

p_value <- coef_table[
    "motif_overlap",
    "Pr(>|z|)"
]

separation_risk <- (
    n_high_overlap == 0 ||
    n_high_not_overlap == 0 ||
    n_low_overlap == 0 ||
    n_low_not_overlap == 0
)

finite_result <- all(
    is.finite(
        c(
            estimate,
            standard_error,
            z_value,
            p_value
        )
    )
)

model_status <- if (!model$converged) {
    "fit_not_converged"
} else if (!finite_result) {
    "fit_nonfinite_coefficient"
} else if (separation_risk) {
    "fit_complete_or_quasi_separation_risk"
} else {
    "fit_ok"
}

result <- cbind(
    base_result,
    data.table(
        model_status = model_status,
        converged = model$converged,
        iterations = model$iter,
        coefficient = estimate,
        standard_error = standard_error,
        z_value = z_value,
        p_value = p_value,
        odds_ratio = exp(estimate),
        ci_lower_95 = exp(
            estimate - 1.96 * standard_error
        ),
        ci_upper_95 = exp(
            estimate + 1.96 * standard_error
        ),
        vif_motif_overlap = unname(
            vif_values["motif_overlap"]
        ),
        vif_splice_distance = unname(
            vif_values["log1p_splice_distance"]
        ),
        vif_tss_distance = unname(
            vif_values[
                "log1p_nearest_gene_tss_distance"
            ]
        ),
        vif_gc_fraction = unname(
            vif_values["gc_fraction"]
        ),
        residual_deviance = deviance(model),
        null_deviance = model$null.deviance,
        AIC = AIC(model),
        warning_text = if (
            length(warning_messages) == 0
        ) {
            NA_character_
        } else {
            paste(
                unique(warning_messages),
                collapse = "; "
            )
        }
    )
)

fwrite(
    result,
    output_file,
    sep = "\t"
)

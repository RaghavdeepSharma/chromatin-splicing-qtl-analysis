#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(logistf)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
    stop(
        "Usage: Rscript 02_fit_firth_model.R ",
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

dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
)

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

# Recreate motif_overlap exactly as in the original logistic script.
if (
    file.exists(overlap_file) &&
    file.info(overlap_file)$size > 0
) {

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

empty_result <- function(
    status,
    warning_text = NA_character_
) {
    cbind(
        base_result,
        data.table(
            firth_status = status,
            converged = FALSE,
            iterations = NA_integer_,
            coefficient = NA_real_,
            standard_error = NA_real_,
            p_value = NA_real_,
            odds_ratio = NA_real_,
            ci_lower_95 = NA_real_,
            ci_upper_95 = NA_real_,
            log_likelihood_penalized = NA_real_,
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

warning_messages <- character()

fit <- tryCatch(
    withCallingHandlers(
        logistf(
            high_pip ~
                motif_overlap +
                log1p_splice_distance +
                log1p_nearest_gene_tss_distance +
                gc_fraction,
            data = dat,
            firth = TRUE,
            pl = TRUE,
            control = logistf.control(
                maxit = 1000,
                maxstep = 5
            ),
            plcontrol = logistpl.control(
                maxit = 1000
            )
        ),
        warning = function(w) {

            warning_messages <<- c(
                warning_messages,
                conditionMessage(w)
            )

            invokeRestart("muffleWarning")
        }
    ),
    error = function(e) {

        warning_messages <<- c(
            warning_messages,
            conditionMessage(e)
        )

        NULL
    }
)

if (is.null(fit)) {

    fwrite(
        empty_result(
            "firth_fit_error",
            paste(
                unique(warning_messages),
                collapse = "; "
            )
        ),
        output_file,
        sep = "\t"
    )

    quit(save = "no", status = 0)
}

if (!"motif_overlap" %in% names(fit$coefficients)) {

    fwrite(
        empty_result(
            "firth_missing_motif_coefficient",
            paste(
                unique(warning_messages),
                collapse = "; "
            )
        ),
        output_file,
        sep = "\t"
    )

    quit(save = "no", status = 0)
}

estimate <- unname(
    fit$coefficients["motif_overlap"]
)

standard_error <- tryCatch(
    {
        variance_matrix <- as.matrix(fit$var)

        if (
            "motif_overlap" %in% rownames(variance_matrix) &&
            "motif_overlap" %in% colnames(variance_matrix)
        ) {
            sqrt(
                variance_matrix[
                    "motif_overlap",
                    "motif_overlap"
                ]
            )
        } else {
            NA_real_
        }
    },
    error = function(e) NA_real_
)

p_value <- tryCatch(
    unname(
        fit$prob["motif_overlap"]
    ),
    error = function(e) NA_real_
)

ci_lower <- tryCatch(
    unname(
        fit$ci.lower["motif_overlap"]
    ),
    error = function(e) NA_real_
)

ci_upper <- tryCatch(
    unname(
        fit$ci.upper["motif_overlap"]
    ),
    error = function(e) NA_real_
)

iterations <- tryCatch(
    as.integer(max(fit$iter, na.rm = TRUE)),
    error = function(e) NA_integer_
)

penalized_loglik <- tryCatch(
    as.numeric(tail(fit$loglik, 1)),
    error = function(e) NA_real_
)

finite_result <- all(
    is.finite(
        c(
            estimate,
            exp(estimate)
        )
    )
)

firth_status <- if (!finite_result) {
    "firth_nonfinite_result"
} else if (
    is.na(ci_lower) ||
    is.na(ci_upper)
) {
    "firth_fit_ci_unavailable"
} else {
    "firth_fit_ok"
}

result <- cbind(
    base_result,
    data.table(
        firth_status = firth_status,
        converged = firth_status == "firth_fit_ok",
        iterations = iterations,
        coefficient = estimate,
        standard_error = standard_error,
        p_value = p_value,
        odds_ratio = exp(estimate),
        ci_lower_95 = exp(ci_lower),
        ci_upper_95 = exp(ci_upper),
        log_likelihood_penalized =
            penalized_loglik,
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

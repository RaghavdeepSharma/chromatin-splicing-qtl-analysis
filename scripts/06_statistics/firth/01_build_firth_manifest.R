#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

root <- Sys.getenv("KERIMOV_BASE")
if (root == "") stop("KERIMOV_BASE is not set.")

input_file <- file.path(
  root,
  "07_logistic_regression/04_full_analysis/combined_results",
  "all_strict_primary_logistic_results.tsv.gz"
)

out_dir <- file.path(
  root,
  "07_logistic_regression/05_firth_comparison/tables"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(input_file)

dt[, firth_eligible :=
     n_overlap > 0 &
     n_overlap < n_total &
     n_high > 0 &
     n_low > 0]

dt[, firth_reason := fifelse(
  n_overlap == 0,
  "no_motif_overlap",
  fifelse(
    n_overlap == n_total,
    "motif_overlap_constant_one",
    fifelse(
      n_high == 0 | n_low == 0,
      "outcome_no_variation",
      fifelse(
        model_status == "fit_complete_or_quasi_separation_risk",
        "separation_candidate",
        fifelse(
          model_status == "fit_ok",
          "ordinary_logistic_fit_ok",
          "other_eligible_status"
        )
      )
    )
  )
)]

manifest <- dt[
  firth_eligible == TRUE,
  .(
    motif_archetype,
    calderon_celltype,
    analysis_set,
    n_total,
    n_high,
    n_low,
    n_overlap,
    n_high_overlap,
    n_high_not_overlap,
    n_low_overlap,
    n_low_not_overlap,
    logistic_model_status = model_status,
    logistic_converged = converged,
    firth_reason
  )
]

setorder(manifest, calderon_celltype, motif_archetype)

fwrite(
  manifest,
  file.path(out_dir, "firth_model_manifest.tsv"),
  sep="\t"
)

summary <- dt[, .(
  n_models = .N,
  n_firth_eligible = sum(firth_eligible),
  n_no_overlap = sum(n_overlap == 0),
  n_logistic_fit_ok = sum(model_status == "fit_ok"),
  n_separation_candidates =
    sum(model_status == "fit_complete_or_quasi_separation_risk"),
  n_high_overlap_zero = sum(n_high_overlap == 0),
  n_high_overlap_1_to_4 =
    sum(n_high_overlap >= 1 & n_high_overlap < 5),
  n_high_overlap_at_least_5 =
    sum(n_high_overlap >= 5)
), by=calderon_celltype]

setorder(summary, calderon_celltype)

fwrite(
  summary,
  file.path(out_dir, "firth_eligibility_summary.tsv"),
  sep="\t"
)

cat("Firth-eligible models:", nrow(manifest), "\n")
cat(
  "Separation candidates:",
  sum(manifest$firth_reason == "separation_candidate"),
  "\n"
)
cat(
  "Ordinary fit-ok models:",
  sum(manifest$firth_reason == "ordinary_logistic_fit_ok"),
  "\n"
)

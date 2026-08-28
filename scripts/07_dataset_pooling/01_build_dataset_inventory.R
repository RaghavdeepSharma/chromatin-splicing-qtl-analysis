#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

base <- Sys.getenv("KERIMOV_BASE")
if (base == "") stop("KERIMOV_BASE is not set.")

available_file <- file.path(
    base,
    "02_manifests/susie_credible_sets_manifest.tsv"
)

pooling_file <- file.path(
    base,
    "02_manifests/motif_sqtl_deduplicated_celltype_manifest.tsv"
)

out_dir <- file.path(
    base,
    "08_dataset_pooling_audit/tables"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(available_file)) {
    stop("Missing available dataset manifest: ", available_file)
}

if (!file.exists(pooling_file)) {
    stop("Missing pooling manifest: ", pooling_file)
}

available <- fread(available_file)
pooling <- fread(pooling_file)

required_available <- c(
    "study_id",
    "dataset_id",
    "study_label",
    "sample_group",
    "tissue_label",
    "condition_label",
    "state_interpretation",
    "sample_size",
    "preliminary_cell_class",
    "calderon_priority",
    "download_tier"
)

missing_available <- setdiff(
    required_available,
    names(available)
)

if (length(missing_available) > 0) {
    stop(
        "Missing columns in available manifest: ",
        paste(missing_available, collapse = ", ")
    )
}

strict_pooling <- pooling[
    analysis_set == "strict_primary"
]

used_long <- strict_pooling[
    ,
    .(
        dataset_id = trimws(
            unlist(strsplit(dataset_ids, ",", fixed = TRUE))
        )
    ),
    by = .(
        analysis_set,
        pooled_calderon_celltype = calderon_celltype,
        accessibility_column,
        n_datasets,
        raw_dataset_variant_gene_rows,
        unique_variant_gene_pairs,
        duplicate_rows_collapsed,
        n_high_pip_pairs,
        n_low_pip_pairs,
        deduplicated_sqtl_bed
    )
]

inventory <- merge(
    available,
    used_long,
    by = "dataset_id",
    all.x = TRUE
)

inventory[, used_in_strict_primary :=
    !is.na(pooled_calderon_celltype)
]

inventory[, inclusion_status := fifelse(
    used_in_strict_primary,
    "used_strict_primary",
    fifelse(
        download_tier == "primary",
        "available_primary_not_used",
        "available_secondary_not_used"
    )
)]

setorder(
    inventory,
    -used_in_strict_primary,
    pooled_calderon_celltype,
    study_label,
    dataset_id
)

fwrite(
    inventory,
    file.path(
        out_dir,
        "available_vs_used_dataset_inventory.tsv"
    ),
    sep = "\t",
    na = "NA"
)

fwrite(
    inventory[used_in_strict_primary == TRUE],
    file.path(
        out_dir,
        "strict_primary_used_dataset_inventory.tsv"
    ),
    sep = "\t",
    na = "NA"
)

fwrite(
    strict_pooling,
    file.path(
        out_dir,
        "strict_primary_pooling_before_after.tsv"
    ),
    sep = "\t",
    na = "NA"
)

summary <- inventory[, .(
    n_available = .N,
    n_used_strict_primary =
        sum(used_in_strict_primary),
    n_not_used =
        sum(!used_in_strict_primary)
)]

fwrite(
    summary,
    file.path(
        out_dir,
        "dataset_inventory_summary.tsv"
    ),
    sep = "\t"
)

cat("Available datasets:", nrow(inventory), "\n")
cat(
    "Used strict-primary datasets:",
    sum(inventory$used_in_strict_primary),
    "\n"
)
cat(
    "Strict-primary pooled cell types:",
    nrow(strict_pooling),
    "\n"
)

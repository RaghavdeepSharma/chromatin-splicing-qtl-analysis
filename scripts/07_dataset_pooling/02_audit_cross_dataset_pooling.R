#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

base <- Sys.getenv("KERIMOV_BASE")
if (base == "") stop("KERIMOV_BASE is not set.")

variant_file <- file.path(
    base,
    "05_processed",
    "leafcutter_dataset_variants_protein_coding_lncRNA.tsv.gz"
)

mapping_file <- file.path(
    base,
    "05_processed",
    "gene_filtered_sqtl_bed_manifest.tsv"
)

out_dir <- file.path(
    base,
    "08_dataset_pooling_audit",
    "tables"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

variants <- fread(
    cmd = paste("gzip -cd", shQuote(variant_file))
)

mapping <- fread(mapping_file)

required <- c(
    "dataset_id",
    "variant",
    "gene_id_unversioned",
    "max_pip",
    "calderon_celltype"
)

missing <- setdiff(required, names(variants))

if (length(missing) > 0) {
    stop(
        "Missing required columns: ",
        paste(missing, collapse = ", ")
    )
}

variants[, max_pip := as.numeric(max_pip)]

strict_ids <- mapping[
    include_strict_primary == "YES",
    dataset_id
]

strict <- variants[
    dataset_id %in% strict_ids
]

strict[, variant_gene_key := paste(
    variant,
    gene_id_unversioned,
    sep = "|"
)]

strict[, source_pip_class := fifelse(
    max_pip > 0.90,
    "high_PIP",
    "low_PIP"
)]

# Optional metadata fields.
optional_fields <- c(
    "study_label",
    "sample_group",
    "condition_label",
    "state_interpretation"
)

available_optional <- optional_fields[
    optional_fields %in% names(strict)
]

# Determine which source row would be retained by the current logic.
setorder(
    strict,
    calderon_celltype,
    variant_gene_key,
    -max_pip,
    dataset_id
)

strict[, retained_by_current_rule :=
    seq_len(.N) == 1L,
    by = .(
        calderon_celltype,
        variant_gene_key
    )
]

source_level <- strict[, c(
    list(
        calderon_celltype =
            calderon_celltype,
        dataset_id =
            dataset_id,
        variant_gene_key =
            variant_gene_key,
        variant =
            variant,
        gene_id_unversioned =
            gene_id_unversioned,
        max_pip =
            max_pip,
        source_pip_class =
            source_pip_class,
        retained_by_current_rule =
            retained_by_current_rule
    ),
    mget(available_optional)
)]

fwrite(
    source_level,
    file.path(
        out_dir,
        "strict_primary_source_level_variant_gene_rows.tsv.gz"
    ),
    sep = "\t",
    na = "NA"
)

key_audit <- strict[, {

    order_index <- order(
        -max_pip,
        dataset_id
    )

    retained_index <- order_index[1]

    output <- list(
        n_source_rows = .N,
        n_source_datasets = uniqueN(dataset_id),
        source_dataset_ids = paste(
            sort(unique(dataset_id)),
            collapse = ","
        ),
        min_pip = min(max_pip),
        max_pip = max(max_pip),
        mean_pip = mean(max_pip),
        pip_range = max(max_pip) - min(max_pip),
        n_high_sources = sum(max_pip > 0.90),
        n_low_sources = sum(max_pip <= 0.90),
        any_high = any(max_pip > 0.90),
        all_high = all(max_pip > 0.90),
        discordant_high_low =
            any(max_pip > 0.90) &
            any(max_pip <= 0.90),
        retained_dataset_id =
            dataset_id[retained_index],
        retained_max_pip =
            max_pip[retained_index],
        retained_pip_class =
            source_pip_class[retained_index]
    )

    for (field in available_optional) {
        output[[paste0("retained_", field)]] <-
            get(field)[retained_index]
    }

    output
},
by = .(
    calderon_celltype,
    variant_gene_key,
    variant,
    gene_id_unversioned
)]

key_audit[, high_supported_by_one_source_only :=
    retained_pip_class == "high_PIP" &
    n_high_sources == 1
]

key_audit[, high_replicated_in_multiple_sources :=
    retained_pip_class == "high_PIP" &
    n_high_sources >= 2
]

fwrite(
    key_audit,
    file.path(
        out_dir,
        "strict_primary_variant_gene_pooling_audit.tsv.gz"
    ),
    sep = "\t",
    na = "NA"
)

summary <- key_audit[, .(
    n_unique_variant_gene_pairs = .N,
    n_pairs_seen_in_multiple_sources =
        sum(n_source_datasets > 1),
    n_discordant_high_low_pairs =
        sum(discordant_high_low),
    n_retained_high_pairs =
        sum(retained_pip_class == "high_PIP"),
    n_retained_high_supported_by_one_source =
        sum(high_supported_by_one_source_only),
    n_retained_high_replicated_across_sources =
        sum(high_replicated_in_multiple_sources),
    n_retained_low_pairs =
        sum(retained_pip_class == "low_PIP")
),
by = calderon_celltype]

setorder(summary, calderon_celltype)

fwrite(
    summary,
    file.path(
        out_dir,
        "strict_primary_pooling_conflict_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
)

dataset_summary <- strict[, .(
    n_rows_before_cross_dataset_pooling = .N,
    n_unique_variant_gene_pairs =
        uniqueN(variant_gene_key),
    n_high_rows =
        sum(max_pip > 0.90),
    n_low_rows =
        sum(max_pip <= 0.90),
    min_pip =
        min(max_pip),
    median_pip =
        median(max_pip),
    max_pip =
        max(max_pip)
),
by = .(
    calderon_celltype,
    dataset_id
)]

setorder(
    dataset_summary,
    calderon_celltype,
    dataset_id
)

fwrite(
    dataset_summary,
    file.path(
        out_dir,
        "strict_primary_dataset_level_before_pooling.tsv"
    ),
    sep = "\t",
    na = "NA"
)

cat("Strict-primary source rows:", nrow(strict), "\n")
cat("Unique pooled variant-gene keys:", nrow(key_audit), "\n")
cat(
    "Pairs represented in multiple sources:",
    sum(key_audit$n_source_datasets > 1),
    "\n"
)
cat(
    "High/Low-discordant repeated pairs:",
    sum(key_audit$discordant_high_low),
    "\n"
)

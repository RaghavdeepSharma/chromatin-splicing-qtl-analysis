#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop(
        "Usage: Rscript 02_add_tss_distance.R ",
        "<variant_table.tsv.gz> <gencode_features.rds> <output.tsv.gz>"
    )
}

variant_file <- args[1]
feature_file <- args[2]
output_file <- args[3]

cat("Reading variant table:\n", variant_file, "\n")

variant <- fread(
    cmd = paste("gzip -dc", shQuote(variant_file))
)

required_variant <- c(
    "gene_id",
    "chromosome",
    "position"
)

missing_variant <- setdiff(
    required_variant,
    names(variant)
)

if (length(missing_variant) > 0) {
    stop(
        "Variant table is missing: ",
        paste(missing_variant, collapse = ", ")
    )
}

cat("Reading GENCODE features:\n", feature_file, "\n")

features <- readRDS(feature_file)

if (!"tss" %in% names(features)) {
    stop("GENCODE feature object does not contain a tss table.")
}

tss <- as.data.table(features$tss)

required_tss <- c(
    "seqnames",
    "strand",
    "transcript_id",
    "gene_id",
    "tss"
)

missing_tss <- setdiff(
    required_tss,
    names(tss)
)

if (length(missing_tss) > 0) {
    stop(
        "TSS table is missing: ",
        paste(missing_tss, collapse = ", ")
    )
}

# Remove GENCODE version suffixes for matching.
variant[, gene_id_clean := sub("\\..*$", "", gene_id)]
tss[, gene_id_clean := sub("\\..*$", "", gene_id)]

variant[, position := as.integer(position)]
tss[, tss := as.integer(tss)]

# Keep only usable transcript TSS records.
tss <- tss[
    !is.na(gene_id_clean) &
    !is.na(seqnames) &
    !is.na(tss) &
    strand %in% c("+", "-")
]

# Collapse duplicate transcript TSS coordinates.
tss_unique <- unique(
    tss[, .(
        gene_id_clean,
        seqnames,
        strand,
        tss
    )]
)

cat("Unique transcript TSS records:", nrow(tss_unique), "\n")

# Create one row per distinct variant location and gene.
variant_lookup <- unique(
    variant[, .(
        gene_id_clean,
        chromosome,
        position
    )]
)

cat("Unique variant-gene positions:", nrow(variant_lookup), "\n")

# Join each variant to all transcript TSSs belonging to its associated gene.
joined <- merge(
    variant_lookup,
    tss_unique,
    by.x = c("gene_id_clean", "chromosome"),
    by.y = c("gene_id_clean", "seqnames"),
    all.x = TRUE,
    allow.cartesian = TRUE,
    sort = FALSE
)

joined[, tss_distance_candidate := abs(position - tss)]

# For each variant-gene pair, choose the nearest transcript TSS.
nearest <- joined[
    !is.na(tss_distance_candidate),
    .SD[which.min(tss_distance_candidate)],
    by = .(
        gene_id_clean,
        chromosome,
        position
    )
]

nearest <- nearest[, .(
    gene_id_clean,
    chromosome,
    position,
    nearest_gene_tss = tss,
    nearest_gene_tss_strand = strand,
    nearest_gene_tss_distance = tss_distance_candidate
)]

cat("Annotated unique variant-gene positions:", nrow(nearest), "\n")

# Preserve original row order.
variant[, original_row_order := .I]

variant <- merge(
    variant,
    nearest,
    by = c(
        "gene_id_clean",
        "chromosome",
        "position"
    ),
    all.x = TRUE,
    sort = FALSE
)

setorder(
    variant,
    original_row_order
)

variant[, original_row_order := NULL]

variant[
    ,
    log1p_nearest_gene_tss_distance :=
        log1p(nearest_gene_tss_distance)
]

cat("Rows after merge:", nrow(variant), "\n")
cat(
    "Rows annotated:",
    sum(!is.na(variant$nearest_gene_tss_distance)),
    "\n"
)
cat(
    "Rows missing:",
    sum(is.na(variant$nearest_gene_tss_distance)),
    "\n"
)

if (any(!is.na(variant$nearest_gene_tss_distance))) {
    cat(
        "Minimum distance:",
        min(
            variant$nearest_gene_tss_distance,
            na.rm = TRUE
        ),
        "\n"
    )

    cat(
        "Median distance:",
        median(
            variant$nearest_gene_tss_distance,
            na.rm = TRUE
        ),
        "\n"
    )

    cat(
        "Maximum distance:",
        max(
            variant$nearest_gene_tss_distance,
            na.rm = TRUE
        ),
        "\n"
    )
}

cat("Writing output:\n", output_file, "\n")

fwrite(
    variant,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    compress = "gzip"
)

cat("Done.\n")

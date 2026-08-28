#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop(
        "Usage: Rscript 01_add_splice_distance.R ",
        "<variant_level.tsv.gz> ",
        "<association_table.tsv.gz> ",
        "<output.tsv.gz>"
    )
}

variant_file <- args[1]
association_file <- args[2]
output_file <- args[3]

cat("Reading variant-level table:\n", variant_file, "\n")

variant_dat <- read.delim(
    gzfile(variant_file),
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

cat("Reading source association table:\n", association_file, "\n")

assoc <- read.delim(
    gzfile(association_file),
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

merge_columns <- c(
    "dataset_id",
    "variant",
    "gene_id"
)

required_assoc <- c(
    merge_columns,
    "position",
    "junction_start",
    "junction_end"
)

missing_variant <- setdiff(
    merge_columns,
    colnames(variant_dat)
)

missing_assoc <- setdiff(
    required_assoc,
    colnames(assoc)
)

if (length(missing_variant) > 0) {
    stop(
        "Variant-level table missing columns: ",
        paste(missing_variant, collapse = ", ")
    )
}

if (length(missing_assoc) > 0) {
    stop(
        "Association table missing columns: ",
        paste(missing_assoc, collapse = ", ")
    )
}

assoc$position <- suppressWarnings(
    as.numeric(assoc$position)
)

assoc$junction_start <- suppressWarnings(
    as.numeric(assoc$junction_start)
)

assoc$junction_end <- suppressWarnings(
    as.numeric(assoc$junction_end)
)

assoc$distance_to_junction_start <- abs(
    assoc$position - assoc$junction_start
)

assoc$distance_to_junction_end <- abs(
    assoc$position - assoc$junction_end
)

assoc$splice_distance <- pmin(
    assoc$distance_to_junction_start,
    assoc$distance_to_junction_end
)

annotation <- assoc[
    !is.na(assoc$splice_distance),
    c(
        "dataset_id",
        "variant",
        "gene_id",
        "splice_distance"
    )
]

# A dataset-variant-gene combination may appear for multiple
# LeafCutter junctions. Keep the nearest splice boundary.
annotation <- aggregate(
    splice_distance ~ dataset_id + variant + gene_id,
    data = annotation,
    FUN = min
)

cat(
    "Unique annotation keys:",
    nrow(annotation),
    "\n"
)

# Preserve original regression-table order
variant_dat$original_row_order <- seq_len(
    nrow(variant_dat)
)

merged <- merge(
    variant_dat,
    annotation,
    by = merge_columns,
    all.x = TRUE,
    sort = FALSE
)

merged <- merged[
    order(merged$original_row_order),
]

merged$original_row_order <- NULL

merged$log1p_splice_distance <- log1p(
    merged$splice_distance
)

cat("Rows before merge:", nrow(variant_dat), "\n")
cat("Rows after merge:", nrow(merged), "\n")

cat(
    "Rows with splice distance:",
    sum(!is.na(merged$splice_distance)),
    "\n"
)

cat(
    "Rows missing splice distance:",
    sum(is.na(merged$splice_distance)),
    "\n"
)

cat(
    "Minimum splice distance:",
    min(merged$splice_distance, na.rm = TRUE),
    "\n"
)

cat(
    "Median splice distance:",
    median(merged$splice_distance, na.rm = TRUE),
    "\n"
)

cat(
    "Maximum splice distance:",
    max(merged$splice_distance, na.rm = TRUE),
    "\n"
)

con <- gzfile(output_file, "wt")

write.table(
    merged,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

close(con)

cat("Created:\n", output_file, "\n")

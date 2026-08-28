#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage:"
    echo "  bash 08_add_gc_content_100bp.sh \\"
    echo "    <variant_table.tsv.gz> <hg38.fa> <output.tsv.gz>"
    exit 1
fi

INPUT="$1"
FASTA="$2"
OUTPUT="$3"

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: input file not found: $INPUT" >&2
    exit 1
fi

if [[ ! -f "$FASTA" ]]; then
    echo "ERROR: FASTA not found: $FASTA" >&2
    exit 1
fi

if [[ ! -f "${FASTA}.fai" ]]; then
    echo "ERROR: FASTA index not found: ${FASTA}.fai" >&2
    echo "Run: samtools faidx $FASTA" >&2
    exit 1
fi

command -v bedtools >/dev/null 2>&1 || {
    echo "ERROR: bedtools is not available." >&2
    exit 1
}

TMPDIR_PATH=$(mktemp -d)

cleanup() {
    rm -rf "$TMPDIR_PATH"
}
trap cleanup EXIT

VARIANTS="$TMPDIR_PATH/unique_variants.tsv"
WINDOWS_RAW="$TMPDIR_PATH/windows.raw.bed"
WINDOWS="$TMPDIR_PATH/windows.clipped.bed"
SEQUENCES="$TMPDIR_PATH/windows.sequences.tsv"
GC_TABLE="$TMPDIR_PATH/gc_content.tsv"
CHROMSIZES="$TMPDIR_PATH/hg38.chrom.sizes"

echo "Input:  $INPUT"
echo "FASTA:  $FASTA"
echo "Output: $OUTPUT"

cut -f1,2 "${FASTA}.fai" > "$CHROMSIZES"

############################################################
# Extract unique chromosome-position combinations.
#
# Input position is 1-based.
# For a ±100 bp window centered on the variant:
#
# BED start = position - 101
# BED end   = position + 100
#
# This produces 201 bp:
# 100 upstream + variant base + 100 downstream.
############################################################

gzip -dc "$INPUT" |
awk -F'\t' 'BEGIN {OFS="\t"}
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "chromosome") chr_col = i
        if ($i == "position") pos_col = i
    }

    if (!chr_col || !pos_col) {
        print "ERROR: chromosome or position column not found." > "/dev/stderr"
        exit 1
    }

    next
}
{
    key = $chr_col OFS $pos_col

    if (!seen[key]++) {
        print $chr_col, $pos_col
    }
}' > "$VARIANTS"

echo "Unique chromosome-position combinations: $(wc -l < "$VARIANTS")"

awk -F'\t' 'BEGIN {OFS="\t"}
{
    chr = $1
    pos = $2

    start = pos - 101
    end = pos + 100

    if (start < 0) {
        start = 0
    }

    name = chr "|" pos

    print chr, start, end, name
}' "$VARIANTS" > "$WINDOWS_RAW"

# Clip windows to chromosome boundaries.
bedtools slop \
    -i "$WINDOWS_RAW" \
    -g "$CHROMSIZES" \
    -b 0 \
    > "$WINDOWS"

############################################################
# Extract sequence.
############################################################

bedtools getfasta \
    -fi "$FASTA" \
    -bed "$WINDOWS" \
    -name \
    -tab \
    > "$SEQUENCES"

############################################################
# Calculate:
#   gc_count
#   valid_base_count
#   window_length
#   n_count
#   gc_fraction
#
# GC denominator uses A/C/G/T bases only.
############################################################

awk -F'\t' 'BEGIN {
    OFS = "\t"
    print "chromosome", "position", "gc_count", \
          "valid_base_count", "window_length", \
          "n_count", "gc_fraction"
}
{
    header = $1
    sequence = toupper($2)

    sub(/::.*/, "", header)

    split(header, id, "|")
    chr = id[1]
    pos = id[2]

    gc = 0
    valid = 0
    n_count = 0

    for (i = 1; i <= length(sequence); i++) {
        base = substr(sequence, i, 1)

        if (base == "G" || base == "C") {
            gc++
            valid++
        } else if (base == "A" || base == "T") {
            valid++
        } else {
            n_count++
        }
    }

    window_length = length(sequence)

    if (valid > 0) {
        gc_fraction = gc / valid
    } else {
        gc_fraction = "NA"
    }

    print chr, pos, gc, valid, window_length, \
          n_count, gc_fraction
}' "$SEQUENCES" > "$GC_TABLE"

############################################################
# Merge GC annotation back onto the regression table using R.
############################################################

Rscript - "$INPUT" "$GC_TABLE" "$OUTPUT" <<'RSCRIPT'
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
gc_file <- args[2]
output_file <- args[3]

cat("Reading regression table...\n")

dat <- fread(
    cmd = paste("gzip -dc", shQuote(input_file))
)

cat("Reading GC table...\n")

gc <- fread(gc_file)

dat[, chromosome := as.character(chromosome)]
dat[, position := as.integer(position)]

gc[, chromosome := as.character(chromosome)]
gc[, position := as.integer(position)]

dat[, original_row_order := .I]

merged <- merge(
    dat,
    gc,
    by = c("chromosome", "position"),
    all.x = TRUE,
    sort = FALSE
)

setorder(merged, original_row_order)
merged[, original_row_order := NULL]

cat("Rows before merge:", nrow(dat), "\n")
cat("Rows after merge:", nrow(merged), "\n")
cat(
    "Rows with GC content:",
    sum(!is.na(merged$gc_fraction)),
    "\n"
)
cat(
    "Rows missing GC content:",
    sum(is.na(merged$gc_fraction)),
    "\n"
)

cat("\nGC fraction summary:\n")
print(summary(merged$gc_fraction))

cat("\nWindow-length summary:\n")
print(summary(merged$window_length))

cat("\nRows containing at least one N base:",
    sum(merged$n_count > 0, na.rm = TRUE),
    "\n"
)

fwrite(
    merged,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    compress = "gzip"
)

cat("\nCreated:", output_file, "\n")
RSCRIPT

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript 02_annotate_motif_features.R <motif.bed.gz> <gencode_features.rds> <outdir>")
}

motif_file <- args[1]
feature_rds <- args[2]
outdir <- args[3]

dir.create(file.path(outdir, "annotated"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

archetype <- basename(motif_file)
archetype <- sub("\\.bed\\.gz$", "", archetype)
archetype <- sub("\\.bed$", "", archetype)

message("annotating motif archetype: ", archetype)
message("reading motif file: ", motif_file)
message("reading prepared GENCODE features: ", feature_rds)

features <- readRDS(feature_rds)

genes <- features$genes
five_utrs <- features$five_utrs
three_utrs <- features$three_utrs
exons <- features$exons
introns <- features$introns
tss_dt <- as.data.table(features$tss)

flip_strand <- function(x) {
  fifelse(x == "+", "-", fifelse(x == "-", "+", x))
}

collapse_unique <- function(x) {
  x <- unique(na.omit(as.character(x)))
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = ";")
}

pick_priority_region <- function(x) {
  x <- unique(unlist(strsplit(x, ";", fixed = TRUE)))
  x <- x[!is.na(x) & x != ""]
  priority <- c("5UTR", "3UTR", "exon", "intron")
  hit <- priority[priority %in% x]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# Read the motif BED file directly from gzip.
motif <- fread(cmd = paste("gzip -dc", shQuote(motif_file)), header = FALSE)

if (ncol(motif) < 6) {
  stop("Expected at least 6 BED columns: chr, start, end, name, score, strand.")
}

# Keep the first 6 BED columns explicit and preserve any extra motif columns.
base_names <- c("chr", "start", "end", "motif_name", "score", "original_strand")
extra_names <- paste0("extra_", seq_len(max(0, ncol(motif) - 6)))
setnames(motif, c(base_names, extra_names))

motif[, input_row_id := .I]
motif[, archetype := archetype]
motif[, tfbs_id := paste(archetype, input_row_id, chr, start, end, motif_name, sep = "|")]

n_input <- nrow(motif)

if (!all(motif$original_strand %in% c("+", "-"))) {
  bad <- unique(motif[!original_strand %in% c("+", "-"), original_strand])
  stop("Found motif rows with non + or - strand values: ", paste(bad, collapse = ", "))
}

# Duplicate each TFBS row to annotate both motif orientations.
# The default row keeps the Vierstra strand.
# The reverse-complement row flips the strand.
default_rows <- copy(motif)
default_rows[, motif_orientation := "default"]
default_rows[, annot_strand := original_strand]

reverse_rows <- copy(motif)
reverse_rows[, motif_orientation := "reverse_complement"]
reverse_rows[, annot_strand := flip_strand(original_strand)]

motif2 <- rbindlist(list(default_rows, reverse_rows), use.names = TRUE, fill = TRUE)
motif2[, output_row_id := .I]

if (nrow(motif2) != 2 * n_input) {
  stop("Expected exactly 2x rows after orientation expansion, but the row count did not match.")
}

# Convert BED coordinates to 1-based GRanges.
# BED start is 0-based and BED end is half-open, so start becomes start + 1.
motif_gr <- GRanges(
  seqnames = motif2$chr,
  ranges = IRanges(start = motif2$start + 1, end = motif2$end),
  strand = motif2$annot_strand
)

mcols(motif_gr)$output_row_id <- motif2$output_row_id

# Initialize annotation columns.
motif2[, region_class := "intergenic"]
motif2[, overlaps_gene_same_strand := FALSE]
motif2[, gene_id := NA_character_]
motif2[, gene_name := NA_character_]
motif2[, gene_type := NA_character_]
motif2[, overlaps_protein_coding_or_lncRNA := FALSE]
motif2[, genic_region_all := NA_character_]
motif2[, genic_region_priority := NA_character_]
motif2[, upstream_tss_distance := NA_integer_]
motif2[, upstream_tss_transcript_id := NA_character_]
motif2[, upstream_tss_gene_id := NA_character_]
motif2[, upstream_tss_gene_name := NA_character_]

# Perform same-strand gene overlaps only.
gene_hits <- findOverlaps(motif_gr, genes, ignore.strand = FALSE)

motif_gene_pairs <- data.table()

if (length(gene_hits) > 0) {
  gene_hit_dt <- data.table(
    output_row_id = queryHits(gene_hits),
    gene_index = subjectHits(gene_hits)
  )

  gene_meta <- as.data.table(as.data.frame(mcols(genes)))
  gene_meta[, gene_index := .I]

  gene_hit_dt <- merge(gene_hit_dt, gene_meta, by = "gene_index", all.x = TRUE)

  motif_gene_pairs <- gene_hit_dt[, .(
    output_row_id,
    gene_id = as.character(gene_id),
    gene_name = as.character(gene_name),
    gene_type = as.character(gene_type)
  )]

  gene_summary <- motif_gene_pairs[, .(
    gene_id = collapse_unique(gene_id),
    gene_name = collapse_unique(gene_name),
    gene_type = collapse_unique(gene_type),
    overlaps_protein_coding_or_lncRNA = any(gene_type %in% c("protein_coding", "lncRNA"))
  ), by = output_row_id]

  motif2[gene_summary$output_row_id, region_class := "intragenic"]
  motif2[gene_summary$output_row_id, overlaps_gene_same_strand := TRUE]
  motif2[gene_summary$output_row_id, gene_id := gene_summary$gene_id]
  motif2[gene_summary$output_row_id, gene_name := gene_summary$gene_name]
  motif2[gene_summary$output_row_id, gene_type := gene_summary$gene_type]
  motif2[gene_summary$output_row_id, overlaps_protein_coding_or_lncRNA := gene_summary$overlaps_protein_coding_or_lncRNA]
}

# Annotate same-strand genic subregions.
add_region_hits <- function(subject_gr, label) {
  if (length(subject_gr) == 0) {
    return(data.table())
  }

  hits <- findOverlaps(motif_gr, subject_gr, ignore.strand = FALSE)

  if (length(hits) == 0) {
    return(data.table())
  }

  data.table(
    output_row_id = queryHits(hits),
    genic_region = label
  )
}

region_hits <- rbindlist(
  list(
    add_region_hits(five_utrs, "5UTR"),
    add_region_hits(three_utrs, "3UTR"),
    add_region_hits(exons, "exon"),
    add_region_hits(introns, "intron")
  ),
  use.names = TRUE,
  fill = TRUE
)

eligible_ids <- motif2[overlaps_protein_coding_or_lncRNA == TRUE, output_row_id]

if (nrow(region_hits) > 0 && length(eligible_ids) > 0) {
  region_hits <- unique(region_hits[output_row_id %in% eligible_ids])

  region_summary <- region_hits[, .(
    genic_region_all = paste(sort(unique(genic_region)), collapse = ";")
  ), by = output_row_id]

  region_summary[, genic_region_priority := vapply(genic_region_all, pick_priority_region, character(1))]

  motif2[region_summary$output_row_id, genic_region_all := region_summary$genic_region_all]
  motif2[region_summary$output_row_id, genic_region_priority := region_summary$genic_region_priority]
}

motif2[region_class == "intergenic", genic_region_all := "intergenic"]
motif2[region_class == "intergenic", genic_region_priority := "intergenic"]

motif2[region_class == "intragenic" & overlaps_protein_coding_or_lncRNA == FALSE, genic_region_all := "other_gene_type"]
motif2[region_class == "intragenic" & overlaps_protein_coding_or_lncRNA == FALSE, genic_region_priority := "other_gene_type"]

motif2[
  region_class == "intragenic" &
    overlaps_protein_coding_or_lncRNA == TRUE &
    is.na(genic_region_priority),
  genic_region_all := "genic_other_or_unresolved"
]

motif2[
  region_class == "intragenic" &
    overlaps_protein_coding_or_lncRNA == TRUE &
    is.na(genic_region_priority),
  genic_region_priority := "genic_other_or_unresolved"
]

# Calculate the nearest upstream TSS in a strand-specific way.
# For plus-strand motifs, upstream means lower coordinate.
# For minus-strand motifs, upstream means higher coordinate.
if (nrow(motif_gene_pairs) > 0) {
  pc_lnc_pairs <- motif_gene_pairs[gene_type %in% c("protein_coding", "lncRNA")]
  pc_lnc_pairs <- unique(pc_lnc_pairs[, .(output_row_id, gene_id)])

  motif_pos <- motif2[, .(
    output_row_id,
    chr,
    motif_start_1based = start + 1,
    motif_end_1based = end,
    annot_strand
  )]

  tss_pc_lnc <- tss_dt[gene_type %in% c("protein_coding", "lncRNA")]

  tss_candidates <- merge(pc_lnc_pairs, motif_pos, by = "output_row_id", allow.cartesian = TRUE)
  tss_candidates <- merge(tss_candidates, tss_pc_lnc, by = "gene_id", allow.cartesian = TRUE)

  tss_candidates <- tss_candidates[
    chr == seqnames &
      annot_strand == strand
  ]

  if (nrow(tss_candidates) > 0) {
    tss_candidates[, upstream_tss_distance := fifelse(
      annot_strand == "+",
      motif_start_1based - tss,
      tss - motif_end_1based
    )]

    tss_candidates <- tss_candidates[upstream_tss_distance >= 0]

    if (nrow(tss_candidates) > 0) {
      setorder(tss_candidates, output_row_id, upstream_tss_distance)
      best_tss <- tss_candidates[, .SD[1], by = output_row_id]

      motif2[best_tss$output_row_id, upstream_tss_distance := as.integer(best_tss$upstream_tss_distance)]
      motif2[best_tss$output_row_id, upstream_tss_transcript_id := best_tss$transcript_id]
      motif2[best_tss$output_row_id, upstream_tss_gene_id := best_tss$gene_id]
      motif2[best_tss$output_row_id, upstream_tss_gene_name := best_tss$gene_name]
    }
  }
}

# Keep TSS distance only where the task asks for it.
motif2[overlaps_protein_coding_or_lncRNA == FALSE, upstream_tss_distance := NA_integer_]
motif2[overlaps_protein_coding_or_lncRNA == FALSE, upstream_tss_transcript_id := NA_character_]
motif2[overlaps_protein_coding_or_lncRNA == FALSE, upstream_tss_gene_id := NA_character_]
motif2[overlaps_protein_coding_or_lncRNA == FALSE, upstream_tss_gene_name := NA_character_]

# Run explicit sanity checks before writing output.
if (motif2[motif_orientation == "default", .N] != n_input) {
  stop("Default orientation row count does not match input row count.")
}

if (motif2[motif_orientation == "reverse_complement", .N] != n_input) {
  stop("Reverse-complement row count does not match input row count.")
}

strand_check <- motif2[, .N, by = .(motif_orientation, original_strand, annot_strand)]

summary_dt <- data.table(
  archetype = archetype,
  input_rows = n_input,
  output_rows = nrow(motif2),
  expected_output_rows = 2 * n_input,
  row_doubling_pass = nrow(motif2) == 2 * n_input,
  default_rows = motif2[motif_orientation == "default", .N],
  reverse_complement_rows = motif2[motif_orientation == "reverse_complement", .N],
  intragenic_rows = motif2[region_class == "intragenic", .N],
  intergenic_rows = motif2[region_class == "intergenic", .N],
  protein_coding_or_lncRNA_rows = motif2[overlaps_protein_coding_or_lncRNA == TRUE, .N],
  rows_with_genic_region = motif2[!is.na(genic_region_priority), .N],
  rows_with_upstream_tss_distance = motif2[!is.na(upstream_tss_distance), .N]
)

out_tsv <- file.path(outdir, "annotated", paste0(archetype, ".general_annotated.tsv"))
out_gz <- paste0(out_tsv, ".gz")
summary_file <- file.path(outdir, "tables", paste0(archetype, ".summary.tsv"))
strand_file <- file.path(outdir, "tables", paste0(archetype, ".strand_orientation_check.tsv"))

fwrite(motif2, out_tsv, sep = "\t")
system2("gzip", c("-f", out_tsv))

fwrite(summary_dt, summary_file, sep = "\t")
fwrite(strand_check, strand_file, sep = "\t")

message("Annotated motif file written: ", out_gz)
message("Summary file written: ", summary_file)
message("Finished motif archetype: ", archetype)

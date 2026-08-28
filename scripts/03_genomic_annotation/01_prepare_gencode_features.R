#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(GenomicFeatures)
  library(txdbmaker)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript 00_prepare_gencode_v47_features.R <gencode_v47_comprehensive.gtf.gz> <outdir>")
}

gtf_file <- args[1]
outdir <- args[2]

dir.create(file.path(outdir, "reference"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

feature_rds <- file.path(outdir, "reference", "gencode_v47_comprehensive_features.rds")

message("reading GENCODE V47 comprehensive from: ", gtf_file)

gtf <- import(gtf_file)

message("Imported ", length(gtf), " total GTF records.")

gtf_types <- data.table(feature_type = as.character(gtf$type))[, .N, by = feature_type][order(-N)]
fwrite(gtf_types, file.path(outdir, "tables", "gencode_v47_feature_type_counts.tsv"), sep = "\t")

mcols_gtf <- as.data.table(as.data.frame(mcols(gtf)))

gene_type_col <- if ("gene_type" %in% names(mcols_gtf)) {
  "gene_type"
} else if ("gene_biotype" %in% names(mcols_gtf)) {
  "gene_biotype"
} else {
  NA_character_
}

if (is.na(gene_type_col)) {
  stop("Could not find gene_type or gene_biotype in the GTF metadata.")
}

if (!"gene_id" %in% names(mcols_gtf)) {
  stop("Could not find gene_id in the GTF metadata.")
}

gene_name_col <- if ("gene_name" %in% names(mcols_gtf)) "gene_name" else NA_character_
tx_id_col <- if ("transcript_id" %in% names(mcols_gtf)) "transcript_id" else NA_character_

message("extracting gene features for same-strand gene overlap annotation.")

genes <- gtf[gtf$type == "gene"]

if (length(genes) == 0) {
  stop("Did not find any gene features in the GTF.")
}

mcols(genes)$gene_id <- as.character(mcols(genes)$gene_id)
mcols(genes)$gene_type <- as.character(mcols(genes)[[gene_type_col]])

if (!is.na(gene_name_col)) {
  mcols(genes)$gene_name <- as.character(mcols(genes)[[gene_name_col]])
} else {
  mcols(genes)$gene_name <- NA_character_
}

message("extracting transcript features for strand-specific upstream TSS calculation.")

transcripts <- gtf[gtf$type == "transcript"]

if (length(transcripts) == 0) {
  stop("Did not find any transcript features in the GTF.")
}

tx_dt <- as.data.table(as.data.frame(transcripts))
tx_meta <- as.data.table(as.data.frame(mcols(transcripts)))

tx_dt[, transcript_id := if (!is.na(tx_id_col)) as.character(tx_meta[[tx_id_col]]) else paste0("tx_", .I)]
tx_dt[, gene_id := as.character(tx_meta$gene_id)]
tx_dt[, gene_type := as.character(tx_meta[[gene_type_col]])]
tx_dt[, gene_name := if (!is.na(gene_name_col)) as.character(tx_meta[[gene_name_col]]) else NA_character_]

tx_dt <- tx_dt[strand %in% c("+", "-")]

# Define the TSS from transcript coordinates.
# For plus-strand transcripts, the TSS is the transcript start.
# For minus-strand transcripts, the TSS is the transcript end.
tx_dt[, tss := ifelse(strand == "+", start, end)]

tss_dt <- tx_dt[, .(
  seqnames = as.character(seqnames),
  strand = as.character(strand),
  transcript_id,
  gene_id,
  gene_name,
  gene_type,
  tss = as.integer(tss)
)]

message("building a TxDb object from the same GENCODE V47 GTF.")

txdb <- txdbmaker::makeTxDbFromGFF(gtf_file, format = "gtf")

message("extracting 5'UTR, 3'UTR, exon, and intron features from TxDb.")

five_utrs <- suppressWarnings(unlist(fiveUTRsByTranscript(txdb, use.names = TRUE), use.names = FALSE))
three_utrs <- suppressWarnings(unlist(threeUTRsByTranscript(txdb, use.names = TRUE), use.names = FALSE))
exons <- suppressWarnings(unlist(exonsBy(txdb, by = "tx", use.names = TRUE), use.names = FALSE))
introns <- suppressWarnings(unlist(intronsByTranscript(txdb, use.names = TRUE), use.names = FALSE))

features <- list(
  gtf_file = gtf_file,
  genes = genes,
  five_utrs = five_utrs,
  three_utrs = three_utrs,
  exons = exons,
  introns = introns,
  tss = tss_dt
)

saveRDS(features, feature_rds)

message("Saved the prepared GENCODE feature object to: ", feature_rds)

feature_summary <- data.table(
  feature = c("genes", "transcripts", "five_utrs", "three_utrs", "exons", "introns", "tss_records"),
  count = c(length(genes), length(transcripts), length(five_utrs), length(three_utrs), length(exons), length(introns), nrow(tss_dt))
)

fwrite(feature_summary, file.path(outdir, "tables", "gencode_v47_prepared_feature_summary.tsv"), sep = "\t")

gene_type_summary <- as.data.table(as.data.frame(mcols(genes)))
gene_type_summary <- gene_type_summary[, .N, by = gene_type][order(-N)]
fwrite(gene_type_summary, file.path(outdir, "tables", "gencode_v47_gene_type_summary.tsv"), sep = "\t")

p1 <- ggplot(feature_summary, aes(x = reorder(feature, count), y = count)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "GENCODE V47 prepared feature counts",
    x = "Feature",
    y = "Count"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(outdir, "plots", "gencode_v47_prepared_feature_counts.png"), p1, width = 7, height = 4, dpi = 300)

top_gene_types <- gene_type_summary[1:min(.N, 20)]

p2 <- ggplot(top_gene_types, aes(x = reorder(gene_type, N), y = N)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top GENCODE V47 gene types",
    x = "Gene type",
    y = "Number of genes"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(outdir, "plots", "gencode_v47_top_gene_types.png"), p2, width = 7, height = 5, dpi = 300)

message("Finished preparing GENCODE V47 features.")

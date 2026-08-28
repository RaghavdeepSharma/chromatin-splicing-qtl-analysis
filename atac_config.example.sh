#!/usr/bin/env bash
# Calderon ATAC-seq exact-style configuration.
# Copy this file, edit paths only, then run: bash scripts/01_atac_seq/00_run_pipeline.sh path/to/config.sh

set -euo pipefail

# Analysis output root. Override BASE or PROJECT_ROOT for your environment.
export BASE="${BASE:-${PROJECT_ROOT:-$PWD}/work/atac}"

# Manifest must be CSV or TSV with a header and at least these columns:
# sample_id,fastq1,fastq2
# Optional but recommended columns:
# biosample_id,donor,cell_type,condition,technical_replicate
# If biosample_id is absent or blank, sample_id is treated as the biological sample.
export MANIFEST="${MANIFEST:-/path/to/calderon_manifest.tsv}"

# Reference inputs.
# Bowtie2 hg19 index basename, e.g. /n/groups/shared_databases/bowtie2_indexes/hg19
export HG19_BOWTIE2_INDEX="${HG19_BOWTIE2_INDEX:-/path/to/hg19/bowtie2/index/base}"

# RefSeq annotation for TSS enrichment.
# Accepted input: BED-like refGene table with chrom, txStart, txEnd, strand columns, or UCSC refGene txt.
# The helper script detects common UCSC refGene column order.
export REFSEQ_ANNOTATION="${REFSEQ_ANNOTATION:-/path/to/refseq/refGene.txt}"

# ENCODE/hg19 blacklist BED.
export HG19_BLACKLIST="${HG19_BLACKLIST:-/path/to/hg19-blacklist.bed}"

# Optional chromosome sizes for bigWig generation.
export HG19_CHROMSIZES="${HG19_CHROMSIZES:-}"

# Conda environment used on O2.
export CONDA_ENV="${CONDA_ENV:-project}"

# Calderon Methods parameters.
export CUTADAPT_MINLEN="${CUTADAPT_MINLEN:-20}"
export CUTADAPT_OVERLAP="${CUTADAPT_OVERLAP:-5}"
export BOWTIE2_MAX_INSERT="${BOWTIE2_MAX_INSERT:-2000}"
export MAPQ_MIN="${MAPQ_MIN:-30}"
export SAMTOOLS_EXCLUDE_FLAGS="${SAMTOOLS_EXCLUDE_FLAGS:-1804}"
export SAMTOOLS_REQUIRE_FLAGS="${SAMTOOLS_REQUIRE_FLAGS:-2}"
export MACS2_GENOME="${MACS2_GENOME:-hs}"
export MAX_PEAK_WIDTH="${MAX_PEAK_WIDTH:-3000}"
export CONSENSUS_MIN_SUPPORT="${CONSENSUS_MIN_SUPPORT:-2}"

# Calderon QC thresholds.
export MIN_FILTERED_READS_PER_BIOSAMPLE="${MIN_FILTERED_READS_PER_BIOSAMPLE:-5000000}"
export MAX_MITO_FRAC="${MAX_MITO_FRAC:-0.25}"
export MAX_BLACKLIST_FRAC="${MAX_BLACKLIST_FRAC:-0.005}"
export MIN_TSS_ENRICHMENT="${MIN_TSS_ENRICHMENT:-4}"
export PCA_OUTLIER_MAD_Z="${PCA_OUTLIER_MAD_Z:-6}"

# Adapter. Leave as the standard Nextera transposase adapter unless you have verified otherwise.
export CUTADAPT_ADAPTER="${CUTADAPT_ADAPTER:-CTGTCTCTTATACACATCT}"

# SLURM defaults. Override at runtime if needed, e.g. export SLURM_PARTITION=medium.
export SLURM_PARTITION="${SLURM_PARTITION:-short}"
export ARRAY_CONCURRENCY="${ARRAY_CONCURRENCY:-50}"

# Counting method. Calderon used nucleoATAC get_count.
# The count script fails if get_count is unavailable unless ALLOW_BEDTOOLS_COUNT_FALLBACK=1.
export COUNT_METHOD="${COUNT_METHOD:-nucleoatac_get_count}"
export ALLOW_BEDTOOLS_COUNT_FALLBACK="${ALLOW_BEDTOOLS_COUNT_FALLBACK:-0}"

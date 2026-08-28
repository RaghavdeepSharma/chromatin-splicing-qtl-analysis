#!/usr/bin/env bash
# Example configuration for the sQTL integration/statistical-analysis modules.
# Copy to a private config file and edit paths for your environment.

export KERIMOV_BASE="${KERIMOV_BASE:-${PROJECT_ROOT:-$PWD}/work/sqtl}"
export KERIMOV_METADATA="$KERIMOV_BASE/01_metadata"
export KERIMOV_MANIFESTS="$KERIMOV_BASE/02_manifests"
export KERIMOV_RAW_SUSIE="$KERIMOV_BASE/03_raw_susie"
export KERIMOV_PROCESSED="$KERIMOV_BASE/05_processed"
export KERIMOV_RESULTS="$KERIMOV_BASE/06_results"
export KERIMOV_LOGS="$KERIMOV_BASE/logs"
export KERIMOV_TMP="$KERIMOV_BASE/tmp"

# External/public resources
export GENCODE_FEATURE_RDS="${GENCODE_FEATURE_RDS:-/path/to/gencode_features.rds}"
export HG38_FASTA="${HG38_FASTA:-/path/to/hg38.fa}"

# Optional environment name used by SLURM wrappers.
export CONDA_ENV="${CONDA_ENV:-project}"

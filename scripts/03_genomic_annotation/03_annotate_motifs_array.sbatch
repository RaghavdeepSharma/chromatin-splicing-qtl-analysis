#!/bin/bash
#SBATCH --job-name=motif_general_annot
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --partition=medium
#SBATCH --array=1-286%25

set -euo pipefail

# Source the environment with nounset disabled because activation may reference unset variables.
set +u
source ~/.bashrc
set -u

conda activate "${CONDA_ENV:-project}"

OUT="${MOTIF_ANNOTATION_ROOT:?Set MOTIF_ANNOTATION_ROOT}"
MOTIF_LIST=$OUT/motif_files.list
FEATURE_RDS=$OUT/reference/gencode_v47_comprehensive_features.rds

MOTIF_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MOTIF_LIST")

echo "SLURM task: ${SLURM_ARRAY_TASK_ID}"
echo "Motif file: $MOTIF_FILE"
echo "Feature RDS: $FEATURE_RDS"
echo "Output directory: $OUT"

if [[ ! -f "$MOTIF_FILE" ]]; then
  echo "ERROR: Motif file not found: $MOTIF_FILE"
  exit 1
fi

if [[ ! -f "$FEATURE_RDS" ]]; then
  echo "ERROR: Prepared GENCODE feature RDS not found: $FEATURE_RDS"
  exit 1
fi

Rscript "$(dirname "$0")/02_annotate_motif_features.R" \
  "$MOTIF_FILE" \
  "$FEATURE_RDS" \
  "$OUT"

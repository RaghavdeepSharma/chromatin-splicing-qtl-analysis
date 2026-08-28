#!/bin/bash
#SBATCH --job-name=q05_tfbs_acc
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-284%8
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -eo pipefail

set +u
source ~/.bashrc
conda activate "${CONDA_ENV:-project}"
set -u

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BASE="${PROJECT_ROOT:?Set PROJECT_ROOT}"
ACC_OUT=$BASE/tfbs_celltype_accessibility_20260701
BAMPE_Q05_OUT=$ACC_OUT/BAMPE_raw_processed_accessibility_q05_significant_peaks

LIST=$BAMPE_Q05_OUT/tables/q05_normal_284_general_annotated_motif_files.list
R_SCRIPT="$SCRIPT_DIR/03_annotate_accessibility_worker.R"
MANIFEST=$BAMPE_Q05_OUT/tables/q05_celltype_accessible_peak_files.hg38.tsv

motif_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")

echo "Task: $SLURM_ARRAY_TASK_ID"
echo "Motif file: $motif_file"

if [[ ! -s "$motif_file" ]]; then
  echo "ERROR: motif annotation file missing: $motif_file" >&2
  exit 1
fi

Rscript "$R_SCRIPT" \
  "$motif_file" \
  "$MANIFEST" \
  "$BAMPE_Q05_OUT"

echo "Done task $SLURM_ARRAY_TASK_ID"

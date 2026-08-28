#!/usr/bin/env bash
#SBATCH --job-name=cal_bigwig_C
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
set -euo pipefail
CONFIG=${1:?Usage: sbatch --array=1-N 10_generate_bigwig.sbatch config.sh}
set +u
source ~/.bashrc
set -u
source "$CONFIG"
conda activate "$CONDA_ENV"
PASS="$BASE/qc/pass_hard_qc_biosamples.txt"
biosample_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PASS")
[[ -n "$biosample_id" ]] || { echo "No passing biosample for task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
mkdir -p "$BASE/bigwig"
if command -v bamCoverage >/dev/null 2>&1; then
  bamCoverage \
    -b "$BASE/bam/merged/${biosample_id}.merged.dedup.bam" \
    -o "$BASE/bigwig/${biosample_id}.RPKM.extendReads.bin10.bw" \
    --binSize 10 \
    --normalizeUsing RPKM \
    --extendReads \
    -p "$SLURM_CPUS_PER_TASK"
else
  echo "ERROR: bamCoverage from deepTools is required for Calderon-style RPKM bigWig tracks." >&2
  exit 1
fi

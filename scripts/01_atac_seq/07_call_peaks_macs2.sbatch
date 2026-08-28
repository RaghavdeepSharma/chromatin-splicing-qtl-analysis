#!/usr/bin/env bash
#SBATCH --job-name=cal_macs2_C
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
set -euo pipefail
CONFIG=${1:?Usage: sbatch --array=1-N 07_call_peaks_macs2.sbatch config.sh}
set +u
source ~/.bashrc
set -u
source "$CONFIG"
conda activate "$CONDA_ENV"
PASS="$BASE/qc/pass_hard_qc_biosamples.txt"
biosample_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PASS")
[[ -n "$biosample_id" ]] || { echo "No passing biosample for task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
mkdir -p "$BASE/peaks/${biosample_id}"
macs2 callpeak \
  -t "$BASE/bam/merged/${biosample_id}.merged.dedup.bam" \
  -n "$biosample_id" \
  --outdir "$BASE/peaks/${biosample_id}" \
  --nomodel \
  --nolambda \
  --keep-dup all \
  --call-summits \
  -g "$MACS2_GENOME" \
  2> "$BASE/peaks/${biosample_id}/${biosample_id}.macs2.log"

#!/usr/bin/env bash
#SBATCH --job-name=cal_trim_C
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
set -euo pipefail
CONFIG=${1:?Usage: sbatch --array=1-N 02_trim_reads.sbatch config.sh}
set +u
source ~/.bashrc
set -u
source "$CONFIG"
conda activate "$CONDA_ENV"
LIBS="$BASE/metadata/libraries.tsv"
row=$(awk -v i="$SLURM_ARRAY_TASK_ID" 'BEGIN{FS=OFS="\t"} NR==i+1 {print}' "$LIBS")
[[ -n "$row" ]] || { echo "No manifest row for task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
IFS=$'\t' read -r sample_id biosample_id fastq1 fastq2 donor cell_type condition techrep <<< "$row"
mkdir -p "$BASE/trimmed" "$BASE/logs"
cutadapt \
  -a "$CUTADAPT_ADAPTER" \
  -A "$CUTADAPT_ADAPTER" \
  --minimum-length "$CUTADAPT_MINLEN" \
  --overlap "$CUTADAPT_OVERLAP" \
  -j "$SLURM_CPUS_PER_TASK" \
  -o "$BASE/trimmed/${sample_id}_R1.trimmed.fastq.gz" \
  -p "$BASE/trimmed/${sample_id}_R2.trimmed.fastq.gz" \
  "$fastq1" "$fastq2" \
  > "$BASE/trimmed/${sample_id}.cutadapt.txt" 2>&1

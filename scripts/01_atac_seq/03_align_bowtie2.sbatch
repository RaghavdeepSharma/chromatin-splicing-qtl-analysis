#!/usr/bin/env bash
#SBATCH --job-name=cal_align_C
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
set -euo pipefail
CONFIG=${1:?Usage: sbatch --array=1-N 03_align_bowtie2.sbatch config.sh}
set +u
source ~/.bashrc
set -u
source "$CONFIG"
conda activate "$CONDA_ENV"
LIBS="$BASE/metadata/libraries.tsv"
row=$(awk -v i="$SLURM_ARRAY_TASK_ID" 'BEGIN{FS=OFS="\t"} NR==i+1 {print}' "$LIBS")
IFS=$'\t' read -r sample_id biosample_id fastq1 fastq2 donor cell_type condition techrep <<< "$row"
mkdir -p "$BASE/bam/sorted"
bowtie2 \
  -x "$HG19_BOWTIE2_INDEX" \
  -1 "$BASE/trimmed/${sample_id}_R1.trimmed.fastq.gz" \
  -2 "$BASE/trimmed/${sample_id}_R2.trimmed.fastq.gz" \
  -X "$BOWTIE2_MAX_INSERT" \
  -p "$SLURM_CPUS_PER_TASK" \
  2> "$BASE/bam/sorted/${sample_id}.bowtie2.log" | \
  samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$BASE/bam/sorted/${sample_id}.sorted.bam" -
samtools index -@ "$SLURM_CPUS_PER_TASK" "$BASE/bam/sorted/${sample_id}.sorted.bam"

#!/usr/bin/env bash
#SBATCH --job-name=cal_dedup_C
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
set -euo pipefail
CONFIG=${1:?Usage: sbatch --array=1-N 05_deduplicate_bam.sbatch config.sh}
set +u
source ~/.bashrc
set -u
source "$CONFIG"
conda activate "$CONDA_ENV"
LIBS="$BASE/metadata/libraries.tsv"
row=$(awk -v i="$SLURM_ARRAY_TASK_ID" 'BEGIN{FS=OFS="\t"} NR==i+1 {print}' "$LIBS")
IFS=$'\t' read -r sample_id biosample_id fastq1 fastq2 donor cell_type condition techrep <<< "$row"
mkdir -p "$BASE/bam/dedup" "$BASE/bam/tmp" "$BASE/qc"
IN="$BASE/bam/filtered/${sample_id}.filtered.bam"
RG="$BASE/bam/tmp/${sample_id}.rg.bam"
OUT="$BASE/bam/dedup/${sample_id}.dedup.bam"
samtools addreplacerg -@ "$SLURM_CPUS_PER_TASK" \
  -r "ID:${sample_id}\tSM:${biosample_id}\tPL:ILLUMINA\tLB:${sample_id}" \
  -o "$RG" "$IN"
picard MarkDuplicates \
  -I "$RG" \
  -O "$OUT" \
  -M "$BASE/qc/${sample_id}.dup_metrics.txt" \
  --REMOVE_DUPLICATES true \
  --VALIDATION_STRINGENCY LENIENT \
  --TMP_DIR "$BASE/bam/tmp"
samtools index -@ "$SLURM_CPUS_PER_TASK" "$OUT"
rm -f "$RG"
DEDUP=$(samtools view -c "$OUT")
awk -v s="$sample_id" -v b="$biosample_id" -v d="$DEDUP" 'BEGIN{OFS="\t"; print s,b,d}' > "$BASE/qc/${sample_id}.dedup_counts.tsv"

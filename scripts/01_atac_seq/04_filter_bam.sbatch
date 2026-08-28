#!/usr/bin/env bash
#SBATCH --job-name=cal_filter_C
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
set -euo pipefail
CONFIG=${1:?Usage: sbatch --array=1-N 04_filter_bam.sbatch config.sh}
set +u
source ~/.bashrc
set -u
source "$CONFIG"
conda activate "$CONDA_ENV"
LIBS="$BASE/metadata/libraries.tsv"
row=$(awk -v i="$SLURM_ARRAY_TASK_ID" 'BEGIN{FS=OFS="\t"} NR==i+1 {print}' "$LIBS")
IFS=$'\t' read -r sample_id biosample_id fastq1 fastq2 donor cell_type condition techrep <<< "$row"
mkdir -p "$BASE/bam/filtered" "$BASE/qc"
IN="$BASE/bam/sorted/${sample_id}.sorted.bam"
OUT="$BASE/bam/filtered/${sample_id}.filtered.bam"
CHROMS=$(samtools view -H "$IN" | awk -F'\t' '/^@SQ/ {for(i=1;i<=NF;i++) if($i ~ /^SN:/){c=substr($i,4); if(c != "chrM" && c != "MT" && c != "M") printf c" "}}')
samtools view -@ "$SLURM_CPUS_PER_TASK" \
  -F "$SAMTOOLS_EXCLUDE_FLAGS" \
  -f "$SAMTOOLS_REQUIRE_FLAGS" \
  -q "$MAPQ_MIN" \
  -b "$IN" $CHROMS \
  -o "$OUT"
samtools index -@ "$SLURM_CPUS_PER_TASK" "$OUT"
# Basic read accounting, including mitochondrial fraction before chrM removal.
TOTAL=$(samtools view -c "$IN")
MITO=$(samtools view -c "$IN" chrM 2>/dev/null || true)
if [[ -z "$MITO" ]]; then MITO=$(samtools view -c "$IN" MT 2>/dev/null || echo 0); fi
FILTERED=$(samtools view -c "$OUT")
awk -v s="$sample_id" -v b="$biosample_id" -v total="$TOTAL" -v mito="$MITO" -v filt="$FILTERED" 'BEGIN{OFS="\t"; print s,b,total,mito,(total>0?mito/total:0),filt}' > "$BASE/qc/${sample_id}.read_counts.tsv"

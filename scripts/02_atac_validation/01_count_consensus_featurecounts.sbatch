#!/usr/bin/env bash
#SBATCH --job-name=bampe_counts
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --partition=medium

set -euo pipefail

set +u
source ~/.bashrc
set -u
conda activate "${CONDA_ENV:-project}"

BASE="${PROJECT_ROOT:?Set PROJECT_ROOT to the project root}"
RUN=$BASE/bampe_peak_validation_20260622
CONS=$RUN/consensus
COUNT=$RUN/counts
mkdir -p "$COUNT"

BED=$CONS/master_consensus_bampe_hg19.bed
SAF=$CONS/master_consensus_bampe_hg19.saf
BAMLIST=$COUNT/bampe_bam_list.txt
OUT=$COUNT/bampe_consensus_featureCounts.txt

echo "Creating SAF annotation..."
awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"}
{
  peak=$1"_"$2"_"$3
  print peak,$1,$2+1,$3,"."
}' "$BED" > "$SAF"

echo "Creating BAM list..."
awk -v base="$BASE" '{print base"/bam/merged/"$1".merged.dedup.bam"}' \
  "$BASE/qc/pass_hard_qc_biosamples.txt" > "$BAMLIST"

echo "Checking missing BAMs..."
while read bam; do
  [[ -s "$bam" ]] || { echo "Missing BAM: $bam" >&2; exit 1; }
done < "$BAMLIST"

echo "Running featureCounts..."
featureCounts \
  -T 8 \
  -p \
  -B \
  -C \
  -F SAF \
  -a "$SAF" \
  -o "$OUT" \
  $(cat "$BAMLIST")

echo "DONE featureCounts"
date

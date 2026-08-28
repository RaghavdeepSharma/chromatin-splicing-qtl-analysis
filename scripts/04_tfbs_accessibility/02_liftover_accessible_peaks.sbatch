#!/bin/bash
#SBATCH --job-name=q05_liftover
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -eo pipefail

# Disable nounset during Conda activation because activation scripts may reference unset variables.
set +u
source ~/.bashrc
conda activate "${CONDA_ENV:-project}"
set -u

BASE="${PROJECT_ROOT:?Set PROJECT_ROOT to the project root}"
ACC_OUT=$BASE/tfbs_celltype_accessibility_20260701
BAMPE_Q05_OUT=$ACC_OUT/BAMPE_raw_processed_accessibility_q05_significant_peaks
CHAIN=$BASE/ref/liftover/hg19ToHg38.over.chain.gz

echo "Checking liftOver:"
which liftOver
liftOver 2>&1 | head -5 || true

echo "Checking chain file:"
ls -lh "$CHAIN"

echo "Checking hg19 manifest:"
ls -lh $BAMPE_Q05_OUT/tables/q05_celltype_accessible_peak_files.hg19.tsv
wc -l $BAMPE_Q05_OUT/tables/q05_celltype_accessible_peak_files.hg19.tsv

mkdir -p \
  $BAMPE_Q05_OUT/celltype_accessible_beds_hg38 \
  $BAMPE_Q05_OUT/tables

echo -e "cell_type\thg19_bed\thg38_bed\tunmapped_bed\tn_hg19\tn_hg38\tn_unmapped" \
  > $BAMPE_Q05_OUT/tables/q05_liftover_qc.tsv

echo -e "cell_type\tbed_file" \
  > $BAMPE_Q05_OUT/tables/q05_celltype_accessible_peak_files.hg38.tsv

while read cell_type hg19_bed; do
  if [[ "$cell_type" == "cell_type" ]]; then
    continue
  fi

  echo "Lifting $cell_type"

  if [[ ! -s "$hg19_bed" ]]; then
    echo "ERROR: missing hg19 BED: $hg19_bed" >&2
    exit 1
  fi

  hg38_bed=$BAMPE_Q05_OUT/celltype_accessible_beds_hg38/${cell_type}.q05.accessible.bampe.hg38.bed
  unmapped=$BAMPE_Q05_OUT/celltype_accessible_beds_hg38/${cell_type}.q05.accessible.bampe.hg19_to_hg38.unmapped.bed

  liftOver "$hg19_bed" "$CHAIN" "$hg38_bed" "$unmapped"

  sort -k1,1 -k2,2n "$hg38_bed" -o "$hg38_bed"

  n_hg19=$(wc -l < "$hg19_bed")
  n_hg38=$(wc -l < "$hg38_bed")
  n_unmapped=$(grep -v "^#" "$unmapped" | wc -l)

  echo -e "${cell_type}\t${hg19_bed}\t${hg38_bed}\t${unmapped}\t${n_hg19}\t${n_hg38}\t${n_unmapped}" \
    >> $BAMPE_Q05_OUT/tables/q05_liftover_qc.tsv

  echo -e "${cell_type}\t${hg38_bed}" \
    >> $BAMPE_Q05_OUT/tables/q05_celltype_accessible_peak_files.hg38.tsv

done < $BAMPE_Q05_OUT/tables/q05_celltype_accessible_peak_files.hg19.tsv

echo "Done liftOver."

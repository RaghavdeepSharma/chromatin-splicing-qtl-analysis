#!/bin/bash
#SBATCH --job-name=q05_acc_beds
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

BASE="${PROJECT_ROOT:?Set PROJECT_ROOT to the project root}"
BAMPE_DIR=$BASE/bampe_peak_validation_20260622
ACC_OUT=$BASE/tfbs_celltype_accessibility_20260701
BAMPE_Q05_OUT=$ACC_OUT/BAMPE_raw_processed_accessibility_q05_significant_peaks

CONSENSUS=$BAMPE_DIR/consensus/master_consensus_bampe_hg19.bed
META=$BAMPE_DIR/differential_atac/bampe_sample_metadata.tsv
PEAK_DIR=$BAMPE_DIR/peaks

QLOG10_MIN=1.30103
MIN_SAMPLES=2

mkdir -p \
  $BAMPE_Q05_OUT/tmp/q05_sample_consensus_hits \
  $BAMPE_Q05_OUT/celltype_accessible_beds_hg19 \
  $BAMPE_Q05_OUT/tables

CONSENSUS_ID_BED=$BAMPE_Q05_OUT/tmp/master_consensus_bampe_hg19.with_peak_id.bed

echo "Making consensus BED with peak IDs..."
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3}' "$CONSENSUS" > "$CONSENSUS_ID_BED"

echo "Building sample-cell type manifest..."
awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1,$3}' "$META" > $BAMPE_Q05_OUT/tables/sample_celltype.tsv

echo -e "cell_type\tsafe_cell_type\tn_samples\tqvalue_cutoff\tqlog10_min\tmin_samples\tn_consensus_peaks\tn_accessible_peaks\tfraction_accessible_peaks\tbed_file_hg19" \
  > $BAMPE_Q05_OUT/tables/q05_celltype_accessible_peak_summary_hg19.tsv

echo -e "cell_type\tbed_file" \
  > $BAMPE_Q05_OUT/tables/q05_celltype_accessible_peak_files.hg19.tsv

N_CONSENSUS=$(wc -l < "$CONSENSUS_ID_BED")

cut -f2 $BAMPE_Q05_OUT/tables/sample_celltype.tsv | sort -u | while read cell_type; do
  safe_cell_type=$(echo "$cell_type" | sed 's/[^A-Za-z0-9_]/_/g')

  echo "Processing cell type: $cell_type"

  sample_list=$BAMPE_Q05_OUT/tmp/${safe_cell_type}.samples.txt
  awk -v ct="$cell_type" 'BEGIN{FS=OFS="\t"} $2==ct {print $1}' \
    $BAMPE_Q05_OUT/tables/sample_celltype.tsv > "$sample_list"

  n_samples=$(wc -l < "$sample_list")

  all_hits=$BAMPE_Q05_OUT/tmp/${safe_cell_type}.q05.consensus_peak_hits.all.txt
  : > "$all_hits"

  while read sample; do
    narrowpeak=$PEAK_DIR/$sample/${sample}_peaks.narrowPeak

    if [[ ! -s "$narrowpeak" ]]; then
      echo "WARNING: missing narrowPeak for sample: $sample" >&2
      continue
    fi

    q05_bed=$BAMPE_Q05_OUT/tmp/${sample}.q05.narrowPeak.bed

    # MACS narrowPeak column 9 is -log10(qvalue).
    # q < 0.05 means -log10(qvalue) >= 1.30103.
    awk -v q="$QLOG10_MIN" 'BEGIN{OFS="\t"} $9 >= q {print $1,$2,$3}' "$narrowpeak" \
      | sort -k1,1 -k2,2n \
      > "$q05_bed"

    bedtools intersect -u \
      -a "$CONSENSUS_ID_BED" \
      -b "$q05_bed" \
      | cut -f4 >> "$all_hits"

    rm -f "$q05_bed"
  done < "$sample_list"

  accessible_ids=$BAMPE_Q05_OUT/tmp/${safe_cell_type}.q05.accessible.consensus_peak_ids.txt

  sort "$all_hits" \
    | uniq -c \
    | awk -v min="$MIN_SAMPLES" '$1 >= min {print $2}' \
    > "$accessible_ids"

  out_bed=$BAMPE_Q05_OUT/celltype_accessible_beds_hg19/${safe_cell_type}.q05.accessible.bampe.hg19.bed

  awk 'BEGIN{FS=OFS="\t"} NR==FNR {keep[$1]=1; next} ($4 in keep) {print $1,$2,$3}' \
    "$accessible_ids" \
    "$CONSENSUS_ID_BED" \
    | sort -k1,1 -k2,2n \
    > "$out_bed"

  n_accessible=$(wc -l < "$out_bed")
  frac=$(awk -v a="$n_accessible" -v n="$N_CONSENSUS" 'BEGIN{if(n>0) print a/n; else print "NA"}')

  echo -e "${cell_type}\t${safe_cell_type}\t${n_samples}\t0.05\t${QLOG10_MIN}\t${MIN_SAMPLES}\t${N_CONSENSUS}\t${n_accessible}\t${frac}\t${out_bed}" \
    >> $BAMPE_Q05_OUT/tables/q05_celltype_accessible_peak_summary_hg19.tsv

  echo -e "${safe_cell_type}\t${out_bed}" \
    >> $BAMPE_Q05_OUT/tables/q05_celltype_accessible_peak_files.hg19.tsv
done

echo "Done creating q<0.05 hg19 cell-type significant peak BEDs."

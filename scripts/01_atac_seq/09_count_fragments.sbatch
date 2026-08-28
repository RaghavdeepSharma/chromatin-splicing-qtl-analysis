#!/usr/bin/env bash
#SBATCH --job-name=cal_counts_C
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
set -euo pipefail
CONFIG=${1:?Usage: sbatch 09_count_fragments.sbatch config.sh}
set +u
source ~/.bashrc
set -u
source "$CONFIG"
conda activate "$CONDA_ENV"
mkdir -p "$BASE/counts"
PEAKS="$BASE/consensus/master_consensus_hg19.calderon_exact.bed"
PASS="$BASE/qc/pass_hard_qc_biosamples.txt"
[[ -s "$PEAKS" ]] || { echo "Missing consensus BED: $PEAKS" >&2; exit 1; }
[[ -s "$PASS" ]] || { echo "Missing pass list: $PASS" >&2; exit 1; }

BAM_LIST="$BASE/counts/bams_for_counting.txt"
: > "$BAM_LIST"
while read -r bid; do
  bam="$BASE/bam/merged/${bid}.merged.dedup.bam"
  [[ -s "$bam" ]] || { echo "Missing BAM for $bid: $bam" >&2; exit 1; }
  echo -e "$bid\t$bam" >> "$BAM_LIST"
done < "$PASS"

if conda run -n nucleoatac_py2 pyatac counts -h >/dev/null 2>&1; then
  echo "Using NucleoATAC/PyATAC counts, matching Calderon-style fragment counting."

  OUTMATRIX="$BASE/counts/atac_fragment_counts.pyatac.tsv"
  TMPDIR="$BASE/counts/pyatac_per_sample"
  mkdir -p "$TMPDIR"

  awk 'BEGIN{OFS="\t"} {print $1":"$2"-"$3, $1, $2, $3}' "$PEAKS" > "$BASE/counts/peak_annotation.tmp.tsv"

  cp "$BASE/counts/peak_annotation.tmp.tsv" "$OUTMATRIX"

  while IFS=$'\t' read -r bid bam; do
    echo "Counting fragments for $bid"

    conda run -n nucleoatac_py2 pyatac counts \
      --bam "$bam" \
      --bed "$PEAKS" \
      --out "$TMPDIR/$bid"

    COUNT_FILE=$(ls "$TMPDIR"/"$bid"* 2>/dev/null | head -1)

    if [[ ! -s "$COUNT_FILE" ]]; then
      echo "ERROR: pyatac counts did not create output for $bid" >&2
      exit 1
    fi

    awk 'BEGIN{OFS="\t"} {print $NF}' "$COUNT_FILE" > "$TMPDIR/$bid.counts.only"

    paste "$OUTMATRIX" "$TMPDIR/$bid.counts.only" > "$TMPDIR/counts.next.tmp"
    mv "$TMPDIR/counts.next.tmp" "$OUTMATRIX"
  done < "$BAM_LIST"

  HEADER=$(cut -f1 "$BAM_LIST" | tr '\n' '\t' | sed 's/\t$//')
  sed -i "1ipeak_id\tchr\tstart\tend\t$HEADER" "$OUTMATRIX"

  echo "Final PyATAC/NucleoATAC count matrix: $OUTMATRIX"

elif [[ "$ALLOW_BEDTOOLS_COUNT_FALLBACK" == "1" ]]; then
  echo "WARN: get_count unavailable; using bedtools coverage fallback. This is NOT Calderon-exact." >&2
  cp "$PEAKS" "$BASE/counts/peaks.tmp.bed"
  awk 'BEGIN{OFS="\t"} {print $1":"$2"-"$3}' "$PEAKS" > "$BASE/counts/peak_ids.tmp.txt"
  paste "$BASE/counts/peak_ids.tmp.txt" "$PEAKS" > "$BASE/counts/atac_fragment_counts.bedtools.tsv"
  while IFS=$'\t' read -r bid bam; do
    bedtools coverage -a "$PEAKS" -b "$bam" -counts | awk '{print $NF}' > "$BASE/counts/${bid}.counts.tmp"
    paste "$BASE/counts/atac_fragment_counts.bedtools.tsv" "$BASE/counts/${bid}.counts.tmp" > "$BASE/counts/counts.next.tmp"
    mv "$BASE/counts/counts.next.tmp" "$BASE/counts/atac_fragment_counts.bedtools.tsv"
  done < "$BAM_LIST"
else
  echo "ERROR: nucleoATAC get_count not found. Calderon used get_count; install nucleoATAC or set ALLOW_BEDTOOLS_COUNT_FALLBACK=1 for non-exact fallback." >&2
  exit 1
fi

#!/usr/bin/env bash
#SBATCH --job-name=cal_consensus_C
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
set -euo pipefail
CONFIG=${1:?Usage: sbatch 08_build_consensus_peaks.sbatch config.sh}
set +u
source ~/.bashrc
set -u
source "$CONFIG"
conda activate "$CONDA_ENV"
mkdir -p "$BASE/consensus"
PASS="$BASE/qc/pass_hard_qc_biosamples.txt"
[[ -s "$PASS" ]] || { echo "No passing biosamples file: $PASS" >&2; exit 1; }

# Use only peaks from samples that passed Calderon hard QC.
: > "$BASE/consensus/narrowpeak_files.used.txt"
while read -r bid; do
  f="$BASE/peaks/${bid}/${bid}_peaks.narrowPeak"
  [[ -s "$f" ]] || { echo "ERROR: missing narrowPeak for passing biosample $bid: $f" >&2; exit 1; }
  echo "$f" >> "$BASE/consensus/narrowpeak_files.used.txt"
done < "$PASS"

# Label each interval with its biosample so support can count unique samples.
awk 'BEGIN{OFS="\t"} {n=split(FILENAME,a,"/"); bid=a[n-1]; print $1,$2,$3,bid}' $(cat "$BASE/consensus/narrowpeak_files.used.txt") \
  | sort -k1,1 -k2,2n > "$BASE/consensus/all_peaks_labeled.sorted.bed"

# Calderon consensus: merge overlapping >=1 bp peaks identified in at least two samples.
bedtools merge \
  -i "$BASE/consensus/all_peaks_labeled.sorted.bed" \
  -c 4 \
  -o distinct,count_distinct \
  > "$BASE/consensus/all_peaks_merged_with_support.bed"

awk -v min_sup="$CONSENSUS_MIN_SUPPORT" 'BEGIN{OFS="\t"} $5>=min_sup {print $1,$2,$3,$4,$5}' \
  "$BASE/consensus/all_peaks_merged_with_support.bed" \
  > "$BASE/consensus/consensus_support_ge${CONSENSUS_MIN_SUPPORT}.bed"

# Calderon final filter: remove non-autosomal and >3 kb peaks.
awk -v maxw="$MAX_PEAK_WIDTH" 'BEGIN{OFS="\t"} $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ && ($3-$2) <= maxw {print $1,$2,$3}' \
  "$BASE/consensus/consensus_support_ge${CONSENSUS_MIN_SUPPORT}.bed" \
  > "$BASE/consensus/master_consensus_hg19.calderon_exact.bed"

# Also keep support metadata for auditing.
awk -v maxw="$MAX_PEAK_WIDTH" 'BEGIN{OFS="\t"} $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ && ($3-$2) <= maxw {print $1,$2,$3,$4,$5}' \
  "$BASE/consensus/consensus_support_ge${CONSENSUS_MIN_SUPPORT}.bed" \
  > "$BASE/consensus/master_consensus_hg19.calderon_exact.with_support.bed"

awk 'BEGIN{n=0} {w=$3-$2; a[++n]=w; sum+=w} END{if(n==0){print "no peaks"; exit 1}; asort(a); med=(n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2); print "n_peaks\t"n"\nmedian_width\t"med"\nmean_width\t"sum/n}' \
  "$BASE/consensus/master_consensus_hg19.calderon_exact.bed" \
  > "$BASE/consensus/master_consensus_summary.tsv"
cat "$BASE/consensus/master_consensus_summary.tsv"

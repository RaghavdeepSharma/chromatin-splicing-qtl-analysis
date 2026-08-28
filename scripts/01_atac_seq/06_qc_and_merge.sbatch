#!/usr/bin/env bash
#SBATCH --job-name=cal_qc_C
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
set -euo pipefail
CONFIG=${1:?Usage: sbatch 06_qc_and_merge.sbatch config.sh}
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
set +u
source ~/.bashrc
set -u
source "$CONFIG"
conda activate "$CONDA_ENV"
mkdir -p "$BASE/qc" "$BASE/bam/merged"

# Merge all available per-library read accounting.
{
  echo -e "sample_id\tbiosample_id\ttotal_aligned_reads\tmitochondrial_reads\tmitochondrial_fraction\tfiltered_reads"
  cat "$BASE"/qc/*.read_counts.tsv 2>/dev/null || true
} > "$BASE/qc/library_read_counts.tsv"
{
  echo -e "sample_id\tbiosample_id\tdedup_reads"
  cat "$BASE"/qc/*.dedup_counts.tsv 2>/dev/null || true
} > "$BASE/qc/library_dedup_counts.tsv"

# Merge technical replicates by biosample_id before downstream peak calling.
awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1,$2}' "$BASE/metadata/libraries.tsv" | sort -k2,2 > "$BASE/metadata/library_to_biosample.tsv"
cut -f2 "$BASE/metadata/library_to_biosample.tsv" | sort -u > "$BASE/metadata/biosample_ids.txt"
while read -r bid; do
  [[ -n "$bid" ]] || continue
  mapfile -t libs < <(awk -v b="$bid" 'BEGIN{FS="\t"} $2==b {print $1}' "$BASE/metadata/library_to_biosample.tsv")
  bams=()
  for s in "${libs[@]}"; do
    [[ -s "$BASE/bam/dedup/${s}.dedup.bam" ]] && bams+=("$BASE/bam/dedup/${s}.dedup.bam")
  done
  if [[ "${#bams[@]}" -eq 0 ]]; then
    echo "WARN: no dedup BAMs for biosample $bid" >&2
    continue
  elif [[ "${#bams[@]}" -eq 1 ]]; then
    ln -sf "${bams[0]}" "$BASE/bam/merged/${bid}.merged.dedup.bam"
    ln -sf "${bams[0]}.bai" "$BASE/bam/merged/${bid}.merged.dedup.bam.bai" 2>/dev/null || true
  else
    samtools merge -@ "$SLURM_CPUS_PER_TASK" -f "$BASE/bam/merged/${bid}.merged.dedup.bam" "${bams[@]}"
    samtools index -@ "$SLURM_CPUS_PER_TASK" "$BASE/bam/merged/${bid}.merged.dedup.bam"
  fi
done < "$BASE/metadata/biosample_ids.txt"

# Compute blacklist fraction on merged deduplicated BAMs.
echo -e "biosample_id\tblacklist_reads\tdedup_reads\tblacklist_fraction" > "$BASE/qc/biosample_blacklist.tsv"
while read -r bid; do
  bam="$BASE/bam/merged/${bid}.merged.dedup.bam"
  [[ -s "$bam" ]] || continue
  total=$(samtools view -c "$bam")
  black=$(bedtools intersect -u -abam "$bam" -b "$HG19_BLACKLIST" 2>/dev/null | samtools view -c - || echo 0)
  awk -v b="$bid" -v bl="$black" -v t="$total" 'BEGIN{OFS="\t"; print b,bl,t,(t>0?bl/t:0)}'
done < "$BASE/metadata/biosample_ids.txt" >> "$BASE/qc/biosample_blacklist.tsv"

# Create RefSeq TSS BED and approximate TSS enrichment relative to +/-2kb flanks.
python "$SCRIPT_DIR/11_tss_enrichment.py" \
  --refseq "$REFSEQ_ANNOTATION" \
  --bam-dir "$BASE/bam/merged" \
  --out "$BASE/qc/biosample_tss_enrichment.tsv" \
  --threads "$SLURM_CPUS_PER_TASK"

# Aggregate QC and make pass/fail table using Calderon thresholds.
python - <<'PY' "$BASE" "$MIN_FILTERED_READS_PER_BIOSAMPLE" "$MAX_MITO_FRAC" "$MAX_BLACKLIST_FRAC" "$MIN_TSS_ENRICHMENT"
import sys, pandas as pd, numpy as np, os
base, min_reads, max_mito, max_black, min_tss = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
libs = pd.read_csv(f'{base}/metadata/libraries.tsv', sep='\t')
reads = pd.read_csv(f'{base}/qc/library_read_counts.tsv', sep='\t')
dedup = pd.read_csv(f'{base}/qc/library_dedup_counts.tsv', sep='\t')
black = pd.read_csv(f'{base}/qc/biosample_blacklist.tsv', sep='\t')
tss = pd.read_csv(f'{base}/qc/biosample_tss_enrichment.tsv', sep='\t')
# Per Calderon wording: donor samples with fewer than 5M reads from all technical replicates that passed quality filters are excluded.
agg = reads.groupby('biosample_id', as_index=False).agg(total_aligned_reads=('total_aligned_reads','sum'), mitochondrial_reads=('mitochondrial_reads','sum'), filtered_reads=('filtered_reads','sum'))
agg['mitochondrial_fraction'] = np.where(agg.total_aligned_reads>0, agg.mitochondrial_reads/agg.total_aligned_reads, 0)
ded = dedup.groupby('biosample_id', as_index=False).agg(dedup_reads=('dedup_reads','sum'))
qc = agg.merge(ded, on='biosample_id', how='outer').merge(black, on='biosample_id', how='left', suffixes=('','_black')).merge(tss, on='biosample_id', how='left')
meta = libs.groupby('biosample_id', as_index=False).agg(donor=('donor','first'), cell_type=('cell_type','first'), condition=('condition','first'), n_libraries=('sample_id','count'))
qc = meta.merge(qc, on='biosample_id', how='left')
qc['fail_reasons']=''
def add_reason(mask, reason):
    qc.loc[mask, 'fail_reasons'] = qc.loc[mask, 'fail_reasons'].astype(str) + reason + ';'
add_reason(qc['dedup_reads'].fillna(0) < min_reads, 'dedup_reads_lt_5M')
add_reason(qc['mitochondrial_fraction'].fillna(1) > max_mito, 'mitochondrial_gt_25pct')
add_reason(qc['blacklist_fraction'].fillna(1) > max_black, 'blacklist_gt_0.5pct')
add_reason(qc['tss_enrichment'].fillna(0) < min_tss, 'tss_lt_4')
qc['pass_hard_qc'] = qc['fail_reasons'].eq('')
qc.to_csv(f'{base}/qc/biosample_qc.tsv', sep='\t', index=False)
qc.loc[qc.pass_hard_qc, ['biosample_id']].to_csv(f'{base}/qc/pass_hard_qc_biosamples.txt', sep='\t', index=False, header=False)
qc.loc[~qc.pass_hard_qc, ['biosample_id','fail_reasons']].to_csv(f'{base}/qc/fail_hard_qc_biosamples.tsv', sep='\t', index=False)
print(qc[['biosample_id','dedup_reads','mitochondrial_fraction','blacklist_fraction','tss_enrichment','pass_hard_qc','fail_reasons']].to_string(index=False))
PY

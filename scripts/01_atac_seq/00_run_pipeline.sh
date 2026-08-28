#!/usr/bin/env bash
# Submit Calderon-exact ATAC-seq preprocessing without hardcoded sample counts.
# Usage: bash scripts/01_atac_seq/00_run_pipeline.sh config/atac_config.example.sh
set -euo pipefail
CONFIG=${1:?Usage: bash 00_run_pipeline.sh config.sh}
CONFIG=$(readlink -f "$CONFIG")
source "$CONFIG"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

bash "$SCRIPT_DIR/01_validate_inputs.sh" "$CONFIG"

N_LIBS=$(awk 'NR>1' "$BASE/metadata/libraries.tsv" | wc -l | tr -d ' ')
[[ "$N_LIBS" -gt 0 ]] || { echo "No libraries found after validation" >&2; exit 1; }
ARRAY_LIBS="1-${N_LIBS}%${ARRAY_CONCURRENCY}"

echo "Submitting Calderon-style preprocessing for $N_LIBS libraries"

default_sbatch=(--partition="$SLURM_PARTITION" --time="12:00:00")
JOB02=$(sbatch "${default_sbatch[@]}" --array="$ARRAY_LIBS" --output="$BASE/logs/trim_%A_%a.out" --error="$BASE/logs/trim_%A_%a.err" --parsable "$SCRIPT_DIR/02_trim_reads.sbatch" "$CONFIG")
JOB03=$(sbatch "${default_sbatch[@]}" --array="$ARRAY_LIBS" --dependency=afterok:$JOB02 --output="$BASE/logs/align_%A_%a.out" --error="$BASE/logs/align_%A_%a.err" --parsable "$SCRIPT_DIR/03_align_bowtie2.sbatch" "$CONFIG")
JOB04=$(sbatch "${default_sbatch[@]}" --array="$ARRAY_LIBS" --dependency=afterok:$JOB03 --output="$BASE/logs/filter_%A_%a.out" --error="$BASE/logs/filter_%A_%a.err" --parsable "$SCRIPT_DIR/04_filter_bam.sbatch" "$CONFIG")
JOB05=$(sbatch "${default_sbatch[@]}" --array="$ARRAY_LIBS" --dependency=afterok:$JOB04 --output="$BASE/logs/dedup_%A_%a.out" --error="$BASE/logs/dedup_%A_%a.err" --parsable "$SCRIPT_DIR/05_deduplicate_bam.sbatch" "$CONFIG")
JOB06=$(sbatch "${default_sbatch[@]}" --dependency=afterok:$JOB05 --output="$BASE/logs/qc_%j.out" --error="$BASE/logs/qc_%j.err" --parsable "$SCRIPT_DIR/06_qc_and_merge.sbatch" "$CONFIG")

# Number of passing biosamples is unknown until QC finishes, so submit a small follow-up launcher.
mkdir -p "$BASE/scripts_runtime"
LAUNCH_PEAKS="$BASE/scripts_runtime/launch_after_qc.sh"
cat > "$LAUNCH_PEAKS" <<LAUNCH
#!/usr/bin/env bash
set -euo pipefail
source "$CONFIG"
SCRIPT_DIR="$SCRIPT_DIR"
N_PASS=\$(wc -l < "\$BASE/qc/pass_hard_qc_biosamples.txt" | tr -d ' ')
if [[ "\$N_PASS" -eq 0 ]]; then echo "No hard-QC passing biosamples" >&2; exit 1; fi
ARRAY_PASS="1-\${N_PASS}%${ARRAY_CONCURRENCY}"
JOB07=\$(sbatch --partition="$SLURM_PARTITION" --time="12:00:00" --array="\$ARRAY_PASS" --output="\$BASE/logs/macs2_%A_%a.out" --error="\$BASE/logs/macs2_%A_%a.err" --parsable "\$SCRIPT_DIR/07_call_peaks_macs2.sbatch" "$CONFIG")
JOB08=\$(sbatch --partition="$SLURM_PARTITION" --time="12:00:00" --dependency=afterok:\$JOB07 --output="\$BASE/logs/consensus_%j.out" --error="\$BASE/logs/consensus_%j.err" --parsable "\$SCRIPT_DIR/08_build_consensus_peaks.sbatch" "$CONFIG")
JOB09=\$(sbatch --partition="$SLURM_PARTITION" --time="12:00:00" --dependency=afterok:\$JOB08 --output="\$BASE/logs/counts_%j.out" --error="\$BASE/logs/counts_%j.err" --parsable "\$SCRIPT_DIR/09_count_fragments.sbatch" "$CONFIG")
JOB10=\$(sbatch --partition="$SLURM_PARTITION" --time="12:00:00" --array="\$ARRAY_PASS" --dependency=afterok:\$JOB08 --output="\$BASE/logs/bigwig_%A_%a.out" --error="\$BASE/logs/bigwig_%A_%a.err" --parsable "\$SCRIPT_DIR/10_generate_bigwig.sbatch" "$CONFIG" || true)
echo "Downstream jobs submitted after QC: peaks=\$JOB07 consensus=\$JOB08 counts=\$JOB09 bigwig=\${JOB10:-not_submitted}"
LAUNCH
chmod +x "$LAUNCH_PEAKS"
JOB_LAUNCH=$(sbatch "${default_sbatch[@]}" --dependency=afterok:$JOB06 --output="$BASE/logs/launcher_%j.out" --error="$BASE/logs/launcher_%j.err" --parsable --wrap="bash '$LAUNCH_PEAKS'")

cat <<SUMMARY
Submitted job chain:
  trim:          $JOB02
  align:         $JOB03
  filter:        $JOB04
  dedup:         $JOB05
  QC/merge:      $JOB06
  downstream launcher after QC: $JOB_LAUNCH

Watch:
  squeue -u $USER
  tail -f $BASE/logs/*.out

Key outputs:
  QC table:          $BASE/qc/biosample_qc.tsv
  Passing samples:   $BASE/qc/pass_hard_qc_biosamples.txt
  Consensus peaks:   $BASE/consensus/master_consensus_hg19.calderon_exact.bed
  Counts:            $BASE/counts/
SUMMARY

#!/usr/bin/env bash
#SBATCH --job-name=sqtl_atac
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

source "${SQTL_CONFIG:?Set SQTL_CONFIG to your private copy of config/sqtl_config.example.sh}"

MANIFEST="$KERIMOV_MANIFESTS/sqtl_atac_overlap_manifest.tsv"

OUTPUT_ROOT="$KERIMOV_PROCESSED/sqtl_atac_overlap"
HIGH_OUTPUT_DIR="$OUTPUT_ROOT/high_PIP"
LOW_OUTPUT_DIR="$OUTPUT_ROOT/low_PIP"
SUMMARY_DIR="$KERIMOV_RESULTS/sqtl_atac_overlap_per_dataset"

mkdir -p \
  "$HIGH_OUTPUT_DIR" \
  "$LOW_OUTPUT_DIR" \
  "$SUMMARY_DIR" \
  "$KERIMOV_TMP"

if ! command -v bedtools >/dev/null 2>&1
then
    module load bedtools >/dev/null 2>&1 || true
fi

if ! command -v bedtools >/dev/null 2>&1
then
    echo "ERROR: bedtools is not available." >&2
    echo "Run: module spider bedtools" >&2
    exit 1
fi

task_row=$((SLURM_ARRAY_TASK_ID + 1))

task_line=$(
  awk -F'\t' -v target="$task_row" '
    BEGIN {
      OFS = "\t"
    }

    NR == target {
      print $1, $2, $3, $4, $5, $6, $7, $8, $9
      exit
    }
  ' "$MANIFEST"
)

if [[ -z "$task_line" ]]
then
    echo "ERROR: No manifest row for array task $SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi

IFS=$'\t' read -r \
  dataset_id \
  calderon_celltype \
  match_quality \
  include_strict \
  include_expanded \
  include_sensitivity \
  high_bed \
  low_bed \
  peak_bed \
  <<< "$task_line"

echo "Dataset:            $dataset_id"
echo "Calderon cell type: $calderon_celltype"
echo "Match quality:      $match_quality"
echo "Peak BED:           $peak_bed"
echo "Bedtools:           $(bedtools --version)"

for required_file in "$high_bed" "$low_bed" "$peak_bed"
do
    if [[ ! -f "$required_file" ]]
    then
        echo "ERROR: Missing file: $required_file" >&2
        exit 1
    fi
done

TASK_TMP="$KERIMOV_TMP/sqtl_atac_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

rm -rf "$TASK_TMP"
mkdir -p "$TASK_TMP"

trap 'rm -rf "$TASK_TMP"' EXIT

HIGH_INPUT="$TASK_TMP/${dataset_id}.high.bed"
LOW_INPUT="$TASK_TMP/${dataset_id}.low.bed"

gzip -cd "$high_bed" > "$HIGH_INPUT"
gzip -cd "$low_bed" > "$LOW_INPUT"

validate_bed() {
    local input_file="$1"
    local expected_columns=14

    if [[ ! -s "$input_file" ]]
    then
        return 0
    fi

    awk -v expected="$expected_columns" '
      NF != expected {
        print "Invalid BED column count on line " NR \
              ": observed " NF ", expected " expected > "/dev/stderr"
        bad = 1
        exit
      }

      $1 !~ /^chr/ {
        print "Invalid chromosome on line " NR > "/dev/stderr"
        bad = 1
        exit
      }

      $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ {
        print "Invalid coordinates on line " NR > "/dev/stderr"
        bad = 1
        exit
      }

      $3 <= $2 {
        print "BED end is not greater than start on line " NR > "/dev/stderr"
        bad = 1
        exit
      }

      END {
        exit(bad == 1)
      }
    ' "$input_file"
}

validate_bed "$HIGH_INPUT"
validate_bed "$LOW_INPUT"

HEADER=$'chromosome\tstart\tend\tname\tscore\tstrand\tdataset_id\tcalderon_celltype\tmax_pip\tpip_class\tgene_id\tmolecular_trait_id\tvariant\trsid\tn_peak_overlaps\taccessibility_status'

annotate_bed() {
    local input_file="$1"
    local output_file="$2"
    local temporary_counts="$3"

    if [[ ! -s "$input_file" ]]
    then
        printf "%s\n" "$HEADER" |
        gzip -c > "$output_file"

        printf "0\t0\t0\tNA\n"
        return
    fi

    bedtools intersect \
      -a "$input_file" \
      -b "$peak_bed" \
      -c \
      > "$temporary_counts"

    {
      printf "%s\n" "$HEADER"

      awk '
        BEGIN {
          OFS = "\t"
        }

        {
          status = ($15 > 0) ? "accessible" : "not_accessible"
          print $0, status
        }
      ' "$temporary_counts"
    } |
    gzip -c > "$output_file"

    total=$(
      awk 'END {print NR + 0}' "$temporary_counts"
    )

    accessible=$(
      awk '
        $15 > 0 {
          n++
        }

        END {
          print n + 0
        }
      ' "$temporary_counts"
    )

    inaccessible=$((total - accessible))

    proportion=$(
      awk -v accessible="$accessible" -v total="$total" '
        BEGIN {
          if (total == 0) {
            print "NA"
          } else {
            printf "%.10f", accessible / total
          }
        }
      '
    )

    printf "%s\t%s\t%s\t%s\n" \
      "$total" \
      "$accessible" \
      "$inaccessible" \
      "$proportion"
}

HIGH_OUTPUT="$HIGH_OUTPUT_DIR/${dataset_id}.high_PIP.atac_annotated.tsv.gz"
LOW_OUTPUT="$LOW_OUTPUT_DIR/${dataset_id}.low_PIP.atac_annotated.tsv.gz"

HIGH_COUNTS="$TASK_TMP/${dataset_id}.high.counts.tsv"
LOW_COUNTS="$TASK_TMP/${dataset_id}.low.counts.tsv"

high_stats=$(
  annotate_bed \
    "$HIGH_INPUT" \
    "$HIGH_OUTPUT" \
    "$HIGH_COUNTS"
)

low_stats=$(
  annotate_bed \
    "$LOW_INPUT" \
    "$LOW_OUTPUT" \
    "$LOW_COUNTS"
)

IFS=$'\t' read -r \
  n_high_total \
  n_high_accessible \
  n_high_inaccessible \
  prop_high_accessible \
  <<< "$high_stats"

IFS=$'\t' read -r \
  n_low_total \
  n_low_accessible \
  n_low_inaccessible \
  prop_low_accessible \
  <<< "$low_stats"

SUMMARY_FILE="$SUMMARY_DIR/${dataset_id}.sqtl_atac_summary.tsv"

printf \
"dataset_id\tcalderon_celltype\tmatch_quality\tinclude_strict_primary\tinclude_expanded_primary\tinclude_sensitivity\tn_high_total\tn_high_accessible\tn_high_not_accessible\tprop_high_accessible\tn_low_total\tn_low_accessible\tn_low_not_accessible\tprop_low_accessible\thigh_annotated_file\tlow_annotated_file\tcalderon_peak_file\n" \
> "$SUMMARY_FILE"

printf \
"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  "$dataset_id" \
  "$calderon_celltype" \
  "$match_quality" \
  "$include_strict" \
  "$include_expanded" \
  "$include_sensitivity" \
  "$n_high_total" \
  "$n_high_accessible" \
  "$n_high_inaccessible" \
  "$prop_high_accessible" \
  "$n_low_total" \
  "$n_low_accessible" \
  "$n_low_inaccessible" \
  "$prop_low_accessible" \
  "$HIGH_OUTPUT" \
  "$LOW_OUTPUT" \
  "$peak_bed" \
  >> "$SUMMARY_FILE"

echo
echo "Overlap completed."
echo "High PIP: $n_high_accessible / $n_high_total accessible"
echo "Low PIP:  $n_low_accessible / $n_low_total accessible"
echo "Summary:  $SUMMARY_FILE"

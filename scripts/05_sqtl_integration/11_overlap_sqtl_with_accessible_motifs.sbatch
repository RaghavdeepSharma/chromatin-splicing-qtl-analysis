#!/usr/bin/env bash
#SBATCH --job-name=sqtl_motif
#SBATCH --partition=short
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

source "${SQTL_CONFIG:?Set SQTL_CONFIG to your private copy of config/sqtl_config.example.sh}"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

MOTIF_MANIFEST="$KERIMOV_MANIFESTS/motif_component_manifest.tsv"

CELLTYPE_MANIFEST="$KERIMOV_MANIFESTS/motif_sqtl_celltype_manifest.tsv"

SQTL_MANIFEST="$KERIMOV_MANIFESTS/motif_sqtl_runtime_manifest.tsv"
GENE_LOOKUP="$KERIMOV_METADATA/gencode_v47_gene_lookup.tsv"

OUTPUT_DIR="$KERIMOV_RESULTS/motif_component_overlap_counts"

mkdir -p \
  "$OUTPUT_DIR" \
  "$KERIMOV_TMP"

if ! command -v bedtools >/dev/null 2>&1; then
    module load bedtools >/dev/null 2>&1 || true
fi

if ! command -v bedtools >/dev/null 2>&1; then
    echo "ERROR: bedtools is not available." >&2
    exit 1
fi

for required_file in \
    "$MOTIF_MANIFEST" \
    "$CELLTYPE_MANIFEST" \
    "$SQTL_MANIFEST"
do
    if [[ ! -f "$required_file" ]]; then
        echo "ERROR: Missing required manifest:" >&2
        echo "$required_file" >&2
        exit 1
    fi
done

task_record=$(
    awk \
        -F'\t' \
        -v task_id="$SLURM_ARRAY_TASK_ID" '
        NR > 1 && $1 == task_id {
            printf "%s\034%s\034%s\034%s\034%s\034%s\n", \
                $1, \
                $2, \
                $3, \
                $4, \
                $5, \
                $6

            exit
        }
    ' "$MOTIF_MANIFEST"
)

if [[ -z "$task_record" ]]; then
    echo \
"ERROR: No motif-manifest row for task $SLURM_ARRAY_TASK_ID." \
        >&2

    exit 1
fi

IFS=$'\034' read -r \
    task_id \
    component_name \
    candidate_archetype \
    is_split_component \
    chromosome_component \
    motif_file \
    <<< "$task_record"

if [[ ! -f "$motif_file" ]]; then
    echo "ERROR: Motif file not found:" >&2
    echo "$motif_file" >&2
    exit 1
fi

TASK_TMP="$KERIMOV_TMP/motif_overlap_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

MOTIF_BED_DIR="$TASK_TMP/motif_by_celltype"

mkdir -p "$MOTIF_BED_DIR"

trap 'rm -rf "$TASK_TMP"' EXIT

ACCESS_SPEC=$(
    awk -F'\t' '
        NR > 1 {
            print $1 "=" $3
        }
    ' "$CELLTYPE_MANIFEST" |
    paste -sd';' -
)

if [[ -z "$ACCESS_SPEC" ]]; then
    echo "ERROR: Could not build accessibility-column mapping." >&2
    exit 1
fi

while IFS= read -r celltype
do
    : > "$MOTIF_BED_DIR/${celltype}.motif_gene.bed"
done < <(
    tail -n +2 "$CELLTYPE_MANIFEST" |
    cut -f1
)

echo "Task ID:              $task_id"
echo "Component:            $component_name"
echo "Candidate archetype:  $candidate_archetype"
echo "Split component:      $is_split_component"
echo "Chromosome component: $chromosome_component"
echo "Motif file:           $motif_file"
echo "Bedtools:             $(bedtools --version)"

# One pass through the processed motif file.
#
# Required motif properties:
#   column 15: intragenic
#   column 16: same-strand gene overlap
#   column 19: individual gene types
#   columns 27–51: matched-cell accessibility
#
# gene_id and gene_type can contain parallel semicolon-separated lists.
# Each qualifying protein-coding or lncRNA gene is expanded into its
# own BED row.

gzip -cd "$motif_file" |
awk \
    -F'\t' \
    -v OFS='\t' \
    -v access_spec="$ACCESS_SPEC" \
    -v output_dir="$MOTIF_BED_DIR" \
    -v archetype="$candidate_archetype" \
    -f "$SCRIPT_DIR/10_extract_accessible_motif_genes.awk" "$GENE_LOOKUP" -

OUTPUT_TMP="$TASK_TMP/component_overlap_counts.tsv"

printf \
"task_id\tcomponent_name\tcandidate_archetype\tis_split_component\tchromosome_component\tanalysis_set\tcalderon_celltype\tn_accessible_intragenic_motif_gene_rows\tn_high_total\tn_high_overlap\tn_high_not_overlap\tn_low_total\tn_low_overlap\tn_low_not_overlap\tmotif_file\tsqtl_bed\n" \
> "$OUTPUT_TMP"

while IFS=$'\034' read -r \
    analysis_set \
    celltype \
    n_high_total \
    n_low_total \
    sqtl_bed
do
    motif_bed="$MOTIF_BED_DIR/${celltype}.motif_gene.bed"

    if [[ ! -f "$sqtl_bed" ]]; then
        echo "ERROR: Missing runtime sQTL BED:" >&2
        echo "$sqtl_bed" >&2
        exit 1
    fi

    n_motif_gene_rows=$(
        awk 'END {print NR + 0}' "$motif_bed"
    )

    if [[ "$n_motif_gene_rows" -eq 0 ]]; then

        n_high_overlap=0
        n_low_overlap=0

    else

        read -r n_high_overlap n_low_overlap < <(
            bedtools intersect \
                -a "$sqtl_bed" \
                -b "$motif_bed" \
                -wa \
                -wb |
            awk -F'\t' '
                # A has 14 columns.
                # B gene_id is its fifth column, therefore column 19
                # in the combined bedtools output.
                $11 == $19 {
                    variant_gene_key = $4

                    if (!(variant_gene_key in seen)) {
                        seen[variant_gene_key] = 1

                        if ($10 == "high_PIP") {
                            high++
                        } else if ($10 == "low_PIP") {
                            low++
                        }
                    }
                }

                END {
                    print high + 0, low + 0
                }
            '
        )
    fi

    n_high_not_overlap=$(
        awk \
            -v total="$n_high_total" \
            -v overlap="$n_high_overlap" \
            'BEGIN {print total - overlap}'
    )

    n_low_not_overlap=$(
        awk \
            -v total="$n_low_total" \
            -v overlap="$n_low_overlap" \
            'BEGIN {print total - overlap}'
    )

    if [[ "$n_high_not_overlap" -lt 0 ]]; then
        echo \
"ERROR: High-PIP overlap exceeds total for $analysis_set / $celltype." \
            >&2

        exit 1
    fi

    if [[ "$n_low_not_overlap" -lt 0 ]]; then
        echo \
"ERROR: Low-PIP overlap exceeds total for $analysis_set / $celltype." \
            >&2

        exit 1
    fi

    printf \
"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$task_id" \
        "$component_name" \
        "$candidate_archetype" \
        "$is_split_component" \
        "$chromosome_component" \
        "$analysis_set" \
        "$celltype" \
        "$n_motif_gene_rows" \
        "$n_high_total" \
        "$n_high_overlap" \
        "$n_high_not_overlap" \
        "$n_low_total" \
        "$n_low_overlap" \
        "$n_low_not_overlap" \
        "$motif_file" \
        "$sqtl_bed" \
        >> "$OUTPUT_TMP"

done < <(
    awk -F'\t' '
        NR > 1 {
            printf "%s\034%s\034%s\034%s\034%s\n", \
                $1, \
                $2, \
                $3, \
                $4, \
                $6
        }
    ' "$SQTL_MANIFEST"
)

n_result_rows=$(
    tail -n +2 "$OUTPUT_TMP" |
    wc -l
)

if [[ "$n_result_rows" -ne 29 ]]; then
    echo \
"ERROR: Expected 29 result rows; found $n_result_rows." \
        >&2

    exit 1
fi

FINAL_OUTPUT="$OUTPUT_DIR/task_${task_id}.${component_name}.motif_overlap_counts.tsv"

mv "$OUTPUT_TMP" "$FINAL_OUTPUT"

echo
echo "Motif-component overlap completed."
echo "Result rows: $n_result_rows"
echo "Output: $FINAL_OUTPUT"

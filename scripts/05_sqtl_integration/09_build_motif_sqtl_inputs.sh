#!/usr/bin/env bash
set -euo pipefail

source "${SQTL_CONFIG:?Set SQTL_CONFIG to your private copy of config/sqtl_config.example.sh}"

SQTL_MANIFEST="$KERIMOV_PROCESSED/gene_filtered_sqtl_bed_manifest.tsv"

MOTIF_DIR="${PROJECT_ROOT:?Set PROJECT_ROOT}/tfbs_celltype_accessibility_20260701/BAMPE_raw_processed_accessibility_q05_significant_peaks/annotated_with_accessibility"

REFERENCE_MOTIF="$MOTIF_DIR/AHR.general_annotated.celltype_accessibility.tsv.gz"

CELLTYPE_BED_DIR="$KERIMOV_PROCESSED/gene_filtered_sqtl_beds/by_celltype"

CELLTYPE_MANIFEST="$KERIMOV_MANIFESTS/motif_sqtl_celltype_manifest.tsv"

MOTIF_MANIFEST="$KERIMOV_MANIFESTS/motif_component_manifest.tsv"

SUMMARY_FILE="$KERIMOV_RESULTS/motif_sqtl_input_preparation_summary.tsv"

TASK_TMP="$KERIMOV_TMP/motif_sqtl_input_preparation_$$"

mkdir -p \
    "$CELLTYPE_BED_DIR" \
    "$KERIMOV_MANIFESTS" \
    "$KERIMOV_RESULTS" \
    "$TASK_TMP"

trap 'rm -rf "$TASK_TMP"' EXIT

if [[ ! -f "$SQTL_MANIFEST" ]]; then
    echo "ERROR: Missing sQTL manifest:" >&2
    echo "$SQTL_MANIFEST" >&2
    exit 1
fi

if [[ ! -d "$MOTIF_DIR" ]]; then
    echo "ERROR: Missing approved motif directory:" >&2
    echo "$MOTIF_DIR" >&2
    exit 1
fi

if [[ ! -f "$REFERENCE_MOTIF" ]]; then
    echo "ERROR: Missing reference motif file:" >&2
    echo "$REFERENCE_MOTIF" >&2
    exit 1
fi

get_manifest_column() {
    local column_name="$1"

    awk -F'\t' -v target="$column_name" '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == target) {
                    print i
                    exit
                }
            }
        }
    ' "$SQTL_MANIFEST"
}

DATASET_COLUMN=$(
    get_manifest_column "dataset_id"
)

CELLTYPE_COLUMN=$(
    get_manifest_column "calderon_celltype"
)

HIGH_BED_COLUMN=$(
    get_manifest_column "high_pip_bed"
)

LOW_BED_COLUMN=$(
    get_manifest_column "low_pip_bed"
)

for required_column in \
    "$DATASET_COLUMN" \
    "$CELLTYPE_COLUMN" \
    "$HIGH_BED_COLUMN" \
    "$LOW_BED_COLUMN"
do
    if [[ -z "$required_column" ]]; then
        echo "ERROR: Required manifest column was not found." >&2
        exit 1
    fi
done

REFERENCE_HEADER=$(
    gzip -cd "$REFERENCE_MOTIF" 2>/dev/null |
    head -n 1 ||
    true
)

if [[ -z "$REFERENCE_HEADER" ]]; then
    echo "ERROR: Could not read reference motif header." >&2
    exit 1
fi

REFERENCE_NCOLS=$(
    awk -F'\t' '{print NF}' <<< "$REFERENCE_HEADER"
)

if [[ "$REFERENCE_NCOLS" -ne 51 ]]; then
    echo "ERROR: Expected 51 motif columns; found $REFERENCE_NCOLS." >&2
    exit 1
fi

printf \
"calderon_celltype\taccessibility_column\taccessibility_column_index\tn_datasets\tdataset_ids\tn_variant_rows\tn_high_pip_rows\tn_low_pip_rows\tcombined_variant_bed\n" \
> "$CELLTYPE_MANIFEST"

mapfile -t CELLTYPES < <(
    awk -F'\t' -v celltype_column="$CELLTYPE_COLUMN" '
        NR > 1 && $celltype_column != "" {
            print $celltype_column
        }
    ' "$SQTL_MANIFEST" |
    sort -u
)

if [[ "${#CELLTYPES[@]}" -eq 0 ]]; then
    echo "ERROR: No Calderon cell types found." >&2
    exit 1
fi

echo
echo "Building combined sQTL BEDs for ${#CELLTYPES[@]} cell types..."
echo

for celltype in "${CELLTYPES[@]}"
do
    accessibility_column="accessible_${celltype}"

    accessibility_index=$(
        awk -F'\t' -v target="$accessibility_column" '
            {
                for (i = 1; i <= NF; i++) {
                    if ($i == target) {
                        print i
                        exit
                    }
                }
            }
        ' <<< "$REFERENCE_HEADER"
    )

    if [[ -z "$accessibility_index" ]]; then
        echo "ERROR: Accessibility column not found:" >&2
        echo "$accessibility_column" >&2
        exit 1
    fi

    input_list="$TASK_TMP/${celltype}.input_files.tsv"
    unsorted_bed="$TASK_TMP/${celltype}.unsorted.bed"
    sorted_bed="$TASK_TMP/${celltype}.sorted.bed"

    : > "$unsorted_bed"

    awk \
        -F'\t' \
        -v OFS='\t' \
        -v target_celltype="$celltype" \
        -v dataset_column="$DATASET_COLUMN" \
        -v celltype_column="$CELLTYPE_COLUMN" \
        -v high_column="$HIGH_BED_COLUMN" \
        -v low_column="$LOW_BED_COLUMN" '
            NR > 1 &&
            $celltype_column == target_celltype {
                print \
                    $dataset_column, \
                    $high_column, \
                    $low_column
            }
        ' "$SQTL_MANIFEST" \
        > "$input_list"

    n_datasets=$(
        awk 'END {print NR + 0}' "$input_list"
    )

    if [[ "$n_datasets" -eq 0 ]]; then
        echo "ERROR: No datasets found for $celltype." >&2
        exit 1
    fi

    dataset_ids=$(
        cut -f1 "$input_list" |
        paste -sd, -
    )

    while IFS=$'\t' read -r dataset_id high_bed low_bed
    do
        for input_bed in "$high_bed" "$low_bed"
        do
            if [[ ! -f "$input_bed" ]]; then
                echo "ERROR: Missing filtered BED:" >&2
                echo "$input_bed" >&2
                exit 1
            fi

            gzip -cd "$input_bed" >> "$unsorted_bed"
        done
    done < "$input_list"

    awk \
        -F'\t' \
        -v expected_celltype="$celltype" '
            NF != 14 {
                print \
                    "ERROR: Invalid column count on line", \
                    NR, \
                    "observed", \
                    NF, \
                    "expected 14" \
                    > "/dev/stderr"

                bad = 1
                exit
            }

            $8 != expected_celltype {
                print \
                    "ERROR: Cell-type mismatch on line", \
                    NR, \
                    "observed", \
                    $8, \
                    "expected", \
                    expected_celltype \
                    > "/dev/stderr"

                bad = 1
                exit
            }

            $10 != "high_PIP" &&
            $10 != "low_PIP" {
                print \
                    "ERROR: Invalid PIP class on line", \
                    NR, \
                    $10 \
                    > "/dev/stderr"

                bad = 1
                exit
            }

            END {
                exit(bad == 1)
            }
        ' "$unsorted_bed"

    LC_ALL=C sort \
        -T "$TASK_TMP" \
        -S 2G \
        -k1,1V \
        -k2,2n \
        -k3,3n \
        -k7,7 \
        "$unsorted_bed" \
        > "$sorted_bed"

    output_bed="$CELLTYPE_BED_DIR/${celltype}.protein_coding_lncRNA.all_PIP.hg38.bed.gz"

    gzip -c "$sorted_bed" > "$output_bed"

    read -r n_rows n_high n_low < <(
        awk -F'\t' '
            {
                total++

                if ($10 == "high_PIP") {
                    high++
                }

                if ($10 == "low_PIP") {
                    low++
                }
            }

            END {
                print \
                    total + 0, \
                    high + 0, \
                    low + 0
            }
        ' "$sorted_bed"
    )

    if [[ $((n_high + n_low)) -ne "$n_rows" ]]; then
        echo "ERROR: High/low counts do not sum for $celltype." >&2
        exit 1
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$celltype" \
        "$accessibility_column" \
        "$accessibility_index" \
        "$n_datasets" \
        "$dataset_ids" \
        "$n_rows" \
        "$n_high" \
        "$n_low" \
        "$output_bed" \
        >> "$CELLTYPE_MANIFEST"

    echo \
"$celltype: $n_datasets datasets; $n_rows variants; $n_high high PIP; $n_low low PIP"
done

printf \
"task_id\tcomponent_name\tcandidate_archetype\tis_split_component\tchromosome_component\tmotif_file\theader_ncols\n" \
> "$MOTIF_MANIFEST"

task_id=0

while IFS= read -r -d '' motif_file
do
    filename=$(
        basename "$motif_file"
    )

    component_name="$filename"
    component_name="${component_name%.general_annotated.celltype_accessibility.tsv.gz}"

    candidate_archetype="$component_name"
    is_split_component="NO"
    chromosome_component=""

    if [[ "$component_name" =~ ^(.+)\.(chr[^.]*)$ ]]; then
        candidate_archetype="${BASH_REMATCH[1]}"
        chromosome_component="${BASH_REMATCH[2]}"
        is_split_component="YES"
    fi

    file_header=$(
        gzip -cd "$motif_file" 2>/dev/null |
        head -n 1 ||
        true
    )

    if [[ -z "$file_header" ]]; then
        echo "ERROR: Empty motif file:" >&2
        echo "$motif_file" >&2
        exit 1
    fi

    file_ncols=$(
        awk -F'\t' '{print NF}' <<< "$file_header"
    )

    if [[ "$file_ncols" -ne 51 ]]; then
        echo "ERROR: Unexpected motif schema in:" >&2
        echo "$motif_file" >&2
        echo "Observed columns: $file_ncols" >&2
        exit 1
    fi

    if [[ "$file_header" != "$REFERENCE_HEADER" ]]; then
        echo "ERROR: Header differs from AHR reference:" >&2
        echo "$motif_file" >&2
        exit 1
    fi

    task_id=$((task_id + 1))

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$task_id" \
        "$component_name" \
        "$candidate_archetype" \
        "$is_split_component" \
        "$chromosome_component" \
        "$motif_file" \
        "$file_ncols" \
        >> "$MOTIF_MANIFEST"

done < <(
    find "$MOTIF_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*.general_annotated.celltype_accessibility.tsv.gz" \
        -print0 |
    sort -z
)

n_celltypes=$(
    tail -n +2 "$CELLTYPE_MANIFEST" |
    wc -l
)

n_dataset_assignments=$(
    awk -F'\t' '
        NR > 1 {
            total += $4
        }

        END {
            print total + 0
        }
    ' "$CELLTYPE_MANIFEST"
)

total_variant_rows=$(
    awk -F'\t' '
        NR > 1 {
            total += $6
        }

        END {
            print total + 0
        }
    ' "$CELLTYPE_MANIFEST"
)

total_high_rows=$(
    awk -F'\t' '
        NR > 1 {
            total += $7
        }

        END {
            print total + 0
        }
    ' "$CELLTYPE_MANIFEST"
)

total_low_rows=$(
    awk -F'\t' '
        NR > 1 {
            total += $8
        }

        END {
            print total + 0
        }
    ' "$CELLTYPE_MANIFEST"
)

n_motif_components=$(
    tail -n +2 "$MOTIF_MANIFEST" |
    wc -l
)

n_candidate_archetypes=$(
    tail -n +2 "$MOTIF_MANIFEST" |
    cut -f3 |
    sort -u |
    wc -l
)

n_split_components=$(
    awk -F'\t' '
        NR > 1 && $4 == "YES" {
            n++
        }

        END {
            print n + 0
        }
    ' "$MOTIF_MANIFEST"
)

if [[ "$n_dataset_assignments" -ne 31 ]]; then
    echo "ERROR: Expected 31 dataset assignments; found $n_dataset_assignments." >&2
    exit 1
fi

if [[ "$total_variant_rows" -ne 522903 ]]; then
    echo "ERROR: Expected 522903 filtered rows; found $total_variant_rows." >&2
    exit 1
fi

if [[ "$total_high_rows" -ne 4102 ]]; then
    echo "ERROR: Expected 4102 high-PIP rows; found $total_high_rows." >&2
    exit 1
fi

if [[ "$total_low_rows" -ne 518801 ]]; then
    echo "ERROR: Expected 518801 low-PIP rows; found $total_low_rows." >&2
    exit 1
fi

printf "metric\tvalue\n" > "$SUMMARY_FILE"

printf "mapped_celltypes\t%s\n" \
    "$n_celltypes" \
    >> "$SUMMARY_FILE"

printf "dataset_assignments\t%s\n" \
    "$n_dataset_assignments" \
    >> "$SUMMARY_FILE"

printf "combined_variant_rows\t%s\n" \
    "$total_variant_rows" \
    >> "$SUMMARY_FILE"

printf "combined_high_pip_rows\t%s\n" \
    "$total_high_rows" \
    >> "$SUMMARY_FILE"

printf "combined_low_pip_rows\t%s\n" \
    "$total_low_rows" \
    >> "$SUMMARY_FILE"

printf "motif_component_files\t%s\n" \
    "$n_motif_components" \
    >> "$SUMMARY_FILE"

printf "candidate_motif_archetypes\t%s\n" \
    "$n_candidate_archetypes" \
    >> "$SUMMARY_FILE"

printf "chromosome_split_components\t%s\n" \
    "$n_split_components" \
    >> "$SUMMARY_FILE"

echo
echo "Motif–sQTL input preparation completed."
echo
column -ts $'\t' "$SUMMARY_FILE"

echo
echo "Cell-type manifest:"
echo "$CELLTYPE_MANIFEST"

echo
echo "Motif-component manifest:"
echo "$MOTIF_MANIFEST"

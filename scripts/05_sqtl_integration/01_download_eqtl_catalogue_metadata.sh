#!/usr/bin/env bash
set -euo pipefail

source "${SQTL_CONFIG:?Set SQTL_CONFIG to your private copy of config/sqtl_config.example.sh}"

mkdir -p "$KERIMOV_METADATA" "$KERIMOV_LOGS"

download_file() {
    local url="$1"
    local output="$2"
    local temporary="${output}.tmp"

    echo
    echo "Downloading:"
    echo "  $url"
    echo "Saving as:"
    echo "  $output"

    rm -f "$temporary"

    curl -L --fail --show-error \
        --retry 3 \
        --retry-delay 5 \
        --connect-timeout 30 \
        -o "$temporary" \
        "$url"

    if [[ ! -s "$temporary" ]]; then
        echo "ERROR: Downloaded file is empty: $temporary" >&2
        exit 1
    fi

    if head -c 200 "$temporary" | grep -Eqi '<html|<!DOCTYPE'; then
        echo "ERROR: Download returned an HTML page: $url" >&2
        exit 1
    fi

    mv "$temporary" "$output"
}

echo "Downloading eQTL Catalogue Release 7 metadata..."

download_file \
"https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/data_tables/dataset_metadata_r7.tsv" \
"$KERIMOV_METADATA/dataset_metadata_r7.tsv"

download_file \
"https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv" \
"$KERIMOV_METADATA/tabix_ftp_paths.tsv"

download_file \
"https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/Columns.md" \
"$KERIMOV_METADATA/Columns.md"

# Stable filename for later scripts
ln -sfn dataset_metadata_r7.tsv \
    "$KERIMOV_METADATA/dataset_metadata.tsv"

echo
echo "Downloaded files:"
ls -lh "$KERIMOV_METADATA"

echo
echo "Line counts:"
wc -l \
    "$KERIMOV_METADATA/dataset_metadata_r7.tsv" \
    "$KERIMOV_METADATA/tabix_ftp_paths.tsv" \
    "$KERIMOV_METADATA/Columns.md"

echo
echo "Metadata columns:"
head -n 1 "$KERIMOV_METADATA/dataset_metadata_r7.tsv" |
    tr '\t' '\n' |
    nl -ba

echo
echo "First five records:"
head -n 6 "$KERIMOV_METADATA/dataset_metadata_r7.tsv"

echo
echo "Immune-related metadata preview:"
grep -Ei \
'monocyte|macrophage|neutrophil|lymphocyte|CD4|CD8|T.cell|B.cell|NK.cell|Treg|regulatory|blood' \
"$KERIMOV_METADATA/dataset_metadata_r7.tsv" |
head -n 30 || true

echo
echo "Metadata download completed successfully."

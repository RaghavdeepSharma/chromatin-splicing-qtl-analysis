BEGIN {
    OFS = "\t"

    n_celltypes = split(access_spec, access_entries, ";")

    for (i = 1; i <= n_celltypes; i++) {
        split(access_entries[i], pieces, "=")

        celltype[i] = pieces[1]
        access_column[i] = pieces[2] + 0
        output_file[i] = output_dir "/" celltype[i] ".motif_gene.bed"
    }
}

# First input file: GENCODE lookup.
FNR == NR {

    if (FNR == 1) {
        next
    }

    gene_id_lookup = $2
    gene_type_lookup[gene_id_lookup] = $4

    next
}

# Second input: motif file from standard input.
FNR == 1 {
    next
}

$15 == "intragenic" && $16 == "TRUE" && $20 == "TRUE" {

    n_gene_ids = split($17, gene_ids, ";")

    for (gene_index = 1; gene_index <= n_gene_ids; gene_index++) {

        gene_id = gene_ids[gene_index]

        gsub(/^[[:space:]]+|[[:space:]]+$/, "", gene_id)
        gsub(/"/, "", gene_id)

        sub(/\.[0-9]+$/, "", gene_id)

        if (gene_id == "") {
            continue
        }

        gene_type = gene_type_lookup[gene_id]

        if (gene_type != "protein_coding" && gene_type != "lncRNA") {
            continue
        }

        for (cell_index = 1; cell_index <= n_celltypes; cell_index++) {

            column_index = access_column[cell_index]

            if ($column_index == "1") {
                print \
                    $1, \
                    $2, \
                    $3, \
                    archetype, \
                    gene_id, \
                    $11 \
                    >> output_file[cell_index]
            }
        }
    }
}

END {
    for (i = 1; i <= n_celltypes; i++) {
        close(output_file[i])
    }
}

#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

base <- Sys.getenv("KERIMOV_BASE")
if (base == "") stop("KERIMOV_BASE is not set.")

src_file <- file.path(
    base,
    "08_dataset_pooling_audit/tables",
    "strict_primary_source_level_variant_gene_rows.tsv.gz"
)

outdir <- file.path(
    base,
    "09_dataset_specific_analysis/corrected_pool"
)

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

x <- fread(src_file)

groups <- list(
    Resting_Monocytes = c(
        "QTD000025",
        "QTD000413",
        "QTD000508"
    ),

    Naive_Tregs = c(
        "QTD000040",
        "QTD000473"
    )
)

summary_list <- list()

for(group_name in names(groups)) {

    ids <- groups[[group_name]]

    d <- x[dataset_id %in% ids]

    cat("\n====================================\n")
    cat(group_name, "\n")
    cat("====================================\n")

    cat("Source rows:", nrow(d), "\n")
    cat("Datasets:", paste(sort(unique(d$dataset_id)), collapse=", "), "\n")

    pooled <- d[
        ,
        .(
            n_sources = uniqueN(dataset_id),

            n_high_sources =
                uniqueN(
                    dataset_id[
                        source_pip_class=="high_PIP"
                    ]
                ),

            n_low_sources =
                uniqueN(
                    dataset_id[
                        source_pip_class=="low_PIP"
                    ]
                ),

            max_pip = max(max_pip),

            variant = first(variant),
            gene_id_unversioned =
                first(gene_id_unversioned)
        ),
        by=variant_gene_key
    ]

    pooled[
        ,
        consensus_status :=
            fifelse(
                n_high_sources > 0 &
                n_low_sources > 0,
                "discordant",

                fifelse(
                    n_high_sources > 0,
                    "high_PIP",
                    "low_PIP"
                )
            )
    ]

    keep <- pooled[
        consensus_status!="discordant"
    ]

    cat("Unique variant-gene pairs:", nrow(pooled), "\n")
    cat(
        "Repeated across sources:",
        pooled[n_sources > 1, .N],
        "\n"
    )
    cat(
        "Discordant High/Low:",
        pooled[consensus_status=="discordant", .N],
        "\n"
    )
    cat(
        "Retained:",
        nrow(keep),
        "\n"
    )
    cat(
        "Retained High:",
        keep[consensus_status=="high_PIP", .N],
        "\n"
    )
    cat(
        "Retained Low:",
        keep[consensus_status=="low_PIP", .N],
        "\n"
    )

    fwrite(
        pooled,
        file.path(
            outdir,
            paste0(
                group_name,
                ".candidate_pool_audit.tsv.gz"
            )
        ),
        sep="\t"
    )

    fwrite(
        keep,
        file.path(
            outdir,
            paste0(
                group_name,
                ".candidate_corrected_pool.tsv.gz"
            )
        ),
        sep="\t"
    )

    summary_list[[group_name]] <- data.table(
        group = group_name,
        source_rows = nrow(d),
        unique_variant_gene_pairs = nrow(pooled),
        repeated_pairs = pooled[n_sources > 1, .N],
        discordant_pairs =
            pooled[consensus_status=="discordant", .N],
        retained_pairs = nrow(keep),
        retained_high =
            keep[consensus_status=="high_PIP", .N],
        retained_low =
            keep[consensus_status=="low_PIP", .N]
    )
}

summary <- rbindlist(summary_list)

fwrite(
    summary,
    file.path(
        outdir,
        "candidate_corrected_pool_summary.tsv"
    ),
    sep="\t"
)

cat("\n\nFINAL SUMMARY\n")
print(summary)

#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

BASE <- Sys.getenv("KERIMOV_BASE")
if (BASE == "") stop("KERIMOV_BASE is not set.")

OLD_DIR <- file.path(
    BASE,
    "07_logistic_regression/05_firth_comparison/firth_model_results"
)

NEW_DIR <- file.path(
    BASE,
    "09_dataset_specific_analysis/corrected_firth/model_results"
)

OUT_DIR <- file.path(
    BASE,
    "09_dataset_specific_analysis/final_corrected_analysis"
)

dir.create(
    file.path(OUT_DIR, "tables"),
    recursive=TRUE,
    showWarnings=FALSE
)

# ============================================================
# Existing Firth results
# ============================================================

old_files <- list.files(
    OLD_DIR,
    pattern="\\.strict_primary\\.firth\\.tsv$",
    full.names=TRUE
)

old <- rbindlist(
    lapply(old_files, fread),
    fill=TRUE
)

cat("Original Firth rows:", nrow(old), "\n")

# These four cell types had only one source dataset,
# so there is no cross-dataset pooling problem to correct.

unchanged_cells <- c(
    "Follicular_T_Helper",
    "Memory_Tregs",
    "Naive_B",
    "Naive_CD8_T"
)

unchanged <- old[
    calderon_celltype %in% unchanged_cells
]

unchanged[
    ,
    `:=`(
        model_context=calderon_celltype,
        analysis_version="existing_single_source",
        pooling_strategy="single_source_no_pooling_correction"
    )
]

# ============================================================
# Corrected Firth results
# ============================================================

new_files <- list.files(
    NEW_DIR,
    pattern="\\.corrected_firth\\.tsv$",
    full.names=TRUE
)

cat("Corrected result files:", length(new_files), "\n")

if(length(new_files) != 286) {
    stop(
        "Expected 286 corrected Firth files, found ",
        length(new_files),
        ". Do not build final results until the array is complete."
    )
}

corrected <- rbindlist(
    lapply(new_files, fread),
    fill=TRUE
)

cat("Corrected Firth rows:", nrow(corrected), "\n")

if(nrow(corrected) != 572) {
    stop(
        "Expected 572 corrected model rows, found ",
        nrow(corrected)
    )
}

corrected[
    ,
    model_context := calderon_celltype
]

# Preserve the biological context while restoring
# the common six-cell-type naming for the final table.

corrected[
    calderon_celltype=="Resting_Monocytes",
    calderon_celltype := "Monocytes"
]

corrected[
    model_context=="Resting_Monocytes",
    `:=`(
        analysis_version="corrected_resting_pool",
        pooling_strategy=
            "resting_sources_only_consensus_discordant_removed"
    )
]

corrected[
    model_context=="Naive_Tregs",
    `:=`(
        analysis_version="corrected_consensus_pool",
        pooling_strategy=
            "compatible_sources_consensus_discordant_removed"
    )
]

# ============================================================
# Final 6-cell analysis
# ============================================================

final <- rbindlist(
    list(
        unchanged,
        corrected
    ),
    fill=TRUE
)

setorder(
    final,
    calderon_celltype,
    motif_archetype
)

cat("\nFinal rows:", nrow(final), "\n")
cat(
    "Final motifs:",
    uniqueN(final$motif_archetype),
    "\n"
)
cat(
    "Final cell types:",
    uniqueN(final$calderon_celltype),
    "\n"
)

# ============================================================
# BH FDR across all valid Firth models
# ============================================================

final[, global_fdr := NA_real_]

valid_global <- which(
    !is.na(final$p_value) &
    is.finite(final$p_value)
)

final[
    valid_global,
    global_fdr :=
        p.adjust(
            p_value,
            method="BH"
        )
]

# Cell-specific BH FDR

final[
    ,
    cell_fdr := {

        ans <- rep(
            NA_real_,
            .N
        )

        good <- which(
            !is.na(p_value) &
            is.finite(p_value)
        )

        if(length(good)>0) {
            ans[good] <-
                p.adjust(
                    p_value[good],
                    method="BH"
                )
        }

        ans
    },
    by=calderon_celltype
]

# ============================================================
# QC
# ============================================================

qc <- final[
    ,
    .(
        n_models=.N,
        n_valid_firth=sum(
            !is.na(p_value)
        ),
        n_global_fdr_05=sum(
            global_fdr < 0.05,
            na.rm=TRUE
        ),
        n_cell_fdr_05=sum(
            cell_fdr < 0.05,
            na.rm=TRUE
        ),
        n_high_overlap_ge_5=sum(
            n_high_overlap >= 5,
            na.rm=TRUE
        )
    ),
    by=calderon_celltype
]

fwrite(
    final,
    file.path(
        OUT_DIR,
        "tables",
        "final_corrected_firth_1716_models.tsv"
    ),
    sep="\t"
)

fwrite(
    qc,
    file.path(
        OUT_DIR,
        "tables",
        "final_corrected_firth_qc.tsv"
    ),
    sep="\t"
)

cat("\n====================================\n")
cat("FINAL QC\n")
cat("====================================\n")

print(qc)

if(nrow(final) != 1716) {
    stop(
        "Expected 1716 final models; found ",
        nrow(final)
    )
}

if(uniqueN(final$motif_archetype) != 286) {
    stop("Expected 286 motifs.")
}

if(uniqueN(final$calderon_celltype) != 6) {
    stop("Expected six final cell types.")
}

cat("\nFINAL 1716-MODEL STRUCTURE PASSED.\n")

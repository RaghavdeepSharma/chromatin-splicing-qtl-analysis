#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

BASE <- Sys.getenv("KERIMOV_BASE")
if (BASE == "") stop("KERIMOV_BASE is not set.")

INFILE <- file.path(
    BASE,
    "09_dataset_specific_analysis/final_corrected_analysis/tables",
    "final_corrected_firth_1716_models.tsv"
)

OUTDIR <- file.path(
    BASE,
    "09_dataset_specific_analysis/final_corrected_analysis/final_result_figures"
)

dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

d <- fread(INFILE)

cat("Rows read:", nrow(d), "\n")

# ============================================================
# Select main results for presentation
#
# Include:
#   - global FDR hits
#   - cell-specific FDR hits
#   - SMAD_1 resting-Monocyte candidate
# ============================================================

key <- d[
    (!is.na(global_fdr) & global_fdr < 0.05) |
    (!is.na(cell_fdr) & cell_fdr < 0.05) |
    (calderon_celltype=="Monocytes" &
     motif_archetype=="SMAD_1")
]

key[
    calderon_celltype=="Monocytes",
    display_cell := "Resting Monocytes"
]

key[
    calderon_celltype!="Monocytes",
    display_cell := calderon_celltype
]

key[
    ,
    candidate := paste(
        motif_archetype,
        display_cell,
        sep=" | "
    )
]

# ============================================================
# FIGURE 1
# HIGH-PIP vs LOW-PIP MOTIF OVERLAP RATE
# ============================================================

rates <- copy(key)

rates[
    ,
    `:=`(
        high_overlap_pct =
            100 * n_high_overlap / n_high,

        low_overlap_pct =
            100 * n_low_overlap / n_low
    )
]

rates_long <- melt(
    rates,
    id.vars=c(
        "candidate",
        "motif_archetype",
        "display_cell",
        "n_high",
        "n_low",
        "n_high_overlap",
        "n_low_overlap"
    ),
    measure.vars=c(
        "high_overlap_pct",
        "low_overlap_pct"
    ),
    variable.name="pip_class",
    value.name="overlap_percent"
)

rates_long[
    ,
    pip_class :=
        fifelse(
            pip_class=="high_overlap_pct",
            "High-PIP",
            "Low-PIP"
        )
]

# Correct numerator/denominator label for each bar
rates_long[
    pip_class=="High-PIP",
    support_label :=
        paste0(
            n_high_overlap,
            "/",
            n_high
        )
]

rates_long[
    pip_class=="Low-PIP",
    support_label :=
        paste0(
            n_low_overlap,
            "/",
            n_low
        )
]

# Preserve useful ordering
rates_long[
    ,
    candidate := factor(
        candidate,
        levels=rev(unique(key$candidate))
    )
]

p1 <- ggplot(
    rates_long,
    aes(
        x=candidate,
        y=overlap_percent,
        fill=pip_class
    )
) +
    geom_col(
        position=position_dodge(width=0.8),
        width=0.7
    ) +
    geom_text(
        aes(label=support_label),
        position=position_dodge(width=0.8),
        vjust=-0.4,
        size=3.5
    ) +
    labs(
        title="Accessible TF motif overlap in High-PIP versus Low-PIP sQTLs",
        subtitle="Labels show motif-overlapping variant-gene pairs / total pairs in each PIP class",
        x=NULL,
        y="Motif-overlap rate (%)",
        fill="sQTL class"
    ) +
    theme_bw(base_size=13) +
    theme(
        axis.text.x=element_text(
            angle=30,
            hjust=1
        ),
        legend.position="top"
    )

ggsave(
    file.path(
        OUTDIR,
        "01_high_vs_low_motif_overlap_rate.png"
    ),
    p1,
    width=10,
    height=6,
    dpi=300
)

ggsave(
    file.path(
        OUTDIR,
        "01_high_vs_low_motif_overlap_rate.pdf"
    ),
    p1,
    width=10,
    height=6
)

# ============================================================
# FIGURE 2
# FIRTH FOREST PLOT
# ============================================================

forest <- key[
    !is.na(odds_ratio) &
    is.finite(odds_ratio) &
    !is.na(ci_lower_95) &
    !is.na(ci_upper_95) &
    ci_lower_95 > 0 &
    ci_upper_95 > 0
]

forest[
    ,
    significance :=
        fifelse(
            !is.na(global_fdr) &
            global_fdr < 0.05,
            "Global FDR < 0.05",
            fifelse(
                !is.na(cell_fdr) &
                cell_fdr < 0.05,
                "Cell-specific FDR < 0.05",
                "Nominal candidate"
            )
        )
]

forest[
    ,
    support :=
        paste0(
            candidate,
            "\nH: ",
            n_high_overlap,
            "/",
            n_high,
            "   L: ",
            n_low_overlap,
            "/",
            n_low
        )
]

forest[
    ,
    support := factor(
        support,
        levels=support[
            order(odds_ratio)
        ]
    )
]

p2 <- ggplot(
    forest,
    aes(
        x=odds_ratio,
        y=support,
        shape=significance
    )
) +
    geom_vline(
        xintercept=1,
        linetype="dashed"
    ) +
    geom_segment(
        aes(
            x=ci_lower_95,
            xend=ci_upper_95,
            y=support,
            yend=support
        ),
        linewidth=0.7
    ) +
    geom_point(
        size=3
    ) +
    scale_x_log10() +
    labs(
        title="Firth regression: accessible motif enrichment among High-PIP sQTLs",
        subtitle="Adjusted odds ratios with 95% confidence intervals",
        x="Adjusted odds ratio (log scale)",
        y=NULL,
        shape="Statistical evidence",
        caption="Model adjusted for splice-site distance, nearest-gene TSS distance, and GC content."
    ) +
    theme_bw(base_size=13) +
    theme(
        legend.position="top"
    )

ggsave(
    file.path(
        OUTDIR,
        "02_firth_forest_plot.png"
    ),
    p2,
    width=10,
    height=6.5,
    dpi=300
)

ggsave(
    file.path(
        OUTDIR,
        "02_firth_forest_plot.pdf"
    ),
    p2,
    width=10,
    height=6.5
)

# ============================================================
# FIGURE 3
# OVERALL FIRTH LANDSCAPE
# ============================================================

land <- d[
    !is.na(p_value) &
    is.finite(p_value) &
    p_value > 0 &
    !is.na(odds_ratio) &
    is.finite(odds_ratio) &
    odds_ratio > 0
]

land[
    ,
    `:=`(
        log2_or = log2(odds_ratio),
        minuslog10_p = -log10(p_value),

        global_sig =
            !is.na(global_fdr) &
            global_fdr < 0.05,

        label = ""
    )
]

land[
    global_sig==TRUE,
    label := motif_archetype
]

land[
    calderon_celltype=="Monocytes" &
    motif_archetype=="SMAD_1",
    label := "SMAD_1"
]

land[
    calderon_celltype=="Monocytes",
    display_cell := "Resting Monocytes"
]

land[
    calderon_celltype!="Monocytes",
    display_cell := calderon_celltype
]

p3 <- ggplot(
    land,
    aes(
        x=log2_or,
        y=minuslog10_p
    )
) +
    geom_point(
        aes(shape=global_sig),
        alpha=0.65,
        size=2
    ) +
    geom_text(
        aes(label=label),
        check_overlap=TRUE,
        hjust=-0.15,
        vjust=-0.3,
        size=3.4
    ) +
    facet_wrap(
        ~display_cell,
        scales="free"
    ) +
    labs(
        title="Overall Firth motif-enrichment landscape",
        subtitle="Each point is one motif-cell model; labeled points highlight global FDR hits and SMAD_1",
        x="log2 adjusted odds ratio",
        y="-log10(p-value)",
        shape="Global FDR < 0.05"
    ) +
    theme_bw(base_size=12) +
    theme(
        legend.position="top"
    )

ggsave(
    file.path(
        OUTDIR,
        "03_overall_firth_landscape.png"
    ),
    p3,
    width=12,
    height=8,
    dpi=300
)

ggsave(
    file.path(
        OUTDIR,
        "03_overall_firth_landscape.pdf"
    ),
    p3,
    width=12,
    height=8
)

# ============================================================
# SAVE TABLE USED IN FIGURES
# ============================================================

fwrite(
    key,
    file.path(
        OUTDIR,
        "key_firth_results_used_in_figures.tsv"
    ),
    sep="\t"
)

cat("\n====================================\n")
cat("DONE\n")
cat("====================================\n")

cat("\nMain candidates plotted:\n")

print(
    key[
        ,
        .(
            motif_archetype,
            display_cell,
            n_high,
            n_high_overlap,
            n_low,
            n_low_overlap,
            odds_ratio,
            ci_lower_95,
            ci_upper_95,
            p_value,
            global_fdr,
            cell_fdr
        )
    ]
)

cat("\nFigures:\n")
cat(
    file.path(
        OUTDIR,
        "01_high_vs_low_motif_overlap_rate.png"
    ),
    "\n"
)

cat(
    file.path(
        OUTDIR,
        "02_firth_forest_plot.png"
    ),
    "\n"
)

cat(
    file.path(
        OUTDIR,
        "03_overall_firth_landscape.png"
    ),
    "\n"
)

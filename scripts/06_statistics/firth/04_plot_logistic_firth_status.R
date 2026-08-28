#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

infile <- paste0(
    "07_logistic_regression/05_firth_comparison/tables/",
    "logistic_vs_firth_status_summary.tsv"
)

outdir <- "07_logistic_regression/05_firth_comparison/plots"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(infile)

dt[, calderon_celltype := gsub('"', "", calderon_celltype)]

# Ordinary logistic categories
dt[, logistic_fit_ok := n_logistic_valid]

dt[, logistic_separation_risk :=
    n_firth_rescued_or_only
]

dt[, logistic_no_variation :=
    n_models -
    logistic_fit_ok -
    logistic_separation_risk
]

# Firth categories
dt[, firth_fit_ok := n_firth_valid]

dt[, firth_no_variation :=
    n_models - firth_fit_ok
]

logistic_long <- melt(
    dt[, .(
        calderon_celltype,
        logistic_fit_ok,
        logistic_separation_risk,
        logistic_no_variation
    )],
    id.vars = "calderon_celltype",
    variable.name = "status",
    value.name = "n_models"
)

logistic_long[, method := "Ordinary logistic"]

logistic_long[, status := fifelse(
    status == "logistic_fit_ok",
    "Fit OK",
    fifelse(
        status == "logistic_separation_risk",
        "Separation risk",
        "No motif variation"
    )
)]

firth_long <- melt(
    dt[, .(
        calderon_celltype,
        firth_fit_ok,
        firth_no_variation
    )],
    id.vars = "calderon_celltype",
    variable.name = "status",
    value.name = "n_models"
)

firth_long[, method := "Firth logistic"]

firth_long[, status := fifelse(
    status == "firth_fit_ok",
    "Fit OK",
    "No motif variation"
)]

plot_dt <- rbindlist(
    list(logistic_long, firth_long),
    use.names = TRUE
)

plot_dt[, celltype_label :=
    gsub("_", " ", calderon_celltype)
]

plot_dt[, celltype_label := factor(
    celltype_label,
    levels = c(
        "Follicular T Helper",
        "Memory Tregs",
        "Monocytes",
        "Naive B",
        "Naive CD8 T",
        "Naive Tregs"
    )
)]

plot_dt[, method := factor(
    method,
    levels = c(
        "Ordinary logistic",
        "Firth logistic"
    )
)]

plot_dt[, status := factor(
    status,
    levels = c(
        "Fit OK",
        "Separation risk",
        "No motif variation"
    )
)]

p <- ggplot(
    plot_dt,
    aes(
        x = method,
        y = n_models,
        fill = status
    )
) +
    geom_col(
        width = 0.7
    ) +
    facet_wrap(
        ~ celltype_label,
        ncol = 3
    ) +
    geom_text(
        aes(
            label = ifelse(n_models > 0, n_models, "")
        ),
        position = position_stack(vjust = 0.5),
        size = 3.5
    ) +
    labs(
        title = "Ordinary logistic versus Firth model status",
        subtitle = paste(
            "Firth rescues separation-risk models,",
            "while models with no motif variation remain untestable"
        ),
        x = NULL,
        y = "Number of motif models",
        fill = "Model status"
    ) +
    theme_bw(base_size = 13) +
    theme(
        axis.text.x = element_text(
            angle = 25,
            hjust = 1
        ),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(
            fill = "grey95"
        ),
        strip.text = element_text(
            face = "bold"
        ),
        plot.title = element_text(
            face = "bold"
        ),
        legend.position = "bottom"
    )

ggsave(
    file.path(
        outdir,
        "09_logistic_vs_firth_status_comparison.png"
    ),
    p,
    width = 12,
    height = 8,
    dpi = 300
)

ggsave(
    file.path(
        outdir,
        "09_logistic_vs_firth_status_comparison.pdf"
    ),
    p,
    width = 12,
    height = 8
)

cat("Generated status-comparison plots.\n")

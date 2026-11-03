################################################################################*
# res_5_lofo.R
# Extrapolation cost to an unseen typology family: CRPS of the leave-one-family-out
#   model relative to the size-matched control; bars: range over CV folds.
################################################################################*


library(data.table)
library(ggplot2)


args     <- commandArgs(trailingOnly = TRUE)
in_file  <- if (length(args) >= 1) args[1] else "data/009_qrf_results/lofo_contrast.csv"
out_stem <- if (length(args) >= 2) args[2] else "figures/main_text/Results/res_5"

## ---- Shared style: identical block in res_1 ... res_5 (edit in all five) ----
base_size  <- 9
fig_width  <- 6.5     # inches; MEE two-column width (~166 mm)
metric_levels <- c("PERMANOVA R2", "Classification Strength", "AucZeta",
                   "ANOSIM R", "Fuzzy Mantel", "PERMANOVA Fuzzy R2")
# plotmath labels: R^2 renders in any locale (no superscript character needed)
metric_labels <- expression(PERMANOVA~R^2, "Classification Strength", AucZeta,
                            "ANOSIM R", "Fuzzy Mantel", PERMANOVA~Fuzzy~R^2)
metric_colours <- setNames(c("#1b9e77", "#d95f02", "#7570b3",
                             "#e7298a", "#66a61e", "#e6ab02"), metric_levels)
eval_levels    <- c("Grouped CV", "Test set")
eval_shapes    <- setNames(c(16, 17), eval_levels)
eval_linetypes <- setNames(c("solid", "dashed"), eval_levels)
family_levels  <- c("crossed_thresholds", "kmeans", "spatial_ward")
family_labels  <- c("Crossed thresholds", "k-means", "Spatially constrained Ward")
family_shapes  <- setNames(c(16, 17, 15), family_levels)
ref_colour     <- "grey50"   # all reference lines: dashed, grey50

taxon_label <- function(x) tools::toTitleCase(tolower(as.character(x)))
eval_label  <- function(x) ifelse(grepl("^test", x, ignore.case = TRUE), "Test set",
                                  ifelse(grepl("cv", x, ignore.case = TRUE), "Grouped CV", x))
read_input  <- function(path, need) {
        d <- data.table::fread(path)
        miss <- setdiff(need, names(d))
        if (length(miss)) stop("Missing columns in ", path, ": ", paste(miss, collapse = ", "))
        if ("metric" %in% names(d)) {
                unknown <- setdiff(unique(d$metric), metric_levels)
                if (length(unknown)) stop("Unknown metric(s) in ", path, ": ", paste(unknown, collapse = ", "))
                d[, metric := factor(metric, levels = metric_levels)]
        }
        if ("taxon" %in% names(d)) d[, taxon := factor(taxon_label(taxon))]
        if ("evaluation" %in% names(d)) d[, evaluation := factor(eval_label(evaluation), levels = eval_levels)]
        d[]
}
scale_colour_metric <- function(...) scale_colour_manual(values = metric_colours, labels = metric_labels,
                                                         drop = FALSE, name = "Metric", ...)
theme_res <- function() {
        theme_bw(base_size = base_size) +
                theme(legend.position = "bottom", legend.box = "vertical",
                      legend.spacing.y = unit(0.1, "lines"), panel.grid.minor = element_blank(),
                      panel.spacing = unit(0.8, "lines"))
}
save_res <- function(p, stem, height) {
        dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
        ggsave(paste0(stem, ".png"), p, width = fig_width, height = height, dpi = 300, bg = "white")
        # ggsave(paste0(stem, ".pdf"), p, width = fig_width, height = height,
               # device = if (capabilities("cairo")) cairo_pdf else "pdf")
        cat("saved:", paste0(stem, c(".png", ".pdf")), sep = "\n")
}
## ---- end of shared style ----

evaluation_shown <- "Grouped CV"   # or "Test set" (no fold range on the test set)

## Data -----------------------------------------------------------------------
lc <- read_input(in_file, c("taxon", "metric", "evaluation", "family", "crps_change",
                            "crps_change_fold_min", "crps_change_fold_max"))
lc <- lc[evaluation == evaluation_shown]
unknown <- setdiff(unique(lc$family), family_levels)
if (length(unknown)) stop("Unknown typology family: ", paste(unknown, collapse = ", "))
lc[, family := factor(family, levels = family_levels)]

# Metrics on the y axis in the shared order (first metric on top), families dodged within a metric
n_met <- length(metric_levels); n_fam <- length(family_levels); dodge_width <- 0.5
lc[, y_pos := (n_met + 1 - as.integer(metric)) -
           (as.integer(family) - (n_fam + 1) / 2) * dodge_width / (n_fam - 1)]

## Figure ---------------------------------------------------------------------
p <- ggplot(lc, aes(y = y_pos, colour = metric, shape = family)) +
        geom_vline(xintercept = 0, linetype = 2, colour = ref_colour) +
        geom_linerange(aes(xmin = 100 * crps_change_fold_min, xmax = 100 * crps_change_fold_max),
                       linewidth = 0.5, alpha = 0.6, na.rm = TRUE) +
        geom_point(aes(x = 100 * crps_change), size = 2) +
        facet_wrap(~ taxon, nrow = 1) +
        scale_colour_metric(guide = "none") +     # metrics are named on the y axis
        scale_shape_manual(values = family_shapes, labels = family_labels, name = "Held-out family") +
        scale_y_continuous(breaks = rev(seq_len(n_met)), labels = metric_labels,
                           expand = expansion(add = 0.5)) +
        scale_x_continuous(breaks = function(lims) pretty(lims, n = 3)) +
        labs(x = "CRPS of LOFO relative to the size-matched control (%)", y = NULL) +
        theme_res()

save_res(p, out_stem, height = 3.8)

### - LOFO supplementary figures: calibration and 90% coverage by held-out family - ###
#
# Inputs (written by qrf_summary.R into output_dir):
#   lofo_calibration.csv   coverage per taxon x metric x evaluation x family x model x nominal level
#   lofo_contrast.csv      LOFO vs. size-matched control per taxon x metric x evaluation x family
# Outputs: one calibration and one coverage figure per evaluation (grouped CV, test set).

# ==========================================================================*
# 01. Setup ---------------------------------------------------------------
# ==========================================================================*

## -- libraries -- 
library(data.table)
library(ggplot2)


## --- Paths ---
output_dir <- "~/projects/pulse/01_wp1/MEE_REVISION_1/data/QRF_RESULTS/tables"
fig_dir    <- "figures/supplement/S7"

## -- Settings --
main_lev <- 0.9   # headline interval level of qrf_fit.R

metric_levels_fig <- c("PERMANOVA R2", "PERMANOVA Fuzzy R2", "Classification Strength",
                       "AucZeta", "ANOSIM R", "Fuzzy Mantel")
metric_colours <- c(`PERMANOVA R2`            = "#1b9e77",
                    `Classification Strength` = "#d95f02",
                    AucZeta                   = "#7570b3",
                    `ANOSIM R`                = "#e7298a",
                    `Fuzzy Mantel`            = "#66a61e",
                    `PERMANOVA Fuzzy R2`      = "#e6ab02")
family_levels <- c("crossed_thresholds", "kmeans", "spatial_ward")
family_labels <- c(crossed_thresholds = "Crossed thresholds", kmeans = "k-means",
                   spatial_ward = "Spatially constrained Ward")
family_shapes <- c(crossed_thresholds = 21, kmeans = 24, spatial_ward = 22)   # fillable shapes
model_linetypes <- c(LOFO = "solid", control = "dashed", full = "dotted")
model_labels    <- c(LOFO = "LOFO", control = "Size-matched control", full = "Final model")


## -- custom functions -- 
prep <- function(d) {
        d[, metric := factor(as.character(metric), levels = metric_levels_fig)]
        d[, family := factor(as.character(family), levels = family_levels)]
        stopifnot(!anyNA(d$metric), !anyNA(d$family))   # catches name mismatches
        d[]
}

theme_lofo <- function() {
        theme_bw() + theme(panel.grid.minor = element_blank(), panel.spacing = unit(0.8, "lines"),
                           legend.position = "bottom", legend.box = "vertical")
}

## --  Calibration by held-out family (coverage error, 
## as in the main calibration figure) -- 
plot_calibration <- function() {
        ggplot(cal, aes(nominal, coverage - nominal, colour = metric, linetype = model,
                        group = interaction(metric, model))) +
                geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
                geom_line(linewidth = 0.45) +
                geom_point(size = 0.7) +
                facet_grid(evaluation + fam_lab ~ taxon, scales = "free_y") +
                scale_colour_manual(values = metric_colours, breaks = metric_levels_fig, name = NULL) +
                scale_linetype_manual(values = model_linetypes, labels = model_labels, name = NULL) +
                scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
                labs(x = "Nominal coverage of the central interval",
                     y = "Coverage error (empirical − nominal)") +
                theme_lofo() +
                theme(strip.text.y = element_text(size = 7))
}

## 90% coverage: size-matched control (open) -> LOFO (filled)
plot_coverage <- function() {
        d   <- copy(lc)
        d[, metric_rev := factor(metric, levels = rev(metric_levels_fig))]   # first metric at the top
        off <- 0.25
        d[, y_pos := as.integer(metric_rev) - (as.integer(family) - (nlevels(family) + 1) / 2) * off]
        n_m <- nlevels(d$metric_rev)
        ggplot(d, aes(y = y_pos, colour = metric, shape = family)) +
                geom_hline(yintercept = seq(1.5, n_m - 0.5, 1), colour = "grey90", linewidth = 0.3) +
                geom_vline(xintercept = main_lev, linetype = 2, colour = "grey50") +
                geom_segment(aes(x = coverage_control, xend = coverage_LOFO, yend = y_pos), linewidth = 0.4) +
                geom_point(aes(x = coverage_control), fill = "white", size = 2.2, stroke = 0.6) +
                geom_point(aes(x = coverage_LOFO, fill = metric), size = 2.2, stroke = 0.6) +
                facet_grid(evaluation ~ taxon) +
                scale_y_continuous(breaks = seq_len(n_m), labels = levels(d$metric_rev),
                                   expand = expansion(add = 0.5)) +
                scale_colour_manual(values = metric_colours, guide = "none") +
                scale_fill_manual(values = metric_colours, guide = "none") +
                scale_shape_manual(values = family_shapes, breaks = family_levels,
                                   labels = unname(family_labels[family_levels]), name = "Held-out family") +
                guides(shape = guide_legend(override.aes = list(colour = "grey20", fill = "grey20", size = 2.6))) +
                labs(x = paste0("Coverage of the ", 100 * main_lev, "% interval"), y = NULL) +
                theme_lofo() +
                theme(panel.grid.major.y = element_blank())
}

# CRPS 
plot_cost <- function(e) {
        d   <- copy(lc[evaluation == e])
        d[, metric_rev := factor(metric, levels = rev(metric_levels_fig))]   # first metric at the top
        off <- 0.25
        d[, y_pos := as.integer(metric_rev) - (as.integer(family) - (nlevels(family) + 1) / 2) * off]
        n_m <- nlevels(d$metric_rev)
        ggplot(d, aes(y = y_pos, colour = metric, shape = family)) +
                geom_hline(yintercept = seq(1.5, n_m - 0.5, 1), colour = "grey90", linewidth = 0.3) +
                geom_vline(xintercept = 0, colour = "grey50") +
                geom_linerange(aes(xmin = 100 * crps_change_fold_min, xmax = 100 * crps_change_fold_max),
                               linewidth = 0.5, alpha = 0.6, na.rm = TRUE) +                # fold range; absent on the test set
                geom_point(aes(x = 100 * crps_change, fill = metric), size = 2.4, stroke = 0.6) +
                facet_wrap(~ taxon, nrow = 1) +
                scale_y_continuous(breaks = seq_len(n_m), labels = levels(d$metric_rev),
                                   expand = expansion(add = 0.5)) +
                scale_colour_manual(values = metric_colours, guide = "none") +
                scale_fill_manual(values = metric_colours, guide = "none") +
                scale_shape_manual(values = family_shapes, breaks = family_levels,
                                   labels = unname(family_labels[family_levels]), name = "Held-out family") +
                guides(shape = guide_legend(override.aes = list(colour = "grey20", fill = "grey20", size = 2.6))) +
                labs(x = "CRPS of LOFO relative to the size-matched control (%)", y = NULL) +
                theme_lofo() +
                theme(panel.grid.major.y = element_blank())
}
# ==========================================================================*
# 02. Load Data -----------------------------------------------------------
# ==========================================================================*


cal <- prep(fread(file.path(output_dir, "lofo_calibration.csv")))
lc  <- prep(fread(file.path(output_dir, "lofo_contrast.csv")))


# ==========================================================================*
# 03. Prepare Data --------------------------------------------------------
# ==========================================================================*


cal[, model := factor(model, levels = names(model_linetypes))]
cal[, fam_lab := factor(family_labels[as.character(family)], levels = family_labels[family_levels])]

#===============================================================================
# 3. Save
#===============================================================================
#
for (e in c("grouped CV", "test set")) {
        tag <- gsub(" ", "_", e)
        ggsave(file.path(fig_dir, paste0("lofo_crps_", tag, ".png")), plot_cost(e),
               width = 10, height = 5.5, dpi = 300, bg = "white")
}
ggsave(file.path(fig_dir, "S7_lofo_calibration.png"), plot_calibration(),
       width = 10, height = 12, dpi = 300, bg = "white")
# ggsave(file.path(fig_dir, "S7_lofo_coverage.png"), plot_coverage(),
#        width = 10, height = 9, dpi = 300, bg = "white")
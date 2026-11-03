### - Supplementary figure: influence of individual data sets and data-set weighting - ###
suppressPackageStartupMessages({
        library(data.table)
        library(ggplot2)
        library(egg)
})

## --- Paths and settings ----
output_dir <- "data/QRF_RESULTS/tables"
tab_dir    <- output_dir                       # save_table() writes here in your summary script
fig_dir    <- "figures/supplement/S7/"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

metric_levels_fig <- c("PERMANOVA R2", "PERMANOVA Fuzzy R2", "Classification Strength",
                       "AucZeta", "ANOSIM R", "Fuzzy Mantel")
metric_colours <- c(`PERMANOVA R2`            = "#1b9e77",
                    `Classification Strength` = "#d95f02",
                    AucZeta                   = "#7570b3",
                    `ANOSIM R`                = "#e7298a",
                    `Fuzzy Mantel`            = "#66a61e",
                    `PERMANOVA Fuzzy R2`      = "#e6ab02")
var_labels <- c(vp_env_pure = "Pure environmental fraction", vp_space_pure = "Pure spatial fraction",
                vp_shared = "Shared fraction", env_r2_after = "Environmental R²",
                n_vars = "No. of variables", n_types = "No. of types", type_evenness = "Type evenness",
                n_taxa_used = "No. of taxa", asw_after = "Average silhouette width",
                env_ch_after = "Calinski–Harabasz index")

## --- Data ----
infl <- fread(file.path(tab_dir, "dataset_influence.csv"))           # one row per LODO refit
isum <- fread(file.path(tab_dir, "dataset_influence_summary.csv"))   # one row per model
wt   <- fread(file.path(tab_dir, "weighting_pooled_vs_per_dataset.csv"))
alst <- as.data.table(readRDS(file.path(output_dir, "..", "summary_tables.rds"))$parts$ale_stability)
# one row per refit x tracked predictor
set_metric <- function(d) {
        d[, metric := factor(as.character(metric), levels = metric_levels_fig)]
        stopifnot(!anyNA(d$metric))                                        # catches name mismatches
        invisible(d)
}
for (d in list(infl, isum, wt, alst)) set_metric(d)

wo <- wt[evaluation == "all data sets, out-of-sample"]

## --- Numbers for the text (check before filling brackets) ----
cat(sprintf("Prediction shift: median %.3f SD, max %.3f SD (%d refits)\n",
            median(infl$pred_shift), max(infl$pred_shift), nrow(infl)))
cat(sprintf("Test CRPS change: median %+.2f%%, max %+.1f%%, |change| > 5%%: %.1f%% of refits\n",
            100 * median(infl$test_crps_change), 100 * max(infl$test_crps_change),
            100 * mean(abs(infl$test_crps_change) > 0.05)))
cat(sprintf("Spearman rho, size vs. shift: %.2f to %.2f\n",
            min(isum$rho_size_shift, na.rm = TRUE), max(isum$rho_size_shift, na.rm = TRUE)))
cat(sprintf("ALE: %d comparisons, %d class changes; curve shift median %.1f%%, max %.1f%%\n",
            nrow(alst), sum(alst$changed), 100 * median(alst$shift_rel), 100 * max(alst$shift_rel)))
cat(sprintf("Equal weights: test CRPS %+.1f%% to %+.1f%% (median %+.1f%%)\n",
            100 * min(isum$eq_crps_change, na.rm = TRUE), 100 * max(isum$eq_crps_change, na.rm = TRUE),
            100 * median(isum$eq_crps_change, na.rm = TRUE)))
wd <- wo[, .(d_skill = skill[weighting == "per data set"] - skill[weighting == "pooled (rows)"],
             d_cov   = coverage[weighting == "per data set"] - coverage[weighting == "pooled (rows)"]),
         by = .(taxon, metric)]
cat(sprintf("Per data set - pooled: skill %+.3f to %+.3f, coverage max |diff| %.3f\n",
            min(wd$d_skill), max(wd$d_skill), max(abs(wd$d_cov))))

## --- Shared theme ----
theme_panel <- function() {
        theme_bw() + theme(panel.grid.minor = element_blank(), panel.spacing = unit(1, "lines"),
                           legend.position = "bottom")
}

## --- (a) Prediction shift vs. size of the removed data set ----
pa <- ggplot(infl, aes(100 * share_train, pred_shift, colour = metric)) +
        geom_point(size = 1.4, alpha = 0.8) +
        facet_wrap(~ taxon, nrow = 1) +
        scale_colour_manual(values = metric_colours, guide = "none") +
        labs(x = "Share of training rows of the removed data set (%)",
             y = "Mean abs. shift of test\nmedian predictions (SD)") +
        theme_panel()

## --- (b) Change in test CRPS per refit; diamond: equal data-set weights ----
pb <- ggplot(infl, aes(100 * test_crps_change, metric, colour = metric)) +
        geom_vline(xintercept = 0, colour = "grey50") +
        geom_vline(xintercept = c(-5, 5), linetype = 2, colour = "grey75") +
        geom_point(position = position_jitter(height = 0.15, width = 0, seed = 1), size = 1.3, alpha = 0.7) +
        geom_point(data = isum[!is.na(eq_crps_change)], aes(x = 100 * eq_crps_change),
                   shape = 23, size = 2.6, fill = "white", colour = "black") +
        facet_wrap(~ taxon, nrow = 1) +
        scale_y_discrete(limits = rev(metric_levels_fig)) +
        scale_colour_manual(values = metric_colours, guide = "none") +
        labs(x = "Change in test-set CRPS (%)", y = NULL) +
        theme_panel()

## --- (c) Shift of the main-effect curves of the tracked predictors ----
alst[, var_lab := fifelse(as.character(variable) %in% names(var_labels),
                          var_labels[as.character(variable)], as.character(variable))]
alst[, var_lab := reorder(factor(var_lab), shift_rel, median)]
pc <- ggplot(alst, aes(100 * shift_rel, var_lab)) +
        geom_boxplot(outlier.shape = NA, colour = "grey40", width = 0.5) +
        geom_point(aes(colour = metric), position = position_jitter(height = 0.15, width = 0, seed = 1),
                   size = 1, alpha = 0.6) +
        scale_colour_manual(values = metric_colours, guide = "none") +
        labs(x = "Shift of the median-ALE curve (% of the full-model amplitude)", y = NULL) +
        theme_panel()

## --- (d) Pooled vs. per-data-set evaluation ----
wl  <- melt(wo, id.vars = c("taxon", "metric", "weighting"), measure.vars = c("skill", "coverage"))
wd2 <- dcast(wl, taxon + metric + variable ~ weighting, value.var = "value")
setnames(wd2, c("per data set", "pooled (rows)"), c("per_ds", "pooled"))
wd2[, variable := factor(variable, levels = c("skill", "coverage"),
                         labels = c("CRPS skill", "90% coverage"))]
pd <- ggplot(wd2, aes(pooled, per_ds, colour = metric, shape = taxon)) +
        geom_abline(linetype = 2, colour = "grey60") +
        geom_point(size = 1.8) +
        facet_wrap(~ variable, scales = "free") +
        scale_colour_manual(values = metric_colours, breaks = metric_levels_fig, name = NULL) +
        scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
        labs(x = "Every row weighted equally", y = "Every data set\nweighted equally") +
        theme_panel()

## --- Combine and save ----
p <- egg::ggarrange(pa, pb, pc, pd, ncol = 1, heights = c(1, 1.4, 1, 1),
                    labels = c("a", "b", "c", "d"), draw = FALSE)

ggsave(file.path(fig_dir, "S7_dataset_influence_weighting.png"), p,
       width = 10, height = 13, dpi = 300, bg = "white")

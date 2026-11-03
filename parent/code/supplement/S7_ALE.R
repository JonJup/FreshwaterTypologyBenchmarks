### - ALE results: main effects of the dominant predictors - ###
#
# Inputs (written by qrf_summary.R into output_dir):
#   ale_summary.csv                 one row per model x predictor x output (class, net, mono, amplitude)
#   ale_additivity.csv              share of prediction variance explained by first-order effects
#   importance_tail_vs_centre.csv   tail-minus-centre contrast of pinball importance
#   summary_tables.rds              $parts$ale_curves: ALE curves of every model
# Outputs: CSV tables in <output_dir>/ale_results, PNG figures in fig_dir; key numbers printed.
# ALE outputs of qrf_fit.R: q0.1, q0.5, q0.9 and width = q0.9 - q0.1.

suppressPackageStartupMessages({
        library(data.table)
        library(ggplot2)
})

## --- Paths and settings ----
output_dir <- "~/projects/pulse/01_wp1/MEE_REVISION_1/data/QRF_RESULTS/tables/"
res_dir    <- file.path(output_dir, "ale_results")
fig_dir    <- "figures/supplement/S7/"
for (d in c(res_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

main_vars    <- c("vp_env_pure", "vp_space_pure", "env_r2_after", "n_vars")
focus_metric <- "PERMANOVA R2"     # metric shown in outcome units in the main-text panel

var_labels <- c(vp_env_pure = "Pure environmental fraction", vp_space_pure = "Pure spatial fraction",
                env_r2_after = "Environmental R²", n_vars = "No. of typology variables")
out_labels <- c(q0.1 = "10th percentile", q0.5 = "Median", q0.9 = "90th percentile",
                width = "Width (90th − 10th percentile)")

metric_levels_fig <- c("PERMANOVA R2", "PERMANOVA Fuzzy R2", "Classification Strength",
                       "AucZeta", "ANOSIM R", "Fuzzy Mantel")
metric_colours <- c(`PERMANOVA R2`            = "#1b9e77",
                    `Classification Strength` = "#d95f02",
                    AucZeta                   = "#7570b3",
                    `ANOSIM R`                = "#e7298a",
                    `Fuzzy Mantel`            = "#66a61e",
                    `PERMANOVA Fuzzy R2`      = "#e6ab02")
taxon_linetypes <- c(Diatoms = "solid", Fish = "dashed", Invertebrates = "dotdash", Macrophytes = "dotted")

## --- Data ----
as_metric <- function(d) {
        d[, metric := factor(as.character(metric), levels = metric_levels_fig)]
        stopifnot(!anyNA(d$metric))                       # catches name mismatches
        d[, taxon := as.character(taxon)]
        d[]
}
ale_sum  <- as_metric(fread(file.path(output_dir, "ale_summary.csv")))
ale_add  <- as_metric(fread(file.path(output_dir, "ale_additivity.csv")))
imp_tail <- as_metric(fread(file.path(output_dir, "importance_tail_vs_centre.csv")))
curves   <- as_metric(as.data.table(readRDS(file.path(output_dir, "..", "summary_tables.rds"))$parts$ale_curves))

ale_sum$taxon <- stringr::str_to_title(ale_sum$taxon)


stopifnot(all(main_vars %in% ale_sum$variable),
          all(c("q0.5", "width") %in% ale_sum$output),
          all(names(taxon_linetypes) %in% ale_sum$taxon))
n_models <- uniqueN(ale_sum[, .(taxon, metric)])
cat("Models:", n_models, "\n\n")

#===============================================================================
# 1. Key results ----
#===============================================================================

## --- Direction and shape of the main effects (median and width) ----
dir_long <- ale_sum[variable %in% main_vars & output %in% c("q0.5", "width"),
                    .(n = .N), by = .(variable, output, class)]
dir_tab  <- dcast(dir_long, variable + output ~ class, value.var = "n", fill = 0)
dir_tab[, n_models := n_models]
fwrite(dir_tab, file.path(res_dir, "ale_direction_counts.csv"))
cat("Class counts of the main effects (models per class):\n"); print(dir_tab)

cat("\nStatements for the text:\n")
for (v in main_vars) for (o in c("q0.5", "width")) {
        d   <- dir_long[variable == v & output == o][order(-n)]
        cat(sprintf("  %-28s %-6s %s\n", var_labels[v], o,
                    paste(sprintf("%s in %d/%d", d$class, d$n, n_models), collapse = ", ")))
}

## --- Tails vs. centre of the importance ----
tail_tab <- imp_tail[variable %in% main_vars,
                     .(median_contrast_pp = median(contrast), min_pp = min(contrast), max_pp = max(contrast),
                       n_more_in_tails = sum(contrast > 0), n = .N), by = variable]
# fwrite(tail_tab, file.path(res_dir, "ale_tail_vs_centre.csv"))
cat("\nTail - centre contrast of the importance share (percentage points):\n"); print(tail_tab)

## --- Additivity ----
add_tab <- ale_add[, .(min = min(r2_main_effects), median = median(r2_main_effects),
                       max = max(r2_main_effects), n_negative = sum(r2_main_effects < 0), n = .N),
                   by = output]
# fwrite(add_tab, file.path(res_dir, "ale_additivity_summary.csv"))
cat("\nVariance of the predictions explained by first-order ALE effects:\n"); print(add_tab)

#===============================================================================
# 2. Figures ----
#===============================================================================
theme_ale <- function() {
        theme_bw() + theme(panel.grid.minor = element_blank(), panel.spacing = unit(0.8, "lines"),
                           strip.placement = "outside", strip.background.x = element_blank(),
                           legend.position = "bottom")
}
prep_curves <- function(d) {
        d <- d[variable %in% main_vars & output %in% c("q0.5", "width")]
        d[, var_lab := factor(var_labels[variable], levels = var_labels[main_vars])]
        d[, out_lab := factor(out_labels[output], levels = out_labels[c("q0.5", "width")])]
        d[]
}

## --- Main text: all models, curves scaled to unit amplitude ----
# Scaling removes the outcome units, so the shape and direction of every model's effect
# share one axis; effect sizes are shown by the importance figure.
cv <- prep_curves(copy(curves))
cv[, amp := max(ale) - min(ale), by = .(taxon, metric, variable, output)]
cv <- cv[amp > 0][, ale_scaled := ale / amp]

p_all <- 
        cv %>%
        filter(out_lab == "Median") %>%
        ggplot(aes(x, ale_scaled, colour = metric, group = interaction(taxon, metric))) +
        geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
        geom_line(linewidth = 0.4, alpha = 0.7) +
        facet_grid(. ~ var_lab, scales = "free_x", switch = "x") +
        scale_colour_manual(values = metric_colours, breaks = metric_levels_fig, name = NULL) +
        labs(x = NULL, y = "Accumulated local effect (scaled to unit amplitude)") +
        theme_ale()
p_all
ggsave(file.path(fig_dir, "S7_ale_main_effects_all_models.png"), p_all,
       width = 10, height = 6, dpi = 300, bg = "white")

# ## --- One metric in outcome units, with group-bootstrap CIs ----
# plot_metric <- function(m) {
#         d   <- prep_curves(curves[metric == m])
#         col <- metric_colours[[m]]
#         ggplot(d, aes(x, ale, group = taxon)) +
#                 geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
#                 geom_ribbon(aes(ymin = lo, ymax = hi), fill = col, alpha = 0.12) +
#                 geom_line(aes(linetype = taxon), colour = col, linewidth = 0.6) +
#                 facet_grid(out_lab ~ var_lab, scales = "free", switch = "x") +
#                 scale_linetype_manual(values = taxon_linetypes, name = NULL) +
#                 labs(x = NULL, y = paste0("Accumulated local effect (", m, ")")) +
#                 theme_ale()
# }
# ggsave(file.path(fig_dir, "S7_ale_main_effects_focus.png"), plot_metric(focus_metric),
#        width = 10, height = 6, dpi = 300, bg = "white")

# Supplement: one figure per metric
# for (m in metric_levels_fig) {
#         ggsave(file.path(fig_dir, paste0("ale_main_effects_", gsub("[^A-Za-z0-9]+", "_", m), ".png")),
#                plot_metric(m), width = 10, height = 6, dpi = 300, bg = "white")
# }

## --- Supplement: additivity ----
# ale_add[, out_lab := factor(out_labels[output], levels = out_labels)]
# p_add <- ggplot(ale_add, aes(r2_main_effects, metric, colour = metric, shape = taxon)) +
#         geom_vline(xintercept = 0, colour = "grey50") +
#         geom_vline(xintercept = 1, linetype = 2, colour = "grey70") +
#         geom_point(size = 2) +
#         facet_wrap(~ out_lab, nrow = 1, scales = "free_x") +
#         scale_y_discrete(limits = rev(metric_levels_fig)) +
#         scale_colour_manual(values = metric_colours, guide = "none") +
#         scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
#         labs(x = "Variance of the predictions explained by first-order ALE effects", y = NULL) +
#         theme_ale()
# ggsave(file.path(fig_dir, "ale_additivity.png"), p_add, width = 11, height = 4, dpi = 300, bg = "white")


#===============================================================================
# 3. Importance across quantile levels: tails vs. centre ----
#     Share of a predictor in a model's total positive pinball importance at each tau
#     (grouped CV). Tails: tau <= 0.1 or >= 0.9; centre: |tau - 0.5| <= 0.1,
#     as in importance_tail_vs_centre.csv.
#===============================================================================
imp <- as.data.table(readRDS(file.path(output_dir,"..", "summary_tables.rds"))$parts$importance)
pb  <- imp[evaluation == "grouped CV" & measure == "pinball"]
pb[, share := 100 * pmax(rel, 0) / sum(pmax(rel, 0)), by = .(taxon, metric, tau)]
pb[!is.finite(share), share := 0]
pb <- as_metric(pb)

## --- Numbers: lower tail, centre and upper tail separately ----
eps  <- 1e-9
prof <- pb[variable %in% main_vars,
           .(lower  = mean(share[tau <= 0.1 + eps]),
             centre = mean(share[abs(tau - 0.5) <= 0.1 + eps]),
             upper  = mean(share[tau >= 0.9 - eps])),
           by = .(taxon, metric, variable)]
prof_tab <- prof[, .(lower = median(lower), centre = median(centre), upper = median(upper),
                     n_lower_gt_centre = sum(lower > centre),
                     n_upper_gt_centre = sum(upper > centre), n = .N), by = variable]
# fwrite(prof_tab, file.path(res_dir, "importance_by_tail.csv"))
cat("\nImportance share (%) in the lower tail, centre and upper tail (median over models):\n")
print(prof_tab, digits = 3)

## --- Figure ----
# d   <- pb[variable %in% main_vars]
# d[, var_lab := factor(var_labels[as.character(variable)], levels = var_labels[main_vars])]
# med <- d[, .(share = median(share)), by = .(var_lab, tau)]
# 
# p_tau <- ggplot(d, aes(tau, share)) +
#         annotate("rect", xmin = -Inf, xmax = 0.1, ymin = -Inf, ymax = Inf, fill = "grey92") +
#         annotate("rect", xmin = 0.9,  xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey92") +
#         annotate("rect", xmin = 0.4,  xmax = 0.6, ymin = -Inf, ymax = Inf, fill = "grey97") +
#         geom_line(aes(colour = metric, group = interaction(taxon, metric)), linewidth = 0.35, alpha = 0.6) +
#         geom_line(data = med, colour = "black", linewidth = 0.9) +
#         facet_wrap(. ~ var_lab, nrow = 1) + # , scales = "free_y"
#         scale_colour_manual(values = metric_colours, breaks = metric_levels_fig, name = NULL) +
#         scale_x_continuous(breaks = c(0.1, 0.5, 0.9)) +
#         labs(x = expression(tau), y = "Share of pinball importance at τ (%)") +
#         theme_ale()
# p_tau

## --- Figure: one line per metric (median over taxa) ----
# d   <- pb[variable %in% main_vars]
# d[, var_lab := factor(var_labels[as.character(variable)], levels = var_labels[main_vars])]
# dm  <- d[, .(share = median(share)), by = .(metric, var_lab, tau)]    # median over taxa
# med <- d[, .(share = median(share)), by = .(var_lab, tau)]            # median over all 24 models (as in the text)
# 
# p_tau <- ggplot(dm, aes(tau, share, colour = metric)) +
#         annotate("rect", xmin = -Inf, xmax = 0.1, ymin = -Inf, ymax = Inf, fill = "grey92") +
#         annotate("rect", xmin = 0.9,  xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey92") +
#         annotate("rect", xmin = 0.4,  xmax = 0.6, ymin = -Inf, ymax = Inf, fill = "grey97") +
#         geom_line(linewidth = 0.6) +
#         geom_line(data = med, aes(tau, share), inherit.aes = FALSE, colour = "black", linewidth = 0.9) +
#         facet_wrap(~ var_lab, nrow = 2) +
#         scale_colour_manual(values = metric_colours, breaks = metric_levels_fig, name = NULL) +
#         scale_x_continuous(breaks = c(0.1, 0.5, 0.9)) +
#         labs(x = expression(tau), y = "Share of pinball importance at τ (%)") +
#         theme_ale()
# p_tau
# 
# ggsave(file.path(fig_dir, "S7_importance_by_tau.png"), p_tau, width = 6, height = 6, dpi = 300, bg = "white")



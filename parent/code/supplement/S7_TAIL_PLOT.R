### - Figure S28: share of pinball importance across quantile levels - ###
#
# Reads the two tables written by the importance-summary script:
#   importance_by_tau.csv   per-tau shares (curves)      <- fwrite(per_tau, ...)
#   importance_by_band.csv  pooled band shares (text)    <- fwrite(tail_imp, ...)
# Coloured lines: median across taxa per metric; black line: median across all
# 24 models. Shading: lower tail (tau <= 0.1), centre (0.4-0.6), upper tail (>= 0.9).


library(data.table)
library(ggplot2)


## --- Settings ----
tau_file  <- "data/supplement/S7_importance_by_tau.csv"
band_file <- "data/supplement/S7_importance_by_band.csv"
out_file  <- "figures/supplement/S7/figure_S28.png"
smooth    <- FALSE        # TRUE: loess-smoothed curves (display only)

key_vars <- c(vp_env_pure   = "Pure environmental fraction",
              vp_space_pure = "Pure spatial fraction",
              env_r2_after  = "Environmental R\u00b2",
              n_vars        = "No. of typology variables")
metric_labels <- c(PERMANOVA_R2            = "PERMANOVA R2",
                   Classification_Strength = "Classification Strength",
                   ANOSIM_R                = "ANOSIM R",
                   PERMANOVA_Fuzzy_R2      = "PERMANOVA Fuzzy R2",
                   AucZeta                 = "AucZeta",
                   Fuzzy_Mantel            = "Fuzzy Mantel")
metric_cols <- c("PERMANOVA R2" = "#1B9E77", "Classification Strength" = "#D95F02",
                 "ANOSIM R" = "#E7298A", "PERMANOVA Fuzzy R2" = "#E6AB02",
                 "AucZeta" = "#7570B3", "Fuzzy Mantel" = "#66A61E")

## --- Data ----
per_tau <- fread(tau_file)
stopifnot(uniqueN(per_tau$run) == 24, all(per_tau$metric %in% names(metric_labels)))

d <- per_tau[variable %in% names(key_vars)]
d[, `:=`(panel  = factor(key_vars[variable], levels = key_vars),
         metric = factor(metric_labels[metric], levels = metric_labels))]
by_metric <- d[, .(share = median(share)), by = .(panel, metric, tau)]
overall   <- d[, .(share = median(share)), by = .(panel, tau)]

## --- Plot ----
bands <- data.table(xmin = c(-Inf, 0.4, 0.9), xmax = c(0.1, 0.6, Inf),
                    fill = c("tail", "centre", "tail"))
line_geom <- function(data, ...) {
        if (smooth) geom_smooth(data = data, method = "loess", formula = y ~ x,
                                span = 0.5, se = FALSE, ...)
        else geom_line(data = data, ...)
}

p <- ggplot() +
        geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
                  inherit.aes = FALSE, show.legend = FALSE) +
        scale_fill_manual(values = c(tail = "grey90", centre = "grey96")) +
        line_geom(by_metric, mapping = aes(x = tau, y = share, colour = metric, group = metric),
                  linewidth = 0.7) +
        line_geom(overall, mapping = aes(x = tau, y = share), colour = "black", linewidth = 1.3) +
        facet_wrap(~ panel, ncol = 2) +
        scale_colour_manual(values = metric_cols, name = NULL) +
        scale_x_continuous(breaks = c(0.1, 0.5, 0.9)) +
        scale_y_continuous(limits = c(0, NA), breaks = seq(0, 100, 20)) +
        labs(x = expression(tau), y = expression("Share of pinball importance at " * tau * " (%)")) +
        guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
        theme_bw(base_size = 11) +
        theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
              strip.background = element_blank(), strip.text = element_text(size = 10),
              legend.position = "bottom", legend.key.width = unit(1.5, "lines"))

ggsave(out_file, p, width = 7, height = 7, dpi = 300, bg = "white")
cat("saved:", out_file, "\n")

## --- Band numbers for S7.6 / 3.3 (same files as the figure) ----
band <- fread(band_file)[variable %in% names(key_vars)]
print(dcast(band[, .(share = round(median(share), 1)), by = .(variable, band)],
            variable ~ band, value.var = "share")[, .(variable, lower, centre, upper)])
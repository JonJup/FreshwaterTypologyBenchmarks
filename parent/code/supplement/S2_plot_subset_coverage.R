#!/usr/bin/env Rscript
## plot_subset_coverage.R
## Figures for subset_coverage.R. Run from the project root.
##   fig_subset_coverage_curves      coverage vs. pressure threshold, all groups;
##                                   dashed = all sites (no subset rules, no filter)
##   fig_subset_coverage_loss        loss of covered area, split into subset rules
##                                   and pressure filter
##   fig_subset_dropout_map_<g>      threshold at which each catchment stops being covered
##   fig_subset_lost_mod_<g>         what limits coverage in catchments covered by all
##                                   sites but not at the strictest threshold

library(data.table)
library(sf)
library(ggplot2)
source("code/functions/coverage_metrics.R")   # lab()

## ---- settings -----------------------------------------------------------------
cov_dir     <- "data/supplement/S1_environmental_coverage/subsets"
fig_dir     <- file.path("figures/supplement/S2/")
path_geom_s <- "E://Arbeit/data/river_network/eu_hydro_dem/processedData/all_eu_hydro_dem_catchments_simplified_250m.rds"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

s   <- fread(file.path(cov_dir, "subset_coverage_summary.csv"))
thr <- s[level != "all", sort(unique(threshold_pct))]

## ---- shared elements ----------------------------------------------------------
theme_set(theme_bw(base_size = 10))
x_thr <- scale_x_reverse(breaks = thr,
                         name = "Pressure threshold (quantile of dist_ref, %)  ->  stricter")
pct   <- scales::label_percent(accuracy = 1)
metric_labs <- c(within_ranges = "Within univariate ranges (MESS >= 0)",
                 covered       = "Covered in joint space (NT2 <= 1)")

save_fig <- function(p, name, w = 7, h = 5) {
        for (ext in c("png", "pdf"))
                ggsave(file.path(fig_dir, paste0(name, ".", ext)), p, width = w, height = h, dpi = 300)
}

## ---- Fig: coverage curves -----------------------------------------------------
long <- melt(s, id.vars = c("group", "level", "threshold_pct"),
             measure.vars = names(metric_labs), variable.name = "metric", value.name = "share")

p1 <- ggplot(long[level != "all"], aes(threshold_pct, share, colour = group)) +
        geom_hline(data = long[level == "all"], aes(yintercept = share, colour = group),
                   linetype = 2) +
        geom_line() +
        geom_point() +
        facet_wrap(~ metric, labeller = as_labeller(metric_labs)) +
        x_thr +
        scale_y_continuous(labels = pct) +
        scale_colour_brewer(palette = "Dark2", name = NULL) +
        labs(y = "Share of European catchment area",
             caption = "Dashed: all sites, without subset rules or pressure filter") +
        theme(legend.position = "bottom")
save_fig(p1, "fig_subset_coverage_curves", w = 8, h = 4.5)

## ---- per group: drop-out map and limiting variables ---------------------------
geom     <- readRDS(path_geom_s)
pc_files <- list.files(cov_dir, "^subset_coverage_per_catchment_.*\\.rds$", full.names = TRUE)

q_lv   <- sprintf("q%03d", sort(thr, decreasing = TRUE))   # loosest -> strictest
q_num  <- sort(thr, decreasing = TRUE)
lost_labs <- sprintf("Lost at %d %%", sort(thr))
drop_cols <- c("Covered at all thresholds" = "#009E73",
               setNames(tail(hcl.colors(length(thr) + 2, "YlOrRd", rev = TRUE), length(thr)),
                        lost_labs),
               "Lost through subset rules" = "#CC79A7",
               "Not covered by all sites"  = "grey75",
               "No data"                   = "white")

for (f in pc_files) {
        g <- sub("^subset_coverage_per_catchment_(.*)\\.rds$", "\\1", basename(f))
        d <- readRDS(f)
        W <- sum(d$w)

        ## levels with too few sites have no class column: treated as not covered
        covered_at <- function(lv) {
                col <- paste0("class_", lv)
                if (col %in% names(d)) d[[col]] %in% "covered" else rep(FALSE, nrow(d))
        }
        base <- covered_at("all")
        fail <- !vapply(q_lv, covered_at, logical(nrow(d)))     # catchments x levels
        first_fail <- max.col(fail, ties.method = "first")
        first_fail[rowSums(fail) == 0] <- NA

        drop <- rep("Covered at all thresholds", nrow(d))
        lost <- !is.na(first_fail)
        drop[lost] <- sprintf("Lost at %d %%", q_num[first_fail[lost]])
        drop[lost & first_fail == 1L] <- "Lost through subset rules"
        drop[!base] <- "Not covered by all sites"

        ## ---- map ----
        i  <- match(geom$ID, d$ID)
        gd <- geom
        gd$drop <- factor(ifelse(is.na(i), "No data", drop[i]), levels = names(drop_cols))
        p3 <- ggplot() +
                geom_sf(data = gd, aes(fill = drop), colour = NA) +
                scale_fill_manual(values = drop_cols, name = NULL, drop = FALSE) +
                guides(fill = guide_legend(ncol = 4)) +
                coord_sf(datum = NA) +
                theme_void(base_size = 10) + theme(legend.position = "bottom")
        save_fig(p3, sprintf("fig_subset_dropout_map_%s", g), h = 8)
}

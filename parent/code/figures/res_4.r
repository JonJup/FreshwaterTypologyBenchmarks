################################################################################*
# res_4_importance.R
# Median share of CRPS permutation importance per metric and predictor across taxa
#   (grouped CV); labels: number of taxa in which the predictor is relevant.
################################################################################*


# setup -------------------------------------------------------------------
library(data.table)
library(ggplot2)


## settings ----------------------------------------------------------------
args     <- commandArgs(trailingOnly = TRUE)
in_file  <- if (length(args) >= 1) args[1] else "data/009_qrf_results/importance_crps_cv.csv"
out_stem <- if (length(args) >= 2) args[2] else "figures/main_text/Results/res_4"

## Shared style: identical block in res_1 ... res_5 (edit in all five) ----
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

## Settings -------------------------------------------------------------------
# Predictors shown, by class (top to bottom). Set show_sampling_classes to TRUE to
# add the sampling-related predictors (negligible in most models).
show_sampling_classes <- FALSE
rows <- list(
        "Environmental separation" = c(asw_after = "Average silhouette width",
                                       env_r2_after = "Environmental R^2",
                                       env_ch_after = "Calinski-Harabasz index"),
        "Taxon richness"           = c(n_taxa_used = "No. of taxa"),
        "Typology structure"       = c(type_evenness = "Type evenness", n_vars = "No. of variables",
                                       n_types = "No. of types"),
        "Variation partitioning"   = c(vp_shared = "Shared fraction", vp_space_pure = "Pure spatial fraction",
                                       vp_env_pure = "Pure environmental fraction"))
if (show_sampling_classes) rows <- c(rows, list(
        "Sampling effort"      = c(samples = "No. of samples"),
        "Spatial extent"       = c(min_distance = "Minimum distance", median_distance = "Median distance"),
        "Taxonomic resolution" = c(species_rank = "Species level (%)", genus_rank = "Genus level (%)",
                                   family_rank = "Family level (%)", higher_rank = "Higher levels (%)"),
        "Sampling season"      = c(mean_doy = "Mean sampling day", circ_sd_days = "Seasonal spread")))
lab_threshold <- 35   # share (%) above which cell labels turn white

## Data -----------------------------------------------------------------------
imp <- read_input(in_file, c("taxon", "metric", "variable", "share", "relevant"))
n_tax <- nlevels(imp$taxon)
lk <- data.table(pclass   = rep(names(rows), lengths(rows)),
                 variable = unlist(lapply(rows, names)),
                 label    = unlist(rows, use.names = FALSE))
# all metric x predictor combinations; a predictor absent from a model counts as share 0, not relevant
grid <- CJ(taxon = levels(imp$taxon), metric = metric_levels, variable = lk$variable)
grid[, `:=`(taxon = factor(taxon, levels = levels(imp$taxon)), metric = factor(metric, levels = metric_levels))]
d <- imp[grid, on = .(taxon, metric, variable)][is.na(share), `:=`(share = 0, relevant = FALSE)]
d <- d[, .(share = median(share), n_rel = sum(relevant %in% TRUE)), by = .(metric, variable)]
d <- lk[d, on = "variable"]
d[, `:=`(pclass   = factor(pclass, levels = names(rows)),
         variable = factor(variable, levels = rev(lk$variable)),
         text     = fifelse(n_rel > 0, paste0(n_rel, "/", n_tax), ""))]
y_labels <- setNames(lk$label, lk$variable)

## Figure ---------------------------------------------------------------------
p <- ggplot(d, aes(x = metric, y = variable, fill = share)) +
        geom_tile(colour = "white", linewidth = 0.4) +
        geom_text(aes(label = text, colour = share > lab_threshold), size = base_size / 3.4, show.legend = FALSE) +
        scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white")) +
        scale_fill_gradient(low = "white", high = "#b2182b", limits = c(0, NA),
                            name = "Median share of\nCRPS importance (%)") +
        scale_x_discrete(labels = metric_labels) +
        scale_y_discrete(labels = function(v) as.expression(lapply(y_labels[as.character(v)], function(s)
                if (grepl("^", s, fixed = TRUE)) str2lang(gsub(" ", "~", s)) else s))) +
        facet_grid(pclass ~ ., scales = "free_y", space = "free_y") +
        labs(x = NULL, y = NULL) +
        theme_res() +
        theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
              strip.text.y = element_text(angle = 0, hjust = 0))

save_res(p, out_stem, height = if (show_sampling_classes) 6.8 else 5)

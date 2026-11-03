################################################################################
# Script Name:        pressure_threshold_sensitivity.R
# Description:        Sensitivity of the scheme definition (define_schemes.R) to
#                     the human-pressure threshold. For each taxon, counts how
#                     many yearly subsets (data set x year) survive when sites
#                     are restricted to dist_ref <= q, where q is the 10 %,
#                     20 %, ..., 100 % quantile of dist_ref across all
#                     catchments in the stressor database.
#
#                     A yearly subset "survives" under exactly the same rules as
#                     in define_schemes.R:
#                       1. the year has >= min_samples unique events,
#                       2. after restricting to the focal months (three
#                          consecutive months with most samples, or the
#                          fallback rules for <= 2 months), >= min_samples
#                          unique events remain.
#                     The schemes themselves are not built; only counted.
#
#                     Outputs
#                       - summary table per taxon x threshold:
#                           n_subsets   number of yearly subsets
#                           n_schemes   number of schemes (subsets x sample sizes)
#                           n_data_sets number of data sets contributing
#                       - four figures (see section 5)
# Notes:              Sites without a dist_ref value (no match in the stressor
#                     database) are excluded at every threshold, including
#                     100 %. Their number is reported per taxon in section 3.
################################################################################

# 1. setup ----
## 1.1 wd, packages, scripts ----
setwd(rstudioapi::getActiveProject())
library(data.table)
library(ggplot2)
source("code/functions/find_max_consequtive_sum.R")

## 1.2 parameters (must match define_schemes.R) ----
min_samples    <- 100
multiplier     <- 1.5
max_iter       <- 10
max_schemes    <- 6
threshold_pct  <- seq(10L, 100L, by = 10L)
log_sequence   <- unique(round(min_samples * multiplier^(0:max_iter)))

fig_dir <- "figures/pressure_thresholds"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

## 1.3 load data ----
files     <- list.files("data/000_biota", full.names = TRUE, pattern = "02_")
bio.names <- sub("_w_environment.rds", "", sub("^02_", "", basename(files)), fixed = FALSE)
bio.list  <- lapply(files, readRDS)

pressures <- readRDS("E://Arbeit/data/my_stressor_db/processedData/PCA_stressors.rds")
p_sub     <- as.data.table(pressures)[, .(ID, dist_ref)]

#- thresholds as quantiles of dist_ref over all catchments (as l10 in define_schemes.R)
thresholds <- data.table(
        threshold_pct = threshold_pct,
        dist_ref_max  = quantile(p_sub$dist_ref, threshold_pct / 100,
                                 na.rm = TRUE, names = FALSE)
)

# 2. functions ----

#' Join pressure data to one taxon's biological data
#'
#' @param bio data.table with one taxon's samples (columns ID, data.set,
#'   eventID, eventYear, eventDate).
#' @return bio with column dist_ref added; list-type eventYear is flattened.
prepare_taxon <- function(bio, missing_passes = T) {
        bio <- copy(as.data.table(bio))
        if (is.list(bio$eventYear)) bio[, eventYear := unlist(eventYear)]
        bio <- p_sub[bio, on = "ID"]
        if (missing_passes) bio[is.na(dist_ref), dist_ref := -Inf]
        bio
}

#' Restrict one year of data to its focal months
#'
#' Replicates the seasonal-window rules of define_schemes.R:
#' > 2 months sampled -> three consecutive months with most samples;
#' 2 months <= 2 apart -> both; otherwise -> the single month with most
#' samples, provided it holds > min_samples events.
#'
#' @param x data.table of one data set and one year.
#' @param min_samples minimum number of unique events.
#' @return x restricted to the focal months (zero rows if none qualify).
select_focal_months <- function(x, min_samples) {
        x           <- copy(x)[, month := month(eventDate)]
        month_table <- table(unique(x, by = "eventID")$month)
        
        if (length(month_table) > 2) {
                focal <- names(find_max_consecutive_sum(month_table)$values)
        } else {
                month_diff  <- diff(as.numeric(names(month_table)))
                consecutive <- length(month_diff) == 1 && month_diff <= 2
                if (consecutive) {
                        focal <- names(month_table)
                } else if (any(month_table > min_samples)) {
                        focal <- names(which.max(month_table))
                } else {
                        return(x[0])
                }
        }
        x[month %in% focal]
}

#' Number of schemes (sample sizes) for a subset of n events
#'
#' Replicates the sample-scheme rules of define_schemes.R: geometric sizes
#' <= n, plus n itself; the second-largest size is dropped if it is within
#' 25 % of n; at most max_schemes sizes are kept.
#'
#' @param n number of unique events in the subset (>= min_samples).
#' @return integer, number of schemes.
count_schemes <- function(n, log_sequence, max_schemes) {
        s <- log_sequence[log_sequence <= n]
        if (max(s) < n) s <- c(s, n)
        if (length(s) > 1 && (n - s[length(s) - 1]) / n < 0.25) {
                s <- s[-(length(s) - 1)]
        }
        as.integer(min(length(s), max_schemes))
}

#' Yearly subsets of one data set that pass all sample-size rules
#'
#' @param ds data.table of one data set, already filtered by pressure.
#' @return data.table with one row per surviving year (eventYear, n_samples,
#'   n_schemes), or NULL if none survive.
count_year_subsets <- function(ds, min_samples, log_sequence, max_schemes) {
        if (nrow(ds) == 0) return(NULL)
        
        years_all_na <- all(is.na(ds$eventYear))
        use_months   <- !years_all_na && !all(is.na(ds$eventDate))
        
        n_per_year <- ds[, .(n = uniqueN(eventID)), by = eventYear][n >= min_samples]
        #- an NA year is only kept if no year information exists at all
        if (!years_all_na) n_per_year <- n_per_year[!is.na(eventYear)]
        if (nrow(n_per_year) == 0) return(NULL)
        
        rbindlist(lapply(n_per_year$eventYear, function(yr) {
                x <- ds[eventYear %in% yr]
                if (use_months) x <- select_focal_months(x, min_samples)
                n <- uniqueN(x$eventID)
                if (n < min_samples) return(NULL)
                data.table(eventYear = yr,
                           n_samples = n,
                           n_schemes = count_schemes(n, log_sequence, max_schemes))
        }))
}

# 3. count surviving subsets ----
subset_table <- rbindlist(lapply(seq_along(bio.list), function(b) {
        bio <- prepare_taxon(bio.list[[b]])
        message(sprintf("%s: %d of %d events without dist_ref (always excluded)",
                        bio.names[b],
                        uniqueN(bio[is.na(dist_ref), eventID]),
                        uniqueN(bio$eventID)))
        
        rbindlist(
                lapply(
                        seq_len(nrow(thresholds)), 
                        function(t) {
                                bio[dist_ref <= thresholds$dist_ref_max[t],
                                    count_year_subsets(.SD, min_samples, log_sequence, max_schemes),
                                    by = data.set
                                    ][, `:=`(taxon = bio.names[b],
                                             threshold_pct = thresholds$threshold_pct[t])
                                      ]
        }), fill = TRUE)
}), fill = TRUE)

# 4. summarise ----
## 4.1 per taxon x threshold (zeros filled in) ----
summary_dt <- subset_table[, .(n_subsets   = .N,
                               n_schemes   = sum(n_schemes),
                               n_data_sets = uniqueN(data.set)),
                           by = .(taxon, threshold_pct)]
# summary_dt <- summary_dt[CJ(taxon = bio.names, threshold_pct = threshold_pct),
#                          on = .(taxon, threshold_pct)]
# count_cols <- c("n_subsets", "n_schemes", "n_data_sets")
# summary_dt[, (count_cols) := lapply(.SD, function(x) fifelse(is.na(x), 0L, x)),
#            .SDcols = count_cols]
summary_dt <- thresholds[summary_dt, on = "threshold_pct"]
summary_dt[, retention := n_subsets / n_subsets[threshold_pct == 100], by = taxon]
setcolorder(summary_dt, c("taxon", "threshold_pct", "dist_ref_max"))
setorder(summary_dt, taxon, threshold_pct)

print(dcast(summary_dt, taxon ~ threshold_pct, value.var = "n_schemes"))

## 4.2 sanity check: 10 % must reproduce define_schemes.R output ----
scheme_files <- paste0("data/000_biota/03_", bio.names, "_scheme.rds")
if (all(file.exists(scheme_files))) {
        original <- data.table(
                taxon    = bio.names,
                original = vapply(scheme_files,
                                  function(f) uniqueN(readRDS(f), by = c("data.set", "eventYear", "samples")),
                                  integer(1), USE.NAMES = FALSE)
        )
        check <- summary_dt[threshold_pct == 100, .(taxon, recount = n_schemes)][original, on = "taxon"]
        print(check)
        if (check[, any(recount != original)]) warning("Recount at 10 % differs from define_schemes.R output.")
}

## 4.3 save ----
saveRDS(list(summary = summary_dt, subsets = subset_table, thresholds = thresholds),
        "data/000_biota/pressure_threshold_counts.rds")

# 5. figures ----
## 5.1 shared plot elements ----
theme_set(theme_bw(base_size = 11))
taxon_colours     <- scale_colour_brewer(palette = "Dark2", name = "Taxon")
current_threshold <- geom_vline(xintercept = 10, linetype = "dashed", colour = "grey40")
x_threshold       <- scale_x_continuous(breaks = threshold_pct,
                                        name = "Pressure threshold (quantile of dist_ref, %)")
save_plot <- function(plot, name, width = 8, height = 5) {
        ggsave(file.path(fig_dir, paste0(name, ".png")), plot,
               width = width, height = height, dpi = 300)
}

## 5.2 absolute counts: subsets, schemes, data sets ----
metric_labels <- c(n_subsets   = "Yearly subsets",
                   n_schemes   = "Schemes (subsets x sample sizes)",
                   n_data_sets = "Contributing data sets")
counts_long <- melt(summary_dt, id.vars = c("taxon", "threshold_pct"),
                    measure.vars = count_cols, variable.name = "metric", value.name = "count")

p_counts <- ggplot(counts_long, aes(threshold_pct, count, colour = taxon)) +
        current_threshold +
        geom_line() +
        geom_point() +
        facet_wrap(~ metric, scales = "free_y", ncol = 1,
                   labeller = as_labeller(metric_labels)) +
        x_threshold + taxon_colours +
        labs(y = "Count",
             title = "Retained subsets across pressure thresholds",
             subtitle = "Dashed line: 10% threshold")
save_plot(p_counts, "01_counts", height = 8)

## 5.3 relative retention (share of subsets available at 100 %) ----
p_retention <- ggplot(summary_dt, aes(threshold_pct, retention, colour = taxon)) +
        current_threshold +
        geom_abline(slope = 1 / 100, intercept = 0, linetype = "dotted", colour = "grey60") +
        geom_line() +
        geom_point() +
        x_threshold + taxon_colours +
        scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1),
                           name = "Share of subsets retained at 100 %") +
        labs(title = "Relative retention of yearly subsets",
             subtitle = "Dotted line: retention proportional to the share of catchments admitted")
save_plot(p_retention, "02_retention")

## 5.4 heatmap: which data sets enter at which threshold ----
heat_dt <- unique(subset_table[, .(taxon, data.set)])[, .(threshold_pct = threshold_pct),
                                                      by = .(taxon, data.set)]
heat_dt <- subset_table[, .(n_subsets = .N), by = .(taxon, data.set, threshold_pct)
][heat_dt, on = .(taxon, data.set, threshold_pct)]

p_heat <- ggplot(heat_dt, aes(threshold_pct, data.set, fill = n_subsets)) +
        geom_tile(colour = "white") +
        geom_text(aes(label = n_subsets), size = 2.5, na.rm = TRUE) +
        facet_grid(taxon ~ ., scales = "free_y", space = "free_y") +
        scale_fill_viridis_c(name = "Yearly\nsubsets", na.value = "grey92") +
        x_threshold +
        labs(y = "Data set",
             title = "Yearly subsets per data set and threshold",
             subtitle = "Grey: data set contributes no subset at this threshold")
save_plot(p_heat, "03_data_set_heatmap",
          height = 2 + 0.18 * uniqueN(heat_dt, by = c("taxon", "data.set")))

## 5.5 subset size: do subsets grow as the threshold is relaxed? ----
p_size <- ggplot(subset_table, aes(threshold_pct, n_samples, group = threshold_pct)) +
        geom_hline(yintercept = min_samples, linetype = "dashed", colour = "grey40") +
        geom_boxplot(outlier.size = 0.8) +
        facet_wrap(~ taxon) +
        scale_y_log10(name = "Samples per yearly subset (log scale)") +
        x_threshold +
        labs(title = "Size of retained yearly subsets",
             subtitle = sprintf("Dashed line: minimum subset size (%d samples)", min_samples))
save_plot(p_size, "04_subset_size")

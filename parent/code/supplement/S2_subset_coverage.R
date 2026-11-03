## subset_coverage.R
## Environmental coverage of European river catchments (EU-Hydro) by the sites that
## actually enter the analysis: the sites of the yearly subsets built with the rules
## of define_schemes.R (code/functions/scheme_rules.R), at pressure thresholds of
## 10 %, 20 %, ..., 100 %. Run from the project root.
##
## Levels
##   all    every site of the taxon in the biota data: no subset rules, no pressure
##          filter. Reference for all losses.
##   qXXX   sites of the yearly subsets after keeping sites with dist_ref <= the
##          XXX % quantile of dist_ref over all catchments (as l10 in define_schemes.R).
##
## Reference sample per group and level: the union of the sites of all yearly
## subsets (one row per site; "all" = union over taxa).
##
## Loss (shares of European catchment area, per metric)
##   <metric>_loss_total     all  -> qXXX
##   <metric>_loss_pressure  q100 -> qXXX   (pressure filter alone)
##   <metric>_loss_total at q100 = loss caused by the subset rules alone.
##   lost / gained: area covered at "all" but not at qXXX, and vice versa. NT2 is not
##   monotone in the reference sample, so net loss of "covered" = lost - gained.
##
## Metrics as in mess_coverage.R (MESS, ExDet NT2). No random-removal null.

library(data.table)
source("code/functions/find_max_consequtive_sum.R")
source("code/functions/scheme_rules.R")
source("code/functions/coverage_metrics.R")

## ---- settings -----------------------------------------------------------------
path_eu    <- "E://Arbeit/data/river_network/eu_hydro_dem/processedData/all_catchments_w_environment.rds"
areas_dir  <- "E://Arbeit/data/river_network/eu_hydro_dem/processedData/areas/"
path_sites <- "data/new_environment.rds"
path_press <- "E://Arbeit/data/my_stressor_db/processedData/PCA_stressors.rds"
biota_dir  <- "data/000_biota"
out_dir    <- "data/environmental_coverage/subsets"

nat_vars <- c("slope", "lake_area", "soil_oc", "soil_pH",
              "inundation_depth", "flooded_area", "spi_mean", "area_calcareous",
              "area_siliceous", "area_sediment", "glacial_area", "bioclim05", "bioclim06",
              "bioclim13", "bioclim14", "Rfactor_max", "Rfactor_avg", "Rfactor_min",
              "roughness", "elevation", "mean_snow_equivalent", "mean_discharge",
              "max_discharge", "min_discharge", "groundwater_table",
              "upstream_catchment_area", "saturated_soil_water_content", "VBM_mean")

min_samples        <- 100                     # as in define_schemes.R
threshold_pct      <- seq(10L, 100L, by = 10L)
na_pressure_passes <- TRUE                    # sites without dist_ref kept at every threshold
min_n              <- length(nat_vars) + 2    # minimum reference size for NT2
class_levels       <- c("covered", "combinatorial", "univariate")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- data ---------------------------------------------------------------------
eu    <- load_catchments(path_eu, areas_dir, nat_vars)
w     <- eu$area
sites <- load_sites(path_sites, nat_vars)
p_sub <- as.data.table(readRDS(path_press))[, .(ID, dist_ref)]

files <- list.files(biota_dir, pattern = "^02_.*_w_environment\\.rds$", full.names = TRUE)
taxa  <- sub("^02_(.*)_w_environment\\.rds$", "\\1", basename(files))
if (length(setdiff(taxa, sites$taxon)))
        stop("Taxa missing from the site table: ",
             paste(setdiff(taxa, sites$taxon), collapse = ", "))

levels_dt <- data.table(
        level         = c("all", sprintf("q%03d", threshold_pct)),
        threshold_pct = c(NA_integer_, threshold_pct),
        dist_ref_max  = c(NA_real_, quantile(p_sub$dist_ref, threshold_pct / 100,
                                             na.rm = TRUE, names = FALSE))
)

## ---- catchments used per taxon and level ---------------------------------------
## Events used at one level; max_dist = NA gives the unfiltered baseline.
level_events <- function(bio, max_dist) {
        if (is.na(max_dist)) return(unique(bio[, .(data.set, eventYear, eventID, ID)]))
        filter_pressure(bio, max_dist, na_pressure_passes)[
                , build_year_subsets(.SD, min_samples), by = data.set]
}

used <- rbindlist(lapply(seq_along(files), function(b) {
        bio <- prepare_taxon(readRDS(files[b]), p_sub)
        rbindlist(lapply(seq_len(nrow(levels_dt)), function(k) {
                ev  <- level_events(bio, levels_dt$dist_ref_max[k])
                has <- "ID" %in% names(ev)   # FALSE if no subset survives
                data.table(
                        group     = taxa[b],
                        level     = levels_dt$level[k],
                        n_subsets = if (k == 1) NA_integer_
                                    else if (has) uniqueN(ev, by = c("data.set", "eventYear"))
                                    else 0L,
                        ids       = list(if (has) unique(ev$ID) else NULL))
        }))
}))
used <- rbind(used,
              used[, .(group = "all", n_subsets = sum(n_subsets),
                       ids = list(unique(unlist(ids)))), by = level],
              use.names = TRUE)

ref_sites <- function(g, ids) {
        s <- if (g == "all") sites else sites[taxon == g]
        unique(s[ID %in% ids], by = "siteID")
}

## ---- coverage per group and level ----------------------------------------------
res <- list()
for (g in unique(used$group)) {
        u        <- used[group == g]
        pc       <- data.table(ID = eu$ID, w = w)
        base_cov <- NULL

        for (k in seq_len(nrow(u))) {
                lv  <- u$level[k]
                ref <- ref_sites(g, u$ids[[k]])
                row <- list(group = g, level = lv, n_subsets = u$n_subsets[k],
                            n_sites = nrow(ref))

                if (nrow(ref) < min_n) {
                        warning(sprintf("%s, %s: %d sites, too few to compute coverage",
                                        g, lv, nrow(ref)))
                } else {
                        cl  <- classify(eu, ref, nat_vars)
                        cov <- cl$class == "covered"
                        if (lv == "all") base_cov <- cov
                        pc[, (paste0(c("class_", "mod_"), lv)) := list(
                                factor(cl$class, levels = class_levels),
                                factor(fifelse(cl$mess < 0, cl$mod, NA_character_),
                                       levels = nat_vars))]
                        row <- c(row, as.list(shares(cl, w)))
                        if (!is.null(base_cov))
                                row <- c(row, lost   = sum(w[base_cov & !cov]) / sum(w),
                                              gained = sum(w[!base_cov & cov]) / sum(w))
                }
                res[[length(res) + 1]] <- as.data.table(row)
                message(sprintf("%s, %s done", g, lv))
        }
        ## per-catchment classes per level, ready to join to the EU-Hydro geometry
        saveRDS(pc, file.path(out_dir, sprintf("subset_coverage_per_catchment_%s.rds", g)))
}

## ---- summary and losses --------------------------------------------------------
summary_tab <- levels_dt[rbindlist(res, fill = TRUE), on = "level"]
for (m in c("within_ranges", "covered")) {
        summary_tab[, (paste0(m, c("_loss_total", "_loss_pressure"))) := {
                x <- get(m)
                list(x[level == "all"] - x,
                     fifelse(level == "all", NA_real_, x[level == "q100"] - x))
        }, by = group]
}
setcolorder(summary_tab, c("group", "level", "threshold_pct", "dist_ref_max"))

fwrite(summary_tab, file.path(out_dir, "subset_coverage_summary.csv"))
saveRDS(list(summary = summary_tab, levels = levels_dt, nat_vars = nat_vars,
             settings = list(min_samples = min_samples,
                             na_pressure_passes = na_pressure_passes)),
        file.path(out_dir, "subset_coverage.rds"))

print(dcast(summary_tab, group ~ level, value.var = "covered"))

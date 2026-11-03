#!/usr/bin/env Rscript
## plot_coverage.R
## Figures and draft results text for the environmental-coverage analysis
## (outputs of mess_coverage.R). One set of outputs per group found in cov_dir.
##
## Figures (per group):
##   fig_coverage_map        coverage class per catchment + sites
##   fig_mod_map             most dissimilar variable where outside the sampled range
##   fig_mod_bar             area share outside the range, by limiting variable
##   fig_stringency          coverage vs. share of least-disturbed sites retained, vs. random removal
##   fig_lost_map            catchments covered by the full data but not by the strictest subset
##   fig_lost_profile        how lost catchments differ from those that stay covered
## Text: coverage_results_text.md (numbers filled in, interpretation in [brackets])

library(data.table)
library(sf)
library(ggplot2)

## ---- settings -----------------------------------------------------------------
cov_dir      <- "data/environmental_coverage"
fig_dir      <- file.path(cov_dir, "figures")
path_eu      <- "E://Arbeit/data/river_network/eu_hydro_dem/processedData/all_catchments_w_environment.rds"
path_sites   <- "data/new_environment.rds"
path_press   <- "E://Arbeit/data/my_stressor_db/processedData/PCA_stressors.rds"
path_geom_s  <- file.path("E://Arbeit/data/river_network/eu_hydro_dem/processedData/all_eu_hydro_dem_catchments_simplified_250m.rds")  # cache for drawing
id_col       <- "ID"


nat_vars <- c("slope", "lake_area", "soil_oc", "soil_pH",
              "inundation_depth", "flooded_area", "spi_mean", "area_calcareous",
              "area_siliceous", "area_sediment", "glacial_area", "bioclim05", "bioclim06",
              "bioclim13", "bioclim14", "Rfactor_max", "Rfactor_avg", "Rfactor_min",
              "roughness", "elevation", "mean_snow_equivalent", "mean_discharge",
              "max_discharge", "min_discharge", "groundwater_table",
              "upstream_catchment_area", "saturated_soil_water_content", "VBM_mean")

## readable labels; unlisted variables keep their column names
var_labels <- c(bioclim05 = "Max. temperature, warmest month",
                bioclim06 = "Min. temperature, coldest month",
                bioclim13 = "Precipitation, wettest month",
                bioclim14 = "Precipitation, driest month")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- helpers ------------------------------------------------------------------
lab <- function(v) unname(ifelse(v %in% names(var_labels), var_labels[v], v))
pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x)
big <- function(x) format(round(x), big.mark = ",", scientific = FALSE)
and_list <- function(x) if (length(x) < 2) x else
        paste(paste(head(x, -1), collapse = ", "), "and", tail(x, 1))

# oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
#         blue = "#0072B2", vermillion = "#D55E00", purple = "#CC79A7")

oi <- c(orange     = "#E69F00",
        sky        = "#56B4E9",
        green      = "#009E73",
        blue       = "#0072B2",
        vermillion = "#D55E00",
        purple     = "#CC79A7",
        # extensions
        yellow     = "#F0E442",  
        indigo     = "#332288",  
        wine       = "#882255",  
        olive      = "#999933")  

class_cols <- c(covered = oi[["green"]], combinatorial = oi[["orange"]],
                univariate = oi[["vermillion"]], nodata = "grey85")
class_labs <- c(covered = "Covered", combinatorial = "Novel combination",
                univariate = "Outside sampled range", nodata = "No data")

theme_map <- theme_void(base_size = 10) + theme(legend.position = "bottom")

save_fig <- function(p, name, w = 7, h = 7) {
        for (ext in c("png", "pdf"))
                ggsave(file.path(fig_dir, paste0(name, ".", ext)), p, width = w, height = h, dpi = 300)
}

## ---- data ---------------------------------------------------------------------
# summary_tab <- fread(file.path(cov_dir, "coverage_summary.csv"))
pc_files    <- list.files(cov_dir, "^coverage_per_catchment_.*\\.rds$", full.names = TRUE)
groups      <- sub("^coverage_per_catchment_(.*)\\.rds$", "\\1", basename(pc_files))
pc          <- setNames(lapply(pc_files, function(f) setDT(readRDS(f))), groups)

eu <- setDT(readRDS(path_eu))[, c(id_col, nat_vars), with = FALSE]
geom <- readRDS(path_geom_s)


nice_names <- 
        c(
 lake_area                      = "Lake area",
 area_calcareous                = "Calcareous area",
 area_sediment                  = "Sediment area",
 slope                          = "Slope",
 soil_oc                        = "Soil organic Carbon",
 bioclim05                      = "Max Temperature of Warmest Month",
 area_siliceous                 = "Siliceous area",
 spi_mean                       = "Stream Power Index",
 glacial_area                   = "Glacial Area",
 inundation_depth               = "Inundation Depth at flood",
 roughness                      = "Rougness",
 mean_snow_equivalent           = "Mean snow equaivalent",
 VBM_mean                       = "Valley Bottom Flatness",
 Rfactor_min                    = "R factor (minimum)",
 elevation                      = "Elevation",
 flooded_area                   = "Flooded area",
 upstream_catchment_area        = "Upstream catchment area",
 bioclim06                      = "Min Temperature of Coldest Month",
 Rfactor_max                    = "R factor (maximum)",
 max_discharge                  = "Discharge (maximum)",
 mean_discharge                 = "Discharge (mean)",
 min_discharge                  = "Discharge (minimum)",
 Rfactor_avg                    = "R factor (mean)",
 bioclim13                      = "Precipitation of Wettest Month",
 groundwater_table              = "Groundwater table",
 soil_pH                        = "Soil pH",
 bioclim14                      = "Precipitation of Driest Month",
 saturated_soil_water_content   = "Saturated soil water content"
)


## same site set as in mess_coverage.R
sites <- readRDS(path_sites)
for (i in seq_along(sites)) sites[[i]]$taxon <- names(sites)[i]
sites <- rbindlist(sites, fill = TRUE)
sites <- sites[complete.cases(sites[, ..nat_vars])]
sites <- sites[ID %in% readRDS(path_press)$ID]

site_pts <- function(g) {
        x <- if (g == "all") unique(sites, by = c("x.coord", "y.coord")) else sites[taxon == g]
        st_as_sf(x, coords = c("x.coord", "y.coord"), crs = 3035)
}

## ---- per group ----------------------------------------------------------------
txt <- c("# Environmental coverage: draft results", "")

for (g in groups) {
        
        if (g == "all") next
        
        d <- pc[[g]]
        
        W <- sum(d$w)
        i <- match(geom[[id_col]], d[[id_col]])   
        gd <- geom
        
        d$mod2 <- nice_names[d$mod]
        
        sh  <- d[, .(share = sum(w) / W), by = class_full]
        shr <- function(cl) sum(sh[class_full == cl, share])
        within  <- d[mess >= 0, sum(w)] / W
        mod_tab <- d[mess < 0, .(share = sum(w) / W), by = mod2][order(-share)]
        
        ## ---- Fig: coverage classes ----
        gd$class <- factor(fifelse(is.na(i), "nodata", d$class_full[i]), levels = names(class_cols))
        p1 <- ggplot() +
                geom_sf(data = gd, aes(fill = class), colour = NA) +
                scale_fill_manual(values = class_cols, labels = class_labs, name = NULL, drop = FALSE) +
                coord_sf(datum = NA) + theme_map + ggtitle(g)
        save_fig(p1, sprintf("fig_coverage_map_%s", g))

}

# writeLines(txt, file.path(cov_dir, "coverage_results_text.md"))
# cat(txt, sep = "\n")

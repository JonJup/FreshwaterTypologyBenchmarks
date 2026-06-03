################################################################################
# Script Name:        combine_env.R
# Description:        Combines 18 environmental predictor tables (one parquet
#                     file per catchment tile per variable) into a single
#                     wide parquet file per catchment tile by merging on the
#                     shared catchment ID column.
# Notes:              - All variable directories must contain the same set of
#                       catchment tiles (same filenames, same count).
#                     - "Hondo" and "Iceland" tiles are excluded project-wide.
#                     - Output goes to data/catchments/ (created if absent).
################################################################################

# setup -----------------------------------------------------------------------
library(fs)
library(magrittr)
library(stringr)
library(arrow)
library(data.table)

# load data -------------------------------------------------------------------
slope                   <- dir_ls("D:/Arbeit/data/river_network/hydrography90m/processedData/slope_w_catchments/")
lakes                   <- dir_ls("D:/Arbeit/data/river_network/hydroLAKES/processedData/fractionLake_w_catchment/")
soil_oc                 <- dir_ls("D:/Arbeit/data/soil/new soil/processedData/oc_w_catchment/")
soil_pH                 <- dir_ls("D:/Arbeit/data/soil/soil_ph_europe/processedData/w_catchments/")
floodpl                 <- dir_ls("D:/Arbeit/data/river_network/flood_area/processedData/w_catchment/")
streamp                 <- dir_ls("D:/Arbeit/data/river_network/hydrography90m/processedData/spi_w_catchment/")
geology                 <- dir_ls("D:/Arbeit/data/LULC/IHME1500_v11/w_catchments/")
glacier                 <- dir_ls("D:/Arbeit/data/LULC/randolph glacier inventory/w_catchments/")
bioclim                 <- dir_ls("D:/Arbeit/data/climate/bioclim_chelsa/processedData/w_catchments/")
erosion                 <- dir_ls("D:/Arbeit/data/misc/Gloreda/processedData/w_catchements/")
roughness               <- dir_ls("D:/Arbeit/data/DEM/DTM_Europe/processedData/roughness_w_catchment/")
elevation               <- dir_ls("D:/Arbeit/data/DEM/DTM_Europe/processedData/elevation_w_catchment/")
snowdepth               <- dir_ls("D:/Arbeit/data/CERRA/processedData/snow_depth_w_catchments/")
discharge               <- dir_ls("D:/Arbeit/data/river_network/discharge copernicus/processedData/w_catchments/")
groundwater             <- dir_ls("D:/Arbeit/data/river_network/BasinATLAS_Data_v10/processedData/groundwater_area_w_catchment/")
upstream_area           <- dir_ls("D:/Arbeit/data/river_network/BasinATLAS_Data_v10/processedData/upstream_area_w_catchment/")
saturated_water_content <- dir_ls("D:/Arbeit/data/river_network/soilHydraulicPropertiesEurope/processedData/soiHydraulicaProperties_w_catchments/")
valley_bottom_flatness  <- dir_ls("D:/Arbeit/data/river_network/valley_bottom_flattness/processedData/valleyBottomFlattness_w_catchments/")

# Collect all variable path vectors in a named list.
# Naming makes downstream error messages and indexing far more readable.
all.var <- list(
        slope                   = slope,
        lakes                   = lakes,
        soil_oc                 = soil_oc,
        soil_pH                 = soil_pH,
        floodpl                 = floodpl,
        streamp                 = streamp,
        geology                 = geology,
        glacier                 = glacier,
        bioclim                 = bioclim,
        erosion                 = erosion,
        roughness               = roughness,
        elevation               = elevation,
        snowdepth               = snowdepth,
        discharge               = discharge,
        groundwater             = groundwater,
        upstream_area           = upstream_area,
        saturated_water_content = saturated_water_content,
        valley_bottom_flatness  = valley_bottom_flatness
)

n_vars <- length(all.var)   # 18 — used throughout to avoid magic numbers

# helper: drop tiles whose filename matches a given pattern ---------------
rmv <- function(x, y) {
        # Remove any path in x whose filename contains the string y.
        if (any(str_detect(x, y))) {
                x <- x[!str_detect(x, y)]
        }
        x
}

# Exclude geographic regions not in scope for this project
all.var %<>% lapply(rmv, "Hondo")
all.var %<>% lapply(rmv, "Iceland")

# Sanity check: every variable must have the same number of tiles.
tile_counts <- sapply(all.var, length)
if (!all(tile_counts == tile_counts[1])) {
        stop(
                "Unequal tile counts across variables:\n",
                paste(names(tile_counts), tile_counts, sep = " = ", collapse = "\n")
        )
}
n_tiles <- tile_counts[1]

# Extract clean catchment names from the first variable's file list.
# The regex captures the bare filename (letters only) before ".parquet".
clean_name_list <- all.var[[1]] %>%
        str_extract("(?i)[a-z]+\\.parquet") %>%   
        str_remove("\\.parquet")


# Ensure output directory exists before the loop writes into it
out_dir <- "data/catchments"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# main loop: one iteration per catchment tile -----------------------------
for (i in seq_len(n_tiles)) {
        
        # --- derive the expected catchment name from variable 1 ----------------
        i.catchment_name <- all.var[[1]][i] %>%
                str_extract("(?i)[a-z]+\\.parquet") %>%
                str_remove("\\.parquet")
        
        message(i, " / ", n_tiles, " — ", i.catchment_name)
        
        # --- alignment check ---------------------------------------------------
        # Verify that every variable's i-th file belongs to the same catchment.
        all_match <- sapply(all.var, function(var) str_detect(var[i], i.catchment_name))
        if (!all(all_match)) {
                mismatched <- names(all.var)[!all_match]
                stop(
                        "Catchment name mismatch at tile ", i, " (", i.catchment_name, ") ",
                        "for variables: ", paste(mismatched, collapse = ", ")
                )
        }
        
        # --- read all 18 parquet files for this tile ---------------------------
        i.f.all <- lapply(all.var, function(var) read_parquet(var[i]))
        
        # --- drop redundant ID2 column if present in any table -----------------
        i.f.all <- lapply(i.f.all, function(dt) {
                # read_parquet returns an Arrow Table or data.frame; setDT() ensures
                # data.table semantics for the := operator used below.
                setDT(dt)
                if ("ID2" %in% names(dt)) dt[, ID2 := NULL]
                dt
        })
        
        # --- merge all tables on shared catchment segment ID -------------------
        # Reduce applies merge() iteratively: merge(merge(f1, f2), f3) …
        # all = TRUE is a full outer join, retaining all IDs even if a predictor
        # is missing for some segments in a given tile.
        i.f.all2 <- Reduce(
                function(...) merge(..., by = "ID", all = TRUE),
                i.f.all
        )
        
        # Drop columns that are not needed in the final output.
        setDT(i.f.all2)
        
        # Remove known redundant columns only if they are actually present,
        # to avoid errors on tiles where these columns were already absent.
        drop_cols <- intersect(c("VBM_min", "spi_min"), names(i.f.all2))
        if (length(drop_cols) > 0) i.f.all2[, (drop_cols) := NULL]
        
        # --- write output ------------------------------------------------------
        write_parquet(
                i.f.all2,
                file.path(out_dir, paste0(i.catchment_name, "_w_variables.parquet"))
        )
        
        # Free memory before the next tile
        rm(i.f.all, i.f.all2)
}

rm(i)
gc()
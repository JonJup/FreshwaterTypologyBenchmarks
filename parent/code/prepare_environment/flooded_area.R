# ################################################################################
# # Script name:  flooded_area.R
# #
# # Purpose:      This script prepares the soil organic carbon data
# #
# # Inputs:       - EU-Hydro DEM Parquet (one per River Basin District); 
# #                       Created with: euhydro_to_parquet.R
#                 - Europe_RP10_filled_depth.tif
# #                 River flood hazard maps for Europe and the Mediterranean Basin region
#                   This file specifies values for a flood with a ten year return period. (RP10).
#                   https://data.jrc.ec.europa.eu/dataset/1d128b6c-a4ee-4858-9e34-6210707f3c81-
#                   Cell values indicate water depth (in m).
#                   Downloaded: 25.11.2024
# #
# # Output:       - flood_area/processedData/w_catchment//*
#                       parquet files with 
#                               flooded_area: fraction of catchment flooded 
#                               inundation_depth: depth of inundation
# ################################################################################


# setup -------------------------------------------------------------------
library(terra)
library(sf)
library(sfarrow)
library(data.table)
library(dplyr)

out.path <- "E:Arbeit/data/river_network/flood_area/processedData/w_catchment/"
dir.create(out.path, showWarnings = FALSE, recursive = T)

vector.files <- fs::dir_ls("E://Arbeit/Data/river_network/eu_hydro_dem/processedData/parquet_catchments/")

# load data ---------------------------------------------------------------
flood <- rast("E://Arbeit/Data/river_network/flood_area/Europe_RP10_filled_depth.tif")

for (i in seq_along(vector.files)){
        
        i.cat <- st_read_parquet(dsn = vector.files[i])
        i.name <- stringr::str_remove(i.cat$ID[1], "1")
        i.cat$ID2 <- 1:nrow(i.cat)
        i.cat2 <- 
                i.cat %>%
                st_drop_geometry %>%
                setDT
        i.ext <- i.cat %>%
                st_bbox() %>% 
                st_as_sfc() %>% 
                st_as_sf(crs = st_crs(i.cat)) %>% 
                st_buffer(10000)
        
        i.raster <- crop(flood, i.ext)
        i.extr <-
                terra::extract(
                        x = i.raster, 
                        y = i.cat, 
                        touch = TRUE
                )
        setDT(i.extr)
        names(i.extr)[2] <- "raster"
        
        i.extr[, flooded_area := sum(!is.na(raster)), by = "ID"]
        i.extr[, flooded_area2 := .N, by = "ID"]
        i.extr[flooded_area != 0, flooded_area3 := flooded_area / flooded_area2]
        i.extr[flooded_area == 0, flooded_area3 := 0]
        i.extr[, inundation_depth := mean(raster, na.rm = T), by = "ID"]
        i.extr[is.nan(inundation_depth), inundation_depth := 0]
        i.extr <- unique(i.extr, by = "ID") 
        
        i.cat2$flooded_area <- i.extr$flooded_area3
        i.cat2$inundation_depth <- i.extr$inundation_depth
        
        i.uri <- paste0(out.path,i.name,".parquet")
        arrow::write_parquet(
                x = i.cat2, 
                sink = i.uri
        )
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
}

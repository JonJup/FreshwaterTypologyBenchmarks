# ################################################################################
# # Script name:  valleyBottomFlatness.R
# #
# # Purpose:      This script summarizes the Valley Bottom Flatnessfor each
#                 EU Hydro DEM catchment.
#                 
# # Inputs:       - EU-Hydro DEM Parquet (one per River Basin District); 
# #                       Created with: euhydro_to_parquet.R.
#               Valley Bottom Flattness:
#                available from https://zenodo.org/records/4495449
#                Continental Europe Digital Terrain Model geomorphometry derivatives at 30 m, 100 m and 250 m
# # Output:       
# ################################################################################

# setup -------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(sfarrow)
library(terra)
library(sf)


# load data ---------------------------------------------------------------
vf <- fs::dir_ls("D://Arbeit/data/river_network/eu_hydro_dem/parquet/")
r <- rast("D://Arbeit/Data/river_network/valley_bottom_flattness/rawData/dtm_vbf_gedi.saga.gis_m_100m_0..0cm_2000..2018_eumap_epsg3035_v0.2.tif")


# loop over catchments ----------------------------------------------------
for (i in 1:length(vf)){
        print(i)
        i.cat  <- st_read_parquet(vf[i])
        i.cat  <- st_transform(i.cat, crs(r))
        i.name <- stringr::str_remove(i.cat$ID[1], "1")
        i.cat2 <- 
                i.cat %>%
                st_drop_geometry %>%
                setDT
        i.ext <- i.cat %>%
                st_bbox() %>% 
                st_as_sfc() %>% 
                st_as_sf(crs = st_crs(i.cat))
        
        i.raster <- crop(r, i.ext)
        stats <- exact_extract(i.raster, i.cat, c("mean", "min", "max"))
        i.cat2$VBM_mean <- stats$mean
        i.cat2$VBM_min  <- stats$min
        i.cat2$VBM_max  <- stats$max

        i.uri <- paste0("E:Arbeit/data/river_network/valley_bottom_flattness/processedData/valleyBottomFlattness_w_catchments/",i.name,".parquet")
        arrow::write_parquet(
                x = i.cat2, 
                sink = i.uri
        )
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
}

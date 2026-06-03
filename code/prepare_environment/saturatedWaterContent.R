# ################################################################################
# # Script name:  saturatedSoilWaterContent.R
# #
# # Purpose:      This script summarizes the Saturated Soil Water content for 
#                 each EU Hydro DEM catchment.
#                 
# # Inputs:       - EU-Hydro DEM Parquet (one per River Basin District); 
# #                       Created with: euhydro_to_parquet.R
#                 - Saturated Soil water content/ field capacity. 
#                       - fc_fao.tif
#                       - Available from: 
#                       https://esdac.jrc.ec.europa.eu/content/maps-indicators-soil-hydraulic-properties-europe
#                       Based on Toth, B., Weynants, M., Nemes, A., Mako, A., 
#                       Bilas, G., Toth, G., 2014. New generation of hydraulic 
#                       pedotransfer functions for Europe. European Journal of 
#                       Soil Science
#                       - 
# # Output:       
# ################################################################################


# setup -------------------------------------------------------------------
library(terra)
library(sf)
library(sfarrow)
library(data.table)
library(dplyr)
library(exactextractr)

# load data ---------------------------------------------------------------
sswc <- rast("E://Arbeit/Data/river_network/Maps of indicators of soil hydraulic properties for Europe/fc_fao.tif")
vector.files <- fs::dir_ls("E://Arbeit/Data/river_network/eu_hydro_dem/parquet")


# loop over catchments ----------------------------------------------------
for (i in 1:length(vector.files)){
          print(i)
        i.cat  <- st_read_parquet(dsn = vector.files[i])
        i.name <- stringr::str_remove(i.cat$ID[1], "1")
        i.cat$ID2 <- 1:nrow(i.cat)
        i.cat <- st_transform(i.cat, crs(sswc))
        i.cat2 <- 
                i.cat %>%
                st_drop_geometry %>%
                setDT
        i.ext <- i.cat %>%
                st_bbox() %>% 
                st_as_sfc() %>% 
                st_as_sf(crs = st_crs(i.cat)) %>% 
                st_buffer(10000)
        
        i.raster <- crop(sswc, i.ext)
        i.extr <- exact_extract(sswc, i.ext, fun = "mean")
        i.extr[, sswc := round(i.extr, 3)]
        i.uri <- paste0("E:Arbeit/Data/river_network/Maps of indicators of soil hydraulic properties for Europe/w_catchments/",i.name,".parquet")
        arrow::write_parquet(
                x = i.cat2, 
                sink = i.uri
        )
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
}

# ################################################################################
# # Script name:  glacier.R
# #
# # Purpose:      This script prepares the glacier data from the Randolph 
#                 Glacier Inventory for the EU Hydro DEM catchments. 
#                 
# # Inputs:       - EU-Hydro DEM Parquet (one per River Basin District); 
# #                       Created with: euhydro_to_parquet.R
#                 - rgi60_global.shp
#                       Randolph Glacier Inventory
#                       Available under: 
#                               https://cds.climate.copernicus.eu/cdsapp#!/dataset/insitu-glaciers-extent?tab=overview
#                       We used the version 6.
# #
# # Output:       - oc_w_catchment/*
#                       parquet files with mean topsoil organic carbon fraction 
#                       in catchment. 
# ################################################################################


# setup -------------------------------------------------------------------
library(sf)
library(sfarrow)
library(data.table)
library(tidyverse)
library(terra)
library(mapview)
library(magrittr)

# load data ---------------------------------------------------------------
vf <- fs::dir_ls("E://Arbeit/Data/river_network/eu_hydro_dem/processedData/parquet_catchments//")
g <- vect("E://Arbeit/Data/randolph glacier inventory/rawData/rgi60_global.shp")
r <- rast("E://Arbeit/Data/river_network/valley_bottom_flattness/rawData/dtm_vbf_gedi.saga.gis_m_100m_0..0cm_2000..2018_eumap_epsg3035_v0.2.tif")


g <- project(g, r)
g <- crop(g, r)
g <- st_as_sf(g)
# g <- rasterize(g, r)


for (i in 1:length(vf)){

        i.v <- st_read_parquet(vf[i])
        i.v %<>% st_transform(crs(g))
        i.n <- stringr::str_remove(i.v$ID[1], "1")
        i.v$ID2 <- 1:nrow(i.v)
        i.v2 <- 
                i.v %>%
                st_drop_geometry %>%
                setDT
        i.v2[, glacial_area := 0]
        i.ext <- i.v %>%
                st_bbox() %>% 
                st_as_sfc() %>% 
                st_as_sf(crs = st_crs(i.v))
        i.g  <- st_filter(g, i.ext)
        i.i  <- st_intersects(i.v, i.g)
        i.i2 <- which(lengths(i.i) != 0)

        for (j in i.i2){
                j.v <- i.v[j,]
                j.v2 <- st_intersection(j.v,i.g)
                i.v2[j, glacial_area := units::drop_units(sum(st_area(j.v2))/st_area(j.v))]
                rm(j.v, j.v2)
        };rm(j)
        i.uri <- paste0("E:Arbeit/Data/randolph glacier inventory/processedData/w_catchments/",i.n,".parquet")
        arrow::write_parquet(
                x = i.v2, 
                sink = i.uri
        )
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
}

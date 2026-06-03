# ################################################################################
# # Script name:  soil_pH.R
# #
# # Purpose:      This script prepares the soil pH data
# #                 https://stac.ecodatacube.eu/ph.h2o_iso.10390.2021.index/collection.json
# #
# # Inputs:       - EU-Hydro DEM Parquet (one per River Basin District); 
# #                       Created with: euhydro_to_parquet.R
#                 - sol_ph.h2o_lucas.iso.10694_m_30m_s0..0cm_2020_eumap_epsg3035_v0.2.tif
# # Output:       - soil_ph_europe/processedData/w_catchments/*
#                       parquet files with mean soi pH 
#                       in catchment. 
# ################################################################################
# 
# setup -------------------------------------------------------------------

library(terra)
library(sf)
library(sfarrow)
library(data.table)
library(dplyr)


out.path <- "E:/Arbeit/Data/soil/soil_ph_europe/processedData/w_catchments/"
dir.create(out.path, showWarnings = FALSE, recursive = T)
vector.files <- fs::dir_ls("E://Arbeit/Data/river_network/eu_hydro_dem/processedData/parquet_catchments")
# data input --------------------------------------------------------------
soil <- rast("E://Arbeit/Data/soil/soil_ph_europe/rawData/Europe soil ph/sol_ph.h2o_lucas.iso.10694_m_30m_s0..0cm_2020_eumap_epsg3035_v0.2.tif")

soil <- terra::aggregate(soil, fact = 2)

for (i in seq_along(files)){
        print(i)
        i.cat <- st_read_parquet(dsn = vector.files[i])
        i.name <- stringr::str_remove(i.cat$ID[1], "1")
        i.cat <- filter(i.cat, !ID %in% dropv)
        i.cat <- st_transform(i.cat, crs = crs(soil))

        i.cat$ID2 <- 1:nrow(i.cat)
        i.cat2 <- 
                i.cat %>%
                st_drop_geometry %>%
                setDT
        i.ext <- i.cat %>%
                st_bbox() %>% 
                st_as_sfc() %>% 
                st_as_sf(crs = st_crs(i.cat))%>%
                st_buffer(30000)
        
        i.raster <- crop(soil, i.ext)
        i.extr <-
                terra::extract(
                        x = i.raster, 
                        y = i.cat,
                        touches = T
        )
        setDT(i.extr)
        names(i.extr)[2] <- "raster"
        i.extr[, soil_pH := mean(raster, na.rm = T), by = "ID"]
        i.extr <- unique(i.extr, by = "ID") 
        i.cat2$soil_pH <- i.extr$soil_pH
        
        ids <- filter(i.cat2, is.na(soil_pH)) %>% pull (ID)
        for (k in seq_along(ids)){
                
                k.cat <- i.cat %>% filter(ID == ids[k])
                k.cat <- st_buffer(k.cat, 2000)
                k.ext <- terra::extract(i.raster, k.cat, touches = TRUE)
                k.ext <- mean(k.ext[,2], na.rm = T)
                i.cat2[ID == ids[k], soil_pH := k.ext]
                rm(list = ls()[grepl(pattern = "^k\\.", x = ls())])
        }
        rm(k)
        rm(ids)
        
        i.uri <- paste0(out.path,i.name,".parquet")
        arrow::write_parquet(
                x = i.cat2, 
                sink = i.uri
        )
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
}

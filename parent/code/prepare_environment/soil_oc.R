# ################################################################################
# # Script name:  soil_oc.R
# #
# # Purpose:      This script prepares the soil organic carbon data
#                 Panagos, P. et al. (2022). European Soil Data Centre 2.0: Soil
#                 data and knowledge in support of the EU policies. European 
#                 Journal of Soil Science, 73(6), e13315. 
#                 DOI: 10.1111/ejss.13315
#                  
#                 Available upon request from https://esdac.jrc.ec.europa.eu/
# #
# # Inputs:       - EU-Hydro DEM Parquet (one per River Basin District); 
# #                       Created with: euhydro_to_parquet.R
# #               - OC_TOP.tif Soil Organic Carbon Raster file
#                       OC_TOP = Topsoil organic carbon content.
#                       H = High ( > 6 %)
#                       M = Medium (2 - 6 %)
#                       L = Low (1 - 2 %)
#                       V = Very low ( < 1 %)
# #
# # Output:       - oc_w_catchment/*
#                       parquet files with mean topsoil organic carbon fraction 
#                       in catchment. 
# ################################################################################


# setup -------------------------------------------------------------------
library(terra)
library(sf)
library(sfarrow)
library(data.table)
library(dplyr)

vector.files <- fs::dir_ls("E://Arbeit/Data/river_network/eu_hydro_dem/parquet")
out.path = "E://Arbeit/data/soil/new soil/processedData/oc_w_catchment/"
dir.create(out.path, showWarnings = FALSE, recursive = T)

# load data ---------------------------------------------------------------
soil <- rast("E://Arbeit/Data/soil/new soil/ESDB-Raster-Library-1k-GeoTIFF-20240507/ESDB-Raster-Library-1k-GeoTIFF-20240507/OC_TOP/OC_TOP.tif")

translation_table <- 
        data.table(
                OC_TOP = c(1,2,3,4,5), 
                OC_TOP_perc = c(6, 4,1.5, 0.5, NA)
        )

for (i in seq_along(vector.files)){

        i.cat <- st_read_parquet(dsn = vector.files[i])
        i.cat <- st_transform(i.cat, crs = crs(soil))
        i.name <- stringr::str_remove(i.cat$ID[1], "1")
        i.cat$ID2 <- 1:nrow(i.cat)
        i.cat2 <- 
                i.cat %>%
                st_drop_geometry %>%
                setDT
        i.ext <- i.cat %>%
                st_bbox() %>% 
                st_as_sfc() %>% 
                st_as_sf(crs = st_crs(i.cat))%>%
                st_buffer(10000)
        
        i.raster <- crop(soil, i.ext)
        i.extr <-
                terra::extract(
                        x = i.raster, 
                        y = i.cat,
                        touches = T
        )
        setDT(i.extr)
        i.extr2 <- translation_table[i.extr, on = "OC_TOP"]
        names(i.extr2)[2] <- "raster"
        i.extr2[, soil_oc:=mean(raster, na.rm = T), by = "ID"]
        i.extr2 <- unique(i.extr2, by = "ID") 
        i.cat2$soil_oc <- i.extr2$soil_oc
        ids <- filter(i.cat2, is.na(soil_oc)) %>% pull (ID)
        for (k in seq_along(ids)){
                
                k.cat <- i.cat %>% filter(ID == ids[k])
                k.cat <- st_buffer(k.cat, 2000)
                k.ext <- terra::extract(i.raster, k.cat, touches = TRUE)
                k.ext <- mean(k.ext$last, na.rm = T)
                i.cat2[ID == ids[k], soil_oc := k.ext]
                rm(list = ls()[grepl(pattern = "^k\\.", x = ls())])
        }
        

        i.uri <- paste0(out.path,i.name,".parquet")
        arrow::write_parquet(
                x = i.cat2, 
                sink = i.uri
        )
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
}

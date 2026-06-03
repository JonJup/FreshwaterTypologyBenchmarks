# ################################################################################
# # Script name:  erosivity.R
# #
# # Purpose:      This script prepares the Erosivity for EU Hydro DEM
#                 catchments. 
#                 
# # Inputs:       - EU-Hydro DEM Parquet (one per River Basin District); 
# #                       Created with: euhydro_to_parquet.R
#                 - Global Rainfall Erosivity (GLOREDA). One Raster (.tiff) file
#                   per month of the year. 
#                   Data available from: 
#                        https://esdac.jrc.ec.europa.eu/content/global-rainfall-erosivity
# #
# # Output:       - 
# ################################################################################


# setup -------------------------------------------------------------------
library(terra)
library(sf)
library(sfarrow)
library(data.table)
library(dplyr)
library(exactextractr)

# load data ---------------------------------------------------------------
rast.files <- fs::dir_ls("E://Arbeit/Data/Gloreda/rawData/", regexp = "tif")
vector.files <- fs::dir_ls("E://Arbeit/Data/river_network/eu_hydro_dem/processedData/parquet_catchments/")


# produce catchment-specific output ---------------------------------------
for (i in 1:length(vector.files)){
        print(paste("Starting vector", i))
        i.cat  <- st_read_parquet(dsn = vector.files[i])
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
        for (j in 1:length(rast.files)){
                print(paste("Starting raster", j, "in vector", i))
                j.rast <- rast(rast.files[j])
                j.raster.name <- stringr::str_remove(names(j.rast), "Rfactor_")
                j.rast <- crop(j.rast, i.ext)
                j.stats <- exact_extract(j.rast, i.cat, c("mean"))
                i.cat2[,j.raster.name] <- j.stats
                rm(list = ls()[grepl(pattern = "^j\\.", x = ls())])
        }
        i.cat2[, Rfactor_avg := rowMeans(.SD, na.rm = TRUE), .SDcols = -c(1,2)]
        i.cat2[, Rfactor_max := apply(.SD, 1, max, na.rm = T), .SDcols = -c(1,2)]
        i.cat2[, Rfactor_min := apply(.SD, 1, min, na.rm = T), .SDcols = -c(1,2)]
        i.cat2[!is.finite(Rfactor_max), Rfactor_max := NA]
        i.cat2[!is.finite(Rfactor_min), Rfactor_min := NA]
        i.cat2 <- i.cat2[, c("ID", "Rfactor_max", "Rfactor_avg", "Rfactor_min")]
        i.uri <- paste0("E:Arbeit/Data/Gloreda/processedData/w_catchements/",i.name,".parquet")
        arrow::write_parquet(
                x = i.cat2, 
                sink = i.uri
        )
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
}

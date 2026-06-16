# ################################################################################
# # Purpose:      This script prepares the BioClimatic Variables from CHELSA for
#                 EU Hydro DEM catchments
#                 
# # Inputs:       - EU-Hydro DEM Parquet (one per River Basin District); 
# #                       Created with: euhydro_to_parquet.R
#                 - CHELSA BIOCLIM variables 5,6,13,14
#                 FILES: 
#                       CHELSA_bio5_1981-2010_V.2.1.tif
#                       CHELSA_bio6_1981-2010_V.2.1.tif
#                       CHELSA_bio13_1981-2010_V.2.1.tif
#                       CHELSA_bio14_1981-2010_V.2.1.tif
#               Files are available under: 
#                       https://www.chelsa-climate.org/datasets/chelsa_bioclim
# #
# # Output:       
# ################################################################################

# library -----------------------------------------------------------------
library(terra)
library(sf)
library(sfarrow)
library(data.table)
library(dplyr)


# load data ---------------------------------------------------------------
bioclim05 <- rast("E://Arbeit/Data/climate/bioclim_chelsa/rawData/V2.1/CHELSA_bio5_1981-2010_V.2.1.tif")
bioclim06 <- rast("E://Arbeit/Data/climate/bioclim_chelsa/rawData/V2.1/CHELSA_bio6_1981-2010_V.2.1.tif")
bioclim13 <- rast("E://Arbeit/Data/climate/bioclim_chelsa/rawData/V2.1/CHELSA_bio13_1981-2010_V.2.1.tif")
bioclim14 <- rast("E://Arbeit/Data/climate/bioclim_chelsa/rawData/V2.1/CHELSA_bio14_1981-2010_V.2.1.tif")
vector.files <- fs::dir_ls("E://Arbeit/Data/river_network/eu_hydro_dem/processedData/parquet_catchments/")

for (i in 2:length(vector.files)){

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
                st_as_sf(crs = st_crs(i.cat))
        
        i.05 <- crop(bioclim05, i.ext)
        i.06 <- crop(bioclim06, i.ext)
        i.13 <- crop(bioclim13, i.ext)
        i.14 <- crop(bioclim14, i.ext)
        i.extr.05 <- terra::extract(x = i.05,y = i.cat)
        i.extr.06 <- terra::extract(x = i.06,y = i.cat)
        i.extr.13 <- terra::extract(x = i.13,y = i.cat)
        i.extr.14 <- terra::extract(x = i.14,y = i.cat)
        
        setDT(i.extr.05)
        setDT(i.extr.06)
        setDT(i.extr.13)
        setDT(i.extr.14)

        names(i.extr.05)[2] <- "bc"
        names(i.extr.06)[2] <- "bc"
        names(i.extr.13)[2] <- "bc"
        names(i.extr.14)[2] <- "bc"

        i.extr.05[, bc05:=mean(bc, na.rm = T), by = "ID"]
        i.extr.06[, bc06:=mean(bc, na.rm = T), by = "ID"]
        i.extr.13[, bc13:=mean(bc, na.rm = T), by = "ID"]
        i.extr.14[, bc14:=mean(bc, na.rm = T), by = "ID"]


        i.extr.05 <- unique(i.extr.05, by = "ID") 
        i.extr.06 <- unique(i.extr.06, by = "ID") 
        i.extr.13 <- unique(i.extr.13, by = "ID") 
        i.extr.14 <- unique(i.extr.14, by = "ID") 
        
        i.cat2$bioclim05 <- i.extr.05$bc05
        i.cat2$bioclim06 <- i.extr.06$bc06
        i.cat2$bioclim13 <- i.extr.13$bc13
        i.cat2$bioclim14 <- i.extr.14$bc14
        
        i.uri <- paste0("E:Arbeit/Data/climate/bioclim_chelsa/processed_data/w_catchments/",i.name,".parquet")
        arrow::write_parquet(
                x = i.cat2, 
                sink = i.uri
        )
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
}

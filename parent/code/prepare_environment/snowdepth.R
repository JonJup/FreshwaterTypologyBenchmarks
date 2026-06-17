# ################################################################################
# # Script name:  snowdepth.R
# #
# # Purpose:      This script prepares the snow depth data for EU Hydro DEM
#                 catchments. 
#                 
# # Inputs:       - EU-Hydro DEM Parquet (one per River Basin District); 
# #                       Created with: euhydro_to_parquet.R
#                 - Snow Depth data from CERRA. 
#                 Available from: https://cds.climate.copernicus.eu/datasets/reanalysis-cerra-land?tab=overview
#                 The variable is called: Snow depth water equivalent 
#                 We used data for February 2021 
#                 The values are provided in kg per square meter 
#                 "Snow depth water equivalent expresses the snow depth in kg of
#                  snow over one square metre. The unit corresponds to 1 mm of 
#                  water equivalent. It is given as the mean for the grid area. 
#                  Given values are instantaneous.
# #
# # Output:       - 
# ################################################################################

# setup -------------------------------------------------------------------
library(terra)
library(sf)
library(sfarrow)
library(data.table)
library(ncdf4)
library(tidyverse)
library(exactextractr)

# load original data ------------------------------------------------------
nc1 <- nc_open("E:/Arbeit/Data/CERRA/rawData/snow_depth_water_equivalent.nc")
nc.sf <- st_as_sf(
        data.frame(
                lon = c(ncvar_get(nc1, "longitude")), 
                lat = c(ncvar_get(nc1, "latitude")),
                sd  = c(ncvar_get(nc1, "sd"))
        ), 
        coords = c("lon", "lat"),
        crs = 4326
)
nc.sf <- dplyr::filter(nc.sf, !is.na(sd))
st_write_parquet(nc.sf, "E:Arbeit/Data/CERRA/processedData/snow_depth_water_equivalent.parquet")

# load processed data -----------------------------------------------------
d1 <- st_read_parquet("E:Arbeit/Data/CERRA/processedData/snow_depth_water_equivalent.parquet")
vector.files <- fs::dir_ls("E://Arbeit/Data/river_network/eu_hydro_dem/processedData/parquet_catchments/")

for (i in 1:length(vector.files)){
        
        print(i)
        
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
        
        i.d1 <- st_filter(d1, i.ext)
        i.d1   <- vect(i.d1)
        i.ext1 <- ext(i.d1)
        r1 <- rast(i.ext1 , resolution = 0.099)
        r1 <- rasterize(i.d1, r1, field = "sd")
        i.stats <- exact_extract(r1, i.cat, fun = "mean") %>% 
                round(2)
        i.cat2$mean_snow_equivalent <- i.stats
        i.uri <- paste0("E:Arbeit/Data/CERRA/processedData/snow_depth_w_catchments/",i.name,".parquet")
        arrow::write_parquet(
                x = i.cat2, 
                sink = i.uri
        )
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
}
rm(i)
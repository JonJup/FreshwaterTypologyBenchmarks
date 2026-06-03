# ################################################################################
# # Script name:  roughness.R
# #
# # Purpose:      This script prepares the stream power index for Eu Hydro DEM
#                 catchments. 
#                 
# # Inputs:       - EU-Hydro DEM Parquet (one per River Basin District); 
# #                       Created with: euhydro_to_parquet.R
#                 - DEM: Continental Europe Digital Terrain Model at 30 m 
#                 resolution based on GEDI, ICESat-2, AW3D, GLO-30, EUDEM, MERIT
#                  DEM and background layers
#                 Available from: 
#                 https://zenodo.org/records/4724549
#                 DOI: https://doi.org/10.5281/zenodo.4724549
# #
# # Output:       - 
# ################################################################################

# setup ---------------------------------------------------------------
library(terra)
library(sf)
library(sfarrow)
library(data.table)
library(dplyr)
library(exactextractr)

# load data ---------------------------------------------------------------
dtm          <- rast("E://Arbeit/Data/DEM/DTM_Europe/rawData/dtm_elev.lowestmode_gedi.eml_mf_30m_0..0cm_2000..2018_eumap_epsg3035_v0.3.tif")
vector.files <- fs::dir_ls("E://Arbeit/Data/river_network/eu_hydro_dem/processedData/parquet_catchments/")

# prepare data ------------------------------------------------------------
dtm.agg  <- terra::aggregate(dtm, fact = 3)
dtm.agg2 <- terra::terrain(x = dtm.agg, v="roughness")

# loop over catchments ----------------------------------------------------
for (i in 1:length(vector.files)){
        
        i.cat     <- st_read_parquet(dsn = vector.files[i])
        i.cat     <- st_transform(i.cat, crs = crs(dtm.agg2))
        i.name    <- stringr::str_remove(i.cat$ID[1], "1")
        i.cat$ID2 <- 1:nrow(i.cat)
        i.cat2 <- 
                i.cat %>%
                st_drop_geometry %>%
                setDT
        i.ext <- i.cat %>%
                st_bbox() %>% 
                st_as_sfc() %>% 
                st_as_sf(crs = st_crs(i.cat))
        
        i.raster <- crop(dtm.agg2, i.ext)
        i.stats <- exact_extract(i.raster, i.cat, fun = "mean")
        i.stats <- round(i.stats, 2)
        i.cat2$roughness <- i.stats
        i.uri <- paste0("E:Arbeit/Data/DEM/DTM_Europe/processedData/roughness_w_catchment/",i.name,".parquet")
        arrow::write_parquet(
                x = i.cat2, 
                sink = i.uri
        )
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
        rm(i)
}

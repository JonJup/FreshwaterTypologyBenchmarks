# ################################################################################
# # Script name:  upstream.R
# #
# # Purpose:      This script summarizes the upstream catchment area for each
#                 EU Hydro DEM catchment.
#                 
# # Inputs:       - EU-Hydro DEM Parquet (one per River Basin District); 
# #                       Created with: euhydro_to_parquet.R
#                 - Groundwater table depth [cm]. Obtained through BasinATLAS
#                 and originally from Fan et al (2013). 
#               Files are available under: 
#                       https://www.hydrosheds.org/products
# #             Fan, Y., Li, H., & Miguez-Macho, G. (2013). Global patterns of 
#               groundwater table depth. Science, 339(6122), 940-943.
#               - Valley Bottom Flattness: available from https://zenodo.org/records/4495449
# # Output:       
# ################################################################################

# setup -------------------------------------------------------------------
library(terra)
library(sf)
library(sfarrow)
library(data.table)
library(dplyr)
library(magrittr)


# load data -------------------------------------------------------------------------
vf <- fs::dir_ls("D://Arbeit/data/river_network/eu_hydro_dem/parquet")
ca       <- vect("D://Arbeit/data/river_network/BasinATLAS_Data_v10/rawData/BasinATLASeurope_manual.gpkg")
r        <- rast("D://Arbeit/data/river_network/valley_bottom_flattness/rawData/dtm_vbf_gedi.saga.gis_m_100m_0..0cm_2000..2018_eumap_epsg3035_v0.2.tif")


# prepare data ------------------------------------------------------------
ca <- project(x = ca, y = r)
ca <- rasterize(x = ca, y = r, field = "UP_AREA", fun = "mean")

rm(r);gc()

# loop over catchments ----------------------------------------------------
out <- c()
for (i in 1:length(vf)){
        i.v <- st_read_parquet(vf[i])
        i.v %<>% st_transform(crs(ca))
        i.n <- stringr::str_remove(i.v$ID[1], "1")
        i.v$ID2 <- 1:nrow(i.v)
        i.v2 <- 
                i.v %>%
                st_drop_geometry %>%
                setDT
        i.ext <- i.v %>%
                st_bbox() %>% 
                st_as_sfc() %>% 
                st_as_sf(crs = st_crs(i.v))
        
        i.r <- crop(ca, i.ext)
        i.extr <- exact_extract(i.r, i.v, fun = "mean")
        i.v2$groundwater_table <- round(i.extr,3)
        i.uri <- paste0("E:Arbeit/Data/river_network/BasinATLAS_Data_v10/processedData/upstream_area_w_catchment/",i.n,".parquet")
        arrow::write_parquet(
                x = i.v2, 
                sink = i.uri
        )
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
        rm(i)
}

saveRDS(out, "D://Arbeit/data/river_network/BasinATLAS_Data_v10/processedData/upstream_hydroDEM_mad.rds")

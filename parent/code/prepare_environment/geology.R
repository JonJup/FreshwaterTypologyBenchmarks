# ################################################################################
# # Script name:  geology.R
# #
# # Purpose:      This script computes the the compositional geology variables
#                 `area_calcareous`, "area_siliceous* and 'area_sediment' for EU
#                  Hydro DEM catchments. 
#                 
# # Inputs:       - EU-Hydro DEM Parquet (one per River Basin District); 
# #                       Created with: euhydro_to_parquet.R
#                 - ihme_1500_litho4changed.shp
#                       A Modified version of the IHME1500 data set.
#                       The original is available here: 
#                               https://www.bgr.bund.de/EN/Themen/Grundwasser/Projekte/Flaechen-Rauminformationen/Ihme1500/ihme1500.html
#                       This modified version was provided by Florian Brogwardt (Pletterbauer).
#                       It can be provided upon request.
#                       Also note that the above link will provide you with v12, while we used v11. 
# #
# # Output:       -  
# ################################################################################


# setup -------------------------------------------------------------------
library(sf)
library(sfarrow)
library(data.table)
library(tidyverse)
library(mapview)
library(terra)
library(magrittr)
library(units)


# load data ---------------------------------------------------------------
vf <- fs::dir_ls("E://Arbeit/Data/river_network/eu_hydro_dem/processedData/parquet_catchments//")
g  <- st_read("E://Arbeit/Data/IHME1500_v11/ihme_1500_litho4changed.shp")

for (i in seq_along(vf)){
#for (i in 2){
        print(paste("i == ", i))
        i.v     <- st_read_parquet(vf[i])
        i.v     %<>% st_transform(st_crs(g))
        i.n     <- stringr::str_remove(i.v$ID[1], "1")
        i.v$ID2 <- 1:nrow(i.v)
        i.v2    <- 
                i.v %>%
                st_drop_geometry %>%
                setDT
        i.v2[, `:=` (area_calcareous = 0,
                     area_siliceous  = 0,
                     area_sediment   = 0)]
        i.ext <- i.v %>%
                st_bbox() %>% 
                st_as_sfc() %>% 
                st_as_sf(crs = st_crs(i.v))%>%
                st_buffer(10000)
        
        i.g  <- st_filter(g, i.ext)
        i.i  <- st_intersects(i.v, i.g)
        i.i2 <- lengths(i.i)
        i.i2 <- which(lengths(i.i) != 0)
        
        for (j in i.i2){
                #if (i == 1 & j < i.i2[5338]) next()
                print(paste(which(i.i2 == j), "of", length(i.i2)))
                j.v  <- i.v[j,]
                j.v2 <- st_intersection(j.v,i.g)
                j.ta <- sum(drop_units(st_area(j.v2)))
                
                j.calc <- j.v2 %>%
                        filter(acid_basic == "calcareous") %>%
                        st_area %>%
                        drop_units %>%
                        {
                                \(x) x / j.ta
                        }() %>% 
                        sum
                
                j.sili <- j.v2 %>%
                        filter(acid_basic == "siliceous") %>%
                        st_area %>%
                        drop_units %>%
                        {
                                \(x) x / j.ta
                        }() %>% sum
                
                j.sedi <- j.v2 %>%
                        filter(acid_basic == "sediments") %>%
                        st_area %>%
                        drop_units %>%
                        {
                                \(x) x / j.ta
                        }()  %>%
                        sum()
                if (
                        max(j.calc, 0) == 0 & 
                        max(j.sedi, 0) == 0 & 
                        max(j.sili, 0) == 0
                        ){
                        i.v2[j, `:=` (
                                area_calcareous = NA,
                                area_siliceous  = NA,
                                area_sediment   = NA
                        )]
                } else {
                        i.v2[j, `:=` (
                                area_calcareous = max(j.calc, 0),
                                area_siliceous  = max(j.sili, 0),
                                area_sediment   = max(j.sedi, 0)
                        )]
                }
                rm(j.v, j.v2, j.calc, j.sili, j.sedi, j.ta)
        };rm(j)
        i.uri <- paste0("E:/Arbeit/Data/IHME1500_v11/processedData/w_catchments/",i.n,".parquet")
        arrow::write_parquet(
                x = i.v2, 
                sink = i.uri
        )
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
}

# ————————————————————————————————————————————————————— #
# ——— Combine environmental data at catchment level ——— # 
# ————————————————————————————————————————————————————— #

# Jonathan Jupke (jonjup@protonmail.com)
# 26.11.2024

# setup -----------------------------------------------------------------------------
setwd(rstudioapi::getActiveProject())

library(fs)
library(magrittr)
library(stringr)
library(arrow)
library(data.table)

# load data -------------------------------------------------------------------------
slope                   <- dir_ls("D://Arbeit/data/river_network/hydrography90m/processedData/slope_w_catchments/")
lakes                   <- dir_ls("D:Arbeit/data/river_network/hydroLAKES/processedData/fractionLake_w_catchment//")
soil_oc                 <- dir_ls("D:Arbeit/data/soil/new soil/processedData/oc_w_catchment//")
soil_pH                 <- dir_ls("D:Arbeit/data/soil/soil_ph_europe/processedData/w_catchments/")
floodpl                 <- dir_ls("D:Arbeit/data/river_network/flood_area/processedData/w_catchment/")
streamp                 <- dir_ls("D:Arbeit/data/river_network/hydrography90m/processedData/spi_w_catchment//")
geology                 <- dir_ls("D:Arbeit/data/LULC/IHME1500_v11/w_catchments/")
glacier                 <- dir_ls("D:Arbeit/data/LULC/randolph glacier inventory/w_catchments/")
bioclim                 <- dir_ls("D:Arbeit/data/climate/bioclim_chelsa/processedData/w_catchments/")
erosion                 <- dir_ls("D:Arbeit/data/misc/Gloreda/processedData/w_catchements/")
roughness               <- dir_ls("D:Arbeit/data/DEM/DTM_Europe/processedData/roughness_w_catchment/")
elevation               <- dir_ls("D:Arbeit/data/DEM/DTM_Europe/processedData/elevation_w_catchment/")
snowdepth               <- dir_ls("D:Arbeit/data/CERRA/processedData/snow_depth_w_catchments/")
discharge               <- dir_ls("D:Arbeit/data/river_network/discharge copernicus/processedData/w_catchments//")
groundwater             <- dir_ls("D:Arbeit/data/river_network/BasinATLAS_Data_v10/processedData/groundwater_area_w_catchment//")
upstream_area           <- dir_ls("D:Arbeit/data/river_network/BasinATLAS_Data_v10/processedData/upstream_area_w_catchment/")
saturated_water_content <- dir_ls("D:Arbeit/data/river_network/soilHydraulicPropertiesEurope/processedData/soiHydraulicaProperties_w_catchments/")
valley_bottom_flatness  <- dir_ls("D:Arbeit/data/river_network/valley_bottom_flattness/processedData/valleyBottomFlattness_w_catchments/")

#dropv <- readRDS("data/dropv_catchments.rds")
all.var <- list(
                slope,
                lakes,
                soil_oc,
                soil_pH,
                floodpl,
                streamp,
                geology,
                glacier,
                bioclim,
                erosion,
                roughness,
                elevation,
                snowdepth,
                discharge,
                groundwater,
                upstream_area,
                saturated_water_content,
                valley_bottom_flatness
                )


# custom function to remove catchments
rmv <- function(x, y){
        if (any(str_detect(x,y))){
                x <- x[!str_detect(x,y)]
        }
        x
}

# remove whole areas not considered in the project
all.var %<>% lapply(rmv, "Hondo")
all.var %<>% lapply(rmv, "Iceland")
all(sapply(all.var, length) == 30)

# #- find missing 
# av13 <- all.var[[13]] %>% str_extract_all("s\\/.*\\.parquet") %>% str_remove("s\\/")
# av14 <- all.var[[14]] %>% str_extract_all("s\\/.*\\.parquet") %>% str_remove("s\\/w_catchments\\/")
# av14[which(!av14 %in% av13)]


clean_name_list <- all.var[[1]] %>% 
        str_extract("(?i)[a-z]*.parquet") %>%
        str_remove("\\.parquet")



# Loop over catchments 
for (i in 1:length(all.var[[1]])) {
        i.catchment_name <-
                all.var[[1]][i] %>%
                str_extract("(?i)[a-z]*.parquet") %>%
                str_remove("\\.parquet")
        print(paste(i, "-", i.catchment_name))
        
        
        #- check that it is the same catchment for everyone
        if (!str_detect(all.var[[1]][i], i.catchment_name) &
            !str_detect(all.var[[2]][i], i.catchment_name) &
            !str_detect(all.var[[3]][i], i.catchment_name) &
            !str_detect(all.var[[4]][i], i.catchment_name) &
            !str_detect(all.var[[5]][i], i.catchment_name) &
            !str_detect(all.var[[6]][i], i.catchment_name) &
            !str_detect(all.var[[7]][i], i.catchment_name) &
            !str_detect(all.var[[8]][i], i.catchment_name) &
            !str_detect(all.var[[9]][i], i.catchment_name) &
            !str_detect(all.var[[10]][i], i.catchment_name) &
            !str_detect(all.var[[11]][i], i.catchment_name) &
            !str_detect(all.var[[12]][i], i.catchment_name) &
            !str_detect(all.var[[13]][i], i.catchment_name) &
            !str_detect(all.var[[14]][i], i.catchment_name) &
            !str_detect(all.var[[15]][i], i.catchment_name) &
            !str_detect(all.var[[16]][i], i.catchment_name) &
            !str_detect(all.var[[17]][i], i.catchment_name) &
            !str_detect(all.var[[18]][i], i.catchment_name)) {
                stop("Differences in name")
        }
        i.f1  <- read_parquet(all.var[[1]][i])
        i.f2  <- read_parquet(all.var[[2]][i])
        i.f3  <- read_parquet(all.var[[3]][i])
        i.f4  <- read_parquet(all.var[[4]][i])
        i.f5  <- read_parquet(all.var[[5]][i])
        i.f6  <- read_parquet(all.var[[6]][i])
        i.f7  <- read_parquet(all.var[[7]][i])
        i.f8  <- read_parquet(all.var[[8]][i])
        i.f9  <- read_parquet(all.var[[9]][i])
        i.f10 <- read_parquet(all.var[[10]][i])
        i.f11 <- read_parquet(all.var[[11]][i])
        i.f12 <- read_parquet(all.var[[12]][i])
        i.f13 <- read_parquet(all.var[[13]][i])
        i.f14 <- read_parquet(all.var[[14]][i])
        i.f15 <- read_parquet(all.var[[15]][i])
        i.f16 <- read_parquet(all.var[[16]][i])
        i.f17 <- read_parquet(all.var[[17]][i])
        i.f18 <- read_parquet(all.var[[18]][i])

        
        i.f.all <- list(
                i.f1,
                i.f2,
                i.f3,
                i.f4,
                i.f5,
                i.f6,
                i.f7,
                i.f8,
                i.f9,
                i.f10,
                i.f11,
                i.f12,
                i.f13,
                i.f14,
                i.f15,
                i.f16,
                i.f17,
                i.f18
        )
        
        for (j in seq_along(i.f.all)) {
                if ("ID2" %in% names(i.f.all[[j]])) {
                        i.f.all[[j]] <- i.f.all[[j]][, ID2 := NULL]
                }
        }
        
        #- thanks to Michael Ohlrogge for this solution
        #- https://stackoverflow.com/questions/13273833/merging-multiple-data-tables
        i.f.all2 = Reduce(function(...)
                merge(..., by = "ID", all = TRUE), i.f.all)
        i.f.all2[, c("VBM_min", "spi_min") := NULL]
        write_parquet(
                i.f.all2,
                paste0(
                        "data/catchments/",
                        i.catchment_name,
                        "_w_variables.parquet"
                )
        )
        
}

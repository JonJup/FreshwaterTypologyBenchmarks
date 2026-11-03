## Prepare data for Fish Example analysis 
## 
## This script extract all necesarry data to run the 'prepare_for_app.R' script
## which can prepare your data for the benchmark generation model shiny app. 


# setup -------------------------------------------------------------------
library(sf)
library(data.table)

# load data ---------------------------------------------------------------
# Available at https://doi.pangaea.de/10.1594/PANGAEA.908578
broadRiverTypes <- st_read("E://Arbeit/data/typology_systems/broad_river_types/rawData/m_river_fec_broad_type.shp")
# Available from MIDFIRE Zenodo 
fish <- readRDS("E://Arbeit/data/biota/MIDFIRE/data/midfire_fish.rds")
fish_w_env <- readRDS("data/000_biota/02_fish_w_environment.rds")

# prepare data ------------------------------------------------------------
fish[, month := lubridate::month(eventDate)]
pressures <- readRDS("E://Arbeit/data/my_stressor_db/processedData/PCA_stressors.rds")
fishSub <- 
        fish[data.set == "fish_france_monitoring_naiades" & 
                     eventYear == 2016 & 
                     month %in% c(7,8,9)
        ]


# community  --------------------------------------------------------------
community <- copy(fishSub)
community[, site := siteID]
community[, taxon := working.taxon]
community[, value := ifelse(organismQuantity >0,1,0)]
community[, date := eventDate]
community <- community[, c("site", "taxon", "value", "date")]



# environment --------------------------------------------------------------
env <- copy(fish_w_env)
env <- env[, slope:VBM_max]
env$site <- fish_w_env$siteID




# sites -------------------------------------------------------------------
sites <- unique(fishSub, by = "siteID")
sites <- data.table(site = sites$siteID, 
                    x = sites$x.coord, 
                    y = sites$y.coord, 
                    crs = sites$EPSG)

# typology as vector 

# save to file  -----------------------------------------------------------
fwrite(community, "data/application_example/brt_fish_community.csv")
fwrite(env, "data/application_example/brt_fish_env.csv")
fwrite(sites, "data/application_example/brt_fish_sites.csv")


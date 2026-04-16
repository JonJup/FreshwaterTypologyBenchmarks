################################################################################
# Script Name:        add_euHydro_to_biota.R
# Description:        Short description of the script
#
# Author:             Jonathan Jupke
# Date Created:       2025-09-18
# Last Modified:      2025-09-18
#
# R Version:          R 4.5.2
# Required Packages:  package1, package2
#
# Notes:              Any notes or assumptions
################################################################################

# setup -------------------------------------------------------------------

setwd(rstudioapi::getActiveProject())

library(sf)
library(sfarrow)
library(dplyr)
library(data.table)

# load data -------------------------------------------------------------------------

# files      <- list.files("data/biota", regexp = "00_")
files <- c(
        "E://data/biota/STREAM/data/stream_diatoms.rds",
        "E://data/biota/STREAM/data/stream_fish.rds",
        "E://data/biota/invertebrate data/00_combine_data/PULSE/pulse_invertebrates.rds",
        "E://data/biota/STREAM/data/stream_macrophytes.rds"
)
# More than one slash after D: causes an error
vector_data <- list.files("E:/data/hydroDEM_parquet/", full.names = T)
bio.names   <- c("diatoms", "fish", "invertebrates", "macrophytes")

for (i in 1:4){
        #if (i %in% c(1,2,4)) next()
        print(paste("STARING", bio.names[i]))
        biota       <- readRDS(files[i])
        
        # prepare sites  ----------------------------------------------------------
        sites <- unique(biota, by = c("siteID"))
        sites <- sites[!is.na(sites$x.coord), ]
        sites <- st_as_sf(sites, coords = c("x.coord", "y.coord"), crs = 3035)
        sites <- st_transform(sites, 4326)
        
        out.ls <- vector(mode = "list", length = length(vector_data))
        
        for (ii in 1:length(vector_data)){
                print(ii)
                i.vec  <- st_read_parquet(vector_data[ii])
                i.vec  <- st_make_valid(i.vec)
                i.join <- st_join(sites, i.vec)
                
                i.join <- i.join[!is.na(i.join$ID), ]
                i.join <- i.join[, c("siteID", "ID")]
                i.join <- st_drop_geometry(i.join)
                i.join <- unique(i.join, by = "siteID")
                
                out.ls[[ii]] <- i.join
                rm(i.vec, i.join, ii)
        }
        
        out.ls <- bind_rows(out.ls)
        sites  <- left_join(sites, out.ls, by = c("siteID"))
        sites  <- select(sites, siteID, ID) %>% st_drop_geometry %>% setDT
        sites <- unique(sites, by = "siteID")
        biota  <- sites[biota, on = "siteID"]
        
        saveRDS(biota, paste0("data/biota/01_", bio.names[i],"_w_catchment_id.rds"))
        rm(out.ls, biota, sites)
        
}; rm(i); gc()
                
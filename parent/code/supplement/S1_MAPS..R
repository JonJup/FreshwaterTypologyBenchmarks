## Data for SM: CHANGES INCURRED THROUGH FILTERING

# setup -------------------------------------------------------------------
library(data.table)
library(tidyverse)
library(magrittr)
library(sf)
library(maptiles)
library(tmap)
library(scales)

# working directory
setwd(rstudioapi::getActiveProject())

# custom function 
source("code/functions/find_max_consequtive_sum.R")

stadia_key <- Sys.getenv("STADIA_KEY")
tile_url <- paste0("https://tiles.stadiamaps.com/tiles/stamen_terrain_background/",
                   "{z}/{x}/{y}.png?api_key=", stadia_key)


# MIDFIRE data ------------------------------------------------------------
path_midfire <- "E://Arbeit/data/biota/MIDFIRE/data/"
pathToProject <- c("~/projects/pulse/01_wp1/MEE_REVISION_1/data/000_biota/")
path_turso <- "//turso-fs.cs.helsinki.fi/turso/home/jonjupke"
# load data ---------------------------------------------------------------
d <- readRDS(file.path(path_midfire, "midfire_diatoms.rds"))
f <- readRDS(file.path(path_midfire, "midfire_fish.rds"))
i <- readRDS("E://Arbeit/data/biota/invertebrate data/00_combine_data/midfire_invertebrates.rds")
m <- readRDS(file.path(path_midfire, "midfire_macrophytes.rds"))

# schemes <- list.files(pathToProject, pattern = "^03", full.names = TRUE)
# schemes %<>% lapply(readRDS)
# schemes2 <- rbindlist(schemes)

# Combine midfire data
d[, taxon := "diatoms"]
f[, taxon := "fish"]
i[, taxon := "invertebrates"]
m[, taxon := "macrophytes"]

# combine all four taxa into one data set 
data <- rbindlist(list(d,f,i,m), use.names=T)
# one row per sampling site  
data2 <- unique(data, by = c("taxon", "eventID"))

# prepare spatial data for all MIDFIRE sites ------------------------------
full <- st_as_sf(data2, coords = c("x.coord", "y.coord"), crs = 3035)
full$taxon <- factor(full$taxon)

pdf("figures/supplement/S1/maps.pdf", width = 8, height = 6, onefile = TRUE)

for (lev in levels(full$taxon)) {
        
        p <- tm_basemap(server = tile_url) +
                tm_shape(full[full$taxon == lev, ], is.main = TRUE) +
                tm_dots(fill_alpha = 0.6) +
                tm_layout(panel.label.bg.color = "white", legend.outside = TRUE) + 
                tm_title(lev)
        print(p)   # print() is needed inside loops
}

dev.off()


data[, max(eventYear,na.rm = T), by = "taxon"]
data[, min(siteID), by = c("taxon", "data.set")][, sum(V1), by = "taxon"]

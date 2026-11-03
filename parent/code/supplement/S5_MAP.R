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
path_turso <- ""
# load data ---------------------------------------------------------------
d <- readRDS(file.path(path_midfire, "midfire_diatoms.rds"))
f <- readRDS(file.path(path_midfire, "midfire_fish.rds"))
i <- readRDS("E://Arbeit/data/biota/invertebrate data/00_combine_data/midfire_invertebrates.rds")
m <- readRDS(file.path(path_midfire, "midfire_macrophytes.rds"))

schemes <- list.files(pathToProject, pattern = "^03", full.names = TRUE)
schemes %<>% lapply(readRDS)
schemes2 <- rbindlist(schemes)

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


# coordinates from fitted models  -----------------------------------------
schemes2 <- schemes2[sample_type == "full"]
taxa <- sort(unique(full$taxon))
site.list <- vector(mode = "list", length = length(taxa))
names(site.list) <- taxa
for (tax in taxa){
        
        tax_shot <- switch(tax, 
                           diatoms = "D", 
                           fish = "F", 
                           invertebrates = "I", 
                           macrophytes = "M"
                           )
        
        full_path = paste0(path_turso,"/pulse",tax_shot,"/data/misc/site_coordinates/")
        full_crd_path <- list.files(full_path, full.names = TRUE)
        scheme_check <- gsub(pattern = full_path, replacement = "", x = full_crd_path)
        scheme_check <- gsub(pattern = "\\.rds" , replacement = "", x = scheme_check)
        full_crd_path <- full_crd_path[which(scheme_check %in% schemes2$scheme_id)]
        site.list[[tax]] <- rbindlist(lapply(full_crd_path, readRDS), use.names = TRUE)
        site.list[[tax]]$taxon <- tax
        
}
sub <- rbindlist(site.list)
sub <- st_as_sf(sub, coords = c("x.coord", "y.coord"), crs = 3035)
taxa <- sort(unique(full$taxon))
hex_long <- do.call(rbind, lapply(taxa, function(tx) {
        h <- hex
        h$taxon  <- tx
        h$n_full <- lengths(st_intersects(h, full[full$taxon == tx, ]))
        h$n_sub  <- lengths(st_intersects(h, sub[sub$taxon == tx, ]))
        h
}))
hex_long <- hex_long[hex_long$n_full > 0, ]
hex_long$prop <- hex_long$n_sub / hex_long$n_full

p <- tm_basemap(server = tile_url) +
        tm_shape(hex_long, is.main = TRUE) +
        tm_polygons(fill = "prop", col = NULL,
                    fill.scale  = tm_scale_continuous(values = "viridis", limits = c(0, 1), outliers.trunc = c(TRUE, TRUE)),
                    fill.legend = tm_legend(title = "Proportion retained",
                                            orientation = "landscape",
                                            position = tm_pos_out("center", "bottom"))
                    ) +
        tm_facets_wrap(by = "taxon", ncol  = 2) +
        tm_layout(panel.label.bg.color = "white", legend.outside = TRUE)

p
tmap_save(tm = p,
          filename = "figures/supplement/S5/S5_MAP.png",
          width = 18, height = 22, units = "cm", dpi = 300)



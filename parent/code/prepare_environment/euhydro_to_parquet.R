################################################################################
# Script name:  euhydro_to_parquet.R
#
# Purpose:      Convert the downloaded EU-Hydro DEM river-network files
#               (GeoPackage, one per River Basin District) into per-basin
#               GeoParquet files, and assemble a single combined layer of all
#               basins.
#
# Inputs:       - EU-Hydro DEM GeoPackages (one per River Basin District)
#                 https://land.copernicus.eu/en/products/eu-hydro/eu-hydro-river-network-database#download
#
# Output:       - parquet/<basin>.parquet             (one GeoParquet per basin)
################################################################################


# setup ------------------------------------------------------------------------
library(sf)
library(sfarrow)    
library(fs)         
library(dplyr)
library(magrittr)   

# All relative paths below resolve against this directory.
base_dir <- "E:/Arbeit/Data/river_network/eu_hydro_dem/"
setwd(base_dir)

# Make sure the parquet output folder exists before the loop writes into it.
dir.create("processedData/", showWarnings = FALSE)
dir.create("processedData/parquet_catchments/", showWarnings = FALSE)

# Each EU-Hydro River Basin District download is a folder whose name contains
# "GPKG" (case-sensitive); list those folders.
# This assumes you stored the downloaded data in a folder called "rawData".
# Adjust as needed. 
files <- dir_ls("rawData/", regexp = "GPKG")

# Pre-allocate a list to collect each basin's geometry for the combined layer.
store_list <- vector("list", length(files))


# convert each basin to GeoParquet ---------------------------------------------
for (i in seq_along(files)) {
        message("basin ", i, "/", length(files))
        
        # Inside the basin folder, pick the GeoPackage holding the drainage network.
        i.file <- dir_info(files[i])$path
        i.file <- i.file[grep("drainage", i.file)]
        
        # A drainage GeoPackage can hold several layers; take the one with the most
        # attribute fields (the full river-network layer rather than helper layers).
        i.layers <- st_layers(i.file)
        i.layers <- i.layers$name[which.max(i.layers$fields)]
        
        # Read that layer and reproject to WGS84 (EPSG:4326).
        i.sf <- st_read(i.file, layer = i.layers) %>%
                st_transform(4326)
        
        # Basin name, used for the output file name. Assumes one basin per drainage
        # layer, so unique() returns a single value.
        i.write.name <- unique(i.sf$NAMEBASIN)
        
        # Keep only a unique feature ID (basin name + running number) and geometry.
        i.sf <- i.sf %>%
                mutate(ID = paste0(NAMEBASIN, row_number())) %>%
                select(ID)
        
        st_write_parquet(i.sf, paste0("processedData/parquet_catchments/", i.write.name, ".parquet"))
        
        rm(list = ls()[grepl("^i\\.", ls())])
}

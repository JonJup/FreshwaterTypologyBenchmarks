################################################################################
# Script name:  slope.R
#
# Purpose:      This script prepares the slope data from hydrograph90m. The 
#               tiles are first merged and then values are extracted and 
#               averaged for each catchment.
#
# Inputs:       - EU-Hydro DEM Parquet (one per River Basin District); 
#                       Created with: euhydro_to_parquet.R
#               - Elevation difference between focal cell and downstream cell from Hydrography90m
#               https://public.igb-berlin.de/index.php/s/agciopgzXjWswF4?path=%2Fr.stream.slope%2Fslope_elv_dw_cel_tiles20d
#               We use the tiles:   
#               "slope_elv_dw_cel_h16v02.tif" 
#               "slope_elv_dw_cel_h16v04.tif" 
#               "slope_elv_dw_cel_h18v00.tif" 
#               "slope_elv_dw_cel_h18v02.tif"
#               "slope_elv_dw_cel_h18v04.tif" 
#               "slope_elv_dw_cel_h20v00.tif" 
#               "slope_elv_dw_cel_h20v02.tif" 
#               "slope_elv_dw_cel_h20v04.tif"
#
# Output:       - 
################################################################################

library(terra)
library(sf)
library(sfarrow)
library(data.table)
library(dplyr)
library(exactextractr)

path_to_slope <- "E:/Arbeit/Data/river_network/hydrography90m/rawData/Elevation difference between focal cell and downstream cell/"
out_file      <- "E:/Arbeit/Data/river_network/hydrography90m/processedData/slope_merged_europe.tif"

# List the input tiles. Filter to .tif and drop the output file, so re-running
# the script does not fold the previous mosaic back into the new one.
raster.files <- fs::dir_ls(path = path_to_slope, regexp = "\\.tif$")
raster.files <- raster.files[basename(raster.files) != basename(out_file)]
if (length(raster.files) == 0) stop("No input tiles found in: ", in_dir)

# merge tiles ------------------------------------------------------------------
# Merge one tile at a time, calling gc() between steps to keep memory in check
#  The tiles are adjacent and do not overlap.
r.all <- rast(raster.files[1])
for (i in seq_along(raster.files)[-1]) {
        message("merging tile ", i, "/", length(raster.files))
        r.all <- merge(r.all, rast(raster.files[i]))
        gc()
}
# write output 1-----------------------------------------------------------------
writeRaster(r.all, filename = out_file, overwrite = TRUE)

# adjust paths to your EU Hydro DEM PARQUET and hydrography90m folders
path_to_parquet <- "E://Arbeit/data/river_network/eu_hydro_dem/processedData/parquet_catchments/"
vector.files <- fs::dir_ls(path_to_parquet)

for (i in 1:length(vector.files)){
        
        i.cat <- st_read_parquet(dsn = vector.files[i])
        i.name <- stringr::str_remove(i.cat$ID[1], "1")
        i.cat$ID2 <- 1:nrow(i.cat)
        i.ext <- i.cat %>%
                st_bbox() %>% 
                st_as_sfc() %>% 
                st_as_sf(crs = st_crs(i.cat))
        
        i.raster <- crop(r.all, i.ext)
        i.stats <- exact_extract(i.raster, i.cat, c("mean"))
        i.cat$slope <- round(i.stats,2)
        i.cat <- st_drop_geometry(i.cat) %>% setDT
        i.uri <- paste0("E:Arbeit/Data/river_network/hydrography90m/processedData/slope_w_catchments/",i.name,".parquet")
        arrow::write_parquet(
                x = i.cat, 
                sink = i.uri
        )
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
}


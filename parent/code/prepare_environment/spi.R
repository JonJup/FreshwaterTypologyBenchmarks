# ################################################################################
# # Script name:  spi.R
# #
# # Purpose:      This script prepares the stream power index for Eu Hydro DEM
#                 catchments. 
#                 
# # Inputs:       - EU-Hydro DEM Parquet (one per River Basin District); 
# #                       Created with: euhydro_to_parquet.R
#                 - Stream Power Index from Hydrograph90m
#                 All files available from: https://public.igb-berlin.de/index.php/s/agciopgzXjWswF4?path=%2Fflow.index%2Fspi_tiles20d 
#                 We used the following files: 
#                 spi_h16v02.tif
#                 spi_h16v04.tif
#                 spi_h18v00.tif
#                 spi_h18v02.tif
#                 spi_h18v04.tif
#                 spi_h20v00.tif
#                 spi_h20v02.tif
#                 spi_h20v04.tif
# #
# # Output:       - oc_w_catchment/*
#                       parquet files with mean topsoil organic carbon fraction 
#                       in catchment. 
# ################################################################################

# setup -------------------------------------------------------------------
library(terra)
library(sf)
library(sfarrow)
library(data.table)
library(dplyr)
library(exactextractr)


# load data ---------------------------------------------------------------
raster.files <- fs::dir_ls(path = "E://Arbeit/Data/river_network/hydrography90m/rawData/spi/", regexp = "\\.tif$")

r.all <- rast(raster.files[1])
for (i in seq_along(raster.files)[-1]) {
        message("merging tile ", i, "/", length(raster.files))
        r.all <- merge(r.all, rast(raster.files[i]))
        gc()
}
# save to file 
terra::writeRaster(r.all, filename = "E://Arbeit/Data/river_network/hydrography90m/processedData/spi_merged_europe.tif")

# spi <- r.all
# rm(r.all)

spi <- rast("E://Arbeit/Data/river_network/hydrography90m/processedData/spi_merged_europe.tif")
vector.files <- fs::dir_ls("E://Arbeit/data/river_network/eu_hydro_dem/processedData/parquet_catchments/")
out.path = "E:Arbeit/data/river_network/hydrography90m/processedData/spi_w_catchment/"
dir.create(out.path, showWarnings = FALSE, recursive = TRUE)

for (i in 1:length(vector.files)){
        print(i)
        i.cat  <- st_read_parquet(dsn = vector.files[i])
        i.name <- stringr::str_remove(i.cat$ID[1], "1")
        i.cat2 <- 
                i.cat %>%
                st_drop_geometry %>%
                setDT
        i.ext <- i.cat %>%
                st_bbox() %>% 
                st_as_sfc() %>% 
                st_as_sf(crs = st_crs(i.cat))
        
        i.raster <- crop(spi, i.ext)
        stats <- exact_extract(i.raster, i.cat, c("mean", "min", "max"))
        i.cat2$spi_mean <- stats$mean   # coverage-weighted
        i.cat2$spi_min  <- stats$min    # unweighted over covered cells
        i.cat2$spi_max  <- stats$max
        i.uri <- paste0(out.path,i.name,".parquet")
        arrow::write_parquet(
                x = i.cat2, 
                sink = i.uri
        )
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
}

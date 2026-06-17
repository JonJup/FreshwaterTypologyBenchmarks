
# setup -------------------------------------------------------------------
library(terra)
library(sf)
library(sfarrow)
library(data.table)
library(dplyr)
library(ncdf4)


# load original data ------------------------------------------------------
nc1 <- nc_open("E:/Arbeit/Data/river_network/discharge copernicus/rdis_tmean_abs_E-HYPEgrid-EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CCLM4-8-17-v1_na_1971-2000_grid5km_v1.nc")
nc2 <- nc_open("E:/Arbeit/Data/river_network/discharge copernicus/rdisyearmax_tmean_abs_E-HYPEgrid-EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CCLM4-8-17-v1_na_1971-2000_grid5km_v1.nc")
nc3 <- nc_open("E:/Arbeit/Data/river_network/discharge copernicus/rdisyearmin_tmean_abs_E-HYPEgrid-EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CCLM4-8-17-v1_na_1971-2000_grid5km_v1.nc")


# transform to sf object class --------------------------------------------
discharge_mean <- st_as_sf(data.frame(lon = c(ncvar_get(nc1, "lon")), lat = c(ncvar_get(nc1, "lat")), mean_discharge = c(ncvar_get(nc1, "rdis_tmean"))), coords = c("lon", "lat"), crs = 4326)
discharge_max  <- st_as_sf(data.frame(lon = c(ncvar_get(nc2, "lon")), lat = c(ncvar_get(nc2, "lat")), max_discharge  = c(ncvar_get(nc2, "rdisyearmax_tmean"))), coords = c("lon", "lat"), crs = 4326)
discharge_min  <- st_as_sf(data.frame(lon = c(ncvar_get(nc3, "lon")), lat = c(ncvar_get(nc3, "lat")), min_discharge  = c(ncvar_get(nc3, "rdisyearmin_tmean"))), coords = c("lon", "lat"), crs = 4326)


# remove missing values  --------------------------------------------------
discharge_mean <- dplyr::filter(discharge_mean, !is.na(mean_discharge))
discharge_max  <- dplyr::filter(discharge_max,  !is.na(max_discharge))
discharge_min  <- dplyr::filter(discharge_min,  !is.na(min_discharge))


# save as parquet ---------------------------------------------------------
st_write_parquet(discharge_mean, "E:Arbeit/Data/river_network/discharge_copernicus/processedData/discharge_mean.parquet")
st_write_parquet(discharge_max,  "E:Arbeit/Data/river_network/discharge_copernicus/processedData/discharge_max.parquet")
st_write_parquet(discharge_min,  "E:Arbeit/Data/river_network/discharge_copernicus/processedData/discharge_min.parquet")

# load data parquet ---------------------------------------------------------------
d1 <- st_read_parquet("E:Arbeit/Data/river_network/discharge copernicus/discharge_mean.parquet")
d2 <- st_read_parquet("E:Arbeit/Data/river_network/discharge copernicus/discharge_max.parquet")
d3 <- st_read_parquet("E:Arbeit/Data/river_network/discharge copernicus/discharge_min.parquet")

# EU Hydro DEM catchment polygons
vector.files <- fs::dir_ls("E://Arbeit/Data/river_network/eu_hydro_dem/processedData/parquet_catchments/")

for (i in 1:length(vector.files)){
        
        i.cat  <- st_read_parquet(dsn = vector.files[i])
        i.name <- stringr::str_remove(i.cat$ID[1], "1")
        i.cat$ID2 <- 1:nrow(i.cat)

        i.ext <- i.cat %>%
                st_bbox() %>% 
                st_as_sfc() %>% 
                st_as_sf(crs = st_crs(i.cat))%>%
                st_buffer(50000)
        
        i.cat2 <- 
                i.cat %>%
                st_drop_geometry %>%
                setDT
        i.d1 <- st_filter(d1, i.ext)
        i.d2 <- st_filter(d2, i.ext)
        i.d3 <- st_filter(d3, i.ext)
        
        i.d1 <- vect(i.d1)
        i.d2 <- vect(i.d2)
        i.d3 <- vect(i.d3)
        
        i.ext1 <- ext(i.d1)
        i.ext2 <- ext(i.d2)
        i.ext3 <- ext(i.d3)
        
        r1 <- rast(i.ext1 , resolution = 0.099)
        r2 <- rast(i.ext2 , resolution = 0.099)
        r3 <- rast(i.ext3 , resolution = 0.099)
        
        r1 <- rasterize(i.d1, r1, field = "mean_discharge")
        r2 <- rasterize(i.d2, r2, field = "max_discharge")
        r3 <- rasterize(i.d3, r3, field = "min_discharge")
       
       i.ext1 <- extract(r1, i.cat, touches = T)
       i.ext2 <- extract(r2, i.cat, touches = T)
       i.ext3 <- extract(r3, i.cat, touches = T)

        setDT(i.ext1)
        setDT(i.ext2)
        setDT(i.ext3)
        
        names(i.ext1)[2] <- "raster"
        names(i.ext2)[2] <- "raster"
        names(i.ext3)[2] <- "raster"
        
        i.ext1[, meanD := mean(raster, na.rm = T), by = "ID"]
        i.ext2[, maxD  := max(raster, na.rm = T), by = "ID"]
        i.ext3[, minD  := min(raster, na.rm = T), by = "ID"]
        
        i.ext1 <- unique(i.ext1, by = "ID")
        i.ext2 <- unique(i.ext2, by = "ID")
        i.ext3 <- unique(i.ext3, by = "ID")
        
        i.cat2$mean_discharge <- i.ext1$meanD
        i.cat2$max_discharge <- i.ext2$maxD
        i.cat2$min_discharge <- i.ext3$minD
        i.uri <- paste0("E:Arbeit/Data/river_network/discharge_copernicus/processedData/w_catchments/",i.name,".parquet")
        arrow::write_parquet(
                x = i.cat2, 
                sink = i.uri
        )
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
}
rm(i)
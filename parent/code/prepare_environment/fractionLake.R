# ################################################################################
# # Script name:  fractionLake.R
# #
# # Purpose:      This script prepares the slope data from hydrograph90m. The 
# #               tiles are first merged and then values are extracted and 
# #               averaged for each catchment.
# #
# # Inputs:       - EU-Hydro DEM Parquet (one per River Basin District); 
# #                       Created with: euhydro_to_parquet.R
# #               - HydroLAKES v1.0 data base. Available at
# #               https://www.hydrosheds.org/products/hydrolakes
# #
# # Output:       - fractionLake_w_catchment/*
#                    - parquet files that have one value per EU-DEM Catchment 
#                    showing the relative area of the catchment covered by lakes. 
# ################################################################################
# 


# setup -------------------------------------------------------------------
library(terra)
library(sf)
library(sfarrow)
library(data.table)
library(dplyr)

path_to_parquet <- "E://Arbeit/data/river_network/eu_hydro_dem/processedData/parquet_catchments/"
vector.files <- fs::dir_ls(path_to_parquet)
out.path = "E://Arbeit/data/river_network/hydroLAKES/processedData/fractionLake_w_catchment/"
dir.create(out.path, showWarnings = FALSE, recursive = T)

# read in data ------------------------------------------------------------
lakes_europe  <- st_read_parquet("E:Arbeit/Data/river_network/hydroLAKES/HydroLAKES_europe.parquet")
invalid_geoms <- st_is_valid(lakes_europe, reason = TRUE)
print(invalid_geoms[invalid_geoms != "Valid Geometry"])

mesima.id <- which(stringr::str_detect(string = vector.files, pattern = "Mesima"))

mkdir("E://Arbeit/data/river_network/hydroLAKES/processedData/fractionLake_w_catchment/")
for (i in seq.along(vector.files)){
        print(i)
        i.cat  <- st_read_parquet(dsn = vector.files[i])
        
        i.name <- stringr::str_remove(i.cat$ID[1], "1")
        i.cat$ID2 <- 1:nrow(i.cat)
        i.cat <- sf::st_make_valid(i.cat)
        i.cat2 <- 
                i.cat %>%
                st_drop_geometry %>%
                setDT
        i.ext <- i.cat %>%
                st_bbox() %>% 
                st_as_sfc() %>% 
                st_as_sf(crs = st_crs(i.cat))
        
        
        #- this created an error because of self-intersection for the Mesima catchment
        if (i == mesima.id){
                i.lakes <- lakes_europe[-which(invalid_geoms != "Valid Geometry"),]
                i.lakes <- st_filter(i.lakes, i.ext)
        } else {
                i.lakes <- st_filter(lakes_europe, i.ext)
        }
        
        # Create spatial index for better performance
        i.polygons2_index <- st_intersects(i.cat, i.lakes, sparse = TRUE)
        i.overlap_areas <- lapply(
                seq_len(nrow(i.cat)), 
                function(i) {
                        intersecting_polys <- i.lakes[i.polygons2_index[[i]], ]
                        if (length(i.polygons2_index[[i]]) > 0) {
                                intersections <- st_intersection(i.cat[i,], intersecting_polys)
                                return(data.frame(
                                        poly2_id = i,
                                        poly1_id = i.polygons2_index[[i]],
                                        area = st_area(intersections)
                                ))
                        }
                        return(NULL)
                }
        )
        #- bind list from lapply into a data.frame
        i.overlap_areas <- do.call(rbind, i.overlap_areas)
        #- turn data.frame into data.table
        setDT(i.overlap_areas)
        #- sum overlapps for each catchment. Only relevant if multiple lakes occur within one catchment.
        i.ola2 <- i.overlap_areas[, sum(area), by = "poly2_id"]
        #- Compute catchment areas for all catchments with lakes to derive relative area of lakes in the next step.
        i.cat.area <- lapply(
                i.ola2$poly2_id, 
                function(i) { 
                        cat_area <- st_area(i.cat[i,])
                        return(cat_area)
                }
        )
        i.ica2 <- unlist(i.cat.area)        
        i.ola2$catchment_area <- i.ica2
        i.ola2[, V1 := units::drop_units(V1)]
        i.ola2[, relative_area := V1/catchment_area]
        i.bind <- data.table(poly2_id = 1:nrow(i.cat2))
        i.ola3 <- merge(i.bind, i.ola2, by = "poly2_id", all.x = T)
        i.ola3[is.na(relative_area), relative_area := 0]
        i.cat2$lake_area <- i.ola3$relative_area
        i.uri <- paste0(out.path,i.name,".parquet")
        arrow::write_parquet(
                x = i.cat2, 
                sink = i.uri
        )
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
}

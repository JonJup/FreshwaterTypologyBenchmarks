library(sfarrow)
library(sf)
library(mapview)

files <- list.files("E:Arbeit/data/river_network/eu_hydro_dem/processedData/parquet_catchments/", full.names = TRUE)

geom <- do.call(rbind, lapply(files, function(f) {
        y <- st_read_parquet(f)
        x <- st_make_valid(y)  
        z <- st_transform(x, 3035)
        aa <- st_simplify(z, preserveTopology = TRUE, dTolerance = 250)
        aa
}))

message(sum(st_is_empty(geom)), " catchments vanished; lower simplify_tol if many")
geom <- geom[!st_is_empty(geom), ]

saveRDS(geom, "E://Arbeit/data/river_network/eu_hydro_dem/processedData/all_eu_hydro_dem_catchments_simplified_250m.rds")
determine_spatial_scale <- function(x){
        
        sites <- unique(x, by = "siteID")
        
        o.sf  <- st_as_sf(sites,
                         coords = c("x.coord", "y.coord"),
                         crs = 3035)
        
        o.sf2 <- o.sf
        o.sf2 <- st_transform(o.sf, crs = 4326)
        o.sf2 <- st_coordinates(o.sf2)
        o.sf <- so.sf2o.sf <- st_distance(o.sf)
        o.sf <- o.sf[lower.tri(o.sf, diag = F)]
        o.sf <- as.numeric(o.sf)
        
        o.sf <- data.table(
                scheme_id =  paste0(bio.names[b], "_", o.scheme.number),
                taxon = o.scheme$taxon,
                data.set  = o.scheme$data.set,
                year = o.scheme$year,
                samples = o.scheme$samples,
                min_distance = min(o.sf),
                mean_distance = mean(o.sf),
                median_distance = median(o.sf),
                max_distance = max(o.sf),
                median_latitude = median(o.sf2[,2]),
                median_longitude = median(o.sf2[,1]),
                max_latitude  = max(o.sf2[,2]),
                max_longitude = max(o.sf2[,1]),
                min_latitude  = min(o.sf2[,2]),
                min_longitude = min(o.sf2[,1]),
                organismQuantityType = o.scheme$organismQuantityType,
                sample_type = o.scheme$sample_type
        )
        
        # o.sf <- list(
        #         min = min(o.sf),
        #         max = max(o.sf),
        #         mean = mean(o.sf),
        #         median = median(o.sf)
        # )
        return(o.sf)
}
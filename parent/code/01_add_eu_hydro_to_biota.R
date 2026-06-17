################################################################################*
# Description:        Spatially joins EU-Hydro catchment IDs to biological
#                     monitoring sites for four organism groups (diatoms, fish,
#                     invertebrates, macrophytes). Each site is matched to its
#                     enclosing catchment polygon via a point-in-polygon join
#                     across all parquet catchment tiles.
# Notes:              - Biota RDS files must follow the same order as bio.names
#                     - Input coordinates are assumed to be in EPSG:3035
#                     - Catchment parquet tiles must contain an "ID" column
#                     - Output is written to data/biota/
################################################################################*

# 1.0 setup -------------------------------------------------------------------
library(sf)
library(sfarrow)   
library(dplyr)
library(data.table)

# 2.0 load data ---------------------------------------------------------------

# Discover all biota RDS files in the MIDFIRE data directory.
# list.files() returns files in alphabetical order, so bio.names below
# must match that order. Verify with: list.files(..., pattern="\\.rds$")

files <- list.files(
                "parent/data/biota/",   
        pattern    = "\\.csv$",                  
        full.names = TRUE
)

# All GeoParquet tiles covering the EU-Hydro catchment polygons
vector_data <- list.files(
        "parent/data/catchments/",
        full.names = TRUE
)

# Organism group labels — must match the alphabetical order of `files`
bio.names <- c("diatoms", "fish", "invertebrates", "macrophytes")

# Sanity check: abort early rather than silently process the wrong files
stopifnot(length(files) == length(bio.names))

# 3.0 main loop ---------------------------------------------------------------
for (i in seq_along(bio.names)) {
        
        message("STARTING ", bio.names[i])
        
        biota <- fread(files[i])
        
        # ========================* 
        ## 3.1 prepare sites ----
        # ========================*
        
        # Reduce to one row per site before the expensive spatial join so we are
        # not doing redundant point-in-polygon tests for repeated visits.
        
        # Coerce to data.table for unique(..., by = "") to work
        if (!is.data.table(biota)) setDT(biota)
        sites <- unique(biota, by = "siteID")          # one row per site
        
        # Drop sites with missing coordinates 
        sites <- sites[!is.na(sites$x.coord) & !is.na(sites$y.coord), ]
        
        # Build an sf object
        sites <- st_as_sf(sites, coords = c("x.coord", "y.coord"), crs = 3035)
        
        # Re-project to WGS84 to match the catchment parquet tiles
        sites <- st_transform(sites, 4326)
        
        # ================================================*
        ## 3.2 spatial join across all catchment tiles ----
        # ================================================*
        
        # Each tile covers a geographic subset.
        # We collect matches from all tiles and then deduplicate.
        out.ls <- vector(mode = "list", length = length(vector_data))
        
        for (ii in seq_along(vector_data)) {
                
                message("  tile ", ii, " / ", length(vector_data))
                
                i.vec <- st_read_parquet(vector_data[ii])
                # repair any invalid geometries in the tile
                i.vec <- st_make_valid(i.vec)   
                
                # Point-in-polygon join: each site row gains the attributes of the
                # catchment polygon it falls inside (left = sites, right = catchments)
                i.join <- st_join(sites, i.vec)
                
                # Keep only sites that actually fell inside a catchment polygon
                i.join <- i.join[!is.na(i.join$ID), ]
                
                # Retain only the two columns needed for the downstream merge
                i.join <- i.join[, c("siteID", "ID")]
                i.join <- st_drop_geometry(i.join)
                
                setDT(i.join)
                i.join <- unique(i.join, by = "siteID")
                
                out.ls[[ii]] <- i.join

                rm(i.vec, i.join)
        }
        
        # ================================================*
        ## 3.3 Combine tile results  -----
        # ================================================*
        out.dt <- rbindlist(out.ls)
        
        # If a site was matched in more than one tile, keep only the first match
        out.dt <- unique(out.dt, by = "siteID")
        
        # Attach catchment ID back to the full biota table 
        # data.table right-join: all biota rows kept, ID added where available
        biota <- out.dt[biota, on = "siteID"]   
        
        # ================================================*
        ## 3.4 Save -----
        # ================================================*
        
        out_dir <- paste0(bio.names[i],"_folder/data/biota/")
        if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
        
        saveRDS(biota, file.path(out_dir, paste0("01_", bio.names[i], "_w_catchment_id.rds")))
        
        rm(out.ls, out.dt, biota, sites)
        message("DONE ", bio.names[i])
}

rm(i)
gc()
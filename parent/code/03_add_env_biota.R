################################################################################*
# Description:        Adds environmental predictors (from per-catchment parquet
#                     files) to each of the four biological datasets. Outputs 
#                     one RDS per taxon group with all predictors joined.
#                     
# Notes:              - Biota files must have been produced by the script
#                       add_eu_Hydro_to_biota.R and contain an "ID" column that 
#                       indicates the EU HYDRO catchment the sampling location 
#                       is in.
#                     - organismQuantity == 0 records are treated as absences
#                       and removed; all retained records receive PA = 1.
################################################################################*

# 1.0 setup -------------------------------------------------------------------
library(arrow)
library(data.table)

# 2.0 load and prepare data ---------------------------------------------------

bio.names <- c("diatoms", "fish", "invertebrates", "macrophytes")

files <- c(
        list.files(paste0(bio.names, "_folder/data/biota/"), full.names = TRUE)
)


# Per-catchment environmental parquet files (one per tile)
catchments <- list.files("parent/data/catchments_w_environment", full.names = TRUE)

# Load all biological datasets into memory
bio.list       <- lapply(files, readRDS)
# ensure data.table throughout
bio.list       <- lapply(bio.list, setDT)   
n_bio_datasets <- length(bio.list)

# Pre-allocate the output accumulator: a list of lists, one per organism group.
out_list <- vector("list", n_bio_datasets)
names(out_list) <- bio.names
for (j in seq_len(n_bio_datasets)) {
        out_list[[j]] <- vector("list", length(catchments))
}



# 3.0 loop over catchment tiles -----------------------------------------------
for (i in seq_along(catchments)) {
        
        message("Catchment tile ", i, " / ", length(catchments))
        
        i.vec <- read_parquet(catchments[i])
        setDT(i.vec)
        
        # Subset each biological dataset to segments present in this tile
        i.data <- lapply(bio.list, function(x) x[ID %in% i.vec$ID])
        
        # Join environmental predictors (all.x = TRUE → left join, keep all bio rows)
        out <- lapply(i.data, function(x) merge(x, i.vec, by = "ID", all.x = TRUE))
        
        # Store this tile's results for each organism group
        for (j in seq_len(n_bio_datasets)) {
                out_list[[j]][[i]] <- out[[j]]
        }
        
        rm(i.vec, i.data, out)
}

# 4.0 Combine per-tile results ------------------------------------------------
# rbindlist(fill = TRUE) handles tiles where some columns may be absent
final_results <- lapply(out_list, function(x) rbindlist(x, fill = TRUE))

# Retain only presence records (organismQuantity != 0) …
final_results <- lapply(final_results, function(x) x[organismQuantity != 0])

# … and label them as presences for downstream modelling
final_results <- lapply(final_results, function(x) x[, PA := 1L])

# 5.0 save to file ------------------------------------------------------------
out_dir <- paste0(bio.names, "_folder/data/biota/")


invisible(lapply(seq_len(n_bio_datasets), function(j) {
        if (!dir.exists(out_dir[j])) dir.create(out_dir[j], recursive = TRUE)
        saveRDS(
                final_results[[j]],
                file.path(out_dir[j], paste0("02_", bio.names[j], "_w_environment.rds"))
        )
        message("Saved: 02_", bio.names[j], "_w_environment.rds")
}))

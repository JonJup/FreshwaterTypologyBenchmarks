################################################################################
# Script Name:        add_env_biota.R
# Description:        Adds environmental predictors (from per-catchment parquet
#                     files) and seven typology classifications to each of the
#                     four biological datasets. Outputs one RDS per organism
#                     group with all predictors and typologies joined.
# Notes:              - Biota files must have been produced by add_euHydro_to_biota.R
#                       (prefix "01_") and contain an "ID" catchment column.
#                     - Typology RDS files must be data.tables with columns
#                       c("ID", <typology_name>) and are read in the fixed order:
#                       FEOW, GLORiC, BGR, BRT, EnZ, HER, IFE.
#                     - organismQuantity == 0 records are treated as absences
#                       and removed; all retained records receive PA = 1.
################################################################################

# setup -------------------------------------------------------------------
library(arrow)
library(data.table)

# load and prepare data ---------------------------------------------------

# Biota files produced by the previous pipeline step (prefix "01_")
files <- list.files("data/biota", pattern = "01_", full.names = TRUE)

# Extract organism-group labels from filenames for use in output naming.
bio.names <- gsub("data/biota/01_|_w_catchment_id\\.rds", "", files)

# Per-catchment environmental parquet files (one per tile)
catchments <- list.files("data/catchments", full.names = TRUE)

# Typology RDS files — order matters; must match the name assignments below
typologies <- list.files("data/typologies", full.names = TRUE)
typologies <- lapply(typologies, readRDS)
# Ensure every element is a data.table (readRDS may return a data.frame)
typologies <- lapply(typologies, setDT)

# Assign a consistent, meaningful second-column name to each typology.
# Each typology file is assumed to have exactly two columns: "ID" and the
# raw typology code column (whatever it was named in the source data).
typo_labels <- c("FEOW", "GLORiC", "BGR", "BRT", "EnZ", "HER", "IFE")
stopifnot(length(typologies) == length(typo_labels))   # fail fast if files added/removed

for (k in seq_along(typologies)) {
        # Guard: confirm we are renaming exactly the non-ID column
        stopifnot(ncol(typologies[[k]]) == 2, "ID" %in% names(typologies[[k]]))
        # Rename the second column (the typology code) to its canonical label
        old_name <- setdiff(names(typologies[[k]]), "ID")
        setnames(typologies[[k]], old_name, typo_labels[k])
}

# Modal typology aggregation ------------------------------------------
# Each catchment segment (ID) can overlap multiple typology polygons, so
# multiple rows per ID may exist. Reduce to one row per ID by taking the
# modal (most frequent) typology class, replacing "no data" sentinel 0 → NA.

mode_func <- function(x) {
        # Returns the modal value of x, ignoring NAs.
        # Returns 0 (later coerced to NA) when x is all-NA.
        ux <- unique(x[!is.na(x)])
        if (length(ux) == 0L) return(0)
        ux[which.max(tabulate(match(x, ux)))]
}

typologies <- lapply(seq_along(typologies), function(k) {
        dt  <- typologies[[k]]
        col <- typo_labels[k]
        # Compute modal class per segment in-place
        dt[, (col) := mode_func(get(col)), by = "ID"]
        dt <- unique(dt, by = "ID")
        # Sentinel 0 means "no data" — replace with NA for downstream modelling
        dt[get(col) == 0, (col) := NA]
        dt
})
names(typologies) <- typo_labels   # carry labels for easier debugging

# Load all biological datasets into memory
bio.list       <- lapply(files, readRDS)
bio.list       <- lapply(bio.list, setDT)   # ensure data.table throughout
n_bio_datasets <- length(bio.list)

# Pre-allocate the output accumulator: a list of lists, one per organism group.
out_list <- vector("list", n_bio_datasets)
names(out_list) <- bio.names
for (j in seq_len(n_bio_datasets)) {
        out_list[[j]] <- vector("list", length(catchments))
}

# loop over catchment tiles -----------------------------------------------
for (i in seq_along(catchments)) {
        
        message("Catchment tile ", i, " / ", length(catchments))
        
        i.vec <- read_parquet(catchments[i])
        setDT(i.vec)
        
        # Subset each biological dataset to segments present in this tile
        i.data <- lapply(bio.list, function(x) x[ID %in% i.vec$ID])
        
        # Subset each typology to segments present in this tile
        i.typology <- lapply(typologies, function(x) x[ID %in% i.vec$ID])
        
        # Join environmental predictors (all.x = TRUE → left join, keep all bio rows)
        out <- lapply(i.data, function(x) merge(x, i.vec, by = "ID", all.x = TRUE))
        
        # Iteratively join each typology onto the already-merged tables
        for (j in seq_along(i.typology)) {
                out <- lapply(out, function(x) merge(x, i.typology[[j]], by = "ID", all.x = TRUE))
        }
        
        # Store this tile's results for each organism group
        for (j in seq_len(n_bio_datasets)) {
                out_list[[j]][[i]] <- out[[j]]
        }
        
        rm(i.vec, i.data, i.typology, out)
}

# Combine per-tile results ------------------------------------------------
# rbindlist(fill = TRUE) handles tiles where some columns may be absent
final_results <- lapply(out_list, function(x) rbindlist(x, fill = TRUE))

# Retain only presence records (organismQuantity > 0) …
final_results <- lapply(final_results, function(x) x[organismQuantity != 0])

# … and label them as presences for downstream modelling
final_results <- lapply(final_results, function(x) x[, PA := 1L])

# save to file ------------------------------------------------------------
out_dir <- "data/biota"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

invisible(lapply(seq_len(n_bio_datasets), function(j) {
        saveRDS(
                final_results[[j]],
                file.path(out_dir, paste0("02_", bio.names[j], "_w_environment.rds"))
        )
        message("Saved: 02_", bio.names[j], "_w_environment.rds")
}))
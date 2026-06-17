################################################################################
# Description:        Evaluates the performance of empirical (k-means) 
#                     typologies on all four taxonomic groups  
################################################################################

# 1. Setup ---------------------------------------------------------------------

library(data.table)
library(tidyverse)
library(sf)
library(mapview)
library(ggplot2)
library(ggh4x)
library(HDInterval)
library(magrittr)
library(flextable)
library(officer)

# Configuration 
base_path <- # path to your folders
# 2. Load data -----------------------------------------------------------------

groups <- list(
        D = list(dir = "diatom_folder", taxon = "diatoms"),
        F = list(dir = "fish_folder", taxon = "fish"),
        I = list(dir = "invertebrate_folder", taxon = "invertebrates"),
        M = list(dir = "macrophyte_folder", taxon = "macrophytes")
)

# Helper Function 
read_dir <- function(path) {
        files <- list.files(path, full.names = TRUE)
        rbindlist(lapply(files, readRDS))
}

# Read all data
data <- schemes <- schem <- ss <- tc <- list()

for (g in names(groups)) {
        print(g)
        grp <- groups[[g]]
        root <- file.path(base_path, grp$dir)
        
        # Evaluations (note: check "Origianls" typo in your M path)
        data[[g]]    <- read_dir(file.path(root, "data/007_evaluations_empirical"))
        # Schemes
        schemes[[g]] <- readRDS(file.path(root, "data/biota", 
                                          paste0("03_", grp$taxon, "_scheme.rds")))
        # Taxonomic resolution
        schem[[g]]   <- read_dir(file.path(root, "data/misc/taxonomic_resolution"))
        # Spatial scale
        ss[[g]]      <- read_dir(file.path(root, "data/misc/spatial_scale"))
        # Taxa counts
        tc[[g]]      <- read_dir(file.path(root, "data/misc/taxa_counts"))
}



# 3. Data Pre-processing  ------------------------------------------------------

# add taxon column
data[[1]][, taxon := "diatoms"]
data[[2]][, taxon := "fish"]
data[[3]][, taxon := "invertebrates"]
data[[4]][, taxon := "macrophytes"]

data <- rbindlist(data, fill =T )

# 4. Data Combination  ------------------------------------------------------
schemes <- rbindlist(schemes)
schem <- rbindlist(schem)
ss <- rbindlist(ss)
tc <- rbindlist(tc)

comb <- tc[ss, on = "scheme_id"]
comb[, c("taxon_group", "data_set", "n_taxa") := NULL]
schem <- schem[, c("taxon", "data.set") := NULL]
schem[, scheme_id := str_remove(scheme_id, "^diatoms_")]
schem[, scheme_id := str_remove(scheme_id, "^fish_")]
schem[, scheme_id := str_remove(scheme_id, "^invertebrates_")]
schem[, scheme_id := str_remove(scheme_id, "^macrophytes_")]
comb <- comb[schem, on = "scheme_id"]
schemes <- schemes[, c("taxon", "data.set", "eventYear", "samples", "sample_type") := NULL]
comb <- comb[schemes, on = "scheme_id"]
comb = rename(comb, "scheme" = "scheme_id")

data <- data[comb, on = "scheme"]

# save to file ------------------------------------------------------------
saveRDS(data, "data/results/results_established_typologies.rds")


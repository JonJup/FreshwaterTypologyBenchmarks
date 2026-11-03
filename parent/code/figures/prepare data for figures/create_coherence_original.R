################################################################################*
# Script Name:        create_coherence_original.R
# Description:        Evaluates the type coherence of all the typologies
#                     build from the original, unmodified environmental variables 
#                     for all four taxonomic groups  
################################################################################*

# ==============================================================================**
# 1. Setup ----
# ==============================================================================**
# 
# 
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


# Helper Function 
read_dir <- function(path) {
        files <- list.files(path, full.names = TRUE)
        rbindlist(lapply(files, readRDS))
}

base_path <- "//turso-fs.cs.helsinki.fi/turso/home/jonjupke"

# ==============================================================================*
# 2. Data Input ----
# ==============================================================================*


groups <- list(
        D = list(dir = "pulseD", taxon = "diatoms"),
        F = list(dir = "pulseF", taxon = "fish"),
        I = list(dir = "pulseI", taxon = "invertebrates"),
        M = list(dir = "pulseM", taxon = "macrophytes")
)

# Read all data
data <- schemes <- schem <- ss <- tc <- list()

for (g in names(groups)) {
        print(g)
        grp <- groups[[g]]
        root <- file.path(base_path, grp$dir)
        
        # Evaluations (note: check "Origianls" typo in your M path)
        data[[g]]    <- read_dir(file.path(root, "data/007_evaluations_originals"))
        # # Schemes
        schemes[[g]] <- readRDS(file.path(root, "data/000_biota",
                                          paste0("03_", grp$taxon, "_scheme.rds")))
        # Taxonomic resolution
        schem[[g]]   <- read_dir(file.path(root, "data/misc/taxonomic_resolution"))
        # Spatial scale
        ss[[g]]      <- read_dir(file.path(root, "data/misc/spatial_scale"))
        # Taxa counts
        tc[[g]]      <- read_dir(file.path(root, "data/misc/taxa_counts"))
}


# ==============================================================================*
# 3. Data Pre-processing ----
# ==============================================================================*

# add taxon 
data[[1]][, taxon := "diatoms"]
data[[2]][, taxon := "fish"]
data[[3]][, taxon := "invertebrates"]
data[[4]][, taxon := "macrophytes"]

data <- rbindlist(data, fill =T )

# combine additional information 
schemes <- rbindlist(schemes)
schem <- rbindlist(schem)
ss <- rbindlist(ss)
tc <- rbindlist(tc)

comb <- tc[ss, on = "scheme_id"]
comb[, c("taxon_group", "data_set", "n_taxa") := NULL]
schem <- schem[, c("taxon", "data.set") := NULL]

comb <- comb[schem, on = "scheme_id"]
schemes <- schemes[, c("taxon", "data.set", "eventYear", "samples", "sample_type") := NULL]
comb <- comb[schemes, on = "scheme_id"]

data <- data[comb, on = "scheme_id"]


# ==============================================================================*
# 4. SAVE TO FILE ------------------------------------------------------------
# ==============================================================================*
setwd(rstudioapi::getActiveProject())
saveRDS(data, "data/figures/coherence_original.rds")


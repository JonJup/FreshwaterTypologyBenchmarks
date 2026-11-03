## ============================================================================
## Compile Data for Figure S23 
## ============================================================================

# ==============================================================================
#  1. Setup
# ==============================================================================
library(data.table)
library(tidyverse)
library(magrittr)
library(mgcv)   
library(patchwork)   # combining plots
library(viridis)     # color scales
library(scales)      # axis formatting
library(pROC)
library(tidytext)

# ==============================================================================
#  2. Load Data 
# ==============================================================================

# Configuration 
base_path <- ""

groups <- list(
        D = list(dir = "pulseD", taxon = "diatoms"),
        F = list(dir = "pulseF", taxon = "fish"),
        I = list(dir = "pulseI", taxon = "invertebrates"),
        M = list(dir = "pulseM", taxon = "macrophytes")
)

# Helper Function 
read_dir <- function(path) {
        files <- list.files(path, full.names = TRUE)
        rbindlist(lapply(files, readRDS))
}

# Read all data
data <- schemes <- schem <- ss <- tc <- list()

for (g in names(groups)) {
        grp <- groups[[g]]
        root <- file.path(base_path, grp$dir)
        
        # Evaluations 
        data[[g]]    <- read_dir(file.path(root, "data/a"))
        # Schemes
        schemes[[g]] <- readRDS(file.path(root, "data/000_biota", 
                                          paste0("03_", grp$taxon, "_scheme.rds")))
        # Taxonomic resolution
        schem[[g]]   <- read_dir(file.path(root, "data/misc/taxonomic_resolution"))
        # Spatial scale
        ss[[g]]      <- read_dir(file.path(root, "data/misc/spatial_scale"))
        # Taxa counts
        tc[[g]]      <- read_dir(file.path(root, "data/misc/taxa_counts"))
}

# ==============================================================================
#  3. Combine Data 
# ==============================================================================

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

tc[, c("taxon_group", "data_set", "scheme_id") := NULL]
tc[, scheme_id := ss$scheme_id]

ss[, c("taxon", "data.set", "sample_type") := NULL]

comb <- tc[ss, on = "scheme_id"]


schem <- schem[, c("taxon", "data.set", "samples") := NULL]
schem[, scheme_id := str_remove(scheme_id, "^diatoms_")]
schem[, scheme_id := str_remove(scheme_id, "^fish_")]
schem[, scheme_id := str_remove(scheme_id, "^invertebrates_")]
schem[, scheme_id := str_remove(scheme_id, "^macrophytes_")]
comb <- comb[schem, on = "scheme_id"]

# schemes <- schemes[, c("taxon", "data.set", "eventYear", "samples", "sample_type") := NULL]
# comb <- comb[schemes, on = "scheme_id"]

data <- data[comb, on = "scheme_id"]


rm(comb, groups, grp, schem, schemes, ss, tc); gc()

data <- data[!is.na(taxon)]
data <- data[!metric %in% c("PERMANOVA F")]
# ==============================================================================
#  4. Subset data 
# ==============================================================================


# limsSep <- c(
#         quantile(data$separation, 0.01, na.rm = T),
#         quantile(data$separation, 0.99, na.rm = T)
# )
# limsCom <- c(
#         quantile(data$compactness, 0.01, na.rm = T),
#         quantile(data$compactness, 0.99, na.rm = T)
# )
# plotData <- data[separation >= limsSep[1] &
#                          separation <= limsSep[2] &
#                          compactness >= limsCom[1] &
#                          compactness <= limsCom[2]]

data2 <- copy(data)

unique(data2$metric)
data2 <- data2[!is.na(metric)]

key_metrics <- c("Classification Strength", "ANOSIM R", "PERMANOVA R2", "AucZeta")
dataKey <- filter(data2, metric %in% key_metrics)
unique(dataKey$metric)
datanoKey <- filter(data2, !metric %in% key_metrics)

# save to file ------------------------------------------------------------
saveRDS(dataKey, "~/projects/pulse/01_wp1/MEE_REVISION_1/data/figures/simulation_results_key_metrics.rds")
saveRDS(datanoKey, "~/projects/pulse/01_wp1/MEE_REVISION_1/data/figures/simulation_results_Nonkey_metrics.rds")



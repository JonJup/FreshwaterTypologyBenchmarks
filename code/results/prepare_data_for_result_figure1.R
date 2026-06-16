## ============================================================================
## Compile Data for the first figure in the result section:
## The distributions of type coherence metrics.
## ============================================================================


#  1. Setup --------------------------------------------------------------------


library(data.table)
library(tidyverse)
library(magrittr)
library(mgcv)
library(patchwork)
library(viridis)
library(scales)
library(pROC)
library(tidytext)


#  2. Load Data ----------------------------------------------------------------


# Configuration
base_path <- "data/"

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
        
        # Evaluations (note: check "Origianls" typo in your M path)
        data[[g]]    <- read_dir(file.path(root, "data/007_evaluations"))
        # Schemes
        schemes[[g]] <- readRDS(file.path(
                root,
                "data/000_biota",
                paste0("03_", grp$taxon, "_scheme.rds")
        ))
        # Taxonomic resolution
        schem[[g]]   <- read_dir(file.path(root, "data/misc/taxonomic_resolution"))
        # Spatial scale
        ss[[g]]      <- read_dir(file.path(root, "data/misc/spatial_scale"))
        # Taxa counts
        tc[[g]]      <- read_dir(file.path(root, "data/misc/taxa_counts"))
}


#  3. Combine Data


# add taxon
data[[1]][, taxon := "diatoms"]
data[[2]][, taxon := "fish"]
data[[3]][, taxon := "invertebrates"]
data[[4]][, taxon := "macrophytes"]
data <- rbindlist(data, fill = T)
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


rm(comb, groups, grp, schem, schemes, ss, tc)
gc()

data <- data[!is.na(taxon)]
data <- data[!metric %in% c("PERMANOVA F")]

data2 <- copy(data)
data2[, metric := fcase(
        metric == "ANOSIM R mean",
        "ANOSIM R<sub>mean</sub>",
        metric == "AucZeta mean",
        "AUCζ<sub>mean</sub>",
        metric == "classification strength",
        "Class. Strength",
        metric == "PERMANOVA R2",
        "PERMANOVA R²",
        metric == "fuzzy_mantel",
        "fuzzy<sub>MANTEL</sub>",
        metric == "PERMANOVA Fuzzy R2",
        "PERMANOVA Fuzzy R²"
)]
unique(data2$metric)
data2 <- data2[!is.na(metric)]

key_metrics <- c("Class. Strength",
                 "ANOSIM R<sub>mean</sub>",
                 "PERMANOVA R²",
                 "AUCζ<sub>mean</sub>")
dataKey <- filter(data2, metric %in% key_metrics)
unique(dataKey$metric)


# save to file ------------------------------------------------------------
saveRDS(dataKey, "data/results/results_simulated_typologies.rds")


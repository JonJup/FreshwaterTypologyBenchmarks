################################################################################
# Script Name:        add_env_biota.R
# Description:        Here, we add environmental data and typologies to the biological data.
#
# Author:             Jonathan Jupke
# Date Created:       2024-11-26
# Last Modified:      2025-09-18
#
# R Version:          R 4.5.1
# Required Packages:  package1, package2
#
# Notes:              Any notes or assumptions
################################################################################

# setup -------------------------------------------------------------------
setwd(rstudioapi::getActiveProject())
library(arrow)
library(data.table)

# load and prepare data -------------------------------------------------------------------------
files      <- 
        list.files("data/biota", pattern = "01_", full.names = T)
bio.names  <- 
        sapply(files, function (x) sub(x = x, pattern = "data/biota/01_", replacement = "")) |>
        sapply(function(x) sub(x = x, pattern = "_w_catchment_id.rds", replacement = ""))
catchments <- list.files("data/catchments", full.names = T) 
typologies <- list.files("data/typologies", full.names = T)
typologies <- lapply(typologies, readRDS)

# Fix typology names 
names(typologies[[1]])[2] <- "FEOW"
names(typologies[[2]])[2] <- "GLORiC"
names(typologies[[3]])[2] <- "BGR"
names(typologies[[4]])[2] <- "BRT"
names(typologies[[5]])[2] <- "EnZ"
names(typologies[[6]])[2] <- "HER"
names(typologies[[7]])[2] <- "IFE"


mode_func <- function(x) {
        ux <- unique(x[!is.na(x)])
        if(length(ux) == 0) return(0)
        ux[which.max(tabulate(match(x, ux)))]
}
typologies[[1]][, FEOW := mode_func(FEOW), by = c("ID")]
typologies[[1]] <- unique(typologies[[1]], by = "ID")
typologies[[1]][FEOW == 0, FEOW := NA]
typologies[[2]][, GLORiC := mode_func(GLORiC), by = c("ID")]
typologies[[2]] <- unique(typologies[[2]], by = "ID")
typologies[[2]][GLORiC == 0, GLORiC := NA]
typologies[[3]][, BGR := mode_func(BGR), by = c("ID")]
typologies[[3]] <- unique(typologies[[3]], by = "ID")
typologies[[3]][BGR == 0, BGR := NA]
typologies[[4]][, BRT := mode_func(BRT), by = c("ID")]
typologies[[4]] <- unique(typologies[[4]], by = "ID")
typologies[[4]][BRT == 0, BRT := NA]
typologies[[5]][, EnZ := mode_func(EnZ), by = c("ID")]
typologies[[5]] <- unique(typologies[[5]], by = "ID")
typologies[[5]][EnZ == 0, EnZ := NA]
typologies[[6]][, HER := mode_func(HER), by = c("ID")]
typologies[[6]] <- unique(typologies[[6]], by = "ID")
typologies[[6]][HER == 0, HER := NA]
typologies[[7]][, IFE := mode_func(IFE), by = c("ID")]
typologies[[7]] <- unique(typologies[[7]], by = "ID")
typologies[[7]][IFE == 0, IFE := NA]

bio.list <- lapply(files, readRDS)
n_bio_datasets <- length(bio.list)
out_list <- vector("list", n_bio_datasets)
names(out_list) <- paste0("out", 1:n_bio_datasets)

# loop over catchments  ---------------------------------------------------
for (i in 1:length(catchments)) {
        
        i.vec      <- read_parquet(catchments[i])
        i.data     <- lapply(bio.list,   function(x) x[ID %in% i.vec$ID])
        i.typology <- lapply(typologies, function(x) x[ID %in% i.vec$ID])
        out <- lapply(i.data, function(x) merge(x, i.vec, by = "ID", all.x = T))
        for (j in 1:length(i.typology)){
                out <- lapply(out, function(x) merge(x, y = i.typology[[j]], by = "ID", all.x = T))
        }
        # Store results for each biological data set
        for (j in 1:n_bio_datasets) {
                out_list[[j]][[i]] <- out[[j]]
        }
        rm(list = ls()[grepl("^i\\.", x = ls())])
        rm(i)
}

# Combine results
final_results <- lapply(out_list, function(x) rbindlist(x, fill = T))

final_results <- lapply(final_results, function(x) x[organismQuantity != 0])
final_results <- lapply(final_results, function(x) x[, PA := 1])

# save to file ----------------------------------------------------------------------
lapply(1:n_bio_datasets,
       function(x)
               saveRDS(
                        final_results[[x]], 
                        paste0("data/biota/02_",bio.names[x],"_w_environment.rds")
               )
       )


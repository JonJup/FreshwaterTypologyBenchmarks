# combine evaluations 

library(data.table)
library(tidyverse)

setwd("E://Arbeit/Projekte/02_ongoing/PULSE/wp1/data/")

# ==============================================================================
#  2. Load Data 
# ==============================================================================

filesD <- list.files("pulseD/data/007_evaluations/", full.names =T)
filesF <- list.files("pulseF/data/007_evaluations/", full.names =T)
filesM <- list.files("pulseI/data/007_evaluations/", full.names =T)
filesI <- list.files("pulseM/data/007_evaluations/", full.names =T)

dataD <- lapply(filesD, readRDS) %>% rbindlist %>%  mutate(taxon = "diatoms")
dataF <- lapply(filesF, readRDS) %>% rbindlist %>%  mutate(taxon = "fish")
dataM <- lapply(filesM, readRDS) %>% rbindlist %>%  mutate(taxon = "invertebrates")
dataI <- lapply(filesI, readRDS) %>% rbindlist %>%  mutate(taxon = "macrophytes")

data <- rbindlist(list(dataD, dataF, dataI, dataM))

rm(dataD, dataI, dataF, dataM, filesD, filesM, filesI, filesF)
gc()
setwd(rstudioapi::getActiveProject())
saveRDS(data, "data/results/combined_data.rds")

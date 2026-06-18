# combine evaluations 


# 1. Setup ----------------------------------------------------------------

library(data.table)
library(tidyverse)


#  2. Load Data ----------------------------------------------------------------


filesD <- list.files("diatoms_folder/data/007_evaluations/", full.names =T)
filesF <- list.files("fish_folder/data/007_evaluations/", full.names =T)
filesI <- list.files("invertebrate_folder/data/007_evaluations/", full.names =T)
filesM <- list.files("macrophytes_folder/data/007_evaluations/", full.names =T)


# 3. Combine data ---------------------------------------------------------


dataD <- lapply(filesD, readRDS) %>% rbindlist %>%  mutate(taxon = "diatoms")
dataF <- lapply(filesF, readRDS) %>% rbindlist %>%  mutate(taxon = "fish")
dataM <- lapply(filesM, readRDS) %>% rbindlist %>%  mutate(taxon = "invertebrates")
dataI <- lapply(filesI, readRDS) %>% rbindlist %>%  mutate(taxon = "macrophytes")

data <- rbindlist(list(dataD, dataF, dataI, dataM))


# 4. Save to file  --------------------------------------------------------

saveRDS(data, "parent/data/results/combined_data.rds")

rm(dataD, dataI, dataF, dataM, filesD, filesM, filesI, filesF)
gc()


# combine evaluations 


# 1. Setup ----------------------------------------------------------------

library(data.table)
library(tidyverse)


data_dir <- ""

#  2. Load Data ----------------------------------------------------------------


filesD <- list.files(file.path(data_dir, "pulseD/data/007_evaluations/"), full.names =T)
filesF <- list.files(file.path(data_dir, "pulseF/data/007_evaluations/"), full.names =T)
filesI <- list.files(file.path(data_dir, "pulseI/data/007_evaluations/"), full.names =T)
filesM <- list.files(file.path(data_dir, "pulseM/data/007_evaluations/"), full.names =T)


# 3. Combine data ---------------------------------------------------------


dataD <- lapply(filesD, readRDS) %>% rbindlist %>% mutate(taxon = "diatoms")
dataF <- lapply(filesF, readRDS) %>% rbindlist %>% mutate(taxon = "fish")
dataI <- lapply(filesM, readRDS) %>% rbindlist %>% mutate(taxon = "invertebrates")
dataM <- lapply(filesI, readRDS) %>% rbindlist %>% mutate(taxon = "macrophytes")

data <- rbindlist(list(dataD, dataF, dataI, dataM))


# 4. Save to file  --------------------------------------------------------

saveRDS(data, "data/results/combined_data.rds")

rm(dataD, dataI, dataF, dataM, filesD, filesM, filesI, filesF)
gc()

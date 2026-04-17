# ==============================================================================
#  Species Specific results of HMSC
# ==============================================================================

# ==============================================================================
#  1. Setup
# ==============================================================================
library(data.table)
library(tidyverse)
library(magrittr)
library(wesanderson)
library(tidytext)


# ==============================================================================
#  2. Load Data 
# ==============================================================================

filesD <- list.files("pulseD/data/005_variation_partitioning/", full.names =T)
filesF <- list.files("pulseF/data/005_variation_partitioning/", full.names =T)
filesM <- list.files("pulseM/data/005_variation_partitioning/", full.names =T)
filesI <- list.files("pulseI/data/005_variation_partitioning/", full.names =T)

dataD <- lapply(filesD, readRDS)
dataF <- lapply(filesF, readRDS)
dataI <- lapply(filesI, readRDS)
dataM <- lapply(filesM, readRDS)

# schemes 
schemeD <- readRDS("pulseD/data/000_biota/03_diatoms_scheme.rds")
schemeF <- readRDS("pulseF/data/000_biota/03_fish_scheme.rds")
schemeI <- readRDS("pulseI/data/000_biota/03_invertebrates_scheme.rds")
schemeM <- readRDS("pulseM/data/000_biota/03_macrophytes_scheme.rds")

# only full schemes 
idD <- schemeD[sample_type == "full", scheme_id]
idF <- schemeF[sample_type == "full", scheme_id]
idI <- schemeI[sample_type == "full", scheme_id]
idM <- schemeM[sample_type == "full", scheme_id]

vpD <- map(dataD, "VP") %>% rbindlist() %>% mutate(group = "diatoms") 
vpF <- map(dataF, "VP") %>% rbindlist() %>% mutate(group = "fish") 
vpI <- map(dataI, "VP") %>% rbindlist() %>% mutate(group = "invertebrates") 
vpM <- map(dataM, "VP") %>% rbindlist() %>% mutate(group = "macrophytes") 

vpD <- vpD[scheme %in% idD]
vpF <- vpF[scheme %in% idF]
vpI <- vpI[scheme %in% idI]
vpM <- vpM[scheme %in% idM]

# most common taxa
mcD <- table(vpD$taxon) %>% sort %>%  tail(20) %>%  names
mcF <- table(vpF$taxon) %>% sort %>%  tail(20) %>%  names
mcI <- table(vpI$taxon) %>% sort %>%  tail(20) %>%  names
mcM <- table(vpM$taxon) %>% sort %>%  tail(20) %>%  names

vpD <- vpD[taxon %in% mcD]
vpF <- vpF[taxon %in% mcF]
vpI <- vpI[taxon %in% mcI]
vpM <- vpM[taxon %in% mcM]



setwd(rstudioapi::getActiveProject())

drme <- vpD[, mean(scaled_values, na.rm = T), by = c("taxon", "driver")] %>% mutate(V1 = round(V1, 2)) %>% dcast(formula=taxon~driver, value.var = "V1") %>% rename(bio_mean = bio, env_mean = env, space_mean = space, stochastic_mean = stochastic)
drsd <- vpD[, sd(scaled_values, na.rm = T), by = c("taxon", "driver")] %>% mutate(V1 = round(V1, 2)) %>% dcast(formula=taxon~driver, value.var = "V1") %>% rename(bio_sd = bio, env_sd = env, space_sd = space, stochastic_sd = stochastic)
r2me <- vpD[, mean(r2, na.rm = T), by = c("taxon")] %>% mutate(V1 = round(V1, 2)) %>% rename(R2_mean = V1)
r2sd <- vpD[, sd(r2, na.rm = T), by = c("taxon")] %>% mutate(V1 = round(V1, 2)) %>% rename(R2_sd = V1)

comb <- drme %>% 
        left_join(drsd) %>% 
        left_join(r2me) %>% 
        left_join(r2sd) %>% 
        mutate(taxon_group = "diatoms")
saveRDS(comb, "data/supplement/taxon_spec_hmsc_diatoms.rds")
drme <- vpF[, mean(scaled_values, na.rm = T), by = c("taxon", "driver")] %>% mutate(V1 = round(V1, 2)) %>% dcast(formula=taxon~driver, value.var = "V1") %>% rename(bio_mean = bio, env_mean = env, space_mean = space, stochastic_mean = stochastic)
drsd <- vpF[, sd(scaled_values, na.rm = T), by = c("taxon", "driver")] %>% mutate(V1 = round(V1, 2)) %>% dcast(formula=taxon~driver, value.var = "V1") %>% rename(bio_sd = bio, env_sd = env, space_sd = space, stochastic_sd = stochastic)
r2me <- vpF[, mean(r2, na.rm = T), by = c("taxon")] %>% mutate(V1 = round(V1, 2)) %>% rename(R2_mean = V1)
r2sd <- vpF[, sd(r2, na.rm = T), by = c("taxon")] %>% mutate(V1 = round(V1, 2)) %>% rename(R2_sd = V1)

comb <- drme %>% 
        left_join(drsd) %>% 
        left_join(r2me) %>% 
        left_join(r2sd) %>% 
        mutate(taxon_group = "fish")
saveRDS(comb, "data/supplement/taxon_spec_hmsc_fish.rds")
drme <- vpI[, mean(scaled_values, na.rm = T), by = c("taxon", "driver")] %>% mutate(V1 = round(V1, 2)) %>% dcast(formula=taxon~driver, value.var = "V1") %>% rename(bio_mean = bio, env_mean = env, space_mean = space, stochastic_mean = stochastic)
drsd <- vpI[, sd(scaled_values, na.rm = T), by = c("taxon", "driver")] %>% mutate(V1 = round(V1, 2)) %>% dcast(formula=taxon~driver, value.var = "V1") %>% rename(bio_sd = bio, env_sd = env, space_sd = space, stochastic_sd = stochastic)
r2me <- vpI[, mean(r2, na.rm = T), by = c("taxon")] %>% mutate(V1 = round(V1, 2)) %>% rename(R2_mean = V1)
r2sd <- vpI[, sd(r2, na.rm = T), by = c("taxon")] %>% mutate(V1 = round(V1, 2)) %>% rename(R2_sd = V1)

comb <- drme %>% 
        left_join(drsd) %>% 
        left_join(r2me) %>% 
        left_join(r2sd) %>% 
        mutate(taxon_group = "invertebrates")
saveRDS(comb, "data/supplement/taxon_spec_hmsc_invertebrates.rds")
drme <- vpM[, mean(scaled_values, na.rm = T), by = c("taxon", "driver")] %>% mutate(V1 = round(V1, 2)) %>% dcast(formula=taxon~driver, value.var = "V1") %>% rename(bio_mean = bio, env_mean = env, space_mean = space, stochastic_mean = stochastic)
drsd <- vpM[, sd(scaled_values, na.rm = T), by = c("taxon", "driver")] %>% mutate(V1 = round(V1, 2)) %>% dcast(formula=taxon~driver, value.var = "V1") %>% rename(bio_sd = bio, env_sd = env, space_sd = space, stochastic_sd = stochastic)
r2me <- vpM[, mean(r2, na.rm = T), by = c("taxon")] %>% mutate(V1 = round(V1, 2)) %>% rename(R2_mean = V1)
r2sd <- vpM[, sd(r2, na.rm = T), by = c("taxon")] %>% mutate(V1 = round(V1, 2)) %>% rename(R2_sd = V1)

comb <- drme %>% 
        left_join(drsd) %>% 
        left_join(r2me) %>% 
        left_join(r2sd) %>% 
        mutate(taxon_group = "macrophytes")
saveRDS(comb, "data/supplement/taxon_spec_hmsc_macrophytes.rds")
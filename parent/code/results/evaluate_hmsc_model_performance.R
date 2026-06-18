
#  Evaluation of Hmsc Model Performance



#  1. Setup -------------------------------------------------------

library(data.table)
library(tidyverse)
library(magrittr)


#  2. Load Data  -------------------------------------------------------


filesD <- list.files("diatoms_folder/data/004_model_fit/", full.names =T)
filesF <- list.files("fish_folder/data/004_model_fit/", full.names =T)
filesI <- list.files("invertebrates_folder/data/004_model_fit/", full.names =T)
filesM <- list.files("macrophytes_folder/data/004_model_fit/", full.names =T)

dataD <- lapply(filesD, readRDS)
dataF <- lapply(filesF, readRDS)
dataI <- lapply(filesI, readRDS)
dataM <- lapply(filesM, readRDS)


#  3. Prepare Data -------------------------------------------------------


dataD %<>% rbindlist() %>% mutate(taxon = "diatoms")
dataF %<>% rbindlist() %>% mutate(taxon = "fish")
dataI %<>% rbindlist() %>% mutate(taxon = "invertebrates")
dataM %<>% rbindlist() %>% mutate(taxon = "macrophytes")

data <- rbindlist(list(dataD, dataF, dataI, dataM))


#  4. Analyses -------------------------------------------------------


# How many passed the psrf
nrow(data)
table(data$taxon)
round(data[psrf_passed == TRUE, .N]/nrow(data),2)
data[psrf_passed == TRUE, .N]
data[AUC_passed == TRUE, .N]

#data[psrf_passed == FALSE, .N, by = "taxon"]
# How many passed the posterior predictive checks?
data[, good_model_n := prevalence_passed + c_score_passed + betaDistr_passed]
data[, good_model := good_model_n >= 2]
data2 <- data[psrf_passed == TRUE]
round(data2[good_model == TRUE, .N]/nrow(data2),2 )
data2[good_model == TRUE, .N]


#  5. Load Detailed Data -------------------------------------------------------

filesD <- list.files("diatoms_folder/data/004_model_fit_detail/", full.names =T)
filesF <- list.files("fish_folder/data/004_model_fit_detail/", full.names =T)
filesM <- list.files("invertebrates_folder/data/004_model_fit_detail/", full.names =T)
filesI <- list.files("macrophytes_folder/data/004_model_fit_detail/", full.names =T)

dataD <- lapply(filesD, readRDS)
dataF <- lapply(filesF, readRDS)
dataI <- lapply(filesI, readRDS)
dataM <- lapply(filesM, readRDS)


#  6. Prepare Data -------------------------------------------------------

dataD %<>% rbindlist() 
dataF %<>% rbindlist() 
dataI %<>% rbindlist() 
dataM %<>% rbindlist() 

data <- rbindlist(list(dataD, dataF, dataI, dataM))


#  7. Analyses -------------------------------------------------------


data[, round(mean(AUC, na.rm = T),2), by = "group"]
data[, round(sd(AUC, na.rm = T)  ,2), by = "group"]
data[, round(mean(TR2, na.rm = T),2), by = "group"]
data[, round(sd(TR2, na.rm = T)  ,2), by = "group"]

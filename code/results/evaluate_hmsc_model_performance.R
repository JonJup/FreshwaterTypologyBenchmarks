# ==============================================================================
#  Evaluation of Hmsc Model Performance
# ==============================================================================

# ==============================================================================
#  1. Setup
# ==============================================================================
library(data.table)
library(tidyverse)
library(magrittr)
setwd("E://Arbeit/Projekte/02_ongoing/PULSE/wp1/data/")

# ==============================================================================
#  2. Load Data 
# ==============================================================================

filesD <- list.files("pulseD/data/004_model_fit/", full.names =T)
filesF <- list.files("pulseF/data/004_model_fit/", full.names =T)
filesI <- list.files("pulseI/data/004_model_fit/", full.names =T)
filesM <- list.files("pulseM/data/004_model_fit/", full.names =T)

dataD <- lapply(filesD, readRDS)
dataF <- lapply(filesF, readRDS)
dataI <- lapply(filesI, readRDS)
dataM <- lapply(filesM, readRDS)

# ==============================================================================
#  3. Prepare Data
# ==============================================================================

dataD %<>% rbindlist() %>% mutate(taxon = "diatoms")
dataF %<>% rbindlist() %>% mutate(taxon = "fish")
dataI %<>% rbindlist() %>% mutate(taxon = "invertebrates")
dataM %<>% rbindlist() %>% mutate(taxon = "macrophytes")

data <- rbindlist(list(dataD, dataF, dataI, dataM))

drop_id <- c(
        paste0("diatoms_0", 329:333), 
        paste0("invertebrates_0", c(374:401, 655:667)), 
        paste0("macrophytes_0", c(244,245))
)
data <- data[!scheme %in% drop_id]


# ==============================================================================
#  4. Analyses
# ==============================================================================

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

# ==============================================================================
#  5. Load Detailed Data 
# ==============================================================================

filesD <- list.files("pulseD/data/004_model_fit_detail/", full.names =T)
filesF <- list.files("pulseF/data/004_model_fit_detail/", full.names =T)
filesM <- list.files("pulseI/data/004_model_fit_detail/", full.names =T)
filesI <- list.files("pulseM/data/004_model_fit_detail/", full.names =T)

dataD <- lapply(filesD, readRDS)
dataF <- lapply(filesF, readRDS)
dataI <- lapply(filesI, readRDS)
dataM <- lapply(filesM, readRDS)

# ==============================================================================
#  3. Prepare Data
# ==============================================================================

dataD %<>% rbindlist() 
dataF %<>% rbindlist() 
dataI %<>% rbindlist() 
dataM %<>% rbindlist() 

data <- rbindlist(list(dataD, dataF, dataI, dataM))

drop_id <- c(
        paste0("diatoms_0", 329:333), 
        paste0("invertebrates_0", c(374:401, 655:667)), 
        paste0("macrophytes_0", c(244,245))
)
data <- data[!scheme_id %in% drop_id]

# ==============================================================================
#  4. Analyses
# ==============================================================================

data[, round(mean(AUC, na.rm = T),2), by = "group"]
data[, round(sd(AUC, na.rm = T)  ,2), by = "group"]
data[, round(mean(TR2, na.rm = T),2), by = "group"]
data[, round(sd(TR2, na.rm = T)  ,2), by = "group"]

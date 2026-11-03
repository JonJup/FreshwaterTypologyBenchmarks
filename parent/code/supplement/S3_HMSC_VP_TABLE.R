# ==============================================================================*
#  Evaluation of Variation Partitioning
# ==============================================================================*

# ==============================================================================*
#  1. Setup ----
# ==============================================================================*
library(data.table)
library(tidyverse)
library(magrittr)
library(wesanderson)
library(tidytext)
library(flextable)
library(officer)

setwd("//turso-fs.cs.helsinki.fi/turso/home/jonjupke/")

# ==============================================================================*
#  2. Load Data ----
# ==============================================================================

filesD <- list.files("pulseD/data/005_variation_partitioning/", full.names =T)
filesF <- list.files("pulseF/data/005_variation_partitioning/", full.names =T)
filesM <- list.files("pulseM/data/005_variation_partitioning/", full.names =T)
filesI <- list.files("pulseI/data/005_variation_partitioning/", full.names =T)

dataD <- lapply(filesD, readRDS)
dataF <- lapply(filesF, readRDS)
dataI <- lapply(filesI, readRDS)
dataM <- lapply(filesM, readRDS)

# ==============================================================================*
#  3. Prepare Data ----
# ==============================================================================*

# ---- Extract variation partitioning results 
vpD <- map(dataD, "VP") %>% rbindlist() %>% mutate(group = "diatoms") 
vpF <- map(dataF, "VP") %>% rbindlist() %>% mutate(group = "fish") 
vpI <- map(dataI, "VP") %>% rbindlist() %>% mutate(group = "invertebrates") 
vpM <- map(dataM, "VP") %>% rbindlist() %>% mutate(group = "macrophytes") 

vpA <- rbindlist(list(vpD, vpF, vpI, vpM))

vpA_summary <- vpA %>%
        group_by(group, driver) %>%
        summarize(mean_value = round(mean(scaled_values, na.rm = TRUE), 2), .groups = 'drop')

vpA_summary <- mutate(
        vpA_summary,
        driver = fcase(
                driver == "bio",
                "co-distributional",
                driver == "env",
                "environmental",
                driver == "space",
                "spatial",
                driver == "stochastic",
                "unexplained"
        )
)

vpA_summary %<>% pivot_wider(id_cols = c(group), names_from = driver,values_from = mean_value)

names(vpA_summary)[1] <- "taxon"
vpA_summary <- vpA_summary[, c(1,3,4,2,5)]


# ==============================================================================*
# save to file ----
# ==============================================================================*

setwd(rstudioapi::getActiveProject())
# saveRDS(vpA_summary, "data/figures/supplement/variationPartitioning_hmsc.rds")

ft <- flextable(head(vpA_summary)) |>
        autofit() |>
        theme_booktabs()

# Save straight to a .docx
save_as_docx(ft, path = "data/supplement/S3_VP_HHMSC_table.docx")
        

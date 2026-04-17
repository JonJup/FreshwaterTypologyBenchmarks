library(flextable)
library(officer)
library(dplyr)

setwd(rstudioapi::getActiveProject())

tsd <- readRDS("data/supplement/taxon_spec_hmsc_diatoms.rds")
tsf <- readRDS("data/supplement/taxon_spec_hmsc_fish.rds")
tsi <- readRDS("data/supplement/taxon_spec_hmsc_invertebrates.rds")
tsm <- readRDS("data/supplement/taxon_spec_hmsc_macrophytes.rds")

tsd <- tsd[, -"taxon_group"]
tsf <- tsf[, -"taxon_group"]
tsi <- tsi[, -"taxon_group"]
tsm <- tsm[, -"taxon_group"]

vars <- c("env", "bio",  "space", "R2")

for(v in vars){
        mean_col <- paste0(v, "_mean")
        sd_col <- paste0(v, "_sd")
        
        tsd[, c(mean_col) := paste0(get(mean_col),"±", get(sd_col))]
        tsf[, c(mean_col) := paste0(get(mean_col),"±", get(sd_col))]
        tsi[, c(mean_col) := paste0(get(mean_col),"±", get(sd_col))]
        tsm[, c(mean_col) := paste0(get(mean_col),"±", get(sd_col))]
        
        tsd[, c(sd_col) := NULL]
        tsf[, c(sd_col) := NULL]
        tsi[, c(sd_col) := NULL]
        tsm[, c(sd_col) := NULL]
}        
        
tsd_t <- flextable(tsd)
tsf_t <- flextable(tsf)
tsi_t <- flextable(tsi)
tsm_t <- flextable(tsm)    
        

# Clean up the tables
tsd_t2 <- tsd_t %>%
        set_header_labels(
                bio_mean = "CoDist",
                env_mean = "Env",
                space_mean = "Space",
                #stochastic_mean = "Stochastic",
                R2_mean = "R²"
        ) %>%  
        autofit() %>%
        set_table_properties(layout = "autofit", width = 1) %>% 
        theme_booktabs() %>%
        set_caption("Fraction of occurrence variability of the twenty most 
                    common diatom taxa by environmental variables (Env), 
                    co-distribution (CoDist), and spatial factors (Space) as 
                    well as Tjur's R² 
                    (R²).") %>%
        font(fontname = "Times New Roman", part = "all") %>%
        fontsize(size = 8, part = "all") %>%
        bold(part = "header")


tsf_t2 <- tsf_t %>%
        set_header_labels(
                bio_mean = "CoDist",
                env_mean = "Env",
                space_mean = "Space",
                #stochastic_mean = "Stochastic",
                R2_mean = "R²"
        ) %>%  
        autofit() %>%
        theme_booktabs() %>%
        set_caption("Fraction of occurrence variability of the twenty most 
                    common fish taxa by environmental variables (Env), 
                    co-distribution (CoDist), and spatial factors (Space), as 
                    well Tjur's R²
                    (R²).") %>%
        font(fontname = "Times New Roman", part = "all") %>%
        fontsize(size = 8, part = "all") %>%
        bold(part = "header")

tsi_t2 <- tsi_t %>%
        set_header_labels(
                bio_mean = "CoDist",
                env_mean = "Env",
                space_mean = "Space",
                #stochastic_mean = "Stochastic",
                R2_mean = "R²"
        ) %>%  
        autofit() %>%
        theme_booktabs() %>%
        set_caption("Fraction of occurrence variability of the twenty most 
                    common invertebrate taxa by co-distribution (CoDist), 
                    environmental variables (Env), and spatial factors (Space), as 
                    well as Tjur's R² 
                    (R²).") %>%
        font(fontname = "Times New Roman", part = "all") %>%
        fontsize(size = 8, part = "all") %>%
        bold(part = "header")

tsm_t2 <- tsm_t %>%
        set_header_labels(
                bio_mean = "CoDist",
                env_mean = "Env",
                space_mean = "Space",
                #stochastic_mean = "Stochastic",
                R2_mean = "R²"
        ) %>%  
        autofit() %>%
        theme_booktabs() %>%
        set_caption("Fraction of occurrence variability of the twenty most 
                    common macrophyte taxa co-distribution (CoDist), 
                    environmental variables (Env), and spatial factors (Space), as 
                    well as Tjur's R² 
                    (R²).") %>%
        font(fontname = "Times New Roman", part = "all") %>%
        fontsize(size = 8, part = "all") %>%
        bold(part = "header")


out <- list(
        "diatoms" = tsd_t2,
        "fish" =  tsf_t2,
        "invertebrates" =  tsi_t2,
        "macrophytes" =  tsm_t2
)
saveRDS(out, "data/supplement/taxon_spec_hmsc_tables.rds")

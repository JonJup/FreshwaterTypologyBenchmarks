# ==============================================================================
# 1. SETUP
# ==============================================================================
library(data.table)
library(tidyverse)
library(flextable)
library(officer)

# ==============================================================================
# 2. LOAD DATA 
# ==============================================================================
data <- readRDS("data/results/results_established_typologies.rds")

# ==============================================================================
# 3. Preparation
# ==============================================================================
data <- data[!is.na(metric)]
data <- data[metric != "PERMANOVA F"]
data <- data[metric != "fuzzy_divergence"]
data <- data[metric != "fuzzy_mantel"]
data <- data[typology != "fuzzy C means"]
data <- data[typology != "kmeans"]

# ==============================================================================
# 4. OVERVIEW TABLE FOR SUPPLEMENTARY MATERIALS
# Create Flextable and save as such 
# ==============================================================================
data[metric == "cs", metric := "classification strength"]

summary_stats <- 
        data %>%  
        filter(!typology %in% c("fuzzy C means", "kmeans")) %>% 
        group_by(taxon, typology, metric) %>% 
        summarise(
                mean = mean(value, na.rm = TRUE),
                sd = sd(value,na.rm = TRUE),
                minimum = range(value, na.rm = TRUE)[1],
                maximum = range(value, na.rm = TRUE)[2]
        )

table_data <- summary_stats %>%
        # Create the mean ± sd text
        mutate(value = sprintf("%.2f ± %.2f", mean, sd)) %>%
        # Select only what we need
        select(taxon, typology, value, metric) %>%
        # Pivot to wide format
        pivot_wider(
                names_from = taxon,
                values_from = value
        )

# Create flextable
ft <- flextable(table_data) %>%
        # Set column names
        set_header_labels(typology = "Typology") %>%
        # Auto-fit columns
        autofit() %>%
        # Add borders
        border_outer(border = fp_border(width = 2)) %>%
        border_inner_h(border = fp_border(width = 1)) %>%
        border_inner_v(border = fp_border(width = 1)) %>%
        # Bold header
        bold(part = "header") %>%
        # Align columns
        align(j = 1, align = "left", part = "all") %>%
        align(j = 2:ncol(table_data), align = "center", part = "all")
saveRDS(ft, "output/tables/original_evaluations.rds")

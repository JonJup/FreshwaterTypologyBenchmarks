### FIGURE 3.1 

# 1. Setup -------------------------------------------------------------------

library(data.table)
library(tidyverse)
library(ggh4x)
library(dplyr)
library(ggtext)
library(patchwork)

zoom_limits <- list(
        "ANOSIM R" = c(-0.1, 0.55),
        "AUCζ" = c(-0.5, 1.2),
        "Class. Strength" = c(-0.05, 0.15),
        "PERMANOVA R²" = c(0, 0.2)
)
setwd(rstudioapi::getActiveProject())

# 2. Load Data ------------------------------------------------------------

data1 <- readRDS("data/results/results_empirical_typologies.rds")
data2 <- readRDS("data/results/results_simulated_typologies.rds")


data3 <- dplyr::filter(data1, typology == "kmeans")
data1 <- dplyr::filter(data1, !typology %in% c("kmeans", NA, "fuzzy C means"))



# 3. Prepare Data ---------------------------------------------------------

combined_data <- bind_rows(
        data2 %>% mutate(dataset_id = "simulated"),
        data3 %>% mutate(dataset_id = "empirical")
)
combined_data <- combined_data[!is.na(metric)]

combined_data[metric == "cs", metric := "Class. Strength"]
combined_data[metric %in% c("AUCζ mean"), metric := "AUCζ"]
combined_data[metric %in% c("ANOSIM R mean"), metric := "ANOSIM R"]
combined_data[metric %in% c("PERMANOVA R2"), metric := "PERMANOVA R²"]
combined_data[metric %in% c("ANOSIM R<sub>mean</sub>"), metric := "ANOSIM R"]
combined_data[metric %in% c("AUCζ<sub>mean</sub>"), metric := "AUCζ"]

keyMetrics <-  c( "ANOSIM R",  "AUCζ", "PERMANOVA R²", "Class. Strength")
combined_data    <- filter(combined_data, metric %in% keyMetrics)

summary_stats <- 
        combined_data %>% 
        group_by(metric, dataset_id) %>%
        summarise(
                mean_value = mean(value, na.rm = TRUE),
                hdi_lower = quantile(value, 0.10, na.rm = TRUE),
                hdi_upper = quantile(value, 0.90, na.rm = TRUE),
                .groups = "drop"
        )


# Define the dodge width so all layers align (0.8 is standard for violins)
pd <- position_dodge(width = 0.8)


# 4. Plot -----------------------------------------------------------------


# Compute outliers per facet × group
outlier_data <- combined_data %>%
        group_by(metric, dataset_id) %>%
        mutate(
                q1 = quantile(value, 0.25),
                q3 = quantile(value, 0.75),
                iqr = q3 - q1,
                is_outlier = value < (q1 - 1.5 * iqr) | value > (q3 + 1.5 * iqr)
        ) %>%
        ungroup() %>%
        filter(is_outlier)

ggplot(combined_data, aes(x = dataset_id, y = value, fill = dataset_id)) +
        geom_hline(yintercept = 0) +
        geom_violin(data = violin_data, alpha = 0.6, trim = FALSE,
                    position = pd, color = NA) +
        
        # Outlier points only
        geom_point(data = outlier_data,
                   aes(color = dataset_id),
                   position = position_jitterdodge(jitter.width = 0.15,
                                                   dodge.width = 0.9),
                   size = 1, alpha = 0.25, shape = 16) +
        
        geom_errorbar(data = summary_stats,
                      aes(y = mean_value, ymin = hdi_lower, ymax = hdi_upper),
                      width = 0.2, linewidth = 0.8, color = "gray30", position = pd) +
        geom_point(data = summary_stats,
                   aes(y = mean_value, group = dataset_id),
                   size = 3, shape = 21, color = "black", fill = "white",
                   stroke = 1, position = pd) +
        facet_wrap(~ metric, scales = "free_y", ncol = 4) +
        scale_fill_manual(values = c("#C97B4A", "#1F4E5C")) +
        scale_color_manual(values = c("#C97B4A", "#1F4E5C")) +
        labs(x = NULL, y = "",
             title = "Distributions of type coherence metrics",
             subtitle = "Comparison of empirical and simulated data sets") +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
              strip.text = element_markdown(size = 15, face = "bold"),
              panel.grid.major.x = element_blank(),
              legend.position = "none",
              plot.subtitle = element_markdown())


# 5. save to file ------------------------------------------------------------


ggsave("output/figures/manuscript/comparative_distribution_coherence_metrics.png")

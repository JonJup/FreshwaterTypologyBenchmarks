### FIGURE RESULTS 1
### OVERALL FIGURE 2 
# ==============================================================================*
## SETUP -----
# ==============================================================================*
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

# ==============================================================================*

# ==============================================================================*
## LOAD DATA -----
# ==============================================================================*


data1 <- readRDS("data/figures/coherence_original.rds")
data2 <- readRDS("data/figures/simulation_results_key_metrics.rds")


# ==============================================================================*
## PREPARE DATA -----
# ==============================================================================*


combined_data <- bind_rows(
        data2 %>% mutate(dataset_id = "modified"),
        data1 %>% mutate(dataset_id = "original")
)


key_metrics <- c(
        "ANOSIM R",
        "AucZeta",
        "PERMANOVA R2",
        "Classification Strength"
)

combined_data    <- dplyr::filter(combined_data, metric %in% key_metrics)

summary_stats <- 
        combined_data %>% 
        group_by(metric, dataset_id) %>%
        summarise(
                mean_value = mean(value, na.rm = TRUE),
                hdi_lower = quantile(value, 0.10, na.rm = TRUE),
                hdi_upper = quantile(value, 0.90, na.rm = TRUE),
                .groups = "drop"
        )


# ==============================================================================*
## PLOT -----
# ==============================================================================*


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

ggplot(combined_data, aes(y = dataset_id, x = value, fill = dataset_id)) +
        geom_hline(yintercept = 0) +
        geom_violin(data = combined_data, alpha = 0.6, trim = FALSE, color = NA) +
        
        # Outlier points only
        geom_point(data = outlier_data,
                   aes(color = dataset_id),
                   position = position_jitterdodge(jitter.width = 0.15,
                                                   dodge.width = 0.9),
                   size = 1, alpha = 0.25, shape = 16) +
        
        geom_errorbar(data = summary_stats,
                      aes(x = mean_value, xmin = hdi_lower, xmax = hdi_upper),
                      width = 0.2, linewidth = 0.8, color = "gray30") +
        geom_point(data = summary_stats,
                   aes(x = mean_value, group = dataset_id),
                   size = 3, shape = 21, color = "black", fill = "white",
                   stroke = 1) +
        facet_wrap(~ metric, scales = "free_x", ncol = 2) +
        scale_fill_manual(name = "Data set", values = c("#C97B4A", "#1F4E5C")) +
        scale_color_manual(name = "Data set", values = c("#C97B4A", "#1F4E5C")) +
        labs(x = NULL, y = NULL) +
        theme_minimal(base_size = 12) +
        geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              strip.text = element_markdown(size = 15, face = "bold"),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              legend.position = "bottom",
              plot.subtitle = element_markdown())


# ==============================================================================*
## SAVE TO FILE -----
# ==============================================================================*

        
ggsave("figures/supplement/S6/S6_CMP_observation_simulation.png")


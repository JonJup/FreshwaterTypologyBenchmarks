# setup -------------------------------------------------------------------
library(data.table)
library(tidyverse)
library(magrittr)
library(broom)
library(HDInterval)
library(tidytext)


# load data ---------------------------------------------------------------
data <- readRDS("data/results/combined_data.rds")


# prep data ---------------------------------------------------------------
selected_data <- data[metric %in% c("PERMANOVA R2", "classification strength", "fuzzy_mantel", "ANOSIM R mean", "AucZeta mean", "PERMANOVA Fuzzy R2")]

selected_data  %<>% mutate(metric = case_when(
        metric == "PERMANOVA R2"  ~ "PERMANOVA R²",
        metric == "ANOSIM R mean" ~ "ANOSIM R",
        metric == "AucZeta mean"  ~ "AUCζ",
        metric == "classification strength"  ~ "Classification Strength",
        metric == "fuzzy_mantel"  ~ "Mantel test",
        metric == "PERMANOVA Fuzzy R2"  ~ "PERMANOVA R² (Fuzzy)"
))

## SPEARMAN WITH FISHER Z
effect_summary <-
        selected_data %>%
        group_by(taxon, metric) %>%
        summarise(
                tidy(cor.test(
                        value, quality, method = "spearman", exact = FALSE
                )),
                n = n(),
                # Capture the sample size for the CI math
                .groups = "drop"
        ) %>%
        mutate(estimate_abs = abs(estimate)) %>% 
        mutate(
                predictor = "quality",
                # Calculate Standard Error of Z
                se_z = 1 / sqrt(n - 3),
                # Calculate 95% CIs and back-transform with tanh()
                conf.low = tanh(atanh(estimate_abs) - 1.96 * se_z),
                conf.high = tanh(atanh(estimate_abs) + 1.96 * se_z)
        )

# plot --------------------------------------------------------------------


effect_summary %>%
        # 1. Reorder 'metric' by 'estimate' within each 'predictor'
        mutate(metric = reorder(metric, estimate_abs, fun = mean)) %>%
        ggplot(aes(x = metric, y = estimate_abs, colour = taxon)) +
        geom_pointrange(
                aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(width = 0.5)
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
        # 2. Use scale_x_reordered to clean up the axis labels
        scale_x_reordered() + 
        coord_flip() +
        labs(x = NULL, y = "Standardized effect size", colour = "Taxon") +
        theme_bw(base_size = 11) +
        theme(
                strip.background = element_rect(fill = "white"),
                strip.text = element_text(face = "bold", size = 16),
                panel.grid.minor = element_blank(),
                legend.position = "bottom",
                legend.text = element_text(size = 12),
                axis.text.x = element_text(size = 12, color = "black"),
                axis.text.y = element_text(size = 8, color = "black"),
                axis.title.x = element_text(size = 14, margin = margin(t = 10)),
                aspect.ratio = 0.6,
                panel.spacing.y=unit(4, "lines")
        )

# save to file ------------------------------------------------------------


ggsave("output/figures/supplement/effect_of_quality_on_metrics.png", dpi = 600, width = 8, height = 7)

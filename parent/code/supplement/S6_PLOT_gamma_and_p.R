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
selected_data <- data[metric %in% c("PERMANOVA R2", 
                                    "Classification Strength",
                                    "Fuzzy Mantel", 
                                    "ANOSIM R", 
                                    "AucZeta",
                                    "PERMANOVA Fuzzy R2")]

## SPEARMAN WITH FISHER Z
effect_summary <-
        selected_data %>%
        group_by(taxon, metric) %>%
        summarise(
                tidy(cor.test(
                        value, gamma_achieved, method = "spearman", exact = FALSE
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
min(e)
# plot --------------------------------------------------------------------


p1 <- effect_summary %>%
        mutate(metric = reorder(metric, estimate)) %>%
        ggplot(aes(x = metric, y = estimate, colour = taxon)) +
        geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
        coord_flip() +
        labs(x = NULL, y = expression("Spearman's" ~ rho ~ "with" ~ gamma), colour = "Taxon") +
        theme_bw(base_size = 11) +
        theme(panel.grid.minor = element_blank(), legend.position = "bottom")

p1 

ggsave("figures/supplement/S6/effect_of_gamma_on_metrics.png",
       dpi = 300, width = 7, height = 4.5)

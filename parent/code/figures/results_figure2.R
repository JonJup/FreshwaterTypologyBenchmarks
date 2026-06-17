### Internal type coherence fuzzy classifications 

# setup -------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(data.table)
library(cowplot)
library(ggh4x)
library(ggrepel)


# load data ---------------------------------------------------------------
data <- readRDS("data/results/combined_data.rds")

key_metrics <- c(
                 "fuzzy_mantel", 
                 "PERMANOVA Fuzzy R2",
                 "PERMANOVA R2")

data2 <- data[metric %in% key_metrics]
all(key_metrics %in% unique(data2$metric))

data2[, metric := case_when(
        metric =="fuzzy_mantel"~ "Mantel test",
        metric =="PERMANOVA Fuzzy R2" ~ "PERMANOVA R2 Fuzzy",
        metric =="PERMANOVA R2"~"PERMANOVA R2 Crisp"
)]

# data2 %>% 
#         filter(metric != "PERMANOVA R2 Crisp") %>% 
#         ggplot(aes(x = quality, y = value)) + 
#         geom_smooth(aes(col = taxon)) + 
#         facet_wrap(.~metric, scales = "free") + 
#         theme_minimal() + 
#         ylab("") + 
#         theme(legend.position = "top")
# 
# ggsave("output/figures/supplement/fuzzy_metric_vs_quality.tiff", dpi = 600, width = 7, height = 5)

## comparison of R2 values 
paired <- dcast(
        data2[metric %in% c("PERMANOVA R2 Crisp", "PERMANOVA R2 Fuzzy")],
        scheme_id ~ metric,
        value.var = "value",
        fun.aggregate = median  # median across iterations / folds if multiple rows per scheme
)

setnames(paired, c("scheme_id",  "Crisp", "Fuzzy"))
paired[, mean_r2 := (Fuzzy + Crisp) / 2]
paired[, delta := Fuzzy - Crisp]
paired[, delta_scaled := delta/mean_r2]
paired[, larger := fifelse(delta > 0, "Fuzzy", "Crisp")]

# --- Bland–Altman style plot -----------------------------------------------

p_ba <- 
        paired %>%  
        filter(!is.na(mean_r2))  %>%  
        ggplot(aes(x = mean_r2, y = delta)) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
        geom_hline(yintercept = mean(paired$delta), linetype = "solid", colour = "black") +
        geom_hline(
                yintercept = mean(paired$delta) + c(-1.96, 1.96) * sd(paired$delta),
                linetype = "dotted", colour = "grey50"
        ) +
        geom_point(aes(colour = larger), alpha = 0.5, size = 1.8) +
        scale_colour_manual(values = c(Crisp = "#C97B4A", Fuzzy = "#1F4E5C"), name = "Larger") +
        labs(
                x = expression(Mean ~ R^2 ~ (Crisp + Fuzzy) / 2),
                y = expression(Delta ~ R^2 ~ (Fuzzy - Crisp))
        ) +
        theme_minimal(base_size = 12)

ggsave("output/figures/manuscript/bland_altman_plot.png", p_ba, width = 7, height = 5, dpi = 300)


# --- 4. Wilcoxon signed-rank test ---------------------------------------------

wt <- wilcox.test(paired$Fuzzy, paired$Crisp, paired = TRUE, conf.int = TRUE)



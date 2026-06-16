# setup -------------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(data.table)


# load data ---------------------------------------------------------------
comp_long <- readRDS("~/projects/pulse/01_wp1/AquaticTypologyBenchmark/data/figures/results_discrimination.rds")


# prepare data ------------------------------------------------------------
taxon_labels <- c(
        "diatoms"        = "Diatoms",
        "fish"           = "Fish",
        "invertebrates"  = "Invertebrates",
        "macrophytes"    = "Macrophytes"
)
comp_long[, taxon_label := taxon_labels[as.character(taxon)]]

criterion_labels <- c(
        "cohens_d"            = "Cohen's *d*",
        "auc"                 = "AUROC",
        "partial_quality_dev" = "Partial *R*²"
)
comp_long[, unique(metric)]

comp_long[, criterion_label := criterion_labels[as.character(criterion)]]
comp_long[, criterion_label := factor(criterion_label, levels = criterion_labels)]

keep_metrics <- unique(comp_long$metric)
plot_data <- comp_long[metric %in% keep_metrics]
plot_data <- plot_data[!metric %in% c("PERMANOVA F")]


plot_data[, metric := case_when(
        metric == "ANOSIM R max"  ~ "ANOSIM R (min)" ,
        metric == "ANOSIM R mean" ~ "ANOSIM R (mean)" ,
        metric == "ANOSIM R min" ~  "ANOSIM R (max)",
        metric == "ANOSIM p max" ~  "ANOSIM p (max)",
        metric == "ANOSIM p min" ~  "ANOSIM p (min)",
        metric == "ANOSIM p mean" ~  "ANOSIM p (mean)",
        metric == "AucZeta max"             ~  "AucZeta (max)",
        metric == "AucZeta mean"            ~  "AucZeta (mean)",
        metric == "AucZeta min"             ~  "AucZeta (min)",
        metric == "PERMANOVA R2"            ~  "PERMANOVA R2",
        metric == "PERMANOVA p"            ~  "PERMANOVA p",
        metric == "classification strength" ~  "Class. Strength",
        metric == "PERMANOVA Fuzzy R2"        ~  "PERMANOVA Fuzzy R2",
        metric == "PERMANOVA Fuzzy p"        ~  "PERMANOVA Fuzzy p",
        metric == "fuzzy_mantel"            ~  "fuzzy Mantel",
        metric == "isa_avg_p"               ~  "IndSpec average p-value",
        metric == "isa_number"              ~  "IndSpec number",
        metric == "isamic"        ~  "ISAMIC",
)]

plot_data <- plot_data[!is.na(metric)]

# Summarise ----------------------------------------------------------------
summary_dt <- plot_data[, .(mean_val = mean(contribution, na.rm = TRUE),
                            sd_val   = sd(contribution, na.rm = TRUE)),
                        by = .(metric, criterion_label)]

# Order metrics by overall mean across criteria
metric_order <- summary_dt[, .(overall = mean(mean_val)), by = metric][order(overall)]
plot_data[,  metric := factor(metric, levels = metric_order$metric)]
summary_dt[, metric := factor(metric, levels = metric_order$metric)]
taxon_cols <- c(
        "Diatoms"        = "#4477AA",
        "Fish"           = "#EE6677",
        "Invertebrates"  = "#228833",
        "Macrophytes"    = "#CCBB44"
)

# plot --------------------------------------------------------------------
p <- ggplot() +
        # geom_point(data = summary_dt,
        #            aes(x = mean_val, y = metric),
        #            size = 3.5, shape = 18, colour = "black") +
        # individual taxon points

        geom_boxplot(data = plot_data, 
                     aes(x=contribution, y = metric)) + 
        geom_point(data = plot_data,
                    aes(x = contribution, y = metric, colour = taxon_label),
                    size = 2.2, alpha = 0.7,
                    #position = position_nudge(y = 0),
                    width = 0.1) +
        # mean ± SD
        # geom_errorbarh(data = summary_dt,
        #                aes(xmin = mean_val - sd_val, xmax = mean_val + sd_val, y = metric),
        #                height = 0.2, colour = "grey30", linewidth = 0.5) +

        facet_wrap(. ~ criterion_label, scales = "free_x", nrow = 1) +
        scale_colour_manual(values = taxon_cols, name = NULL) +
        labs(x = NULL, y = NULL) +
        theme_minimal(base_size = 11) +
        theme(
                strip.text         = ggtext::element_markdown(face = "bold", size = 11),
                panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3),
                panel.grid.minor   = element_blank(),
                panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.3),
                axis.text.y        = element_text(size = 10),
                legend.position    = "bottom",
                plot.margin        = margin(10, 15, 10, 5)
        ) + 
        geom_vline(data = filter(summary_dt, criterion_label=="Cohen's *d*"), aes(xintercept = 0))+
        geom_vline(data = filter(summary_dt, criterion_label=="AUROC"), aes(xintercept = 0.5))+
        geom_vline(data = filter(summary_dt, criterion_label=="Partial *R*²"), aes(xintercept = 0))
p
ggsave("~/projects/pulse/01_wp1/AquaticTypologyBenchmark/output/figures/supplement/results_discrimination.png", width = 12, height = 8, dpi = 300)

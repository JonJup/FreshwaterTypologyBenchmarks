
# 1. setup ---------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(ggtext)
library(ggh4x)

# 2. load data ---------------------------------------------------------------------

dt2 <- readRDS("data/results/qrf_variable_importance.rds")


# 3. Define plotting Function ---------------------------------------------

plot_qrf_importance <- function(data, target_metrics = NULL, output_filename, 
                                facet_cols = 1, axis.text.y.var = 12, 
                                out.height = 8, strip_size = 12, out.width = 6) {
        
        # 1. Filter data (if target_metrics is NULL, keep all data)
        if (!is.null(target_metrics)) {
                dt_filtered <- data %>% filter(metric %in% target_metrics)
        } else {
                dt_filtered <- data
        }
        
        # 2. Calculate ranks and scale
        dt3 <- dt_filtered %>% 
                group_by(taxon, metric) %>%  
                mutate(
                        importance_rank = rank(-importance),
                        n_vars          = n(),
                        rank_scaled     = 1 - (importance_rank - 1) / (n_vars - 1),
                        direction_cat   = case_when(
                                direction_cor >  0.1 ~ "positive",
                                direction_cor < -0.1 ~ "negative",
                                TRUE                 ~ "ambiguous"
                        ),
                        fill_value      = case_when(
                                direction_cat == "positive"  ~  rank_scaled,
                                direction_cat == "negative"  ~ -rank_scaled,
                                direction_cat == "ambiguous" ~  0
                        )
                ) %>% 
                ungroup()
        
        # 3. Apply custom labels
        setDT(dt3)
        dt3[metric == "classification_strength", metric := "Class. Strength"]
        # (Keeping your duplicate overrides as written)
        # ANOSIM 
        dt3[metric == "ANOSIM_R_mean",           metric := "ANOSIM R<sub>mean</sub>"]
        # AUC Zeta variants
        dt3[metric == "AucZeta_mean",            metric := "AUC&zeta;<sub>mean</sub>"]
        # PERMANOVA variants
        dt3[metric == "PERMANOVA_R2",            metric := "PERMANOVA R<sup>2</sup>"]
        dt3[metric == "PERMANOVA_Fuzzy_R2",      metric := "PERMANOVA Fuzzy R<sup>2</sup>"]
        # Classification Strength
        dt3[metric == "classification_strength", metric := "Class. Strength"]
        dt3[metric == "fuzzy_mantel",            metric := "Fuzzy Mantel"]
        
        # Just in case the input data already contains the formatted versions
        dt3[metric == "AUCζ<sub>mean</sub>",     metric := "AUC&zeta;<sub>mean</sub>"]
        dt3[metric == "PERMANOVA R²",            metric := "PERMANOVA R<sup>2</sup>"]
        
        
        dt3[variable == "n_types", variable := "# Types"]
        dt3[variable == "env_asw", variable := "ASW Environment"]
        dt3[variable == "env", variable := "V(env)"]
        dt3[variable == "bio", variable := "V(codist)"]
        dt3[variable == "stochastic", variable := "V(stoch)"]
        dt3[variable == "min_longitude", variable := "Longitude (min)"]
        dt3[variable == "min_latitude", variable := "Latitude (min)"]
        dt3[variable == "median_longitude", variable := "Longitude (med)"]
        dt3[variable == "median_latitude", variable := "Latitude (med)"]
        dt3[variable == "max_longitude", variable := "Longitude (max)"]
        dt3[variable == "max_latitude", variable := "Latitude (max)"]
        dt3[variable == "min_distance", variable := "Distance (min)"]
        dt3[variable == "mean_distance", variable := "Distance (mean)"]
        dt3[variable == "median_distance", variable := "Distance (med)"]
        dt3[variable == "max_distance", variable := "Distance (max)"]
        dt3[variable == "n_taxa", variable := "# Taxa"]
        dt3[variable == "species_rank", variable := "Species Rank"]
        dt3[variable == "genus_rank", variable := "Genus Rank"]
        dt3[variable == "family_rank", variable := "Family Rank"]
        dt3[variable == "samples", variable := "# Samples"]
        dt3[variable == "space", variable := "V(space)"]
        
        # 4. Extract ordering variables
        keyVar <- dt3 %>% 
                filter(direction_cat != "ambiguous") %>% 
                group_by(variable) %>%  
                mutate(sumIm = sum(rank_scaled)) %>% 
                ungroup() %>%  
                arrange(sumIm) %>% 
                setDT() %>% 
                unique(by = "variable") %>% 
                pull(variable)
        
        # 5. Build Plot
        p <- dt3 %>%
                filter(variable %in% keyVar) %>%
                mutate(variable = factor(variable, levels = keyVar)) %>%
                ggplot(aes(x = taxon, y = variable, fill = fill_value)) +
                geom_tile(colour = "white", linewidth = 0.4) +
                facet_wrap( ~ metric, ncol = facet_cols, scales = "free_x") +
                scale_fill_gradient2(
                        low      = "steelblue",
                        mid      = "grey90",
                        high     = "firebrick",
                        midpoint = 0,
                        name     = "",
                        breaks   = c(-1, 0, 1),
                        limits   = c(-1, 1),
                        labels   = c("Important and \nnegative", "Ambiguous or \nunimportant", "Important and \npositive")
                ) +
                labs(x = NULL, y = NULL) +
                theme_minimal(base_size = 11) +
                facetted_pos_scales(
                        x = list(
                                metric == "ANOSIM R<sub>mean</sub>"       ~ scale_x_discrete(), # X-axis labels visible
                                metric == "AUC&zeta;<sub>mean</sub>"      ~ scale_x_discrete(guide = "none"),
                                metric == "Class. Strength"               ~ scale_x_discrete(guide = "none"),
                                metric == "Fuzzy Mantel"                  ~ scale_x_discrete(guide = "none"),
                                metric == "PERMANOVA Fuzzy R<sup>2</sup>" ~ scale_x_discrete(guide = "none"),
                                metric == "PERMANOVA R<sup>2</sup>"       ~ scale_x_discrete(guide = "none")
                        )
                ) + 
                theme(
                        strip.text       = element_markdown(size = strip_size), 
                        #axis.text.x      = element_text(angle = 0, hjust = 1, size = 9),
                       # aspect.ratio     = 2,
                        axis.text.y      = element_markdown(size = axis.text.y.var),
                        panel.grid       = element_blank(),
                        panel.spacing    = unit(0.1, "lines"),
                        #legend.position  = c(.6, -0.15),
                       legend.position = "bottom",
                        legend.key.width = unit(2, "cm"),
                        legend.direction = "horizontal"
                )
        
        # 6. Save Plot
        ggsave(output_filename, plot = p, dpi = 600, width = out.width, height = out.height, bg = "white")
        # Return the plot object in case you want to view it in the RStudio pane
        return(p) 
}



# 4. Create plot -------------------------------------------------------------




# 1. Plot for the SELECTED metrics
selected_metrics <- c("ANOSIM_R_mean", "AucZeta_mean", "classification_strength", "PERMANOVA_R2")

plot_qrf_importance(
        data = dt2, 
        target_metrics = selected_metrics, 
        output_filename = "output/figures/manuscript/results_figure3.png",
        facet_cols = 2,
        out.height = 10,
        out.width = 8
        )


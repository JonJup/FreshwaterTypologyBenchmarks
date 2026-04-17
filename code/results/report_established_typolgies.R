# ==============================================================================
# 1. SETUP
# ==============================================================================
library(data.table)
library(tidyverse)
library(glue)


# ==============================================================================
# 2. LOAD DATA 
# ==============================================================================
#data <- readRDS("data/figures/orginal_results_key_metrics.rds")
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
# 4. Write Report 
# ==============================================================================
# 1. Define a helper function to format the LaTeX math syntax automatically
# Adjust d1 and d2 where needed to control decimal places for mean and sd
fmt <- function(mean_val, sd_val, d1 = 3, d2 = 3) {
        sprintf("$%s\\pm%s$", round(mean_val, d1), round(sd_val, d2))
}

# 2. Pre-calculate the summaries
# PERMANOVA
perm_all <- data[metric == "PERMANOVA R2", .(m = mean(value, na.rm=TRUE), s = sd(value, na.rm=TRUE))]
perm_tax <- data[metric == "PERMANOVA R2", .(m = mean(value, na.rm=TRUE), s = sd(value, na.rm=TRUE)), by = taxon]
perm_typ <- data[metric == "PERMANOVA R2", .(m = mean(value, na.rm=TRUE), s = sd(value, na.rm=TRUE)), by = typology]

# ANOSIM
anosim_all <- data[metric == "ANOSIM R mean", .(m = mean(value, na.rm=TRUE), s = sd(value, na.rm=TRUE))]
# Order taxa highest to lowest for the sequenced sentence
anosim_tax <- data[metric == "ANOSIM R mean", .(m = mean(value, na.rm=TRUE), s = sd(value, na.rm=TRUE)), by = taxon][order(-m)] 
anosim_typ <- data[metric == "ANOSIM R mean", .(m = mean(value, na.rm=TRUE), s = sd(value, na.rm=TRUE)), by = typology]

# CS
cs_all <- data[metric == "cs", .(m = mean(value, na.rm=TRUE), s = sd(value, na.rm=TRUE))]
cs_tax <- data[metric == "cs", .(m = mean(value, na.rm=TRUE), s = sd(value, na.rm=TRUE)), by = taxon]
cs_typ <- data[metric == "cs", .(m = mean(value, na.rm=TRUE), s = sd(value, na.rm=TRUE)), by = typology]

# AUCζ mean
auc_all <- data[metric == "AUCζ mean", .(m = mean(value, na.rm=TRUE), s = sd(value, na.rm=TRUE))]
auc_tax <- data[metric == "AUCζ mean", .(m = mean(value, na.rm=TRUE), s = sd(value, na.rm=TRUE)), by = taxon]
auc_typ <- data[metric == "AUCζ mean", .(m = mean(value, na.rm=TRUE), s = sd(value, na.rm=TRUE)), by = typology]

# 3. Generate the text using glue
# glue requires double brackets {{ }} to output a literal single bracket { }.
report_text <- glue("
For PERMANOVA R², the overall mean was {fmt(perm_all$m, perm_all$s, 2, 2)}.
The highest values were observed for {perm_tax[which.max(m)]$taxon} {fmt(perm_tax[which.max(m)]$m, perm_tax[which.max(m)]$s)} across established typologies and for the {perm_typ[which.max(m)]$typology} {fmt(perm_typ[which.max(m)]$m, perm_typ[which.max(m)]$s)} across taxonomic groups. 
We observed the lowest PERMANOVA R² for {perm_tax[which.min(m)]$taxon} {fmt(perm_tax[which.min(m)]$m, perm_tax[which.min(m)]$s)} and {perm_typ[which.min(m)]$typology} {fmt(perm_typ[which.min(m)]$m, perm_typ[which.min(m)]$s)}.
The standard deviation exceeded the mean for most combinations, indicating strongly skewed distributions and that a substantial fraction of combinations produced values at or below the no-coherence baseline.

Across taxonomic groups and established typologies, the mean ANOSIM R$_{{mean}}$ was {fmt(anosim_all$m, anosim_all$s, 2, 2)} (Figure \\@ref(fig:res1)). 
{anosim_tax$taxon[1]} showed the highest type coherence ({fmt(anosim_tax$m[1], anosim_tax$s[1], 2, 2)}), followed by {anosim_tax$taxon[2]} ({fmt(anosim_tax$m[2], anosim_tax$s[2], 2, 2)}), {anosim_tax$taxon[3]} ({fmt(anosim_tax$m[3], anosim_tax$s[3], 2, 2)}), and {anosim_tax$taxon[4]} ({fmt(anosim_tax$m[4], anosim_tax$s[4], 2, 2)}). 
The {anosim_typ[which.max(m)]$typology} types had the highest type coherence across taxa ({fmt(anosim_typ[which.max(m)]$m, anosim_typ[which.max(m)]$s, 2, 2)}) and the {anosim_typ[which.min(m)]$typology} the lowest ({fmt(anosim_typ[which.min(m)]$m, anosim_typ[which.min(m)]$s, 2, 2)}).

The overall mean classification strength was {fmt(cs_all$m, cs_all$s)}. 
{cs_tax[which.max(m)]$taxon} communities showed the highest coherence with a mean of {fmt(cs_tax[which.max(m)]$m, cs_tax[which.max(m)]$s, 2, 2)}, while {cs_tax[which.min(m)]$taxon} showed the lowest ({fmt(cs_tax[which.min(m)]$m, cs_tax[which.min(m)]$s)}). 
Across taxa, {cs_typ[which.max(m)]$typology} had the highest mean classification strength ({fmt(cs_typ[which.max(m)]$m, cs_typ[which.max(m)]$s, 4, 3)}) and {cs_typ[which.min(m)]$typology} the lowest ({fmt(cs_typ[which.min(m)]$m, cs_typ[which.min(m)]$s)}).

For  AUCζ$_{{mean}}$, the overall mean was {fmt(auc_all$m, auc_all$s, 2, 2)}.
The highest values were observed for {auc_tax[which.max(m)]$taxon} {fmt(auc_tax[which.max(m)]$m, auc_tax[which.max(m)]$s)} across established typologies and for the {auc_typ[which.max(m)]$typology} {fmt(auc_typ[which.max(m)]$m, auc_typ[which.max(m)]$s)} across taxonomic groups. 
We observed the lowest AUCζ$_{{mean}}$ for {auc_tax[which.min(m)]$taxon} {fmt(auc_tax[which.min(m)]$m, auc_tax[which.min(m)]$s)} and {auc_typ[which.min(m)]$typology} {fmt(auc_typ[which.min(m)]$m, auc_typ[which.min(m)]$s)}.
")

# Print the cleanly compiled text
cat(report_text)

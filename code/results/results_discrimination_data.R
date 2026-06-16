
#  1. Setup ----------------------------------------------------------------

library(data.table)
library(tidyverse)
library(magrittr)
library(mgcv)   
library(patchwork)
library(viridis)
library(scales)
library(pROC)
library(tidytext)



#  2. Load Data ----------------------------------------------------------------


base_path <- "data/"

groups <- list(
        D = list(dir = "pulseD", taxon = "diatoms"),
        F = list(dir = "pulseF", taxon = "fish"),
        I = list(dir = "pulseI", taxon = "invertebrates"),
        M = list(dir = "pulseM", taxon = "macrophytes")
)

# Helper Function 
read_dir <- function(path) {
        files <- list.files(path, full.names = TRUE)
        rbindlist(lapply(files, readRDS))
}

# Read all data
data <- schemes <- schem <- ss <- tc <- list()

for (g in names(groups)) {
       
        grp <- groups[[g]]
        root <- file.path(base_path, grp$dir)
        
        # Evaluations 
        data[[g]]    <- read_dir(file.path(root, "data/007_evaluations"))
        # Schemes
        schemes[[g]] <- readRDS(file.path(root, "data/000_biota", 
                                          paste0("03_", grp$taxon, "_scheme.rds")))
        # Taxonomic resolution
        schem[[g]]   <- read_dir(file.path(root, "data/misc/taxonomic_resolution"))
        # Spatial scale
        ss[[g]]      <- read_dir(file.path(root, "data/misc/spatial_scale"))
        # Taxa counts
        tc[[g]]      <- read_dir(file.path(root, "data/misc/taxa_counts"))
}


#  3. Combine Data 


# add taxon 
data[[1]][, taxon := "diatoms"]
data[[2]][, taxon := "fish"]
data[[3]][, taxon := "invertebrates"]
data[[4]][, taxon := "macrophytes"]
data <- rbindlist(data, fill =T )
# combine additional information 
schemes <- rbindlist(schemes)
schem <- rbindlist(schem)
ss <- rbindlist(ss)
tc <- rbindlist(tc)

tc[, c("taxon_group", "data_set", "scheme_id") := NULL]
tc[, scheme_id := ss$scheme_id]

ss[, c("taxon", "data.set", "sample_type") := NULL]

comb <- tc[ss, on = "scheme_id"]


schem <- schem[, c("taxon", "data.set", "samples") := NULL]
schem[, scheme_id := str_remove(scheme_id, "^diatoms_")]
schem[, scheme_id := str_remove(scheme_id, "^fish_")]
schem[, scheme_id := str_remove(scheme_id, "^invertebrates_")]
schem[, scheme_id := str_remove(scheme_id, "^macrophytes_")]
comb <- comb[schem, on = "scheme_id"]

data <- data[comb, on = "scheme_id"]
rm(comb, groups, grp, schem, schemes, ss, tc); gc()

data <- data[!is.na(taxon)]



## 4. PREPARE DATA -------------------------------------------------------------

## 4.1 -------
# # -- Create a summary table per unique typology scheme 
# # Each scheme_id has one separation, one compactness, etc.
# # Multiple metrics per scheme. Aggregate if needed.
# scheme_summary <- unique(data[, .(scheme_id, n_types, quality,
#                                   env_asw, importance, fuzzy_npe, taxon,
#                                   samples)])
# cat("\n=== Unique typology schemes ===\n")
# print(nrow(scheme_summary))

# -- Group variables by role --------------------------------------------------
quality_vars   <- c("quality")
# env_asw, importance, fuzzy_npe are themselves evaluation metrics, 
# not predictors — keep them separate unless they describe the input typology

structure_vars <- c("n_types", "variables")
sample_vars    <- c("samples", "n_taxa")
spatial_vars   <- c("median_distance", "max_distance",
                    "median_latitude", "median_longitude")
resolution_vars <- c("species_rank", "genus_rank", "family_rank", "higher_rank")

# -- Check correlations among predictors first --------------------------------
pred_vars <- c(quality_vars, structure_vars, sample_vars, 
               spatial_vars, resolution_vars, "year")

# Use one row per scheme to avoid inflated correlations
scheme_level <- unique(data[, ..pred_vars])

cor_mat <- cor(scheme_level, use = "pairwise.complete.obs")
cat("High correlations (|r| > 0.7):\n")
high_cor <- which(abs(cor_mat) > 0.7 & upper.tri(cor_mat), arr.ind = TRUE)
for (i in seq_len(nrow(high_cor))) {
        cat(sprintf("  %s — %s: %.2f\n",
                    rownames(cor_mat)[high_cor[i, 1]],
                    colnames(cor_mat)[high_cor[i, 2]],
                    cor_mat[high_cor[i, 1], high_cor[i, 2]]))
}

# drop highly correlated variables
high_cor <- c("max_distance", "species_rank", "higher_rank")
pred_vars <- pred_vars[-which(pred_vars %in% high_cor)]
scheme_level <- unique(data[, ..pred_vars])
rm(high_cor)


## 5. FIT GAMS

# -- Fit full GAM per metric x taxon -----------------------------------------
# This gives you the PARTIAL effect of separation/compactness
# after controlling for everything else

full_results <- list()
partial_r2_list <- list()

for (tx in unique(data$taxon)) {
        for (met in unique(data$metric)) {
                sub <- data[taxon == tx & metric == met]
                if (nrow(sub) < 50) next
                
                # Subsample if large
                if (nrow(sub) > 50000) {
                        sub <- sub[sample(.N, 50000)]
                }
                
                full_fit <- tryCatch(
                        bam(value ~ s(quality, k = 5) +
                                    s(samples, k = 5) +
                                    s(n_taxa, k = 5) +
                                    s(median_distance, k = 5) +
                                    s(median_latitude, k = 5) +
                                    s(n_types, k = 3) +
                                    s(variables, k = 3),
                            data = sub, method = "fREML", discrete = TRUE, nthreads = 4),
                        error = function(e) NULL
                )
                
                reduced_fit <- tryCatch(
                        bam(value ~ s(samples, k = 5) +
                                    s(n_taxa, k = 5) +
                                    s(median_distance, k = 5) +
                                    s(median_latitude, k = 5) +
                                    s(n_types, k = 3) +
                                    s(variables, k = 3),
                            data = sub, method = "fREML", discrete = TRUE, nthreads = 4),
                        error = function(e) NULL
                )
                
                if (!is.null(full_fit) && !is.null(reduced_fit)) {
                        full_dev  <- summary(full_fit)$dev.expl
                        red_dev   <- summary(reduced_fit)$dev.expl
                        partial_dev <- full_dev - red_dev  # partial R² of quality vars
                        
                        full_results[[paste(tx, met)]] <- data.table(
                                taxon = tx,
                                metric = met,
                                full_dev_expl = full_dev,
                                confounders_dev_expl = red_dev,
                                partial_quality_dev = partial_dev,
                                n = nrow(sub)
                        )
                }
        }
}

model_comparison <- rbindlist(full_results)



## 6. WHICH METRIC BEST DISCRIMINATES GOOD vs BAD TYPOLOGIES?


# -- 6a. Define good vs bad typologies ----------------------------------------
# Good = high separation AND high compactness
# Bad  = low separation AND low compactness
# Use median split or quantile-based thresholds

data[, q_q := ecdf(quality)(quality), by = taxon]

data[, quality_class := fcase(
        q_q >= 0.75, "good",
        q_q <= 0.25, "bad",
        default = "intermediate"
)]

# Quality class distribution
print(data[, .N, by = .(taxon, quality_class)][order(taxon, quality_class)])

# We need to reverse the values in IndValP because it is the only variable where large values are bad. 
reverse_var <- c("isa_avg_p","ANOSIM p min","within-type dissimilarity", 
                 "ANOSIM p mean", "ANOSIM p max", "PERMANOVA p", 
                 "PERMANOVA Fuzzy p")
data[metric %in% reverse_var, value := 1-value]

# -- 6b. Cohen's d: standardized effect size good vs bad ---------------------
cohens_d <- data[quality_class %in% c("good", "bad"),
                 {
                         good_vals <- value[quality_class == "good"]
                         bad_vals  <- value[quality_class == "bad"]
                         
                         n1 <- length(good_vals)
                         n2 <- length(bad_vals)
                         
                         m1 <- mean(good_vals, na.rm = TRUE)
                         m2 <- mean(bad_vals, na.rm = TRUE)
                         
                         s_pooled <- sqrt(((n1 - 1) * var(good_vals, na.rm = TRUE) +
                                                   (n2 - 1) * var(bad_vals, na.rm = TRUE)) /
                                                  (n1 + n2 - 2))
                         
                         d <- (m1 - m2) / s_pooled
                         
                         .(cohens_d = d,
                           mean_good = m1,
                           mean_bad = m2,
                           s_pooled = s_pooled,
                           n_good = n1,
                           n_bad = n2)
                 },
                 by = .(taxon, metric)]

# -- 6c. ROC-AUC: discriminatory power as a classifier -----------------------
auc_results <- data[quality_class %in% c("good", "bad"),
                    {
                            label <- as.numeric(quality_class == "good")
                            roc_obj <- 
                                    #tryCatch(
                                    roc(label, value, quiet = TRUE)
                            # error = function(e) NULL
                            #)
                            if (!is.null(roc_obj)) {
                                    .(auc = as.numeric(auc(roc_obj)))
                            } else {
                                    .(auc = NA_real_)
                            }
                    },
                    by = .(taxon, metric)]

# -- 6d. Rank metrics by overall discriminatory power -------------------------
discrimination_summary <- merge(
        cohens_d[, .(taxon, metric, cohens_d)],
        auc_results,
        by = c("taxon", "metric")
)

# Merge in partial R² from the GAM analysis
discrimination_summary <- merge(
        discrimination_summary,
        model_comparison[, .(taxon, metric, partial_quality_dev)],
        by = c("taxon", "metric"),
        all.x = TRUE
)

# remove questionable cases 
discrimination_summary = discrimination_summary[! (taxon == "macrophytes" & metric == "AucZeta mean")]
discrimination_summary = discrimination_summary[! (taxon == "inverterbrates" & metric %in% c("PERMANOVA R2","AucZeta mean"))]

discrimination_summary[, cohens_d := abs(cohens_d)]
discrimination_summary[, `:=`(
        rank_d   = rank(-abs(cohens_d)),
        rank_auc = rank(-auc),
        rank_r2  = rank(-partial_quality_dev)
), by = taxon]
discrimination_summary[, composite := rank_d+ rank_auc+ rank_r2]
comp_long <- melt(discrimination_summary, 
                  id.vars = c("taxon", "metric", "composite"),
                  measure.vars = c("cohens_d", "auc", "partial_quality_dev"),
                  variable.name = "criterion", value.name = "contribution")



# # Rank across all three criteria
# discrimination_summary[, `:=`(
#         rank_d   = rank(-abs(cohens_d)),
#         rank_auc = rank(-auc),
#         rank_r2  = rank(-partial_quality_dev)
# ), by = taxon]

# # mean rank across criteria 
# discrimination_summary[, mean_rank := (rank_d + rank_auc + rank_r2) / 3]
# discrimination_summary[, composite := (scaled_d + scaled_auc + scaled_r2) / 3,
#                        env = list(
#                                scaled_d   = quote(rank(-abs(cohens_d)) / .N),
#                                scaled_auc = quote(rank(-auc) / .N),
#                                scaled_r2  = quote(rank(-partial_quality_dev) / .N)
#                        )]
# comp_long <- melt(discrimination_summary, 
#                   id.vars = c("taxon", "metric", "composite"),
#                   measure.vars = c("cohens_d", "auc", "partial_quality_dev"),
#                   variable.name = "criterion", value.name = "contribution")


saveRDS(comp_long, "~/projects/pulse/01_wp1/AquaticTypologyBenchmark/data/figures/results_discrimination.rds")



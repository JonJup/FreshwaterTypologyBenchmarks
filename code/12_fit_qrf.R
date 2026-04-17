### - Quantile Random Forest Fitting - ###

#===============================================================================
# 1. setup ----
#===============================================================================
args <- commandArgs(trailingOnly = TRUE)

# Default values
metric_name <- NULL 
output_dir  <- NULL

# Parse arguments
for (i in seq_along(args)) {
        if (args[i] == "--metric") {
                metric_name <- args[i + 1]
        } else if (args[i] == "--output") {
                output_dir <- args[i + 1]
        }
}

# Sanity check
if (is.null(metric_name) || is.null(output_dir)) {
        stop("ERROR: Missing --metric or --output argument.")
}

cat("\n========================================\n")
cat("Processing Metric:", metric_name, "\n")
cat("Output directory:", output_dir, "\n")
cat("========================================\n\n")

suppressPackageStartupMessages({
        library(data.table)
        library(tidymodels)
        library(doParallel)
        library(bundle)
        library(finetune)
})

#===============================================================================
# 2. Load data ----
#===============================================================================
cat("Loading evaluation data...\n")
fileList <- list.files("data/007_evaluations/", full.names = TRUE)
evl <- lapply(fileList, readRDS)
evl <- rbindlist(evl)

cat("Loading variation partitioning data...\n")
vp <- list.files("data/005_variation_partitioning/", full.names = TRUE)
vp <- lapply(vp, readRDS)
vp2 <- list()
for (i in seq_along(vp)) {
        x <- vp[[i]]$VP
        y <- x[, mean(scaled_values), by = c("driver")]
        y <- dcast(y, formula = . ~ driver, value.var = "V1")
        y$. <- unique(x$scheme)
        names(y)[1] <- "scheme_id"
        vp2[[i]] <- y
}
vp2 <- rbindlist(vp2)

cat("Loading spatial scale data...\n")
ss <- list.files("data/misc/spatial_scale/", full.names = TRUE)
ss <- lapply(ss, readRDS)
ss <- rbindlist(ss, fill = T)

cat("Loading taxonomic resolution data...\n")
tx <- list.files("data/misc/taxonomic_resolution/", full.names = TRUE)
tx <- lapply(tx, readRDS)
tx <- rbindlist(tx)
if (length(gregexpr("_", tx$scheme_id[1])[[1]]) == 2) {
        tx$scheme_id <- sub("^([A-Za-z]+_)\\1", "\\1", tx$scheme_id)
}


cat("Loading number of taxa...\n")
tn <- list.files("data/misc/taxa_counts/", full.names = TRUE)
tn <- lapply(tn, readRDS)
tn <- rbindlist(tn)
tn[, c("taxon_group", "data_set") := NULL]
tn[, scheme_id := sub("^([A-Za-z0-9]+_)\\1", "\\1", scheme_id)]


#===============================================================================
# --- 3. Prepare data ----
#===============================================================================
cat("Merging datasets...\n")
evl2 <- evl[ss, on = "scheme_id"]
evl2 <- evl2[vp2, on = "scheme_id"]
evl2 <- evl2[tx, on = "scheme_id"]
evl2 <- evl2[tn, on = "scheme_id"]
evl2 <- evl2[scheme_id %in% evl$scheme_id]

# Remove duplicate columns
if (any(grepl("^i\\.", names(evl2)))) {
        drop.id <- names(evl2)[grep("^i\\.", names(evl2))]
        evl2[, (drop.id) := NULL]
}
stopifnot("Row count changed unexpectedly after merge" = 
                  nrow(evl2) == nrow(evl))

# Filter to specific metric
cat("Filtering to metric:", metric_name, "\n")
j.d <- evl2[metric == metric_name, ]

# Check if we have enough data
if (nrow(j.d) < 10) {
        stop(paste("Not enough data for metric:", metric_name, "(only", nrow(j.d), "rows)"))
}

cat("Data contains", nrow(j.d), "rows\n")

#===============================================================================
# 4. Define predictors and model ----
#===============================================================================

predictor_vec <- c(
        "env_asw", "fuzzy_npe",
        "samples", "min_distance", "mean_distance", 
        "median_distance", "max_distance", "median_latitude", "median_longitude",
        "max_latitude", "max_longitude", "min_latitude", "min_longitude",
        "bio", "env", "space", "stochastic", "species_rank",
        "genus_rank", "family_rank", "higher_rank", "n_taxa"
)

if (!metric_name %in% c("fuzzy_divergence", "fuzzy_mantel")){
        predictor_vec <- predictor_vec[!predictor_vec %in% c("fuzzy_npe")]
}

#===============================================================================
# 5. Split data and create folds ----
#===============================================================================

set.seed(1)
j.d <- na.omit(j.d[, .SD, .SDcols = c("value", predictor_vec, "scheme_id")])
j.split <- initial_split(data = j.d)
j.train <- training(j.split)
j.test  <- testing(j.split)
j.folds <- vfold_cv(j.train, v = 5)

cat("Training set:", nrow(j.train), "rows\n")
cat("Test set:", nrow(j.test), "rows\n")

#===============================================================================
# 6. Create recipe and preprocess ----
#===============================================================================
cat("\nCreating recipe...\n")

rf_recipe <- 
        recipe(value ~ ., data = j.train) %>%
        update_role(everything(), new_role = "unused") %>% # Set all to unused first
        update_role(all_of(predictor_vec), new_role = "predictor") %>%
        update_role(value, new_role = "outcome") %>%
        step_zv(all_predictors()) %>%
        step_naomit(all_predictors(), value)


# Prep the recipe
prepped_recipe <- prep(rf_recipe, training = j.train)

# Get final predictors after preprocessing
j.train.processed <- bake(prepped_recipe, new_data = j.train)
final_predictors <- prepped_recipe %>%
        summary() %>%
        filter(role == "predictor") %>%
        pull(variable)

if (length(final_predictors) == 0) {
        stop(paste("No predictors remain after preprocessing for metric:", metric_name))
}

cat("Final predictors after preprocessing:", length(final_predictors), "\n")
cat(paste(final_predictors, collapse = ", "), "\n\n")

# parallel backend 
all_cores <- max(parallel::detectCores(logical = FALSE), 2L, na.rm = TRUE)
cl <- makePSOCKcluster(all_cores - 1)
on.exit(stopCluster(cl), add = TRUE)
registerDoParallel(cl)
cat("Parallel backend registered with", all_cores - 1, "cores.\n")

#===============================================================================
# 7. Define model for TUNING (Fast version) ----
#===============================================================================
# to reduce computational cost, we just set a high number for the 'trees' para-
# meter. Tuning this parameter would likely result in only marginal gains.
cat("Stting up lightweight tuning model ... \n")
rf_tune_spec <- 
        rand_forest(
                mode = "regression",
                trees = 1000,  
                mtry = tune(),
                min_n = tune()
        ) %>%
        set_engine(
                "ranger",
                # when running with do Parallel this prevents over-subscribing 
                # cores
                num.threads = 1
        )

#===============================================================================
# 8. Create tuning workflow ----
#===============================================================================

rf_tune_wf <- 
        workflow() %>%
        add_model(rf_tune_spec) %>%
        add_recipe(rf_recipe)

#===============================================================================
# 9. Create tuning grid ----
#===============================================================================

cat("Creating Latin Hypercube grid...\n")

rf_grid <- grid_space_filling(
        mtry(range = c(1L, length(final_predictors))),
        min_n(),
        size = 20

)

#===============================================================================
# 10. Tune hyperparameters ----
#===============================================================================
cat("\nTuning hyperparameters...\n")
set.seed(1)
j.tune <- 
        tune_race_anova(
                object = rf_tune_wf,
                resamples = j.folds,
                grid = rf_grid
)

#===============================================================================
# 11. Select best model ----
#===============================================================================
cat("Selecting best model...\n")
j.best <- select_best(j.tune, metric = "rmse")
cat("Best parameters:\n")
print(j.best)


#stopCluster(cl)
registerDoSEQ()
cat("Parallel backend stopped. Starting final fit with internal ranger 
    threads.\n")

#===============================================================================
# 12. Finalize and fit heavy model ----
#===============================================================================
rf_final_spec <- 
        rand_forest(
                mode = "regression",
                trees = 1000,
                mtry = j.best$mtry,
                min_n = j.best$min_n
        ) %>% 
        set_engine(
                "ranger",
                quantreg = TRUE,
                num.threads = all_cores,
                importance = "permutation"
        )

# Update workflow with the heavy model
j.wf_final <- workflow() %>%
        add_recipe(rf_recipe) %>%
        add_model(rf_final_spec)

# Fit on the full training set and evaluate on test set
j.final.fit <- last_fit(j.wf_final, j.split)

#===============================================================================
# 13. Extract model components ----
#===============================================================================

cat("\nExtracting model components...\n")
j.rf           <- extract_workflow(j.final.fit)
j.ranger_model <- extract_fit_engine(j.rf)

cat("Computing variable importance and direction...\n")
# 1. Get built-in permutation importance from ranger
var_imp <- j.ranger_model$variable.importance
# 2. Calculate direction (correlation with response)
# We use j.train.processed which was baked in Section 6
X <- j.train.processed[, final_predictors]
y <- j.train.processed$value
cors <- cor(X, y, use = "pairwise.complete.obs")
# 3. Combine into a single tibble
importance_df <- tibble(
        variable      = names(var_imp),
        importance    = as.numeric(var_imp),
        direction_cor = as.numeric(cors[names(var_imp), 1])
) %>%
        arrange(desc(importance)) # Optional: sort by most important

print(head(importance_df))

#===============================================================================
# 14. Save outputs ----
#===============================================================================

cat("\nSaving outputs...\n")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
}

# Clean metric name for filename
metric_clean <- gsub(" ", "_", metric_name)

# Bundle and save model
j.bundle   <- bundle(j.rf)
model_file <- file.path(output_dir, paste0(metric_clean, ".rds"))
saveRDS(j.bundle, model_file)
cat("Model saved to:", model_file, "\n")

# Save model performance metrics
metrics <- collect_metrics(j.final.fit)
metrics_file <- file.path(output_dir, paste0(metric_clean, "_metrics.rds"))
saveRDS(metrics, metrics_file)
cat("Metrics saved to:", metrics_file, "\n")

imp_file <- file.path(output_dir, paste0(metric_clean, "_importance.rds"))
saveRDS(importance_df, imp_file)
cat("Variable importance saved to:", imp_file, "\n")


cat("\n========================================\n")
cat("Processing complete for metric:", metric_name, "\n")
cat("========================================\n\n")
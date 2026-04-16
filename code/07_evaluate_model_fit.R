1################################################################################
# Script Name:        check_model_fit_large_models.R
# Description:        Checks MCMC convergence (PSRF) and predictive performance 
#                     (AUC, Beta diversity, C-score, Prevalence) for HMSC models.
#
# Author:             Jonathan Jupke
################################################################################

cat("R script starting \n")

# 1. Setup & Arguments =========================================================
args <- commandArgs(trailingOnly = TRUE)

# Initialize defaults
iter_id <- NULL
input_file <- NULL
output_file <- NULL

# Parse arguments
for (i in seq_along(args)) {
        if (args[i] == "--iter_id") {
                iter_id <- args[i + 1]
        } else if (args[i] == "--input") {
                input_file <- args[i + 1]
        } else if (args[i] == "--output") {
                output_file <- args[i + 1]
        }
}

if (is.null(input_file) || is.null(output_file)) {
        stop("ERROR: Missing --input or --output argument.")
}

cat("\n========================================\n")
cat("Processing Iteration:", iter_id, "\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n")
cat("========================================\n\n")

# 1.1 Load Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
        library(Hmsc)
        library(data.table)
        library(jsonify)
        library(vegan)
})

# 1.2 Load Custom Functions ----------------------------------------------------
# Use file.path for safer path construction
base_func_path <- "../pulse/code/functions"
cat("Loading custom functions from:", base_func_path, "\n")

source(file.path(base_func_path, "gelman_check.R"))
source(file.path(base_func_path, "c_score.R"))
source(file.path(base_func_path, "spec_env_test.r"))
source(file.path(base_func_path, "site_occupancy.R"))

# =======================================================
# 1.3 Load Data 
# =======================================================
cat("Loading input data...\n")

# Load Unfitted Model
unfittedModel <- readRDS(input_file)
cat("Unfitted Model successfully read.\n")

# Extract Model Name
# Remove path prefix
model.name <- basename(input_file) 
# Remove extension
model.name <- tools::file_path_sans_ext(model.name)
cat("Target Model Name:", model.name, "\n")

# Locate Fitted Models
fitted_files_dir <- "data/003_fitted_hmsc_models/"
all_fitted_files <- list.files(fitted_files_dir, full.names = TRUE)

# Find specific fitted model files (Using fixed=TRUE ensures Model_1 doesn't match Model_10)
target_fitted_indices <- grep(model.name, all_fitted_files, fixed = TRUE)

if (length(target_fitted_indices) == 0) stop("No fitted model files found for: ", model.name)
target_fitted_files <- all_fitted_files[target_fitted_indices]

cat("Found fitted files:\n")
print(target_fitted_files)

# Read fitted models (Assuming list structure where [[1]] is the JSON)
# We read them into a temp list first to avoid variable confusion
loaded_json_list <- lapply(target_fitted_files, function(x) readRDS(x)[[1]])

# Locate Spatial Data
spatial_dir <- "data/misc/spatial_scale/"
all_spatial_files <- list.files(spatial_dir, full.names = TRUE)
target_spatial_index <- grep(model.name, all_spatial_files, fixed = TRUE)

if (length(target_spatial_index) == 0) stop("No spatial file found for: ", model.name)
spatial_data <- readRDS(all_spatial_files[target_spatial_index])
cat("Spatial data successfully read.\n")

# =======================================================
# 1.4 Unpack JSON 
# =======================================================

cat("Unpacking JSON models...\n")

fit_model_objects <- lapply(loaded_json_list, from_json, buffer_size = 524288000)

# Assuming we strictly want the first two chains (or files)
postList <- list(fit_model_objects[[1]][[1]], fit_model_objects[[2]][[1]])

cat("Models unpacked. Importing posterior...\n")

# 2. Check Convergence (PSRF) ==================================================
mcmc.samples <- 2150
fitTF <- importPosteriorFromHPC(
        m = unfittedModel,
        postList = postList,
        nSamples = mcmc.samples,
        thin = 150,
        transient = 10000
)
cat("Posterior imported. Converting to Coda...\n")
# Converting from Hmsc to CODA object class. CODA is used for PSRF computation.
postCoda <- convertToCodaObject(fitTF)
# Using custom function to compute Potential scale reduction factor. 
psrfs    <- gelman_check(postCoda)

# Calculate Exceedance
# Threshold: 1.1
exceedence.rate <- sum(psrfs$taxon$psrf >= 1.1) / nrow(psrfs$taxon)

# Initialize Record Table
iterRecord <- data.table(
        scheme = spatial_data$scheme_id,
        psrf_exceedence = exceedence.rate,
        psrf_passed = exceedence.rate < 0.1, # TRUE if Good
        AUC_failure = 0,
        AUC_passed = NA,
        beta_ses_score = NA_real_,
        betaDistr_passed = NA,
        c_score_z = NA_real_,
        c_score_passed = NA,
        prevalence_flagged = NA_integer_,
        prevalence_passed = NA
)

# =======================================================
#  3. Model Fit Checks
# =======================================================

if (iterRecord$psrf_passed) {
        cat(model.name, ": PSRF Good (Exceedance:", round(exceedence.rate, 3), "). Checking AUC ...\n")
        
        # --- 3.1 Check model fit  ----
        i.preds <- computePredictedValues(fitTF)
        i.MF    <- evaluateModelFit(hM = fitTF, predY = i.preds)
        
        AUC_failure <- sum(i.MF$AUC < 0.75, na.rm = TRUE)
        AUC_failure_fraction <- AUC_failure/sum(!is.na(i.MF$AUC))
        iterRecord[, AUC_failure := AUC_failure_fraction]
        # check that less than 25% of models have a AUC below 0.75
        iterRecord[, AUC_passed := AUC_failure < 0.25]
                
        if (!iterRecord$AUC_passed) {
                cat(model.name, ": Bad AUC:", round(AUC_failure_fraction, 3), ")\n")
        } else {
                cat(model.name, ": Good AUC. Proceeding to PPC...\n")
                
                # --- 3.2 Posterior Predictive Checks ----
                
                # - 3.2.1 Predict Communities - 
                # We want 100 samples total for prediction to save resources
                n_pred_samples <- 100
                chain.samples <- sample(seq_len(mcmc.samples), size = n_pred_samples, replace = FALSE)
                
                pooled.posterior.samples <- poolMcmcChains(fitTF$postList)
                
                predCom <- predict(
                        object = fitTF,
                        post = pooled.posterior.samples[chain.samples],
                        XData = fitTF$XData,
                        expected = FALSE
                )
        
                # -------------------------------------------------------------------------
                # METRIC 1: Beta Diversity (Community Dissimilarity)
                # -------------------------------------------------------------------------
                
                Y_obs_clean <- fitTF$Y
                # Sorenson for binary
                obs_dist <- as.vector(
                        vegdist(
                                Y_obs_clean, 
                                method = "bray", 
                                na.rm = TRUE
                                )
                        ) 
                obs_quantiles <- 
                        quantile(
                                obs_dist, 
                                probs = c(0.25, 0.50, 0.75), 
                                na.rm = TRUE
                                )
                sim_quantiles <- matrix(NA, 
                                        nrow = length(chain.samples), 
                                        ncol = 3
                                        )
                for (i in seq_along(chain.samples)) {
                        # Prepare Simulated Data
                        Y_sim_clean        <- predCom[[i]]
                        sim_dist           <- as.vector(vegdist(Y_sim_clean, method = "bray"))
                        sim_quantiles[i, ] <- quantile(sim_dist, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
                }
                # Create a matrix with 2 rows (2.5% and 97.5% percentiles) and 
                # three columns (0.25, 0.50, 0.75 distance quantiles)

                # Calculate mean and standard deviation of the simulations
                sim_means <- colMeans(sim_quantiles, na.rm = TRUE)
                sim_sds   <- apply(sim_quantiles, 2, sd, na.rm = TRUE)
                
                # Calculate Standardized Effect Size (SES) for each quantile
                # (Observed - Mean) / SD
                ses_values <- (obs_quantiles - sim_means) / sim_sds
                
                # We take the mean of the absolute SES values to get an "Average Deviation"
                beta_deviation_score <- mean(abs(ses_values), na.rm = TRUE)
                iterRecord[, beta_ses_score := beta_deviation_score]
                
                # Compare Average SES to Threshold
                # < 2.0 : KEEP (The deviation is within 2 standard deviations, acceptable for JSDMs)
                # >= 2.0 : FLAG/EXCLUDE (The model structure is significantly biased)
                iterRecord[,  betaDistr_passed := (beta_deviation_score < 2)]
                
                # -------------------------------------------------------------------------
                # METRIC 2: Species Co-occurrence (C-score)
                # -------------------------------------------------------------------------
                c_score_res <- ppc_cooccurrence(obs_comm = fitTF$Y, sim_comm_list = predCom)
                iterRecord[, c_score_z := c_score_res$z_score]
                iterRecord[, c_score_passed := !c_score_res$flag]
                
                # -------------------------------------------------------------------------
                # METRIC 3: Species Prevalence Distribution
                # -------------------------------------------------------------------------
                prev_res <- ppc_prevalence(obs_com = fitTF$Y, sim_comm_list = predCom)
                iterRecord[, prevalence_flagged := sum(prev_res$summary_stats$flag)]
                iterRecord[, prevalence_passed := !prev_res$overall_flag]

                # Closes if (!iterRecord$AUC_passed) if clause 
                } 
        # else from if (iterRecord$psrf_passed) if clause 
        } else { 
                
        cat(model.name, ": PSRF Bad. Skipping predictions.\n")
}

# 4. Save Output ===============================================================
cat("Saving results to:", output_file, "\n")
saveRDS(iterRecord, output_file)
cat("Done.\n")
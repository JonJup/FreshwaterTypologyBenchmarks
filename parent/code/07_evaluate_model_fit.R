################################################################################*
# Description:  Checks MCMC convergence (PSRF) and predictive performance
#               (AUC, beta-diversity SES, C-score, species prevalence) for a
#               single HMSC model fitted on HPC, and writes a one-row QC record.
#               
# Notes:        Intended to be called via Slurm array jobs:
#               Run via Rscript with --iter_id / --input / --output arguments.
#               Expects the unfitted model (--input), the matching fitted JSON
#               models in data/003_fitted_hmsc_models/, and the matching spatial
#               object in data/misc/spatial_scale/. 
#               Assumes a 2-chain fit.
################################################################################*

cat("R script starting \n")

# 1. Setup & Arguments =========================================================
# CLI args passed by the HPC launcher
args <- commandArgs(trailingOnly = TRUE)          

# Initialize defaults (so missing args are detectable below)
iter_id <- NULL
input_file <- NULL
output_file <- NULL
output_file2 <- NULL

# Parse arguments: each flag is followed by its value (args[i + 1])
for (i in seq_along(args)) {
        if (args[i] == "--iter_id") {
                iter_id <- args[i + 1]
        } else if (args[i] == "--input") {
                input_file <- args[i + 1]
        } else if (args[i] == "--output") {
                output_file <- args[i + 1]
        } else if (args[i] == "--output2") {
                output_file2 <- args[i + 1]
        }
}
output_file2 <- paste0(output_file2)
# Hard-fail early if the two mandatory paths are absent
if (is.null(input_file) || is.null(output_file)) {
        stop("ERROR: Missing --input or --output argument.")
}

cat("\n========================================\n")
cat("Processing Iteration:", iter_id, "\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n")
cat("Output file2 :", output_file2, "\n")
cat("========================================\n\n")

# ================================================*
## 1.1 Load Libraries ----------------------------
# ================================================*

suppressPackageStartupMessages({
        library(Hmsc)        # JSDM fitting / posterior import / model-fit metrics
        library(data.table)  # iterRecord container + := updates
        library(jsonify)     # from_json() to unpack the HPC-fitted chains
        library(vegan)       # vegdist() for community dissimilarities
})

# ================================================*
## 1.2 Load Custom Functions -----------------------
# ================================================*

# Use file.path for safer (OS-independent) path construction
base_func_path <- "../R"
cat("Loading custom functions from:", base_func_path, "\n")

source(file.path(base_func_path, "gelman_check.R"))    # gelman_check(): PSRF table
source(file.path(base_func_path, "c_score.R"))         # ppc_cooccurrence(): C-score PPC
source(file.path(base_func_path, "spec_env_test.r"))   # species-environment helpers
source(file.path(base_func_path, "site_occupancy.R"))  # ppc_prevalence(): prevalence PPC

# =======================================================*
## 1.3 Load Data ----
# =======================================================*
cat("Loading input data...\n")

# Load Unfitted Model (the Hmsc skeleton that the posterior is imported into)
unfittedModel <- readRDS(input_file)
cat("Unfitted Model successfully read.\n")

# Extract Model Name from the input filename
model.name <- basename(input_file)                       # drop the directory prefix
model.name <- tools::file_path_sans_ext(model.name)      # drop the file extension
cat("Target Model Name:", model.name, "\n")

# Locate Fitted Models on disk
fitted_files_dir <- "data/003_fitted_hmsc_models/"
all_fitted_files <- list.files(fitted_files_dir, full.names = TRUE)

# Find fitted files belonging to this model.
target_fitted_indices <- grep(model.name, all_fitted_files, fixed = TRUE)

if (length(target_fitted_indices) == 0) stop("No fitted model files found for: ", model.name)
target_fitted_files <- all_fitted_files[target_fitted_indices]

cat("Found fitted files:\n")
print(target_fitted_files)

# Read fitted models into a temporary list ([[1]] holds the JSON string per file)
loaded_json_list <- lapply(target_fitted_files, function(x) readRDS(x)[[1]])

# Locate Spatial Data
spatial_dir <- "data/misc/spatial_scale/"
all_spatial_files <- list.files(spatial_dir, full.names = TRUE)
target_spatial_index <- grep(model.name, all_spatial_files, fixed = TRUE)

if (length(target_spatial_index) == 0) stop("No spatial file found for: ", model.name)
spatial_data <- readRDS(all_spatial_files[target_spatial_index])
cat("Spatial data successfully read.\n")

# =======================================================*
## 1.4 Unpack JSON ----
# =======================================================*

cat("Unpacking JSON models...\n")

# Deserialize each chain's JSON (large buffer needed for big HMSC posteriors)
fit_model_objects <- lapply(loaded_json_list, from_json, buffer_size = 524288000)

# Build the postList Hmsc expects. This hardcodes exactly TWO chains: element
# [[1]] of each fitted object is the chain's posterior.
postList <- list(fit_model_objects[[1]][[1]], fit_model_objects[[2]][[1]])

cat("Models unpacked. Importing posterior...\n")

# 2. Check Convergence (PSRF) --------------------------------------------------

mcmc.samples <- 2150          # recorded posterior samples PER chain
fitTF <- importPosteriorFromHPC(
        m = unfittedModel,
        postList = postList,
        nSamples = mcmc.samples,
        thin = 150,
        transient = 10000
)
cat("Posterior imported. Converting to Coda...\n")
# Convert from Hmsc to CODA object class; CODA is used for PSRF computation.
postCoda <- convertToCodaObject(fitTF)
# Custom function computing the Potential Scale Reduction Factor per taxon.
psrfs    <- gelman_check(postCoda)

# Exceedance = fraction of taxa whose PSRF is at/above the 1.1 convergence threshold
exceedence.rate <- sum(psrfs$taxon$psrf >= 1.1) / nrow(psrfs$taxon)

# Initialize the one-row QC record. Metrics default to NA and are filled only if
# the relevant check is reached (PSRF gate -> AUC gate -> PPC metrics).
iterRecord <- data.table(
        scheme = spatial_data$scheme_id,
        psrf_exceedence = exceedence.rate,
        psrf_passed = exceedence.rate < 0.1,   # TRUE = converged (<10% of taxa exceed 1.1)
        AUC_failure = 0,
        AUC_passed = NA,
        beta_ses_score = NA_real_,
        betaDistr_passed = NA,
        c_score_z = NA_real_,
        c_score_passed = NA,
        prevalence_flagged = NA_integer_,
        prevalence_passed = NA
)


#  3. Check Model Fit -----------------------------------------------------------


if (iterRecord$psrf_passed) {
        cat(model.name, ": PSRF Good (Exceedance:", round(exceedence.rate, 3), "). Checking AUC ...\n")
        
        # =======================================================*
        # --- 3.1 Check AUC  ----
        # =======================================================*
        i.preds <- computePredictedValues(fitTF)                 # posterior predictive Y
        i.MF    <- evaluateModelFit(hM = fitTF, predY = i.preds) # per-species fit metrics
        
        # Save for later inspections 
        out <- data.table(
                scheme_id = spatial_data$scheme_id,
                group = gsub("_.*", "", spatial_data$scheme_id),
                taxon = colnames(fitTF$Y),
                AUC = i.MF$AUC,
                TR2 = i.MF$TjurR2,
                RMSE = i.MF$RMSE
        )
        saveRDS(out, output_file2)
        
        
        AUC_failure <- sum(i.MF$AUC < 0.75, na.rm = TRUE)            # count of poorly-fit species
        AUC_failure_fraction <- AUC_failure / sum(!is.na(i.MF$AUC))  # as a fraction of scored species
        iterRecord[, AUC_failure := AUC_failure_fraction]
        # PASS if fewer than 25% of species fall below AUC 0.75.
        # (Compare the FRACTION explicitly; the bare `AUC_failure` is a count.)
        iterRecord[, AUC_passed := AUC_failure_fraction < 0.25]

        if (!iterRecord$AUC_passed) {
                cat(model.name, ": Bad AUC:", round(AUC_failure_fraction, 3), "\n")
        } else {
                cat(model.name, ": Good AUC. Proceeding to PPC...\n")
                
                # =======================================================*
                ## 3.2 Posterior Predictive Checks ----
                # =======================================================*
        
                ###  3.2.1 Predict Communities ----
                # Use 100 posterior draws for prediction to save resources.
                n_pred_samples <- 100

                # Pool both chains FIRST, then draw sample indices from the full
                # pooled posterior. (Sampling must span all pooled draws, not just
                # the first chain's worth of indices.)
                pooled.posterior.samples <- poolMcmcChains(fitTF$postList)
                chain.samples <- sample(seq_along(pooled.posterior.samples),
                                        size = n_pred_samples, replace = FALSE)

                predCom <- predict(
                        object = fitTF,
                        post = pooled.posterior.samples[chain.samples],
                        XData = fitTF$XData,
                        expected = FALSE          # draw 0/1 occurrences, not probabilities
                )

                
                ### 3.2.2 Beta Diversity (Community Dissimilarity) ----
                

                Y_obs_clean <- fitTF$Y
                # Bray-Curtis dissimilarity. With presence/absence (0/1) data this is
                # equivalent to the Sorensen index.
                obs_dist <- as.vector(
                        vegdist(
                                Y_obs_clean,
                                method = "bray",
                                na.rm = TRUE
                                )
                        )
                # Reference quartiles of the observed dissimilarity distribution
                obs_quantiles <-
                        quantile(
                                obs_dist,
                                probs = c(0.25, 0.50, 0.75),
                                na.rm = TRUE
                                )
                # One row per simulated community, columns = the three quartiles
                sim_quantiles <- matrix(NA,
                                        nrow = length(chain.samples),
                                        ncol = 3
                                        )
                for (i in seq_along(chain.samples)) {
                        # Dissimilarity quartiles of each simulated community
                        Y_sim_clean        <- predCom[[i]]
                        sim_dist           <- as.vector(vegdist(Y_sim_clean, method = "bray"))
                        sim_quantiles[i, ] <- quantile(sim_dist, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
                }
                # sim_quantiles: length(chain.samples) rows x 3 columns
                # (the 0.25, 0.50 and 0.75 dissimilarity quantiles per simulation).

                # Column-wise mean and SD of the simulated quantiles (the null reference)
                sim_means <- colMeans(sim_quantiles, na.rm = TRUE)
                sim_sds   <- apply(sim_quantiles, 2, sd, na.rm = TRUE)

                # Standardized Effect Size (SES) per quantile: (Observed - Mean) / SD
                ses_values <- (obs_quantiles - sim_means) / sim_sds

                # Average absolute SES = mean deviation of observed from simulated structure
                beta_deviation_score <- mean(abs(ses_values), na.rm = TRUE)
                iterRecord[, beta_ses_score := beta_deviation_score]

                # Threshold on the average SES:
                #   < 2.0  : KEEP   (within ~2 SD; acceptable for JSDMs)
                #   >= 2.0 : FLAG   (model structure significantly biased)
                iterRecord[,  betaDistr_passed := (beta_deviation_score < 2)]

               
                ### 3.3.3 Species Co-occurrence (C-score) ----
               
                c_score_res <- ppc_cooccurrence(obs_comm = fitTF$Y, sim_comm_list = predCom)
                iterRecord[, c_score_z := c_score_res$z_score]
                iterRecord[, c_score_passed := !c_score_res$flag]   # passed = not flagged

                
                ### 3.3.4 Species Prevalence Distribution ----
                
                prev_res <- ppc_prevalence(obs_com = fitTF$Y, sim_comm_list = predCom)
                iterRecord[, prevalence_flagged := sum(prev_res$summary_stats$flag)]  # # flagged species
                iterRecord[, prevalence_passed := !prev_res$overall_flag]

        }   # closes the AUC else-branch (the PPC block)

} else {
        # PSRF gate failed: skip all predictive checks, leave metrics as NA.
        cat(model.name, ": PSRF Bad. Skipping predictions.\n")
}

# 4. Save Output ---------------------------------------------------------------
cat("Saving results to:", output_file, "\n")
saveRDS(iterRecord, output_file)
cat("Done.\n")
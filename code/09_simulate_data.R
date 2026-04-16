################################################################################
# Script Name:        diagnostic_simulation_pipeline.R
#
# Purpose:            This script generates synthetic community-composition data
#                     that mimic realistic freshwater biomonitoring samples. It
#                     does so by:
#
#                       1. Loading a fitted HMSC model.
#                       2. Creating artificial environmental gradients by
#                          clustering observed environmental conditions (k-means),
#                          then systematically stretching / compressing those
#                          clusters along their discriminant axes (LDA-based
#                          manipulation) to produce environments with known,
#                          controlled typological quality.
#                       3. Validating those synthetic environments against the
#                          real-world catchment envelope (Isolation Forest +
#                          Mahalanobis distance) to ensure they remain plausible.
#                       4. Feeding validated synthetic environments into the
#                          fitted HMSC model to predict community composition.
#                       5. Performing fuzzy classification (FCM) on the validated
#                          environments and computing Normalised Partitioning
#                          Entropy (NPE) as a measure of classification crispness.
#
#                     These simulated datasets have are used for benchmarking
#                     typology systems. The "quality_factor" parameter controls
#                     how distinct the types are — low values compress clusters
#                     (overlapping types), high values separate them
#                     (crisp types).
#
#                     In addition to production outputs (predictions, cluster
#                     assignments, fuzzy memberships), this version writes a
#                     diagnostic table recording every parameter combination
#                     tried, including those that failed filters, so we can
#                     understand where and why simulations are rejected.
#
# Author:             Jonathan Jupke
# Date Created:       2026-03-16
#
# R Version:          R 4.5.2
# Required Packages:  data.table, kernlab, cluster, vegclust, Hmsc, isotree,
#                     jsonify, sf, arrow, dplyr, MASS 
#
# Execution:          Called from SLURM array jobs. Each array task processes
#                     one HMSC model (one biological group × typology scheme).
#
# Output files:
#   1. <output_file>                — Production list: predictions, clusters,
#                                     fuzzy memberships, ASW, NPE, metadata.
#   2. data/misc/simulation_diagnostics/<output_file> — Diagnostic data.table:
#                                     one row per parameter combo per iteration,
#                                     with filter pass/fail status and metrics.
################################################################################


# ==============================================================================
# 1. SETUP AND CONFIGURATION
# ==============================================================================
# Parse command-line arguments, load libraries, source helper functions, and
# set global simulation parameters.
# ==============================================================================

## 1.1 Arguments Parsing -------------------------------------------------------
# The script expects three named arguments from the command line:
#   --iter_id  : An identifier for this SLURM array task (informational only;
#                not used in computation, just printed for log traceability).
#   --input    : Path to the *unfitted* HMSC model .rds file. The model name
#                is extracted from this path and used to locate all associated
#                files (fitted model, evaluation, variation partitioning, etc.).
#   --output   : Path where the production output list will be saved.

args <- commandArgs(trailingOnly = TRUE)

iter_id     <- NULL
input_file  <- NULL
output_file <- NULL

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
cat("DIAGNOSTIC + PRODUCTION MODE\n")
cat("Processing Iteration:", iter_id, "\n")
cat("Input file:  ", input_file, "\n")
cat("Output file: ", output_file, "\n")
cat("========================================\n\n")

## 1.2 Load Libraries ----------------------------------------------------------

suppressPackageStartupMessages({
        library(data.table)   
        library(kernlab)      
        library(cluster)      
        library(vegclust)     
        library(Hmsc)         
        library(isotree)      
        library(jsonify)      
        library(sf)           
        library(arrow)        
        library(dplyr)        
})

## 1.3 Load Helper Scripts -----------------------------------------------------
# Two custom functions sourced from the PULSE project:
#
# normalized_partitioning_entropy.R:
#   Defines NPE(), which takes a vegclust object and returns Normalised
#   Partitioning Entropy — a [0,1] measure of how fuzzy vs. crisp the
#   classification is. NPE = 0 means perfectly crisp (each site belongs to
#   exactly one cluster); NPE = 1 means maximum fuzziness (uniform membership
#   across all clusters).
#
# adjust_clusters_lda.R:
#   Defines adjust_clusters_lda(), the core manipulation function. Given
#   cluster centroids, observations, assignments, and a quality_factor, it:
#     1. Fits an LDA to find the discriminant axes separating clusters.
#     2. Projects observations into discriminant space.
#     3. Scales the projections by quality_factor (>1 = more separated,
#        <1 = more overlapping).
#     4. Back-projects to the original variable space.
#   Returns the modified observations plus effective separation/dispersion/
#   quality metrics that reflect what the manipulation actually achieved
#   (which may differ from the requested quality_factor due to data geometry).

cat("Loading helper functions ...\n")
source("../pulse/code/functions/normalized_partitioning_entropy.R")
source("../pulse/code/functions/adjust_clusters_lda.R")
# Helper: map quality values to bin indices
get_bin <- function(q) {
        findInterval(q, quality_bins, rightmost.closed = TRUE)
}
# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
# N_TARGET_SUCCESSES: How many validated simulation environments we need.
#   Each "success" is one environment that passed both the Isolation Forest
#   and Mahalanobis plausibility checks, and was used to generate HMSC
#   predictions. The loop runs until we accumulate this many, or hit the
#   iteration cap.
#
# MAX_ITERATIONS: Safety cap to prevent infinite loops when most parameter
#   combinations fail validation (e.g., for models with unusual environmental
#   spaces). Each iteration tries one random variable subset × one random k,
#   then tests five quality_factor values within that setup.
# ---------------------------------------------------------------------------
N_TARGET_SUCCESSES <- 60
MAX_ITERATIONS     <- 3000

# Master seed: set once here so that the data-loading phase (which involves
# random operations like factor-level assignment) is reproducible. The
# simulation loop uses per-iteration seeds (see Section 4) to decouple
# loop randomness from any upstream code changes.
set.seed(1)


# ==============================================================================
# 2. DATA LOADING
# ==============================================================================
# Load the HMSC model (unfitted structure + fitted posterior chains), model
# evaluation metrics, typology scheme metadata, environmental data, and
# catchment-level data used for outlier detection baselines.
# ==============================================================================

## 2.1 Identify Model and Files ------------------------------------------------
# The model_name (e.g., "diatoms_0001") is extracted from the input path
# and used as a pattern to find all associated files across the project's
# directory structure.

cat("Identifying files...\n")
model_name <- gsub(pattern = "data/001_unfitted_hmsc_models/", "", input_file)
model_name <- gsub(pattern = "\\.rds", "", model_name)

path_eval           <- list.files("data/004_model_fit/", full.names = TRUE, pattern = model_name)
path_fitted_files   <- list.files("data/003_fitted_hmsc_models/", full.names = TRUE, pattern = model_name)
path_vp             <- list.files("data/005_variation_partitioning/", full.names = TRUE, pattern = model_name)
path_unscaled_env   <- list.files("data/misc/unscaled_environments/", full.names = TRUE, pattern = model_name)
path_unfitted_model <- list.files("data/001_unfitted_hmsc_models/", full.names = TRUE, pattern = model_name)
path_catchments     <- list.files("../pulse/data/catchments/", full.names = TRUE)

## 2.2 Quality Check -----------------------------------------------------------
# Before investing computation in simulation, verify the HMSC model is
# trustworthy. The evaluation file contains binary pass/fail flags for five
# diagnostics:
#   - psrf_passed:       Potential Scale Reduction Factor (Gelman-Rubin) ≈ 1
#                         → chains have converged.
#   - AUC_passed:        Area Under ROC Curve is acceptable → model has
#                         predictive discrimination.
#   - betaDistr_passed:  Beta (species response) distributions are well-behaved.
#   - c_score_passed:    C-score (checkerboard metric) → residual species
#                         co-occurrence patterns are captured.
#   - prevalence_passed: Model handles rare/common species appropriately.
#
# We require at least 4/5 to proceed. This prevents wasting HPC time on
# models whose predictions would be unreliable.

eval_data <- readRDS(path_eval)
eval_data[, good_model := psrf_passed + AUC_passed + betaDistr_passed + c_score_passed + prevalence_passed]
eval_data[is.na(good_model), good_model := 0]

if (eval_data$good_model < 4) {
        cat("Skipping: Model quality insufficient (score < 4)\n")
        quit(save = "no", status = 0)
}

cat("Model quality score:", eval_data$good_model, "/ 5 — proceeding.\n")

## 2.3 Load HMSC Models --------------------------------------------------------
# Two objects are needed:
#   1. unfitted_model: The HMSC model structure (formula, priors, random
#      effects specification, training data). This is the "skeleton" that
#      importPosteriorFromHPC() will populate with fitted parameters.
#   2. post_list: The MCMC posterior samples, which were fitted on HPC and
#      exported as JSON. We have two chains, each stored as a separate file.
#      The from_json() call deserializes them back into R list structures.
#      buffer_size is set high (500 MB) because these posteriors can be large
#      (hundreds of taxa × thousands of MCMC samples × multiple parameters).

cat("Reading unfitted and fitted models...\n")
unfitted_model <- readRDS(path_unfitted_model)
cat("Unfitted model loaded.\n")


loaded_json_list  <- lapply(path_fitted_files, function(x) readRDS(x)[[1]])
fit_model_objects <- lapply(loaded_json_list, from_json, buffer_size = 524288000)
post_list         <- list(fit_model_objects[[1]][[1]], fit_model_objects[[2]][[1]])
cat("Fitted model loaded.\n")

## 2.4 Load Environmental & Spatial Data ---------------------------------------
# schemes: A data.table mapping each model (scheme_id) to its associated
#   typology metadata — which Environmental Zones (EnZ) it covers, how many
#   samples it has, etc. This determines which catchment data to use for
#   the outlier detection baseline.
#
# env_val: Environmental validation data (loaded but not directly used in
#   this script; kept for compatibility with downstream code).
#
# eu_hydro_enz: A lookup table linking EU-Hydro river segment IDs to their
#   Environmental Zone (EnZ) classification. This is joined with catchment
#   data to assign each catchment to its EnZ.

cat("Reading auxiliary files ...\n")

schemes_file   <- list.files("data/000_biota", pattern = "03_", full.names = TRUE)
schemes        <- readRDS(schemes_file)
current_scheme <- schemes[scheme_id == model_name]
cat("Schemes loaded.\n")

env_val       <- readRDS("../pulse/data/environmental_validation.rds")
eu_hydro_enz  <- readRDS("../pulse/data/eu_hydro_dem_w_enz.rds")
cat("Validation data loaded.\n")

# Catchment processing:
# Each parquet file contains environmental conditions for river catchments in
# one region. We join each with the EnZ lookup, keep only rows with valid EnZ
# assignments, then combine everything and split by zone. This gives us
# zone-specific environmental distributions that serve as the reference for
# plausibility checking — simulated environments should look like they could
# have come from these real catchments.

cat("Processing catchments ... ")
catchment_list <- list()
for (i in seq_along(path_catchments)){
        cat_data    <- read_parquet(path_catchments[i])
        joined_data <- eu_hydro_enz[cat_data, on = "ID"]
        joined_data$EnZ_name <- factor(joined_data$EnZ_name)
        joined_data <- joined_data[!is.na(EnZ_name)]
        joined_data <- na.omit(joined_data)
        catchment_list[[i]] <- joined_data
};  rm(i)

all_catchments          <- rbindlist(catchment_list, use.names = TRUE)
catchment_split_by_zone <- split(all_catchments, by = "EnZ_name")
rm(path_catchments, eu_hydro_enz, catchment_list, all_catchments)
cat("done.\n")

## 2.5 Load Unscaled Environment -----------------------------------------------
# unscaled_env: The environmental conditions at the HMSC model's training
# sites, in their original (unscaled) measurement units. The HMSC model
# internally works with scaled (standardised) variables, but the outlier
# detection checks (ISO forest, Mahalanobis) operate on unscaled values
# because the catchment reference data is also unscaled.
#
# This distinction matters: we manipulate clusters in *both* spaces:
#   - Unscaled space: for validation (do the synthetic environments look
#     like real catchments?)
#   - Scaled space: for HMSC prediction (the model expects standardised inputs)

unscaled_env <- readRDS(path_unscaled_env)


# ==============================================================================
# 3. PRE-PROCESSING: OUTLIER DETECTION BASELINE
# ==============================================================================
# Before the simulation loop, establish reference distributions from the real
# catchment data. These are used in the loop to check whether each synthetic
# environment is "plausible" — i.e., whether it falls within the range of
# conditions that actually exist in the landscape.
#
# Two complementary methods:
#   1. Isolation Forest: A non-parametric outlier detector. It builds random
#      trees that partition the data; outliers are isolated in fewer splits.
#      We fit it on the catchment data, then score simulated environments.
#   2. Mahalanobis Distance: A parametric measure of how far a point is from
#      the multivariate mean, accounting for variable correlations. Points
#      beyond a chi-squared threshold are considered outliers.
# ==============================================================================

## 3.1 Isolation Forest on Catchment Data --------------------------------------
# Subset catchments to only the Environmental Zones relevant to this model,
# then train the isolation forest on those catchments.
#
# Key parameters:
#   ntrees = 1000:                Enough trees for stable anomaly scores.
#   sample_size = 256:            Standard subsample size per tree (as per
#                                 the original Isolation Forest paper).
#   ndim = ncol(catchment_matrix): Use all variables for splitting.
#   prob_pick_pooled_gain = 0:    Disable gain-based split selection; use
#                                 purely random splits (classic iForest).

target_zones      <- unlist(current_scheme$catchments)
target_catchments <- rbindlist(catchment_split_by_zone[target_zones])
target_catchments[, c("ID", "EnZ_name") := NULL]
catchment_matrix  <- as.matrix(as.data.frame(target_catchments))
catchment_matrix  <- catchment_matrix[, c(names(unscaled_env))]

isolation_forest <- isolation.forest(
        catchment_matrix,
        ntrees            = 1000,
        sample_size       = 256,
        ndim              = ncol(catchment_matrix),
        prob_pick_pooled_gain = 0
)

# Score the original (real) training sites to establish a baseline distribution.
# Simulated environments will be compared against this baseline using Z-scores.
iso_scores <- predict(isolation_forest, unscaled_env)
score_mean <- mean(iso_scores)
score_sd   <- sd(iso_scores)

# Two thresholds for isolation forests:
# iso_threshold (0.10): The fraction of sites allowed to fail — used to compute
#   the minimum number of sites that must pass the Z-score filter for an entire
#   simulated environment to be accepted. With 0.10, we require ≥ 90% of sites
#   to pass.
# iso_threshold2 (1): The per-site Z-score ceiling. Sites with anomaly scores
#   more than 1 standard deviations above the mean are flagged as implausible
#   and removed before the count-based filter is applied.
iso_threshold  <- 0.15
iso_threshold2 <- 1.5

cat(sprintf("Isolation Forest baseline: mean = %.4f, sd = %.4f\n", score_mean, score_sd))

## 3.2 Mahalanobis Distance Baseline -------------------------------------------
# Compute the covariance matrix and mean vector of the catchment data. These
# define the multivariate normal envelope against which simulated environments
# are checked.
#
# The threshold is the 95th percentile of the chi-squared distribution with
# degrees of freedom = number of environmental variables. Points beyond this
# threshold are considered multivariate outliers.

env_cov_mat    <- cov(catchment_matrix)
env_mean       <- colMeans(catchment_matrix)
maha_threshold <- qchisq(0.95, df = ncol(catchment_matrix))

cat(sprintf("Mahalanobis threshold (chi-sq 0.95, df=%d): %.2f\n",
            ncol(catchment_matrix), maha_threshold))


# ==============================================================================
# 4. DIAGNOSTIC SIMULATION LOOP
# ==============================================================================
# This is the heart of the script. Each iteration of the outer loop:
#
#   1. Randomly selects a subset of environmental variables.
#   2. Randomly selects k (number of clusters) and runs k-means.
#   3. Checks whether the resulting clusters have clear discriminant
#      structure (LDA variance check).
#   4. For each of several quality_factor values (0.5 to 1.5):
#      a. Manipulates the cluster structure via LDA-based adjustment.
#      b. Validates the resulting environment against the catchment envelope.
#      c. If valid: builds the HMSC prediction environment, computes ASW,
#         runs fuzzy classification, predicts community composition.
#
# The loop continues until N_TARGET_SUCCESSES environments pass all checks,
# or MAX_ITERATIONS is reached.
#
# DIAGNOSTIC PHILOSOPHY:
# Every parameter combination is recorded — pass or fail — with the specific
# failure reason. This lets us understand the simulation's "acceptance rate"
# and diagnose whether filters are too strict/lenient for particular models.
# ==============================================================================

## 4.0 Prepare HMSC Model and Posterior ----------------------------------------
# Reconstruct the full fitted HMSC model by merging the unfitted structure
# with the HPC-exported posterior chains.
#
# importPosteriorFromHPC() parameters:
#   nSamples  = 2150:  Total posterior samples per chain (after thinning).
#   thin      = 150:   Keep every 150th MCMC iteration.
#   transient = 10000: Discard the first 10,000 iterations as burn-in.
#
# The resulting hmsc_model object contains everything needed for prediction:
#   - $XData:    The original training environmental data (a data.frame).
#   - $postList: List of posterior samples (beta coefficients, etc.).
#   - $samples:  Number of posterior samples available.

cat("Importing posterior from HPC object...\n")

vp_data <- readRDS(path_vp)
# vp_data contains variation partitioning results — specifically,
# vp_data$importance is a named vector giving each environmental variable's
# contribution to explained variance. Used later (Section 4.10) to record
# how much of the model's explanatory power is captured by each simulation's
# variable subset.

hmsc_model <- importPosteriorFromHPC(
        m         = unfitted_model,
        postList  = post_list,
        nSamples  = 2150,
        thin      = 150,
        transient = 10000
)

cat("Posterior imported.\n")

# Extract the environmental data matrix from the model.
# hmsc_XData is a data.frame with columns for:
#   - Environmental variables (e.g., temperature, slope, geology fractions)
#   - MEM variables (Moran's Eigenvector Maps — spatial eigenvectors that
#     capture spatial autocorrelation in the random effects)
#   - A "." column (intercept placeholder in some HMSC versions)
#
# For clustering, we use only the non-spatial, non-intercept variables.
hmsc_XData <- hmsc_model$XData

all_vars <- colnames(hmsc_XData)[!grepl("MEM", colnames(hmsc_XData))]
if (any(colnames(hmsc_XData) == ".")) {
        all_vars <- all_vars[-which(colnames(hmsc_XData) == ".")]
}
all_vars <- sort(all_vars)

# var_selection_counter: Tracks how many times each variable has been part
# of a *successful* simulation. Variables that have already contributed to
# many successes get down-weighted in future selections (via 1/counter
# probability weights), encouraging diversity across simulations. All
# variables start at 1 (not 0, to avoid division by zero).
var_selection_counter <- rep(1, length(all_vars))
names(var_selection_counter) <- all_vars

# geo_vars: The three geology composition variables (calcareous, sediment,
# siliceous) must always appear together because they are compositional —
# they sum to ~1. Including only one or two would create artificial
# constraints in the simulated environments.
geo_vars <- c("area_calcareous", "area_sediment", "area_siliceous")


# Bin edges and quota tracking (initialize before loop)
quality_bins   <- seq(0.5, 1.5, by = 0.10)  # 0.50, 0.75, 1.00, 1.25, 1.50
n_bins         <- length(quality_bins) - 1    # 4 bins
per_bin_quota  <- ceiling(N_TARGET_SUCCESSES / n_bins)
bin_success_counts <- rep(0, n_bins)
names(bin_success_counts) <- paste0("[", quality_bins[-length(quality_bins)],
                                    ",", quality_bins[-1], ")")


## 4.1 Diagnostic Storage ------------------------------------------------------
# diag_records: Collects one data.table row per parameter combination per
# iteration. Failed combinations get NA for metrics computed after their
# failure point, plus a final_status string indicating where they failed.
# This list is rbindlist()'d at the end into a single diagnostic table.
diag_records <- list()
record_id    <- 0L

## 4.1b Production Storage -----------------------------------------------------
# These vectors/lists accumulate results from successful simulations only.
# They are indexed in parallel — res_sim_output[[i]], res_n_clusters[i],
# res_asw[i], etc. all refer to the same simulation.

res_sim_output    <- list()  # HMSC predictions (list of posterior-draw lists)
res_cluster_ids   <- list()  # Hard cluster assignments (integer vectors)
res_fuzzy_ids     <- list()  # vegclust FCM objects (contain membership matrices)
res_npe           <- c()     # Normalised Partitioning Entropy values
res_n_clusters    <- c()     # Number of clusters (k) for each simulation
res_n_variables   <- c()     # Number of env variables selected
res_n_valid  <- c()          # Number of valid environments comming from this iteration
res_importance    <- c()     # Summed variable importance from variation partitioning
res_asw           <- c()     # Average Silhouette Width (full-space, post-manipulation)
res_contr_points  <- c()     # Effective dispersion (within-cluster spread metric)
res_contr_centro  <- c()     # Effective separation (between-centroid distance metric)


success_count <- 0

cat("\nStarting diagnostic loop...\n")
cat(sprintf("Target: %d successes, max %d iterations\n", N_TARGET_SUCCESSES, MAX_ITERATIONS))
cat("---\n")

# =========================================================================
# BEGIN MAIN LOOP
# =========================================================================

for (q in 1:MAX_ITERATIONS) {
        
        # Per-iteration seed: ensures that each iteration's random draws
        # (variable selection, k, chain sampling) are deterministic and
        # independent of how many iterations ran before. This means adding
        # code before the loop (e.g., new data loading steps) won't cascade
        # into changed simulation results.
        set.seed(1000 + q)
        
        if (success_count >= N_TARGET_SUCCESSES) break
        
        #cat(sprintf("\rIteration q=%d : Successes %d/%d", q, success_count, N_TARGET_SUCCESSES))
        cat(sprintf("\rIteration q=%d : Successes %d/%d  Bins: %s",
                    q, success_count, N_TARGET_SUCCESSES,
                    paste(bin_success_counts, collapse = "/")))
        flush.console()
        # ------------------------------------------------------------------
        # 4.2 Select Random Variable Subset
        # ------------------------------------------------------------------
        # Each simulation uses a random subset of environmental variables for
        # clustering. This creates diversity in the typological "dimension" —
        # some simulations define types by temperature + geology, others by
        # slope + land use, etc.
        #
        # prob_weights = 1/counter: Variables that have appeared in many
        # successful simulations are less likely to be selected, pushing the
        # simulation toward underrepresented variable combinations.
        #
        # n_vars_selected is drawn uniformly from [3, n_all_vars - 1].
        # Minimum of 3 ensures LDA has enough dimensions for meaningful
        # discriminant axes; maximum of n-1 avoids using all variables
        # (which would leave no "missing" variables for the full-space ASW
        # comparison to be interesting).
        
        prob_weights    <- 1 / var_selection_counter
        n_vars_selected <- sample(x = 3:(length(all_vars) - 1), size = 1)
        selected_vars   <- sample(all_vars, size = n_vars_selected, replace = FALSE, prob = prob_weights)
        
        # Force geology co-occurrence: if any of the three geology variables
        # was selected, add the others. This may increase n_vars_selected
        # beyond the originally drawn value.
        if (any(geo_vars %in% selected_vars) && !all(geo_vars %in% selected_vars)) {
                vars_to_add   <- setdiff(geo_vars, selected_vars)
                vars_to_add   <- intersect(vars_to_add, all_vars)
                selected_vars <- unique(c(selected_vars, vars_to_add))
                selected_vars <- intersect(selected_vars, all_vars)
                n_vars_selected <- length(selected_vars)
        }
        
        # ------------------------------------------------------------------
        # 4.3 Cluster Environments (K-Means)
        # ------------------------------------------------------------------
        # Cluster the HMSC training sites in the selected variable subspace.
        # This creates the "ground truth" typology that we will then
        # manipulate to vary type distinctness.
        #
        # k is drawn uniformly from [2, 10]. nstart = 100 gives k-means many
        # random starts to find a good solution. iter.max = 100 allows
        # convergence for difficult configurations.
        #
        # Minimum cluster size check: clusters with < 3 points cannot support
        # LDA (need at least k observations per class for a non-degenerate
        # within-class covariance estimate). If any cluster is too small,
        # skip this iteration entirely.
        
        sim_env_data <- copy(hmsc_XData)
        setDT(sim_env_data)
        if ("." %in% names(sim_env_data)) sim_env_data[, "." := NULL]
        
        cluster_data <- as.matrix(sim_env_data[, ..selected_vars])
        k <- sample(2:10, 1)
        km_fit <- kmeans(cluster_data, centers = k, nstart = 100, iter.max = 100)
        
        if (any(km_fit$size < 3)) next
        
        # MIN_LD1_VARIANCE: Adaptive threshold for the LDA discrimination
        # check (Filter 0, Section 4.4b). With more clusters, it's natural
        # for variance to be spread across more discriminant axes, so we
        # relax the threshold by 5 percentage points per additional cluster.
        # k=2: threshold = 0.60 (one LD axis must explain ≥ 60%)
        # k=5: threshold = 0.45
        # k=10: threshold = 0.20
        MIN_LD1_VARIANCE <- 0.7
        
        final_clusters   <- km_fit$cluster
        n_clusters_optim <- k
        
        # ----- Compute baseline ASW (before any cluster manipulation) -----
        # ASW (Average Silhouette Width) measures how well-separated the
        # clusters are. Range: [-1, 1]. Higher = better separation.
        # We compute it in two spaces:
        #
        # Subspace ASW: Using only the selected_vars (the variables that
        #   defined the clustering). This is optimistic because k-means
        #   directly optimises separation in this space.
        #
        # Full-space ASW: Using all non-MEM variables. This is the more
        #   honest measure — it asks whether the types are distinguishable
        #   when you consider ALL environmental conditions, not just the
        #   ones used for clustering. A typology that looks great in its
        #   own subspace but terrible in full space is not very useful.
        
        asw_orig_sub <- {
                sil_orig <- silhouette(final_clusters, dist(cluster_data))
                mean(sil_orig[, 3])
        }
        
        full_data_orig <- copy(sim_env_data)
        full_cols      <- setdiff(names(full_data_orig), grep("MEM", names(full_data_orig), value = TRUE))
        full_data_orig <- full_data_orig[, ..full_cols]
        asw_orig_full  <- {
                sil_full <- silhouette(final_clusters, dist(as.matrix(full_data_orig)))
                mean(sil_full[, 3])
        }
        rm(full_data_orig, full_cols, sil_full)
        
        # ------------------------------------------------------------------
        # 4.4 Compute Centroids in Unscaled Space
        # ------------------------------------------------------------------
        # The LDA-based cluster manipulation operates in unscaled space
        # (for the validation pass) and later in scaled space (for HMSC
        # prediction). Here we set up the unscaled version.
        #
        # unscaled_dt:            Copy of unscaled_env with cluster assignments added.
        # unscaled_cols:          Column names excluding "type" — used for subsetting
        #                         throughout this iteration. Defined once, never mutated.
        # centroids_matrix:       k × p matrix of cluster centroids (means of selected
        #                         variables within each cluster), in unscaled units.
        # unscaled_cluster_data:  n × p matrix of observations (selected variables only).
        # unscaled_missing_vars:  The remaining columns (variables NOT selected for
        #                         clustering + any non-env columns). These are carried
        #                         through unchanged when we reconstruct full environments.
        
        unscaled_dt    <- copy(unscaled_env)
        unscaled_cols  <- setdiff(names(unscaled_dt), "type")
        setDT(unscaled_dt)
        unscaled_dt[, type := final_clusters]
        
        centroids_dt <- unscaled_dt[, lapply(.SD, mean), .SDcols = selected_vars, by = "type"]
        setorderv(centroids_dt, "type")
        centroids_dt[, type := NULL]
        centroids_matrix <- as.matrix(as.data.frame(centroids_dt))
        
        unscaled_cluster_data <- as.matrix(as.data.frame(unscaled_dt)[, selected_vars])
        unscaled_missing_vars <- unscaled_dt[, .SD, .SDcols = !selected_vars]
        unscaled_missing_vars[, type := NULL]
        
        # ------------------------------------------------------------------
        # 4.4b Early skip: constant-within-group variables
        # ------------------------------------------------------------------
        # If any selected variable has zero within-group variance for any
        # cluster, LDA will fail. We check this with the pooled within-class 
        # covariance. This matrix must be invertabile for LDA to work.  
        # compute W the pooled within-class covariance. 
                # split() divides the data into a list of sub-matrices, one per
                # cluster. cov gives the pxp covariance matrix. Because 
                # covariance. Multiply by nrow(x) - 1 to remove the Bessel 
                # Correction. Then we sum across all clusters. 
        W <- Reduce("+",
                    lapply(
                                split(
                                        as.data.frame(unscaled_cluster_data),
                                        final_clusters
                                        ), 
                                function(x) cov(x) * (nrow(x) - 1)
                                )
                    )
        
        is_invertible <- tryCatch({
                solve(W)
                TRUE
        }, error = function(e) {
                FALSE
        })
        if (!is_invertible) next
        rm(W, is_invertible)
        has_constant <- FALSE
        for (j in 1:ncol(unscaled_cluster_data)) {
                gv <- tapply(unscaled_cluster_data[, j], final_clusters, var)
                if (any(gv < 1e-8, na.rm = TRUE)) {
                        has_constant <- TRUE
                        break
                }
        }
        if (has_constant) next
        
        # ------------------------------------------------------------------
        # 4.4c FILTER 0: Discriminant Variance Check (iteration-level)
        # ------------------------------------------------------------------
        # Before trying any quality_factor values, check whether the k-means
        # clustering has clear discriminant structure. If the first Linear
        # Discriminant axis (LD1) explains less than MIN_LD1_VARIANCE of the
        # total between-class variance, the clusters are poorly separated in
        # discriminant space — manipulating them along LDA axes won't produce
        # meaningful variation in type distinctness. Skip the iteration.
        #
        # This is an efficiency filter: it avoids wasting time on cluster
        # configurations that are inherently fuzzy regardless of quality_factor.
        
        lda_check <- MASS::lda(
                x        = unscaled_cluster_data,
                grouping = factor(final_clusters)
        )
        ld_eigenvalues   <- lda_check$svd^2
        ld_var_explained <- ld_eigenvalues / max(sum(ld_eigenvalues), 1e-10)
        ld1_variance     <- ld_var_explained[1]
        
        if (ld1_variance < MIN_LD1_VARIANCE) {
                record_id <- record_id + 1L
                diag_records[[record_id]] <- data.table(
                        iteration          = q,
                        param_id           = NA_integer_,
                        n_vars             = n_vars_selected,
                        n_clusters         = n_clusters_optim,
                        asw_orig_sub       = NA_real_,
                        asw_orig_full      = NA_real_,
                        asw_after_sub      = NA_real_,
                        asw_after_full     = NA_real_,
                        separation         = NA_real_,
                        dispersion         = NA_real_,
                        quality_effective  = NA_real_,
                        quality_original   = NA_real_,
                        n_sites_pass_z     = NA_integer_,
                        n_sites_required   = NA_integer_,
                        avg_z_score        = NA_real_,
                        passed_isoforest   = NA,
                        maha_fail_count    = NA_integer_,
                        passed_mahalanobis = NA,
                        ld1_variance       = ld1_variance,
                        final_status       = "FAIL_low_discrimination"
                )
                next()
        }
        
        # ------------------------------------------------------------------
        # 4.5 Generate Parameter Combinations
        # ------------------------------------------------------------------
        # For each iteration that passes Filter 0, we test a range of
        # quality_factor values from 0.5 (strong compression → overlapping
        # types) to 1.5 (strong separation → crisp types) in steps of 0.25.
        #
        # adjust_clusters_lda() may not achieve the exact requested factor
        # due to data geometry constraints, so it returns effective_separation,
        # effective_dispersion, and effective_quality alongside the modified
        # observations.
        #
        # param_id: A stable integer identifier for each parameter combination.
        # This ID is carried through all downstream data structures so that
        # after filtering rows from params, we can still unambiguously link
        # back to the correct entry in temp_sim_envs.
        
        # how many quality_factor values to try per iteration
        n_param_draws <- 25  
        params <- data.frame(
                quality_original = round(runif(n_param_draws, min = 0.5, max = 1.5),2)
        )  %>%  unique
        params$param_id <- seq_len(nrow(params))
        
        # Build temp_sim_envs: a list of manipulated observation matrices,
        # keyed by param_id. Each entry is an n × p matrix (same dimensions
        # as unscaled_cluster_data) with observations shifted along the
        # discriminant axes.
        temp_sim_envs <- list()
        lda_failed <- FALSE
        for (i in 1:nrow(params)) {
                adj_result <- adjust_clusters_lda(
                        centroids           = centroids_matrix,
                        observations        = unscaled_cluster_data,
                        cluster_assignments = final_clusters,
                        quality_factor      = params$quality_original[i]
                )
                if (is.null(adj_result)) {
                        lda_failed <- TRUE
                        break
                }
                
                params$separation[i]                <- adj_result$effective_separation
                params$dispersion[i]                <- adj_result$effective_dispersion
                params$quality_effective[i]         <- adj_result$effective_quality
                temp_sim_envs[[params$param_id[i]]] <- adj_result$observations
        }
        rm(i, adj_result)
        
        if (lda_failed) next
        
        # Quality control: if the monotonic relationship between requested
        # and achieved quality factors breaks down (correlation < 0.9), the
        # LDA manipulation isn't behaving as expected for this data geometry.
        # Skip the iteration rather than produce misleading simulations.
        cor_val <- suppressWarnings(cor(params$quality_original, params$quality_effective))
        if (is.na(cor_val) || cor_val < 0.9) {
                rm(cor_val)
                next
        }
        rm(cor_val)
        
        # Filter: remove parameter combinations where dispersion (within-cluster
        # spread) has grown too large relative to separation (between-centroid
        # distance). If diff_rel = (dispersion - separation) / separation ≥ 1,
        # the clusters have become so dispersed that they overlap heavily
        # regardless of the centroid positions. This is a geometric sanity check.
        params <- mutate(params, diff = dispersion - separation)
        params <- mutate(params, diff_rel = diff / separation)
        params <- dplyr::filter(params, diff_rel < 1)
        if (nrow(params) == 0) next
        params$bin <- get_bin(params$quality_original)
        params <- params[bin_success_counts[params$bin] < per_bin_quota, ]
        if (nrow(params) == 0) next
        # ------------------------------------------------------------------
        # 4.6 FILTER 1: Isolation Forest
        # ------------------------------------------------------------------
        # For each surviving parameter combination, reconstruct the full
        # (all-variable) unscaled environment by combining the manipulated
        # selected variables with the unchanged non-selected variables.
        # Then score every site with the isolation forest.
        #
        # The filtering is two-stage:
        #   Stage A (per-site): Remove sites with Z-score > iso_threshold2.
        #     These individual sites are implausible outliers.
        #   Stage B (per-environment): After removing outlier sites, count how
        #     many remain. If fewer than (1 - iso_threshold) × n_samples = 90%
        #     survive, the entire simulated environment is rejected — too many
        #     of its sites are outside the real-world envelope.
        
        validation_results <- list()
        
        for (i in 1:nrow(params)) {
                # Reconstruct full environment: manipulated selected vars +
                # unchanged non-selected vars
                current_data <- cbind(
                        temp_sim_envs[[params$param_id[i]]],
                        unscaled_missing_vars
                )
                # Remove MEM columns (spatial eigenvectors are not part of the
                # environmental plausibility check)
                cols_to_keep <- setdiff(
                        names(current_data),
                        grep("MEM", names(current_data), value = TRUE)
                )
                current_data <- current_data[, ..cols_to_keep]
                # Ensure column order matches the isolation forest's training data
                current_data <- current_data[, ..unscaled_cols]
                
                current_scores <- predict(isolation_forest, current_data)
                
                iter_ratings <- data.table(
                        scores          = current_scores,
                        scores_centered = current_scores - score_mean,
                        separation      = params$separation[i],
                        dispersion      = params$dispersion[i],
                        quality_effective = params$quality_effective[i],
                        quality_original  = params$quality_original[i]
                )
                iter_ratings[, scores_z := scores_centered / score_sd]
                # Use param_id (stable across filtering) rather than the loop
                # counter i (which is a row index into the filtered params)
                iter_ratings[, id := params$param_id[i]]
                validation_results[[i]] <- iter_ratings
                rm(current_data, current_scores, iter_ratings)
        }; rm(i)
        validation_results <- rbindlist(validation_results)
        
        # Stage A: Per-site Z-score filter
        n_sites_before_z   <- nrow(validation_results)
        validation_results <- validation_results[scores_z <= iso_threshold2]
        n_sites_after_z    <- nrow(validation_results)
        
        # Aggregate per simulation ID: how many sites survived, and what's
        # their average anomaly score?
        validation_results[, id_count := .N, by = "id"]
        validation_results[, avg_z := mean(scores_z), by = "id"]
        
        # Stage B: Count threshold — at least 90% of original sites must survive
        threshold_count <- floor(current_scheme$samples - iso_threshold * current_scheme$samples)
        
        # Snapshot for diagnostics: record pass/fail for ALL param IDs before
        # we filter the data.table down to survivors only.
        id_summary_pre <- unique(validation_results[, .(id, id_count, avg_z,
                                                        separation, dispersion,
                                                        quality_effective, quality_original)])
        id_summary_pre[, passed_isoforest := id_count >= threshold_count]
        
        # Keep only environments where enough sites passed
        validation_results <- validation_results[id_count >= threshold_count]
        validation_results <- unique(validation_results, by = "id")
        
        if (nrow(validation_results) == 0) {
                # All parameter combos failed the ISO filter. Record each one
                # with final_status = "FAIL_isoforest" and move on.
                for (row_i in 1:nrow(id_summary_pre)) {
                        record_id <- record_id + 1L
                        diag_records[[record_id]] <- data.table(
                                iteration          = q,
                                param_id           = id_summary_pre$id[row_i],
                                n_vars             = n_vars_selected,
                                n_clusters         = n_clusters_optim,
                                asw_orig_sub       = asw_orig_sub,
                                asw_orig_full      = asw_orig_full,
                                asw_after_sub      = NA_real_,
                                asw_after_full     = NA_real_,
                                separation         = NA_real_,
                                dispersion         = NA_real_,
                                quality_effective  = NA_real_,
                                quality_original   = NA_real_,
                                n_sites_pass_z     = id_summary_pre$id_count[row_i],
                                n_sites_required   = threshold_count,
                                avg_z_score        = id_summary_pre$avg_z[row_i],
                                passed_isoforest   = FALSE,
                                maha_fail_count    = NA_integer_,
                                passed_mahalanobis = NA,
                                ld1_variance       = ld1_variance,
                                final_status       = "FAIL_isoforest"
                        )
                }
                next()
        }
        
        # ------------------------------------------------------------------
        # 4.7 FILTER 2: Mahalanobis Distance
        # ------------------------------------------------------------------
        # For environments that passed the Isolation Forest, apply a second
        # parametric plausibility check: Mahalanobis distance from the
        # catchment centroid.
        #
        # For each simulated environment, count how many sites exceed the
        # chi-squared threshold. If more than 10% of sites are Mahalanobis
        # outliers, the environment is rejected.
        #
        # The cor_mat / cor_vec computation detects near-perfect correlations
        # (> 0.98) in the catchment covariance matrix. If present, a small
        # ridge (1e-6 × I) is added to the covariance matrix to prevent
        # numerical instability in the matrix inversion required by
        # mahalanobis().
        
        fail_counts <- numeric(nrow(validation_results))
        
        for (i in 1:nrow(validation_results)) {
                sim_id <- validation_results$id[i]
                
                current_data <- cbind(temp_sim_envs[[sim_id]], unscaled_missing_vars)
                cols_to_keep <- setdiff(names(current_data), grep("MEM", names(current_data), value = TRUE))
                current_data <- current_data[, ..cols_to_keep]
                current_data <- current_data[, ..unscaled_cols]
                
                cor_mat <- cov2cor(env_cov_mat)
                cor_vec <- cor_mat[lower.tri(cor_mat, diag = FALSE)]
                
                if (any(abs(cor_vec) > 0.98)) {
                        temp_cov_mat <- env_cov_mat + diag(1e-6, nrow(env_cov_mat))
                        m_dist <- mahalanobis(current_data, center = env_mean, cov = temp_cov_mat)
                } else {
                        m_dist <- mahalanobis(current_data, center = env_mean, cov = env_cov_mat)
                }
                
                fail_counts[i] <- sum(m_dist > maha_threshold)
                rm(current_data, m_dist, cor_mat, cor_vec)
        }; rm(i)
        
        validation_results$failsum <- fail_counts
        rm(fail_counts)
        
        # Accept if fewer than 10% of sites are Mahalanobis outliers
        validation_results[, passed_mahalanobis := failsum < (0.10 * current_scheme$samples)]
        
        # ------------------------------------------------------------------
        # 4.8 Record Diagnostics for ALL Parameter Combos This Iteration
        # ------------------------------------------------------------------
        # At this point we have three groups:
        #   1. IDs that failed ISO filter (in id_summary_pre but not validation_results)
        #   2. IDs that passed ISO but failed Mahalanobis (in validation_results, passed_mahalanobis = FALSE)
        #   3. IDs that passed both (PASS)
        # Record all three for the diagnostic table.
        
        # --- Group 1: ISO failures ---
        failed_iso_ids <- setdiff(id_summary_pre$id, validation_results$id)
        for (fid in failed_iso_ids) {
                row_data <- id_summary_pre[id == fid]
                record_id <- record_id + 1L
                diag_records[[record_id]] <- data.table(
                        iteration          = q,
                        param_id           = fid,
                        n_vars             = n_vars_selected,
                        n_clusters         = n_clusters_optim,
                        asw_orig_sub       = asw_orig_sub,
                        asw_orig_full      = asw_orig_full,
                        asw_after_sub      = NA_real_,
                        asw_after_full     = NA_real_,
                        quality_effective  = row_data$quality_effective,
                        quality_original   = row_data$quality_original,
                        separation         = row_data$separation,
                        dispersion         = row_data$dispersion,
                        n_sites_pass_z     = row_data$id_count,
                        n_sites_required   = threshold_count,
                        avg_z_score        = row_data$avg_z,
                        passed_isoforest   = FALSE,
                        maha_fail_count    = NA_integer_,
                        passed_mahalanobis = NA,
                        ld1_variance       = ld1_variance,
                        final_status       = "FAIL_isoforest"
                )
        }
        
        # --- Groups 2 & 3: Build scaled environments for ISO-pass combos ---
        # For parameter combinations that passed the ISO filter, we now need
        # to create the *scaled* (standardised) version of the manipulated
        # environment, because the HMSC model expects scaled inputs.
        #
        # The manipulation is re-run in scaled space with the *effective*
        # quality_factor (not the original requested one), ensuring the
        # scaled environment has the same relative cluster structure as the
        # validated unscaled one.
        
        base_scaled_data <- copy(hmsc_XData)
        setDT(base_scaled_data)
        base_scaled_data[, type := final_clusters]
        
        centroids_scaled     <- base_scaled_data[, lapply(.SD, mean), .SDcols = selected_vars, by = "type"]
        setorderv(centroids_scaled, "type")
        centroids_scaled[, type := NULL]
        centroids_scaled_mat <- as.matrix(as.data.frame(centroids_scaled))
        
        scaled_cluster_data  <- as.matrix(as.data.frame(base_scaled_data)[, selected_vars])
        scaled_missing_vars  <- base_scaled_data[, .SD, .SDcols = !selected_vars]
        scaled_missing_vars[, type := NULL]
        
        # validated_envs: Accumulates full prediction-ready environments
        # (with intercept column) for combinations that pass both filters.
        validated_envs     <- list()
        
        # Identify Maha-pass combos and select up to 5
        valid_indices <- which(validation_results$passed_mahalanobis == 1)
        if (length(valid_indices) == 0) next

        
        valid_param_ids <- validation_results$id[valid_indices]
        valid_bins <- params$bin[match(valid_param_ids, params$param_id)]
        bin_deficits <- pmax(per_bin_quota - bin_success_counts, 0)
        
        # For each bin, pick at most one — the one whose quality_original is
        # closest to the bin center (most representative)
        keep <- integer(0)
        for (b in seq_len(n_bins)) {
                candidates <- which(valid_bins == b & bin_deficits[b] > 0)
                if (length(candidates) == 0) next
                
                candidate_ids <- valid_param_ids[candidates]
                candidate_q <- params$quality_original[match(candidate_ids, params$param_id)]
                bin_center <- (quality_bins[b] + quality_bins[b + 1]) / 2
                best <- candidates[which.min(abs(candidate_q - bin_center))]
                keep <- c(keep, best)
        }
        
        if (length(keep) == 0) next
        valid_indices <- valid_indices[keep]
        
        # Now subset validation_results to ONLY the selected 5
        validation_results <- validation_results[valid_indices]
        
        
        asw_after_sub_vec  <- numeric(nrow(validation_results))
        asw_after_full_vec <- numeric(nrow(validation_results))
        
        for (row_i in 1:nrow(validation_results)) {
                vr          <- validation_results[row_i]
                passed_maha <- vr$passed_mahalanobis
                
                asw_after_sub_val  <- NA_real_
                asw_after_full_val <- NA_real_
                
                if (passed_maha) {
                        # Re-run LDA manipulation in scaled space
                        # NOTE: We index params using match() because params has
                        # been filtered (rows removed by dplyr::filter), so
                        # positional indexing would be wrong. match()
                        # finds the row where param_id == vr$id.
                        result <- adjust_clusters_lda(
                                centroids           = centroids_scaled_mat,
                                observations        = scaled_cluster_data,
                                cluster_assignments = final_clusters,
                                quality_factor      = params$quality_effective[match(vr$id, params$param_id)]
                        )
                        # LDA failed on scaled data — treat as failed combo
                        if (is.null(result)) {
                                passed_maha <- FALSE
                        }
                }
                if (passed_maha){
                        
                        # Post-manipulation ASW (subspace)
                        sil_sub <- silhouette(final_clusters, dist(result$observations))
                        asw_after_sub_val <- mean(sil_sub[, 3])
                        
                        # Post-manipulation ASW (full space, excluding MEM)
                        new_env_diag <- cbind(result$observations, scaled_missing_vars)
                        mem_cols     <- grep("MEM", names(new_env_diag), value = TRUE)
                        if (length(mem_cols) > 0) new_env_diag[, (mem_cols) := NULL]
                        dist_mat           <- stats::dist(new_env_diag)
                        sil_full           <- silhouette(final_clusters, dist = dist_mat)
                        asw_after_full_val <- mean(sil_full[, 3])
                        
                        # Construct HMSC prediction matrix:
                        # Combine manipulated selected vars + unchanged non-selected vars,
                        # reorder to match model's expected column order, convert to matrix,
                        # and prepend an intercept column.
                        new_env_full <- cbind(result$observations, scaled_missing_vars)
                        col_order    <- names(hmsc_model$XData)
                        dt_env       <- data.frame(new_env_full)
                        setDT(dt_env)
                        dt_env  <- dt_env[, ..col_order]
                        mat_env <- as.matrix(data.frame(dt_env))
                        mat_env <- cbind(rep(1, nrow(mat_env)), mat_env)
                        colnames(mat_env)[1] <- "(Intercept)"
                        
                        validated_envs[[length(validated_envs) + 1]] <- mat_env
                        
                        rm(result, new_env_diag, new_env_full, dist_mat,
                           sil_sub, sil_full, dt_env, mat_env)
                }
                
                asw_after_sub_vec[row_i]  <- asw_after_sub_val
                asw_after_full_vec[row_i] <- asw_after_full_val
                
                # Record diagnostic row for this parameter combo
                status <- if (passed_maha) "PASS" else "FAIL_mahalanobis"
                
                record_id <- record_id + 1L
                diag_records[[record_id]] <- data.table(
                        iteration          = q,
                        param_id           = vr$id,
                        n_vars             = n_vars_selected,
                        n_clusters         = n_clusters_optim,
                        asw_orig_sub       = asw_orig_sub,
                        asw_orig_full      = asw_orig_full,
                        asw_after_sub      = asw_after_sub_val,
                        asw_after_full     = asw_after_full_val,
                        separation         = params$separation[match(vr$id, params$param_id)],
                        dispersion         = params$dispersion[match(vr$id, params$param_id)],
                        quality_original   = params$quality_original[match(vr$id, params$param_id)],
                        quality_effective  = params$quality_effective[match(vr$id, params$param_id)],
                        n_sites_pass_z     = vr$id_count,
                        n_sites_required   = threshold_count,
                        avg_z_score        = vr$avg_z,
                        passed_isoforest   = TRUE,
                        maha_fail_count    = vr$failsum,
                        passed_mahalanobis = passed_maha,
                        ld1_variance       = ld1_variance,
                        final_status       = status
                )
        }
        
        rm(base_scaled_data, centroids_scaled, centroids_scaled_mat,
           scaled_cluster_data, scaled_missing_vars, col_order)
        
        # ------------------------------------------------------------------
        # 4.9 Update Counters
        # ------------------------------------------------------------------
        # Identify which parameter combos passed both filters.
        n_valid <- length(validated_envs)
        
        success_count <- success_count + n_valid
        
        for (vid in validation_results$id) {
                b <- params$bin[match(vid, params$param_id)]
                bin_success_counts[b] <- bin_success_counts[b] + 1
        }
        
        
        # Deprioritise variables that contributed to successful simulations,
        # encouraging future iterations to explore different variable subsets.
        var_selection_counter[selected_vars] <- var_selection_counter[selected_vars] + 1
        
        # ------------------------------------------------------------------
        # 4.10 Production: Store Metadata
        # ------------------------------------------------------------------
        # Record metadata for each validated simulation. All production storage
        # vectors are indexed in parallel.
        validated_param_ids <- validation_results$id
        
        
        res_n_valid      <- append(res_n_valid, n_valid)
        res_n_clusters   <- append(res_n_clusters, rep(n_clusters_optim, n_valid))
        res_n_variables  <- append(res_n_variables, rep(n_vars_selected, n_valid))
        res_contr_points <- append(res_contr_points,
                                   params$dispersion[match(validated_param_ids, params$param_id)])
        res_contr_centro <- append(res_contr_centro,
                                   params$separation[match(validated_param_ids, params$param_id)])
        
        current_importance <- sum(vp_data$importance[selected_vars])
        res_importance     <- append(res_importance, rep(current_importance, n_valid))
        
        asw_scores <- asw_after_full_vec
        res_asw    <- append(res_asw, asw_scores)
        
        for (k in 1:n_valid) {
                res_cluster_ids[[length(res_cluster_ids) + 1]] <- final_clusters
        }
        
        # ------------------------------------------------------------------
        # 4.11 Production: Fuzzy Classification & NPE
        # ------------------------------------------------------------------
        # Run Fuzzy C-Means (FCM) on each validated environment using only
        # the selected variables. The fuzziness parameter m = 1.5 (moderate
        # fuzziness; m = 1 → hard clustering, m → ∞ → uniform membership).
        #
        # mobileCenters = n_clusters_optim tells vegclust to initialise with
        # the specified number of mobile (optimisable) cluster centers.
        #
        # NPE is then computed from each FCM result. This gives us a measure
        # of how crisp the fuzzy classification is for each simulated
        # environment — environments with well-separated types should have
        # low NPE, environments with overlapping types should have high NPE.
        
        fuzzy_results <- lapply(1:n_valid, function(x) {
                vegclust(
                        validated_envs[[x]][, which(colnames(validated_envs[[x]]) %in% selected_vars)],
                        mobileCenters = n_clusters_optim,
                        method        = "FCM",
                        m             = 1.5
                )
        })
        
        res_fuzzy_ids <- append(res_fuzzy_ids, fuzzy_results)
        
        current_npe <- sapply(fuzzy_results, NPE)
        res_npe     <- append(res_npe, current_npe)
        
        # ------------------------------------------------------------------
        # 4.12 Production: HMSC Community Predictions
        # ------------------------------------------------------------------
        # For each validated environment, predict species occurrence
        # probabilities using the fitted HMSC model.
        #
        # We subsample 5 random posterior draws from the last 200 samples of
        # the pooled chains. Using a small subsample (rather than all ~4000
        # pooled samples) is a computational compromise — predict() is
        # expensive (matrix multiplication for every species × every posterior
        # sample × every site), and we're generating 60 simulated environments.
        #
        # poolMcmcChains() concatenates the two MCMC chains into a single list
        # of posterior samples. chain_indices selects 5 samples near the end
        # of the chain (where mixing is best).
        #
        # predict() returns a list of matrices, one per posterior sample. Each
        # matrix is n_sites × n_species with occurrence probabilities.
        
        chain_indices    <- sample(x = (hmsc_model$samples - 200):hmsc_model$samples,
                                   size = 5, replace = FALSE)
        pooled_posterior <- poolMcmcChains(hmsc_model$postList)
        selected_samples <- pooled_posterior[chain_indices]
        
        predictions <- lapply(validated_envs, function(env_matrix) {
                predict(
                        object = hmsc_model,
                        post   = selected_samples,
                        X      = env_matrix
                )
        })
        
        # Store as a single list element per iteration. Downstream code uses
        # unlist(, recursive = FALSE) to flatten, then rapply() to walk into
        # the nested structure for distance matrix computation.
        res_sim_output[[length(res_sim_output) + 1]] <- predictions
        
        # Cleanup: remove all iteration-specific objects to free memory.
        rm(unscaled_dt, unscaled_cols, unscaled_missing_vars,
           centroids_dt, centroids_matrix, unscaled_cluster_data,
           params, temp_sim_envs, validation_results, id_summary_pre,
           validated_envs, predictions, fuzzy_results,
           validated_param_ids, asw_scores)
}

cat(sprintf("\n\nLoop complete. Total successes: %d\n", success_count))


# ==============================================================================
# 5. COMPILE AND SAVE DIAGNOSTIC REPORT
# ==============================================================================
# Bind all diagnostic records into a single data.table and save.
# fill = TRUE handles the fact that different failure paths produce different
# column sets (e.g., FAIL_low_discrimination records have ld1_variance but
# no maha_fail_count; PASS records have both).

if (length(diag_records) == 0) {
        cat("No records generated — all iterations failed before parameter generation.\n")
        stop("ERROR: No diagnostic records to report.")
}

diag_dt <- rbindlist(diag_records, fill = TRUE)

# ==============================================================================
# 6. SAVE PRODUCTION OUTPUTS
# ==============================================================================
# Flatten the nested prediction list and cap at N_TARGET_SUCCESSES.
# The output list follows the format expected by downstream analysis scripts
# (e.g., the distance matrix computation pipeline).
#
# Key fields:
#   data:                     List of HMSC predictions (each element is a list
#                             of posterior-draw matrices).
#   number_of_clusters:       k for each simulation.
#   contraction_points:       Effective dispersion (within-cluster spread).
#   contraction_centroids:    Effective separation (between-centroid distance).
#   variable_importance:      Sum of VP importance for selected variables.
#   asw:                      Full-space ASW after cluster manipulation.
#   npe:                      Normalised Partitioning Entropy from FCM.
#   hard_cluster_assignment:  Integer vector of cluster labels per site.
#   fuzzy_cluster_assignment: vegclust objects (contain membership matrices).

if (length(res_sim_output) == 0) {
        cat("\nAll simulations failed — no production output.\n")
        stop("ERROR: All simulations failed")
}

cat("\nPreparing production output...\n")

flat_sim_output   <- unlist(res_sim_output, recursive = FALSE)
final_index_limit <- min(length(flat_sim_output), N_TARGET_SUCCESSES)
idx_range         <- 1:final_index_limit

# Adjust res_n_valid so cumulative sums match idx_range
cumul <- cumsum(res_n_valid)
# How many full iterations fit?
last_full <- max(which(cumul <= final_index_limit), 0)

if (last_full < length(res_n_valid)) {
        # Truncate: keep full iterations + partial last one
        remainder <- final_index_limit - ifelse(last_full > 0, cumul[last_full], 0)
        res_n_valid <- c(res_n_valid[seq_len(last_full)], remainder)
}

# Mark truncated PASS rows in diagnostics
pass_rows <- which(diag_dt$final_status == "PASS")
if (length(pass_rows) > final_index_limit) {
        truncated <- pass_rows[(final_index_limit + 1):length(pass_rows)]
        diag_dt[truncated, final_status := "PASS_truncated"]
}

if (!dir.exists("data/misc/simulation_diagnostics")) {
        dir.create("data/misc/simulation_diagnostics")
}
diag_output_file <- paste0("data/misc/simulation_diagnostics/", model_name)
saveRDS(diag_dt, diag_output_file)
cat(sprintf("\nDiagnostic table saved to: %s\n", diag_output_file))
cat(sprintf("  Rows: %d | Columns: %s\n", nrow(diag_dt), paste(names(diag_dt), collapse = ", ")))


out_list <- list(
        data                     = flat_sim_output[idx_range],
        number_of_clusters       = res_n_clusters[idx_range],
        number_of_variables      = res_n_variables[idx_range],
        contraction_points       = res_contr_points[idx_range],
        contraction_centroids    = res_contr_centro[idx_range],
        variable_importance      = res_importance[idx_range],
        n_valid_per_iteration    = res_n_valid[idx_range],
        asw                      = res_asw[idx_range],
        npe                      = res_npe[idx_range],
        hard_cluster_assignment  = res_cluster_ids[idx_range],
        fuzzy_cluster_assignment = res_fuzzy_ids[idx_range]
)

saveRDS(out_list, output_file)
cat(sprintf("Production results saved to: %s\n", output_file))
cat(sprintf("  Simulations: %d | Clusters: %s | NPE range: [%.3f, %.3f]\n",
            final_index_limit,
            paste(unique(res_n_clusters[idx_range]), collapse = ","),
            min(res_npe[idx_range]), max(res_npe[idx_range])))
cat("\nDone.\n")
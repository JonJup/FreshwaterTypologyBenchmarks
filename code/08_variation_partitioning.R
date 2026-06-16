################################################################################
# Script Name: 08_variation_partitioning.R
# Description: For a single fitted HMSC model, import the posterior produced on
#              the HPC, evaluate model fit (Tjur R2), and partition each taxon's
#              explained variance into environmental, spatial (MEM) and biotic
#              (random-effect) components, plus an unexplained "stochastic"
#              fraction. Also derives a normalised per-variable importance score.
#              Writes one .rds (VP table + importance + model-fit object).
# Notes:       Run as a SLURM array job, one model per task:
#                Rscript 08_variation_partitioning.R --iter_id <ID> \
#                        --input <unfitted_model.rds> --output <out.rds>
#              Models that fail PSRF or AUC, or score < 4 across the five
#              posterior-predictive checks, are skipped (clean exit, no output).
################################################################################

# =======================================================
# 1. setup ----
# =======================================================
args <- commandArgs(trailingOnly = TRUE)

# Default values (optional)
iter_id <- NULL
input_file <- NULL
output_file <- NULL

# Simple manual parsing of the named --flag value pairs passed by the launcher
for (i in seq_along(args)) {
        if (args[i] == "--iter_id") {
                iter_id <- args[i + 1]
        } else if (args[i] == "--input") {
                input_file <- args[i + 1]
        } else if (args[i] == "--output") {
                output_file <- args[i + 1]
        }
}
# Sanity check: both --input and --output are required to do anything useful
if (is.null(input_file) || is.null(output_file)) {
        stop("ERROR: Missing --input or --output argument.")
}
cat("\n========================================\n")
cat("Processing Iteration:", iter_id, "\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n")
cat("========================================\n\n")

suppressPackageStartupMessages({
        library(Hmsc)
        library(data.table)
        library(jsonify)
})

# --- 1.4 Load Data ---
cat("Loading input data...\n")
# Read your input data

# Derive the bare model name (e.g. "diatoms_0001") from the unfitted-model path.
# It is reused as the key to locate every associated file for this model.
model.name <- gsub(pattern="data/001_unfitted_hmsc_models/","", input_file)
model.name <- gsub(pattern="\\.rds", "", model.name)
cat("model name is", model.name, "\n")

# Anchor the look-up pattern to the START of the file name so a name like
# "diatoms_0001" does not also match "diatoms_00011" (substring collision).
file_pattern        <- paste0("^", model.name)
# The fitted model is stored as several chain files -> a vector is expected here.
target_fitted_files <- list.files("data/003_fitted_hmsc_models/", full.names = T, pattern = file_pattern)
# The model-fit / evaluation object is a single file for this model.
eval        <- list.files("data/004_model_fit/", full.names = T, pattern = file_pattern)
cat("eval is", eval, "\n")

# Guard: the evaluation look-up must resolve to exactly one file. A bare
# readRDS() on a length-0 (not found) or length-2 (ambiguous) vector would
# otherwise fail with a cryptic error; stop loudly with a clear message instead.
if (length(eval) != 1L) {
        stop(sprintf("Expected exactly one model-fit file for '%s', found %d.",
                     model.name, length(eval)))
}

cat("eval read ...")
eval          <- readRDS(eval)
cat("successful", "\n")
# Check if PSRF (chain-convergence) test failed -> abort: model is untrustworthy
if (!eval$psrf_passed){

        cat("PSRF failed \n")
        stop("PSRF failed")

}
# Check if AUC (predictive discrimination) test failed -> abort
if (!eval$AUC_passed){
        
        cat("AUC failed \n")
        stop("AUC failed")
        
}

# Overall quality score = number of the five posterior-predictive checks passed
eval[, good_model := psrf_passed + AUC_passed + betaDistr_passed + c_score_passed + prevalence_passed]
# A failed (NA) check contributes 0 rather than propagating NA into the score
eval[is.na(good_model), good_model := 0]
cat(eval$good_model, "\n")

# Stop for models that have failed two or three posterior predictive checks.
# (PSRF + AUC already passed above, so < 4 means < 2 of the remaining checks.)
if (eval$good_model < 4) {
        cat("Skipping: Model quality insufficient (score < 4)\n")
        quit(save = "no", status = 0) 
}

# Read fitted models (Assuming list structure where [[1]] is the JSON)
# We read them into a temp list first to avoid variable confusion
loaded_json_list <- lapply(target_fitted_files, function(x) readRDS(x)[[1]])
# Read unfitted model 
unfittedModel    <- readRDS(input_file)

# 1.4 Unpack JSON --------------------------------------------------------------
cat("Unpacking JSON models...\n")
# Using a large buffer size so the (potentially big) JSON chains deserialise
fit_model_objects <- lapply(loaded_json_list, from_json, buffer_size = 524288000)

# Assuming we strictly want the first two chains (or files); both must exist.
postList <- list(fit_model_objects[[1]][[1]], fit_model_objects[[2]][[1]])

cat("Models unpacked. Importing posterior...\n")

# Re-attach the HPC-sampled posterior to the unfitted Hmsc object so that the
# usual Hmsc post-processing functions can operate on it.
fitTF = importPosteriorFromHPC(
                m = unfittedModel,  # the Hmsc object containing the unfitted model
                postList = postList,
                nSamples = 2150,
                thin = 150,
                transient = 10000)

# ==============================================================================
# x. COMPUTE AND EVALUATE PREDICTIONS ============================================
# ==============================================================================

# Compute predicted values from model 
i.preds     <- computePredictedValues(fitTF)
# Evaluate the model fit 
i.MF        <- evaluateModelFit(hM = fitTF,
                                predY = i.preds)
# Extract numeric vector of Tjur R2 values for each species 
R2 <- i.MF$TjurR2


# How many spatial predictors (Moran's Eigenvector Maps) did the model use?
i.n.spatial <- sum(grepl(x = colnames(unfittedModel$XData), pattern = "MEM"))

# ==============================================================================
# x. Variation Partitioning ====================================================
# ==============================================================================
# NOTE: the env/(space)/bio labelling below assumes computeVariancePartitioning
# returns `vals` with one row per fixed-effect group PLUS exactly one random
# effect (the biotic component). The recycled driver labels rely on that layout.

if (i.n.spatial == 0){
        # No spatial predictors: a single "env" group covers all fixed effects.
        i.VP    <- computeVariancePartitioning(
                fitTF,
                group = c(rep(
                        1, ncol(unfittedModel$XData) - i.n.spatial
                )),
                groupnames = c("env")
        )
        # vals has 2 rows/taxon (env, random) -> label them env & bio per taxon.
        i.VP2   <- data.table(
                taxon = rep(colnames(i.VP$vals), each = 2),
                driver = rep(c("env",  "bio"), times = ncol(i.VP$vals)),
                value  = c(i.VP$vals),
                r2     = rep(R2, each = 2)
        )
        # Add an explicit zero "space" component so all models share a schema.
        i.VP22 <- data.table(
                taxon = colnames(i.VP$vals),
                driver = "space",
                value  = 0,
                r2     =R2
        )
        i.VP2 <- rbindlist(list(i.VP2, i.VP22))
} else {
        # With MEMs present: group 1 = environment, group 2 = space (the MEMs).
        i.VP    <- computeVariancePartitioning(
                fitTF,
                group = c(rep(
                        1, ncol(unfittedModel$XData) - i.n.spatial
                ), rep(2, i.n.spatial)),
                groupnames = c("env", "space")
        )
        # vals has 3 rows/taxon (env, space, random) -> env, space, bio per taxon.
        i.VP2   <- data.table(
                taxon = rep(colnames(i.VP$vals), each = 3),
                driver = rep(c("env", "space", "bio"), times = ncol(i.VP$vals)),
                value  = c(i.VP$vals),
                r2     = rep(R2, each = 3)
        )
}
# replace negative scores with zeros (Tjur R2 can dip slightly below 0)
i.VP2[r2 < 0, r2 := 0]
# scale values of env, space, and bio by R2 to make room for stochasticity
i.VP2[, scaled_values := value * r2]
# add stochasticity = 1-R2 (the share of variance the model leaves unexplained)
i.VP3 <- data.table(
        taxon = colnames(i.VP$vals),
        driver = "stochastic",
        value = 0,
        r2 = R2
)
i.VP3[r2 < 0, r2 := 0]
i.VP3[, value := 1 - r2]
i.VP3[, scaled_values := 1 - r2]
# join stochasticity values to env, space, bio (per taxon the four now sum to 1)
i.VP4 <- rbindlist(list(i.VP2, i.VP3))
setorderv(i.VP4, "taxon")
i.VP4$scheme <- model.name

# new VP to determine relative predictor importance (each predictor on its own)
i.VP5 <- computeVariancePartitioning(fitTF)
# total variance explained by each predictor, summed across taxa
i.VP6 <- rowSums(i.VP5$vals)

# drop Moran's eigenvectors and random terms (keep only "real" env predictors)
i.VP6 <- i.VP6[!grepl(x = names(i.VP6), pattern = "MEM|Random")]
# normalize variable importance so the remaining predictors sum to 1
i.VP6 <- i.VP6 / sum(i.VP6)

# =======================================================
# x. Export Results ----
# =======================================================

saveRDS(
        list(
                VP = i.VP4, 
                importance = i.VP6, 
                MF = i.MF
                ), 
        output_file
        )
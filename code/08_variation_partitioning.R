### - Variation Partitioning - ###

# =======================================================
# 1. setup ----
# =======================================================
args <- commandArgs(trailingOnly = TRUE)

# Default values (optional)
iter_id <- NULL
input_file <- NULL
output_file <- NULL

# Simple manual parsing
for (i in seq_along(args)) {
        if (args[i] == "--iter_id") {
                iter_id <- args[i + 1]
        } else if (args[i] == "--input") {
                input_file <- args[i + 1]
        } else if (args[i] == "--output") {
                output_file <- args[i + 1]
        }
}
# Sanity check
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

model.name <- gsub(pattern="data/001_unfitted_hmsc_models/","", input_file)
model.name <- gsub(pattern="\\.rds", "", model.name)
cat("model name is", model.name, "\n")
target_fitted_files <- list.files("data/003_fitted_hmsc_models/", full.names = T, pattern = model.name)
eval        <- list.files("data/004_model_fit/", full.names = T, pattern = model.name)
cat("eval is", eval, "\n")

cat("eval read ...")
eval          <- readRDS(eval)
cat("successful", "\n")
# Check if PSRF test failed
if (!eval$psrf_passed){

        cat("PSRF failed \n")
        stop("PSRF failed")

}
if (!eval$AUC_passed){
        
        cat("AUC failed \n")
        stop("AUC failed")
        
}

eval[, good_model := psrf_passed + AUC_passed + betaDistr_passed + c_score_passed + prevalence_passed]
eval[is.na(good_model), good_model := 0]
cat(eval$good_model, "\n")

# Stop for models that have failed two or three posterior predictive checks
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
# Using a large buffer size 
fit_model_objects <- lapply(loaded_json_list, from_json, buffer_size = 524288000)

# Assuming we strictly want the first two chains (or files)
postList <- list(fit_model_objects[[1]][[1]], fit_model_objects[[2]][[1]])

cat("Models unpacked. Importing posterior...\n")

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

if (i.n.spatial == 0){
        i.VP    <- computeVariancePartitioning(
                fitTF,
                group = c(rep(
                        1, ncol(unfittedModel$XData) - i.n.spatial
                )),
                groupnames = c("env")
        )
        i.VP2   <- data.table(
                taxon = rep(colnames(i.VP$vals), each = 2),
                driver = rep(c("env",  "bio"), times = ncol(i.VP$vals)),
                value  = c(i.VP$vals),
                r2     = rep(R2, each = 2)
        )
        i.VP22 <- data.table(
                taxon = colnames(i.VP$vals),
                driver = "space",
                value  = 0,
                r2     =R2
        )
        i.VP2 <- rbindlist(list(i.VP2, i.VP22))
} else {
        i.VP    <- computeVariancePartitioning(
                fitTF,
                group = c(rep(
                        1, ncol(unfittedModel$XData) - i.n.spatial
                ), rep(2, i.n.spatial)),
                groupnames = c("env", "space")
        )
        i.VP2   <- data.table(
                taxon = rep(colnames(i.VP$vals), each = 3),
                driver = rep(c("env", "space", "bio"), times = ncol(i.VP$vals)),
                value  = c(i.VP$vals),
                r2     = rep(R2, each = 3)
        )
}
# replace negative scores with zeros
i.VP2[r2 < 0, r2 := 0]
# scale values of env, space, and bio by R2 to make room for stochasticity
i.VP2[, scaled_values := value * r2]
# add stochasticity = 1-R2
i.VP3 <- data.table(
        taxon = colnames(i.VP$vals),
        driver = "stochastic",
        value = 0,
        r2 = R2
)
i.VP3[r2 < 0, r2 := 0]
i.VP3[, value := 1 - r2]
i.VP3[, scaled_values := 1 - r2]
# join stochasticity values to env, space, bio
i.VP4 <- rbindlist(list(i.VP2, i.VP3))
setorderv(i.VP4, "taxon")
i.VP4$scheme <- model.name

# new VP to determine relative predictor importance
i.VP5 <- computeVariancePartitioning(fitTF)
i.VP6 <- rowSums(i.VP5$vals)

# drop Moran's eigenvectors and random terms 
i.VP6 <- i.VP6[!grepl(x = names(i.VP6), pattern = "MEM|Random")]
# normalize variable importance
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


################################################################################
# Simulation of manipulated typologies and HMSC community realisations
#
# For one fitted HMSC model (one taxon x data set x sampling scheme), the
# script repeatedly
#   1. draws a classification family (k-means, spatially constrained Ward, or
#      crossed thresholds), balanced so that each family contributes about the
#      same number of retained environments,
#   2. draws a random subset S of environmental variables,
#   3. classifies the sites into types based on S with the drawn family,
#   4. makes the types more or less distinct by multiplying the between-type
#      component of the variables in S by gamma (shift_types()),
#   5. screens each manipulated environment for plausibility against the
#      reference catchments of the scheme's Environmental Zones
#      (Isolation Forest and Mahalanobis distance),
#   6. retains at most one plausible environment per gamma bin that has not
#      yet reached its quota, and
#   7. computes a fuzzy classification (fuzzy c-means) of every retained
#      environment from the same information as the hard classification, and
#   8. predicts N_POSTERIOR_DRAWS presence-absence realisations for every
#      retained environment with the fitted HMSC model.
#
# Terminology (same in manuscript)
#   manipulated environment : one environment per gamma value in an iteration
#   plausible environment   : passed both plausibility checks
#   retained environment    : plausible, chosen for an open gamma bin, simulated
#   successful iteration    : an iteration with >= 1 retained environment
#   classification family   : the method that produced the typology
#
# Units
#   raw units   : the unscaled environment and the reference catchments.
#                 Physical bounds and plausibility checks work here.
#   model units : hmsc_XData = scale(raw) on the model's own sites (see
#                 05_build_hmsc_hpc_models.R). Classification, descriptors
#                 and HMSC prediction work here.
#
#   shift_types() is invariant to linear rescaling, so the manipulation is done
#   once in raw units (where the bounds are defined) and converted to model
#   units with exactly the scaling used for XData. The environment that passes
#   the plausibility checks is therefore exactly the one HMSC predicts from.
#
# Typology descriptors
#   All descriptors are computed with typology_descriptors() in one reference
#   space, all environmental variables of the model in model units, so they do
#   not depend on the classification family and can be computed identically for
#   real typologies. "_orig" = original environment, "_after" = manipulated.
#
# Outputs
#   --output                             production list (section 6)
#   DIAG_DIR/iterations/<model>.rds      one row per iteration
#   DIAG_DIR/environments/<model>.rds    one row per manipulated environment
#   DIAG_DIR/models/<model>.rds          one-row model summary
#
# Inputs besides the model files
#   data/misc/site_coordinates/<model>.rds : x.coord, y.coord in model site
#       order (needed for the spatially constrained family; see
#       extract_site_coordinates.R)
#
# Author: Jonathan Jupke
################################################################################


# ==============================================================================*
# 1. SETUP
# ==============================================================================*

t_start <- Sys.time()

## 1.1 Command-line arguments ---------------------------------------------------
#   --iter_id   SLURM array task id (only echoed to the log)
#   --input     path to the *unfitted* HMSC model (.rds)
#   --output    path for the production output list (.rds)

args    <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
        i <- match(flag, args)
        if (is.na(i) || i == length(args)) NULL else args[i + 1L]
}
iter_id     <- get_arg("--iter_id")
input_file  <- get_arg("--input")
output_file <- get_arg("--output")
if (is.null(input_file) || is.null(output_file)) stop("Missing --input or --output argument.")

cat("\n========================================\n")
cat("Task:  ", if (is.null(iter_id)) "NA" else iter_id, "\n")
cat("Input: ", input_file, "\n")
cat("Output:", output_file, "\n")
cat("========================================\n\n")

## 1.2 Libraries ----------------------------------------------------------------
# MASS is used only through MASS::lda() in ld1_dominance() and is not attached.
suppressPackageStartupMessages({
        library(data.table)
        library(arrow)
        library(jsonify)
        library(Hmsc)
        library(isotree)
        library(cluster)
        library(vegclust)   # fuzzy c-means
        library(vegan)
})

## 1.3 Parameters ---------------------------------------------------------------

# ---- Simulation design -------------------------------------------------------
N_TARGET_RETAINED <- 100L  # retained environments per model
GAMMA_BINS     <- c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)  # edges of the gamma bins
N_GAMMA_DRAWS  <- 25     # gamma values drawn per open bin and iteration
MIN_TYPE_SIZE  <- 3      # typologies with a smaller type are discarded (k-means, Ward)
                         # or merged (crossed thresholds)
SITES_PER_TYPE <- 15     # k_max = round(n_sites / SITES_PER_TYPE), capped at MAX_TYPES
MAX_TYPES      <- 10
MAX_ITERATIONS <- 3000   # hard cap on iterations
STALL_LIMIT    <- 100    # stop after this many consecutive iterations without a retained environment

# ---- Classification families -------------------------------------------------
# In each iteration one family is drawn with probability proportional to its
# deficit relative to an equal share of N_TARGET_RETAINED (see draw_family()).
# A family is no longer drawn after FAMILY_PATIENCE consecutive unsuccessful
# iterations, so a family that cannot produce plausible environments for this
# model does not stall the run.
FAMILIES         <- c("kmeans", "spatial_ward", "crossed_thresholds")
FAMILY_PATIENCE  <- 50
CROSS_N_VARS     <- 2:3   # variables crossed in a threshold typology (before adding geology)
CROSS_CLASSES    <- 2:3   # classes per crossed variable

# ---- Descriptors and fuzzy classification ------------------------------------
KNN_K <- 3L               # neighbours for the classifiability; <= MIN_TYPE_SIZE
FCM_M      <- 1.5         # fuzziness exponent of the fuzzy c-means
FCM_NSTART <- 10L         # random starts of the fuzzy c-means (k-means uses 100; see notes)

# ---- HMSC --------------------------------------------------------------------
MCMC_SAMPLES      <- 2150
MCMC_THIN         <- 150
MCMC_TRANSIENT    <- 10000
N_POSTERIOR_DRAWS <- 5L    # posterior draws (= realisations) per retained environment
N_POST_TAIL       <- 200L  # draws are taken from the last N_POST_TAIL samples of each chain

# ---- Plausibility checks -----------------------------------------------------
# A site is anomalous if its Isolation Forest score exceeds the REF_QUANTILE
# quantile of the reference catchments' scores, or if its squared Mahalanobis
# distance to the reference catchments exceeds the REF_QUANTILE quantile of a
# chi-square distribution. The tolerated share of anomalous sites is anchored
# to the observed sites: TOL_MULTIPLIER x their own anomalous share, bounded to
# TOL_RANGE. An environment is plausible if it stays below the tolerance in
# both checks.
REF_QUANTILE   <- 0.95
TOL_MULTIPLIER <- 1.5
TOL_RANGE      <- c(0.1, 0.5)
RIDGE_RCOND    <- 1e-10  # ridge the reference covariance if its correlation matrix is this ill-conditioned
RIDGE_LAMBDA   <- 1e-6   # ridge size, relative to each variable's variance (unit-free)

# ---- Variables ---------------------------------------------------------------
# A variable is excluded from the candidate set if the SD across the model's
# sites is below NEAR_CONST_RATIO x its SD across the reference catchments,
# i.e. the sites cover too little of the regional variability to be clustered
# on.
NEAR_CONST_RATIO <- 0.05

# Compositional geology shares are always selected together.
GEO_VARS <- c("area_calcareous", "area_sediment", "area_siliceous")

# Physical bounds in raw units, enforced by clipping after the manipulation.
# Variables that are not listed are unbounded.
# !! VERIFY against the units of your data. The script stops in section 2 if
# !! the observed sites or the reference catchments violate a declared bound,
# !! which catches wrong units (e.g. pH stored as pH x 10, shares in %).
PHYSICAL_LOWER <- c(
        slope = 0, lake_area = 0, soil_oc = 0, soil_pH = 0, inundation_depth = 0,
        flooded_area = 0, area_calcareous = 0, area_siliceous = 0, area_sediment = 0,
        glacial_area = 0, bioclim13 = 0, bioclim14 = 0, Rfactor_max = 0,
        Rfactor_avg = 0, Rfactor_min = 0, roughness = 0, mean_snow_equivalent = 0,
        mean_discharge = 0, max_discharge = 0, min_discharge = 0,
        upstream_catchment_area = 0, saturated_soil_water_content = 0
)
PHYSICAL_UPPER <- c(
        soil_pH = 140, area_calcareous = 1, area_siliceous = 1, area_sediment = 1
)


# ---- MEMs for the end-user predictors ----------------------------------------
# Same procedure as in 05_build_hmsc_hpc_models.R, applied to each realisation.
MEM_AUTOCOR  <- "non-null"  # as in script 05
MEM_ALPHA    <- 0.005       # Holm-adjusted p threshold
MEM_MAX      <- 5L          # at most this many MEMs
MEM_MIN_PREV <- 0.05        # taxa below this prevalence are not used for selection

# ---- Paths -------------------------------------------------------------------
DIAG_DIR   <- "data/misc/simulation_diagnostics/"
HELPER_DIR <- "../pulse/code/functions/"

# Seed for the loading phase. Every iteration re-seeds with set.seed(1000 + q),
# so an iteration's random draws do not depend on how many draws came before.
set.seed(1)

## 1.4 Helper files -------------------------------------------------------------
# shift_types.R    : shift_types(), clip_to_bounds(), achieved_gamma(),
#                    ld1_dominance()
# typology_tools.R : spatial_ward(), crossed_thresholds(), geo_coordinates(),
#                    fuzzy_classification(), typology_descriptors() and its components

source(file.path(HELPER_DIR, "shift_types.R"))
source(file.path(HELPER_DIR, "typology_tools.R"))

## 1.5 Helper functions ---------------------------------------------------------

# Sample one element of x. Unlike sample(x, 1), this does not switch to
# sampling from 1:x when x has length one.
resample1 <- function(x) x[sample.int(length(x), 1L)]

# Tolerated share of anomalous sites, anchored to the observed sites.
anchored_tolerance <- function(share_observed) {
        min(TOL_RANGE[2], max(TOL_RANGE[1], TOL_MULTIPLIER * share_observed))
}

# Stop if real data violate a declared physical bound
check_physical_bounds <- function(M, lower, upper, label) {
        msg <- character(0)
        for (v in intersect(names(lower), colnames(M))) {
                n_bad <- sum(M[, v] < lower[[v]] - 1e-9)
                if (n_bad > 0) msg <- c(msg, sprintf("%s < %g in %d rows (min = %g)",
                                                     v, lower[[v]], n_bad, min(M[, v])))
        }
        for (v in intersect(names(upper), colnames(M))) {
                n_bad <- sum(M[, v] > upper[[v]] + 1e-9)
                if (n_bad > 0) msg <- c(msg, sprintf("%s > %g in %d rows (max = %g)",
                                                     v, upper[[v]], n_bad, max(M[, v])))
        }
        if (length(msg)) {
                stop(sprintf("Declared physical bounds are violated by the %s:\n  %s\nCorrect PHYSICAL_LOWER / PHYSICAL_UPPER in section 1.3.",
                             label, paste(msg, collapse = "\n  ")), call. = FALSE)
        }
        invisible(TRUE)
}

# Draw the classification family of an iteration. Weights are the deficits
# relative to an equal share of N_TARGET_RETAINED, so under-represented
# families are drawn more often (analogous to the variable weights). Families
# with FAMILY_PATIENCE consecutive unsuccessful iterations are skipped; if all
# are skipped, or all have met their share, the eligible families are drawn
# with equal probability.
draw_family <- function() {
        eligible <- family_streak < FAMILY_PATIENCE
        if (!any(eligible)) eligible[] <- TRUE
        w <- pmax(0, N_TARGET_RETAINED / length(FAMILIES) - family_counts) * eligible
        if (sum(w) <= 0) w <- as.numeric(eligible)
        sample(FAMILIES, 1L, prob = w)
}

# Descriptors that depend on the environment; stored with suffix _orig (per
# iteration) and _after (per retained environment). n_types and type_evenness
# depend only on the labels and are stored once.
ENV_DESCRIPTORS <- c("env_r2", "env_ch", "disp_mean", "disp_F", "knn_kappa", "asw")

# select_mems(): MEM selection as in script 05. For every taxon with prevalence
# >= MEM_MIN_PREV, fit a binomial GLM on the typology variables, correlate the
# deviance residuals with every candidate MEM, Holm-adjust, and count p < MEM_ALPHA.
# The MEM_MAX most frequently significant MEMs are kept. If none is significant
# for any taxon, the MEM with the smallest adjusted p-value is forced (MEM1 if no
# taxon was testable), so the fractions can always be computed.
# glm.fit() with a model matrix built once replaces glm(y ~ .), because this runs
# n_taxa times per realisation; the fit and the deviance residuals are identical.
select_mems <- function(Y, E, MEM) {
        n   <- nrow(Y)
        Xg  <- cbind(1, as.matrix(E))
        fam <- binomial()
        counts   <- setNames(numeric(ncol(MEM)), colnames(MEM))
        min_p    <- setNames(rep(1, ncol(MEM)), colnames(MEM))
        n_tested <- 0L
        for (j in seq_len(ncol(Y))) {
                y <- Y[, j]
                if (var(y) == 0 || mean(y) < MEM_MIN_PREV) next
                fit <- tryCatch(suppressWarnings(glm.fit(Xg, y, family = fam)),
                                error = function(e) NULL)
                if (is.null(fit)) next
                mu  <- fit$fitted.values
                res <- sign(y - mu) * sqrt(pmax(fam$dev.resids(y, mu, 1), 0))
                if (sd(res) == 0) next
                # Pearson correlation test, vectorised over MEMs (= cor.test())
                r <- drop(cor(MEM, res))
                p <- 2 * pt(-abs(r * sqrt((n - 2) / pmax(1 - r^2, 1e-12))), df = n - 2)
                p <- p.adjust(p, method = "holm")
                counts   <- counts + (p < MEM_ALPHA)
                min_p    <- pmin(min_p, p)
                n_tested <- n_tested + 1L
        }
        if (any(counts > 0)) {
                top    <- names(sort(counts[counts > 0], decreasing = TRUE))   # as in script 05
                mems   <- head(top, MEM_MAX)
                forced <- FALSE
        } else {
                mems   <- names(which.min(min_p))
                forced <- TRUE
        }
        list(mems = mems, forced = forced, n_taxa_tested = n_tested)
}

# community_predictors(): computed per realisation. MEMs are selected from the
# coordinate-derived candidates MEM with select_mems() on this realisation.
#   vp_status: ok | ok_forced (no significant MEM, one forced) |
#              df_exhausted (no MEM fits next to E) | too_few_sites | error
community_predictors <- function(Y, E, MEM) {
        out <- data.table(n_sites_used = NA_integer_, n_taxa_used = NA_integer_,
                          n_mem_selected = 0L, mem_forced = NA, mems_selected = NA_character_,
                          vp_status = NA_character_,
                          dbrda_r2 = NA_real_, vp_env_pure = NA_real_, vp_shared = NA_real_,
                          vp_space_pure = NA_real_, vp_residual = NA_real_)
        Y  <- (as.matrix(Y) > 0) * 1
        sp <- colSums(Y) > 0
        si <- rowSums(Y[, sp, drop = FALSE]) > 0
        Y  <- Y[si, sp, drop = FALSE]
        out$n_sites_used <- nrow(Y)
        out$n_taxa_used  <- ncol(Y)
        if (nrow(Y) <= ncol(E) + 2L || ncol(Y) < 2L) { out$vp_status <- "too_few_sites"; return(out) }

        E   <- as.matrix(E[si, , drop = FALSE])
        MEM <- MEM[si, , drop = FALSE]
        Ed  <- as.data.frame(E)
        d   <- vegdist(Y, method = "jaccard", binary = TRUE)

        # Selection on this realisation; keep enough residual df for [a + b + c]:
        # rank(E) + n_mem <= n - 2
        sel  <- select_mems(Y, E, MEM)
        mems <- head(sel$mems, max(0L, nrow(Y) - 2L - qr(E)$rank))
        out$mem_forced     <- sel$forced
        out$n_mem_selected <- length(mems)
        out$mems_selected  <- paste(mems, collapse = ";")

        res <- tryCatch(
          suppressWarnings({
            if (length(mems) == 0L) {
              list(
                r2     = RsquareAdj(dbrda(
                  d ~ ., data = Ed, sqrt.dist = TRUE
                ))$adj.r.squared,
                ind    = rep(NA_real_, 4),
                status = "df_exhausted"
              )
            } else {
              vp <- varpart(d, Ed, as.data.frame(MEM[, mems, drop = FALSE]), sqrt.dist = TRUE)
              list(
                r2     = vp$part$fract$Adj.R.squared[1],
                # [a + b]
                ind    = vp$part$indfract$Adj.R.squared,
                # [a], [b], [c], [d]
                status = if (sel$forced)
                  "ok_forced"
                else
                  "ok"
              )
            }
          }),
          error = function(e)
            list(error = conditionMessage(e))
        )
        
        if (!is.null(res$error)) {out$vp_status <- paste("error:", res$error); return(out) }
      

        out$vp_status     <- res$status
        out$dbrda_r2      <- res$r2
        out$vp_env_pure   <- res$ind[1]
        out$vp_shared     <- res$ind[2]
        out$vp_space_pure <- res$ind[3]
        out$vp_residual   <- res$ind[4]
        out
}

# ---- Diagnostics -------------------------------------------------------------
# iter_log : one record per iteration (family, typology descriptors, status)
# env_log  : one table per iteration that reached the manipulation step, with
#            one row per manipulated environment
iter_log <- list()
env_log  <- list()

# Family bookkeeping (used by draw_family() and log_iteration()).
family_counts <- setNames(integer(length(FAMILIES)), FAMILIES)   # retained environments
family_streak <- setNames(integer(length(FAMILIES)), FAMILIES)   # consecutive unsuccessful iterations

ITER_STATUSES <- c("FAIL_too_few_types",    # the typology has a single type
                   "FAIL_small_type",       # a type has fewer than MIN_TYPE_SIZE sites
                   "FAIL_isoforest_all",    # no manipulated environment passed the Isolation Forest
                   "FAIL_mahalanobis_all",  # every one that did failed the Mahalanobis check
                   "OK")                    # successful iteration

# Record an iteration and update the family's streak of unsuccessful iterations.
log_iteration <- function(rec, status, env_diag = NULL) {
        rec$status <- status
        if (!is.null(env_diag)) {
                rec$n_gamma          <- nrow(env_diag)
                rec$n_pass_isoforest <- sum(env_diag$passed_isoforest)
                rec$n_plausible      <- sum(env_diag$plausible)
                rec$n_retained       <- sum(env_diag$retained)
                env_log[[length(env_log) + 1L]] <<- env_diag
        }
        family_streak[rec$family] <<- if (status == "OK") 0L else family_streak[rec$family] + 1L
        iter_log[[length(iter_log) + 1L]] <<- as.data.table(rec)
        invisible(NULL)
}

# finalize(): write all three diagnostic files. Called on every exit path, so
# skipped and failed models are also accounted for downstream. Objects that do
# not exist yet (early exits) are reported as NA.
finalize <- function(status) {
        g  <- function(nm, default = NA) get0(nm, envir = globalenv(), inherits = FALSE, ifnotfound = default)
        it <- if (length(iter_log)) rbindlist(iter_log, fill = TRUE) else data.table()
        en <- if (length(env_log))  rbindlist(env_log,  fill = TRUE) else data.table()
        if (nrow(it)) it[, `:=`(model_name = model_name, taxon = taxon)]
        if (nrow(en)) en[, `:=`(model_name = model_name, taxon = taxon)]

        vp <- g("vp_data", NULL)
        vp_fraction <- function(nm) {
                v <- tryCatch(vp$VP[driver == nm, scaled_values], error = function(e) NULL)
                if (is.null(v)) NA_real_ else mean(unlist(v), na.rm = TRUE)
        }
        n_status <- function(s) if (nrow(it)) sum(it$status == s) else 0L
        n_env    <- function(col) if (nrow(en)) sum(en[[col]]) else 0L

        summary_row <- data.table(
                model_name              = model_name,
                taxon                   = taxon,
                run_status              = status,
                stop_reason             = g("stop_reason"),
                model_quality_score     = g("model_quality_score"),
                n_sites                 = g("n_sites"),
                n_taxa                  = g("n_taxa"),
                n_mem                   = g("n_mem"),
                n_candidate_vars        = length(g("all_vars", character(0))),
                near_const_vars         = paste(g("near_const", character(0)), collapse = ";"),
                enz                     = g("enz_string"),
                n_enz                   = g("n_enz"),
                n_reference_catchments  = g("n_reference_catchments"),
                coords_lonlat           = g("coords_lonlat"),
                vp_environmental        = vp_fraction("env"),
                vp_spatial              = vp_fraction("space"),
                vp_codistribution       = vp_fraction("bio"),
                vp_unexplained          = vp_fraction("stochastic"),
                iso_threshold           = unname(g("ISO_THRESHOLD")),
                iso_tol_share           = g("ISO_TOL_SHARE"),
                maha_threshold          = g("MAHA_THRESHOLD"),
                maha_tol_share          = g("MAHA_TOL_SHARE"),
                maha_ridge_used         = g("maha_ridge_used"),
                n_iterations            = nrow(it),
                n_successful_iterations = n_status("OK"),
                n_env_generated         = nrow(en),
                n_env_pass_isoforest    = n_env("passed_isoforest"),
                n_env_plausible         = n_env("plausible"),
                n_env_retained          = n_env("retained"),
                runtime_min             = as.numeric(difftime(Sys.time(), t_start, units = "mins"))
        )
        for (f in FAMILIES) {
                summary_row[, (paste0("n_env_retained_", f)) := family_counts[[f]]]
                summary_row[, (paste0("n_iter_", f))         := if (nrow(it)) sum(it$family == f) else 0L]
        }
        for (s in ITER_STATUSES) summary_row[, (paste0("n_iter_", tolower(s))) := n_status(s)]

        for (d in c("iterations", "environments", "models")) {
                dir.create(file.path(DIAG_DIR, d), recursive = TRUE, showWarnings = FALSE)
        }
        saveRDS(it,          file.path(DIAG_DIR, "iterations",   paste0(model_name, ".rds")))
        saveRDS(en,          file.path(DIAG_DIR, "environments", paste0(model_name, ".rds")))
        saveRDS(summary_row, file.path(DIAG_DIR, "models",       paste0(model_name, ".rds")))
        cat(sprintf("\nDiagnostics saved (%s): %d iterations, %d manipulated environments.\n",
                    status, nrow(it), nrow(en)))
        if (nrow(en)) { cat("Environment funnel:\n"); print(en[, .N, by = status][order(-N)]) }
        invisible(summary_row)
}


# ==============================================================================*
# 2. DATA LOADING
# ==============================================================================*


## 2.1 Model identity and files ---------------------------------------------------
model_name <- sub("\\.rds$", "", basename(input_file))
taxon      <- sub("_[^_]*$", "", model_name)

# Files of this model contain model_name delimited by non-alphanumeric
# characters, so "diatoms_12" does not also match "diatoms_120".
model_pattern    <- paste0("(^|[^[:alnum:]])", model_name, "([^[:alnum:]]|$)")
find_model_files <- function(dir) sort(list.files(dir, pattern = model_pattern, full.names = TRUE))

path_fitted_files <- find_model_files("data/003_fitted_hmsc_models/")
path_eval         <- find_model_files("data/004_model_fit/")
path_vp           <- find_model_files("data/005_variation_partitioning/")
path_unscaled_env <- find_model_files("data/misc/unscaled_environments/")
path_coords       <- find_model_files("data/misc/site_coordinates/")
stopifnot(length(path_eval) == 1L, length(path_fitted_files) == 2L,
          length(path_vp) == 1L, length(path_unscaled_env) == 1L)
if (length(path_coords) != 1L) {
        stop("No site coordinates for ", model_name, " in data/misc/site_coordinates/ ",
             "(needed for spatially constrained Ward; run extract_site_coordinates.R).")
}

## 2.2 Model quality ---------------------------------------------------------------
# Five binary diagnostics (PSRF, AUC, beta distribution, C-score, prevalence);
# at least 4 of 5 required. A missing diagnostic sets the score to 0.

eval_data <- as.data.table(readRDS(path_eval))
stopifnot(nrow(eval_data) == 1L)
model_quality_score <- with(eval_data, psrf_passed + AUC_passed + betaDistr_passed +
                                    c_score_passed + prevalence_passed)
if (is.na(model_quality_score)) model_quality_score <- 0

if (model_quality_score < 4) {
        cat("Skipping: model quality score", model_quality_score, "< 4.\n")
        stop_reason <- "model_quality"
        finalize("SKIPPED_model_quality")
        quit(save = "no", status = 0)
}
cat("Model quality score:", model_quality_score, "/ 5\n")

## 2.3 HMSC models -------------------------------------------------------------------

unfitted_model <- readRDS(input_file)
n_taxa  <- ncol(unfitted_model$Y)
n_sites <- nrow(unfitted_model$Y)
post_list <- lapply(path_fitted_files, function(f) {
        from_json(readRDS(f)[[1]], buffer_size = 524288000)[[1]]
})
cat("HMSC models loaded:", n_sites, "sites,", n_taxa, "taxa.\n")

## 2.4 Sampling scheme ---------------------------------------------------------------
# Each scheme file is filtered before binding, so schemes of other taxa with
# differently typed columns cannot interfere.

scheme_files   <- list.files("data/000_biota", pattern = "03_", full.names = TRUE)
current_scheme <- rbindlist(lapply(scheme_files, function(f) {
        s <- as.data.table(readRDS(f))
        s[as.character(scheme_id) == model_name]
}), fill = TRUE)
stopifnot(nrow(current_scheme) == 1L)

target_zones <- unique(unlist(current_scheme$catchments))
enz_string   <- paste(sort(target_zones), collapse = ";")
n_enz        <- length(target_zones)

## 2.5 Environment in raw and model units -------------------------------------------
# The model-definition script saves the environment after imputation and VIF
# filtering and then builds XData as scale(environment) plus the selected MEMs.
# Here we recover that scaling and verify it, because the conversion
# raw -> model units must reproduce XData exactly.

hmsc_XData <- unfitted_model$XData
mem_cols   <- grep("^MEM", names(hmsc_XData), value = TRUE)
env_cols   <- setdiff(names(hmsc_XData), mem_cols)       # environmental variables, XData order
n_mem      <- length(mem_cols)

unscaled_env <- as.data.table(readRDS(path_unscaled_env))
if (!setequal(names(unscaled_env), env_cols)) {
        stop("Variables of the unscaled environment and of XData differ: ",
             paste(union(setdiff(names(unscaled_env), env_cols),
                         setdiff(env_cols, names(unscaled_env))), collapse = ", "))
}
X_raw_all <- as.matrix(unscaled_env[, env_cols, with = FALSE])   # sites x variables, raw units
if (nrow(X_raw_all) != n_sites) {
        stop(sprintf("Unscaled environment has %d rows, the model %d sites.", nrow(X_raw_all), n_sites))
}
if (current_scheme$samples != n_sites) {
        warning(sprintf("schemes$samples (%d) != number of model sites (%d).", current_scheme$samples, n_sites))
}

SCALE_MU <- colMeans(X_raw_all)
SCALE_SD <- apply(X_raw_all, 2, sd)
to_model_units <- function(M) {
        v <- colnames(M)
        sweep(sweep(M, 2, SCALE_MU[v], "-"), 2, SCALE_SD[v], "/")
}
X_model_all   <- as.matrix(hmsc_XData[, env_cols])           # sites x variables, model units
scaling_error <- max(abs(to_model_units(X_raw_all) - X_model_all))
if (scaling_error > 1e-6) {
        stop(sprintf("XData is not scale(unscaled environment) (max. deviation %.3g).", scaling_error))
}

check_physical_bounds(X_raw_all, PHYSICAL_LOWER, PHYSICAL_UPPER, "observed sites")

## 2.6 Site coordinates ------------------------------------------------------------------
# Used only by the spatially constrained family (hard and fuzzy). Projected
# coordinates are used as they are; longitude/latitude are replaced by a planar
# embedding of the great-circle distances (geo_coordinates()).

site_coords <- as.matrix(as.data.frame(readRDS(path_coords))[, c("x.coord", "y.coord")])
if (nrow(site_coords) != n_sites || anyNA(site_coords)) {
        stop(sprintf("Site coordinates: %d rows (model: %d sites) or missing values.", nrow(site_coords), n_sites))
}
coords_lonlat <- all(abs(site_coords[, 1]) <= 180) && all(abs(site_coords[, 2]) <= 90)
GEO_XY        <- geo_coordinates(site_coords)
cat("Site coordinates:", if (coords_lonlat) "longitude/latitude (planar embedding of great-circle distances)." else "projected.", "\n")

## 2.7 Reference catchments ------------------------------------------------------------
# All catchments in the scheme's Environmental Zones. Only the variables of
# this model are kept, so a missing value in an unused variable does not drop a
# catchment.

eu_hydro_enz <- as.data.table(readRDS("../pulse/data/eu_hydro_dem_w_enz.rds"))
eu_hydro_enz <- eu_hydro_enz[EnZ_name %in% target_zones]

cat("Reading reference catchments ... ")
ref_dt <- rbindlist(lapply(list.files("../pulse/data/catchments/", full.names = TRUE), function(f) {
        joined <- eu_hydro_enz[as.data.table(read_parquet(f)), on = "ID", nomatch = NULL]
        missing_cols <- setdiff(env_cols, names(joined))
        if (length(missing_cols)) stop("Reference catchments lack: ", paste(missing_cols, collapse = ", "))
        joined[, env_cols, with = FALSE]
}), use.names = TRUE)
ref_matrix <- as.matrix(na.omit(ref_dt))
n_reference_catchments <- nrow(ref_matrix)
rm(eu_hydro_enz, ref_dt)
cat(n_reference_catchments, "catchments.\n")

check_physical_bounds(ref_matrix, PHYSICAL_LOWER, PHYSICAL_UPPER, "reference catchments")

## 2.8 Near-constant variables -------------------------------------------------------
near_const <- env_cols[apply(X_raw_all, 2, sd) < NEAR_CONST_RATIO * apply(ref_matrix, 2, sd)]
if (length(near_const)) cat("Near-constant variables excluded from S:", paste(near_const, collapse = ", "), "\n")


# ==============================================================================*
# 3. PLAUSIBILITY BASELINES (raw units, full environment)
# ==============================================================================*

## 3.1 Isolation Forest ---------------------------------------------------------------
# Extended isolation forest (hyperplane splits on all variables) fitted to the
# reference catchments.

isolation_forest <- isolation.forest(
        ref_matrix,
        ntrees                = 1000,
        sample_size           = min(256L, n_reference_catchments),
        ndim                  = ncol(ref_matrix),
        prob_pick_pooled_gain = 0
)
ISO_THRESHOLD <- quantile(predict(isolation_forest, ref_matrix), REF_QUANTILE)
ISO_TOL_SHARE <- anchored_tolerance(mean(predict(isolation_forest, X_raw_all) > ISO_THRESHOLD))

## 3.2 Mahalanobis distance ------------------------------------------------------------
# Squared distance of each site to the mean and covariance of the reference
# catchments. If the covariance is (near-)singular, a ridge proportional to each
# variable's variance is added; this equals a ridge on the correlation matrix
# and is therefore unit-free. The condition check also uses the correlation
# matrix, because the condition number of a covariance matrix depends on units.

REF_MEAN        <- colMeans(ref_matrix)
REF_COV         <- cov(ref_matrix)
maha_ridge_used <- rcond(cov2cor(REF_COV)) < RIDGE_RCOND
if (maha_ridge_used) REF_COV <- REF_COV + diag(RIDGE_LAMBDA * diag(REF_COV), nrow = ncol(REF_COV))
REF_COV_INV     <- solve(REF_COV)
MAHA_THRESHOLD  <- qchisq(REF_QUANTILE, df = ncol(ref_matrix))
MAHA_TOL_SHARE  <- anchored_tolerance(
        mean(mahalanobis(X_raw_all, REF_MEAN, REF_COV_INV, inverted = TRUE) > MAHA_THRESHOLD)
)

cat(sprintf("Isolation Forest: threshold %.3f, tolerated share %.2f\n", ISO_THRESHOLD, ISO_TOL_SHARE))
cat(sprintf("Mahalanobis:      threshold %.2f (df = %d), tolerated share %.2f%s\n",
            MAHA_THRESHOLD, ncol(ref_matrix), MAHA_TOL_SHARE, if (maha_ridge_used) ", ridge used" else ""))


# ==============================================================================*
# 4. HMSC POSTERIOR AND CANDIDATE VARIABLES
# ==============================================================================*

## 4.1 Posterior -----------------------------------------------------------------------
vp_data <- readRDS(path_vp)   # vp_data$importance: per-variable VP contribution

hmsc_model <- importPosteriorFromHPC(
        m         = unfitted_model,
        postList  = post_list,
        nSamples  = MCMC_SAMPLES,
        thin      = MCMC_THIN,
        transient = MCMC_TRANSIENT
)
rm(post_list, unfitted_model)

# Pooled posterior; chains are concatenated in order (chain 1 = 1:n_per_chain, ...).
hmsc_n_samples   <- hmsc_model$samples        # samples per chain
pooled_posterior <- poolMcmcChains(hmsc_model$postList)
n_chains         <- length(hmsc_model$postList)
stopifnot(length(pooled_posterior) == n_chains * hmsc_n_samples)

# Candidate draws: the last N_POST_TAIL samples of every chain.
post_tail_idx <- unlist(lapply(seq_len(n_chains) - 1L, function(ch) {
        ch * hmsc_n_samples + seq(hmsc_n_samples - N_POST_TAIL + 1L, hmsc_n_samples)
}))

## 4.2 Candidate variables -------------------------------------------------------------
# S always leaves at least one environmental variable unmanipulated, so it is
# drawn with size 3 .. P - 1 (crossed thresholds: CROSS_N_VARS, at most P - 1)
# and needs P >= 4 candidates.

all_vars <- sort(setdiff(env_cols, near_const))
geo_vars <- intersect(GEO_VARS, all_vars)
if (length(all_vars) < 4L) {
        cat("Skipping: only", length(all_vars), "candidate variables.\n")
        stop_reason <- "too_few_candidate_variables"
        finalize("SKIPPED_too_few_variables")
        quit(save = "no", status = 0)
}
N_VARS_RANGE     <- 3:(length(all_vars) - 1L)
CROSS_VARS_RANGE <- CROSS_N_VARS[CROSS_N_VARS <= length(all_vars) - 1L]
K_MAX            <- max(2L, min(round(n_sites / SITES_PER_TYPE), MAX_TYPES))

# Down-weighting of variables already used in successful iterations: the
# weight of variable j is 1 / (1 + number of successful iterations including j).
var_selection_counter <- setNames(rep(1, length(all_vars)), all_vars)

## 4.3 Fixed model-unit objects ----------------------------------------------------------
# Prediction template: all XData columns (environment + MEMs) in XData order.
# Per retained environment only the columns in S are replaced.
X_pred_template <- as.matrix(hmsc_XData)
stopifnot(identical(c("(Intercept)", colnames(X_pred_template)), hmsc_model$covNames))

# Reference space of the descriptors for the original environment. The
# distances and nearest neighbours do not depend on the typology, so they are
# computed once and reused for every iteration.
FULL_DMAT_ORIG <- as.matrix(dist(X_model_all))
NN_ORIG        <- knn_index(FULL_DMAT_ORIG, KNN_K)

# Candidate MEMs for the end-user predictors, from the site coordinates as in
# script 05. They depend only on the coordinates; which of them enter the
# variation partitioning is decided per realisation (select_mems()).
mem_matrix <- as.matrix(adespatial::dbmem(site_coords, MEM.autocor = MEM_AUTOCOR, store.listw = FALSE))
cat(sprintf("Candidate variables: %d | k in 2..%d | HMSC MEMs: %d | candidate MEMs: %d",
            length(all_vars), K_MAX, n_mem, ncol(mem_matrix)))

# ==============================================================================*
# 5. SIMULATION LOOP
# ==============================================================================*

# Bin quotas: N_TARGET_RETAINED split as evenly as possible over the bins
# (100 over 6 bins: 17, 17, 17, 17, 16, 16).
n_bins     <- length(GAMMA_BINS) - 1L
BIN_QUOTA  <- N_TARGET_RETAINED %/% n_bins + as.integer(seq_len(n_bins) <= N_TARGET_RETAINED %% n_bins)
bin_labels <- paste0("[", head(GAMMA_BINS, -1), ",", GAMMA_BINS[-1], ")")
bin_counts <- setNames(integer(n_bins), bin_labels)
names(BIN_QUOTA) <- bin_labels
stopifnot(sum(BIN_QUOTA) == N_TARGET_RETAINED)

last_success_iter <- 0L
stop_reason       <- "max_iterations"

# Production storage, one element per retained environment.
realisations  <- list()   # list of N_POSTERIOR_DRAWS presence-absence matrices
env_meta      <- list()   # one-row data.table of descriptors
hard_clusters  <- list()  # type labels
fuzzy_clusters <- list()  # vegclust objects (fuzzy c-means); cluster order arbitrary
comm_pred      <- list()  # data.table, one row per realisation

cat(sprintf("\nTarget: %d retained environments (bins %s), max. %d iterations.\n",
            N_TARGET_RETAINED, paste(BIN_QUOTA, collapse = "/"), MAX_ITERATIONS))

for (q in seq_len(MAX_ITERATIONS)) {

        ## 5.0 Stopping rules and seed -------------------------------------------------
        if (all(bin_counts >= BIN_QUOTA))               { stop_reason <- "target_reached"; break }
        if (q - 1L - last_success_iter >= STALL_LIMIT) { stop_reason <- "stalled";        break }

        set.seed(1000 + q)
        cat(sprintf("\rIteration %d | retained %d/%d | bins %s | families %s", q, sum(bin_counts),
                    N_TARGET_RETAINED, paste(bin_counts, collapse = "/"), paste(family_counts, collapse = "/")))
        flush.console()

        ## 5.1 Classification family and variable subset S --------------------------------
        # Weighted sampling without replacement (successive sampling, R's sample()
        # with prob); the size is uniform on N_VARS_RANGE (crossed thresholds:
        # CROSS_VARS_RANGE). If one geology share is drawn, the others are added so
        # the composition is shifted as a whole.
        family   <- draw_family()
        iter_rec <- list(iteration = q, family = family, alpha = NA_real_)

        size_range    <- if (family == "crossed_thresholds") CROSS_VARS_RANGE else N_VARS_RANGE
        selected_vars <- sample(all_vars, size = resample1(size_range), replace = FALSE,
                                prob = 1 / var_selection_counter)
        if (any(geo_vars %in% selected_vars)) selected_vars <- union(selected_vars, geo_vars)

        iter_rec$n_vars        <- length(selected_vars)
        iter_rec$selected_vars <- paste(sort(selected_vars), collapse = ";")
        iter_rec$vp_importance <- sum(vp_data$importance[selected_vars])

        X_sel_model <- X_model_all[, selected_vars, drop = FALSE]   # S, model units
        X_sel_raw   <- X_raw_all[,   selected_vars, drop = FALSE]   # S, raw units

        ## 5.2 Typology ------------------------------------------------------------------
        # kmeans             : k-means on standardised S, k drawn from 2..K_MAX
        # spatial_ward       : Ward on the combined environmental (S) and geographic
        #                      dissimilarity, alpha ~ U(0, 1), k drawn from 2..K_MAX
        # crossed_thresholds : each variable in S split independently by 1-D k-means
        #                      (2-3 classes), classes crossed; small cells merged, so k
        #                      results from the data (at most K_MAX)
        if (family == "kmeans") {
                cl <- kmeans(X_sel_model, centers = resample1(2:K_MAX), nstart = 100, iter.max = 100)$cluster
        } else if (family == "spatial_ward") {
                iter_rec$alpha <- runif(1)
                cl <- spatial_ward(X_sel_model, GEO_XY, iter_rec$alpha, resample1(2:K_MAX))
        } else {
                cr <- crossed_thresholds(X_sel_model, CROSS_CLASSES, MIN_TYPE_SIZE, K_MAX)
                cl <- cr$cl
                iter_rec$n_classes <- paste(cr$n_classes, collapse = ";")
                iter_rec$n_cells   <- cr$n_cells
        }

        type_sizes <- tabulate(cl)
        k          <- length(type_sizes)
        iter_rec$n_types       <- k
        iter_rec$min_type_size <- min(type_sizes)
        if (k < 2L)                         { log_iteration(iter_rec, "FAIL_too_few_types"); next }
        if (min(type_sizes) < MIN_TYPE_SIZE) { log_iteration(iter_rec, "FAIL_small_type");   next }

        ## 5.3 Descriptors of the original typology ---------------------------------------
        # All in the reference space (all environmental variables, model units),
        # except asw_orig_sub (in S), which only documents the manipulation.
        # ld1_dominance and type_evenness do not change with gamma.
        ld     <- ld1_dominance(X_sel_model, cl)
        d_orig <- typology_descriptors(X_model_all, cl, Dm = FULL_DMAT_ORIG, nn = NN_ORIG)
        iter_rec$ld1_share       <- ld[["ld1_share"]]
        iter_rec$ld1_dominance   <- ld[["ld1_dominance"]]
        iter_rec$n_ld_axes       <- ld[["n_ld_axes"]]
        iter_rec$asw_orig_sub    <- mean(silhouette(cl, dist(X_sel_model))[, 3])
        iter_rec$type_evenness   <- d_orig[["type_evenness"]]
        for (nm in ENV_DESCRIPTORS) iter_rec[[paste0(nm, "_orig")]] <- d_orig[[nm]]

        ## 5.4 gamma values ------------------------------------------------------------------
        # N_GAMMA_DRAWS values per bin that has not reached its quota. runif() never
        # returns the bin edges, so every value lies strictly inside its bin.
        open_bins    <- which(bin_counts < BIN_QUOTA)
        gamma_values <- unlist(lapply(open_bins, function(b) runif(N_GAMMA_DRAWS, GAMMA_BINS[b], GAMMA_BINS[b + 1])))
        gamma_bin    <- rep(open_bins, each = N_GAMMA_DRAWS)
        n_gamma      <- length(gamma_values)

        ## 5.5 Manipulated environments (raw units) -------------------------------------------
        # For each gamma: shift the types in S, clip to physical bounds, and insert
        # the result into the full raw environment (variables outside S unchanged).
        # achieved_gamma() is evaluated in model units so that all variables carry
        # equal weight; it equals gamma unless values were clipped.
        full_envs_raw  <- vector("list", n_gamma)
        gamma_achieved <- numeric(n_gamma)
        n_clipped      <- integer(n_gamma)
        for (i in seq_len(n_gamma)) {
                clipped <- clip_to_bounds(shift_types(X_sel_raw, cl, gamma_values[i]),
                                          PHYSICAL_LOWER, PHYSICAL_UPPER)
                env     <- X_raw_all
                env[, selected_vars] <- clipped$X
                full_envs_raw[[i]]   <- env
                gamma_achieved[i]    <- achieved_gamma(X_sel_model, to_model_units(clipped$X), cl)
                n_clipped[i]         <- clipped$n_clipped
        }

        ## 5.6 Plausibility (raw units, full environment) ----------------------------------------
        # All manipulated environments are scored in one batch. Column j of the
        # score matrices holds the n_sites values of environment j.
        stacked  <- do.call(rbind, full_envs_raw)
        iso_mat  <- matrix(predict(isolation_forest, stacked), nrow = n_sites)
        maha_mat <- matrix(mahalanobis(stacked, REF_MEAN, REF_COV_INV, inverted = TRUE), nrow = n_sites)
        rm(stacked)

        env_diag <- data.table(
                iteration        = q,
                family           = family,
                gamma_id         = seq_len(n_gamma),
                gamma            = gamma_values,
                gamma_bin        = gamma_bin,
                gamma_achieved   = gamma_achieved,
                share_clipped    = n_clipped / length(X_sel_raw),
                iso_n_exceed     = colSums(iso_mat > ISO_THRESHOLD),
                iso_score_mean   = colMeans(iso_mat),
                iso_score_max    = apply(iso_mat, 2, max),
                maha_n_exceed    = colSums(maha_mat > MAHA_THRESHOLD)
        )
        env_diag[, passed_isoforest   := iso_n_exceed  < ISO_TOL_SHARE  * n_sites]
        env_diag[, passed_mahalanobis := maha_n_exceed < MAHA_TOL_SHARE * n_sites]
        env_diag[, plausible          := passed_isoforest & passed_mahalanobis]
        env_diag[, `:=`(retained = FALSE, production_index = NA_integer_)]

        ## 5.7 Retention -------------------------------------------------------------------------
        # In every open bin, the plausible environment whose gamma is closest to the
        # bin centre is retained. gamma is only drawn for open bins, so an iteration
        # with at least one plausible environment is always successful.
        retained_ids <- integer(0)
        for (b in open_bins) {
                cand <- which(env_diag$plausible & env_diag$gamma_bin == b)
                if (length(cand) == 0) next
                centre       <- mean(GAMMA_BINS[c(b, b + 1)])
                retained_ids <- c(retained_ids, cand[which.min(abs(gamma_values[cand] - centre))])
        }
        # Parentheses make data.table use the local vector, not the column 'retained'.
        env_diag[(retained_ids), retained := TRUE]
        env_diag[, status := fifelse(retained, "RETAINED",
                             fifelse(plausible, "PLAUSIBLE_not_retained",
                             fifelse(passed_isoforest, "FAIL_mahalanobis", "FAIL_isoforest")))]

        if (length(retained_ids) == 0) {
                log_iteration(iter_rec, if (!any(env_diag$passed_isoforest)) "FAIL_isoforest_all"
                                        else "FAIL_mahalanobis_all", env_diag)
                next
        }

        ## 5.8 Successful iteration: counters ------------------------------------------------------
        last_success_iter <- q
        var_selection_counter[selected_vars] <- var_selection_counter[selected_vars] + 1
        bin_counts[gamma_bin[retained_ids]]  <- bin_counts[gamma_bin[retained_ids]] + 1L
        family_counts[family]                <- family_counts[family] + length(retained_ids)

        # The same posterior draws are used for all environments of an iteration.
        chain_indices    <- sort(sample(post_tail_idx, size = N_POSTERIOR_DRAWS, replace = FALSE))
        selected_samples <- pooled_posterior[chain_indices]

        ## 5.9 Retained environments: model units, descriptors, predictions ------------------------
        for (i in retained_ids) {
                prod_idx <- length(realisations) + 1L
                env_diag[i, production_index := prod_idx]

                # Exactly the environment that passed the plausibility checks, in model units.
                sel_model <- to_model_units(full_envs_raw[[i]][, selected_vars, drop = FALSE])
                X_pred    <- X_pred_template
                X_pred[, selected_vars] <- sel_model

                # Descriptors of the manipulated environment, same reference space as _orig.
                d_after <- typology_descriptors(X_pred[, env_cols, drop = FALSE], cl, k_nn = KNN_K)
                asw_sub <- mean(silhouette(cl, dist(sel_model))[, 3])

                # Fuzzy classification of the manipulated environment: fuzzy c-means with
                # k clusters and random starts, given the same information as the hard
                # classification of this family (S; for spatial Ward also the coordinates
                # and alpha) and nothing derived from the hard result.
                fcm <- if (family == "spatial_ward") {
                        fuzzy_classification(sel_model, k, FCM_M, FCM_NSTART, GEO_XY, iter_rec$alpha)
                } else {
                        fuzzy_classification(sel_model, k, FCM_M, FCM_NSTART)
                }

                # HMSC prediction. expected = FALSE (default), so every posterior draw
                # yields one realised presence-absence matrix.
                preds <- predict(object = hmsc_model, post = selected_samples,
                                 X = cbind("(Intercept)" = 1, X_pred))
                names(preds) <- paste0("E", prod_idx, "_P", seq_along(preds))

                # End-user predictors per realisation, on the manipulated variables
                # in S (db-RDA is invariant to their units) and the unchanged MEMs.
                comm_pred[[prod_idx]] <- rbindlist(lapply(seq_along(preds), function(p) {
                        cp <- community_predictors(preds[[p]], sel_model, mem_matrix)
                        cp[, `:=`(production_index = prod_idx,
                                  realization      = names(preds)[p],
                                  posterior_index  = chain_indices[p],
                                  chain            = (chain_indices[p] - 1L) %/% hmsc_n_samples + 1L)]
                        cp
                }))

                realisations[[prod_idx]]   <- preds
                hard_clusters[[prod_idx]]  <- cl
                fuzzy_clusters[[prod_idx]] <- fcm
                env_meta[[prod_idx]] <- as.data.table(c(
                        list(production_index = prod_idx,
                             iteration        = q,
                             family           = family,
                             alpha            = iter_rec$alpha,
                             gamma_id         = i,
                             gamma            = gamma_values[i],
                             gamma_bin        = gamma_bin[i],
                             gamma_achieved   = gamma_achieved[i],
                             share_clipped    = env_diag$share_clipped[i],
                             n_types          = k,
                             n_vars           = iter_rec$n_vars,
                             selected_vars    = iter_rec$selected_vars,
                             vp_importance    = iter_rec$vp_importance,
                             type_evenness    = iter_rec$type_evenness,
                             ld1_share        = iter_rec$ld1_share,
                             ld1_dominance    = iter_rec$ld1_dominance,
                             n_ld_axes        = iter_rec$n_ld_axes,
                             asw_orig_sub     = iter_rec$asw_orig_sub,
                             asw_sub_after    = asw_sub,
                             iso_n_exceed     = env_diag$iso_n_exceed[i],
                             maha_n_exceed    = env_diag$maha_n_exceed[i]),
                        iter_rec[paste0(ENV_DESCRIPTORS, "_orig")],
                        as.list(setNames(d_after[ENV_DESCRIPTORS], paste0(ENV_DESCRIPTORS, "_after")))
                ))
        }

        ## 5.10 Close the iteration ----------------------------------------------------------------
        log_iteration(iter_rec, "OK", env_diag)
}

cat(sprintf("\n\nLoop finished (%s): %d retained environments from %d iterations.\n",
            stop_reason, length(realisations), length(iter_log)))


# ==============================================================================*
# 6. PRODUCTION OUTPUT
# ==============================================================================*
# Every element below is indexed by production_index (one per retained
# environment). realization names inside `data` ("E<index>_P<draw>") link to
# community_predictors$realization.

if (length(realisations) == 0L) {
        finalize("NO_RETAINED_ENV")
        stop("No retained environments.")
}

comm_pred_dt <- rbindlist(comm_pred)
setcolorder(comm_pred_dt, c("production_index", "realization", "posterior_index", "chain"))

out_list <- list(
        data                    = realisations,             # PA matrices per retained environment
        environments            = rbindlist(env_meta),      # descriptors per retained environment
        hard_cluster_assignment  = hard_clusters,
        fuzzy_cluster_assignment = fuzzy_clusters,          # vegclust objects; cluster order arbitrary
        community_predictors     = comm_pred_dt,            # one row per realisation
        n_mem                   = n_mem,
        stop_reason             = stop_reason
)
saveRDS(out_list, output_file)

cat(sprintf("Production results saved to: %s\n", output_file))
cat(sprintf("  Retained: %d | per family: %s | types: %s\n",
            length(realisations),
            paste(names(family_counts), family_counts, sep = " = ", collapse = ", "),
            paste(sort(unique(out_list$environments$n_types)), collapse = ",")))
cat(sprintf("  db-RDA adj. R2 (median over realisations): %.3f | NA: %d of %d\n",
            median(comm_pred_dt$dbrda_r2, na.rm = TRUE),
            sum(is.na(comm_pred_dt$dbrda_r2)), nrow(comm_pred_dt)))
cat("  VP status:\n"); print(comm_pred_dt[, .N, by = vp_status][order(-N)])
finalize("OK")
cat("\nDone.\n")
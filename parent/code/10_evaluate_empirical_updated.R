################################################################################*
# Script Name:        evaluate_empirical_updated.R
# Description:        Evaluate coherence between typologies of the ORIGINAL
#                     (unmanipulated) environment and the OBSERVED communities.
#                     Counterpart of evaluate_typologies.R (array version).
#
# Author:             Jonathan Jupke
# R Version:          R 4.5.1
# Required Packages:  data.table, vegan, parallelDist, vegclust, cluster
#
# Usage:              Rscript evaluate_originals.R --iter_id <ID> --input <file> --output <file>
#                     --input is the SIMULATION OUTPUT of one model, i.e. the
#                     same file evaluate_typologies.R reads, so both scripts
#                     can run from the same array list.
#
# Typologies
#   In the simulation, the hard classification `cl` of an iteration is
#   computed on the unmanipulated variables S before shift_types() is applied.
#   hard_cluster_assignment therefore holds typologies of the original
#   environment. All retained environments of an iteration share the same
#   `cl`, so one typology per successful iteration is evaluated here against
#   the observed Y.
#
#   The stored fuzzy classifications were fitted to the MANIPULATED
#   environment. They are refitted here on the original S with the same
#   settings and information as in the simulation (S; for spatial Ward also
#   the site coordinates and alpha).
#
# Metrics
#   Identical to evaluate_typologies.R: Classification Strength, ANOSIM R,
#   AUCzeta (+ per order), PERMANOVA R2, fuzzy PERMANOVA R2, fuzzy Mantel.
#   Binary Jaccard on presence-absence. Only test statistics are kept, so all
#   permutation counts are 0.
################################################################################*


# 1 SETUP -----

## 1.1 Parse command line arguments ----
args    <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
        i <- match(flag, args)
        if (is.na(i) || i == length(args)) NULL else args[i + 1L]
}
iter_id     <- get_arg("--iter_id")
input_file  <- get_arg("--input")
output_file <- get_arg("--output")
if (is.null(input_file) || is.null(output_file)) stop("ERROR: Missing --input or --output argument.")

cat("\n========================================\n")
cat("Processing Iteration:", iter_id, "\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n")
cat("========================================\n\n")

## 1.2 libraries ----
suppressPackageStartupMessages({
        library(data.table)
        library(vegan)
        library(parallelDist)
        library(vegclust)   
        library(cluster)
})

## 1.3 custom scripts ----
HELPER_DIR <- "../pulse/code/functions/"
source(file.path(HELPER_DIR, "render_table.R"))    
source(file.path(HELPER_DIR, "calculate_auc.R"))
source(file.path(HELPER_DIR, "exact_zeta.R"))
source(file.path(HELPER_DIR, "typology_tools.R"))  

## 1.4 Parameters --------
zeta_max   <- 10
FCM_M      <- 1.5
FCM_NSTART <- 10L
THREADS    <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))


# 2 LOAD DATA -----

## 2.1 Simulation output ----
cat("Loading simulation output... ")
if (!file.exists(input_file)) stop("Simulation output not found: ", input_file)
sim        <- readRDS(input_file)
model_name <- sub("\\.rds$", "", basename(input_file))
cat("✅\n")

cat("Loading spatial scale data... ")
ss_path <- file.path("data/misc/spatial_scale/", basename(input_file))
if (!file.exists(ss_path)) stop("Spatial scale file not found: ", ss_path)
i.ss <- readRDS(ss_path)
cat("✅\n")

# Files of this model contain model_name delimited by non-alphanumeric
# characters (as in the simulation script), so "diatoms_12" does not also
# match "diatoms_120".
model_pattern    <- paste0("(^|[^[:alnum:]])", model_name, "([^[:alnum:]]|$)")
find_model_files <- function(dir) sort(list.files(dir, pattern = model_pattern, full.names = TRUE))

## 2.2 Observed community and original environment ----
cat("Loading unfitted HMSC model... ")
path_model <- find_model_files("data/001_unfitted_hmsc_models/")
if (length(path_model) != 1L) stop("Expected one unfitted model for ", model_name, ", found ", length(path_model))
unfitted_model <- readRDS(path_model)
cat("✅\n")

Y <- as.matrix(unfitted_model$Y)
if (anyNA(Y)) stop("Y contains missing values.")
if (!all(Y %in% c(0, 1))) warning("Y is not presence-absence. Will be binarized.")
Y       <- (Y > 0) * 1
n_sites <- nrow(Y)
n_taxa  <- ncol(Y)

if (n_taxa < 2L) {
        cat("WARNING: fewer than two taxa. Aborting.\n")
        saveRDS(NULL, output_file)
        quit(save = "no", status = 0)
}
n_empty <- sum(rowSums(Y) == 0)
if (n_empty) warning(n_empty, " sites without any taxon; pairs of empty sites have undefined Jaccard distances.")

# Environment in model units, exactly the X_model_all of the simulation
XData       <- unfitted_model$XData
env_cols    <- grep("^MEM", names(XData), value = TRUE, invert = TRUE)
X_model_all <- as.matrix(XData[, env_cols])
rm(unfitted_model, XData)


# 3 TYPOLOGIES OF THE ORIGINAL ENVIRONMENT -----

## 3.1 One typology per successful iteration ----
env_dt <- as.data.table(sim$environments)
if (nrow(env_dt) == 0L) stop("Simulation output contains no retained environments.")
stopifnot(nrow(env_dt) == length(sim$hard_cluster_assignment),
          identical(env_dt$production_index, seq_len(nrow(env_dt))))

# All retained environments of an iteration must share one typology
same_cl <- env_dt[, .(ok = all(vapply(
        production_index,
        function(p) identical(sim$hard_cluster_assignment[[p]],
                              sim$hard_cluster_assignment[[production_index[1]]]),
        logical(1)))), by = iteration]
if (!all(same_cl$ok)) stop("Iterations with differing hard classifications: ",
                           paste(same_cl[ok == FALSE, iteration], collapse = ","))

meta_cols <- c("production_index", "family", "alpha", "n_types", "n_vars", "selected_vars",
               "type_evenness", grep("_orig", names(env_dt), value = TRUE))
typ <- env_dt[, c(.SD[1], list(n_retained = .N)), by = iteration, .SDcols = meta_cols]
setorder(typ, iteration)
typ[, ID := .I]
N <- nrow(typ)

hard_of <- function(t) sim$hard_cluster_assignment[[typ$production_index[t]]]

n_ok <- vapply(seq_len(N), function(t) length(hard_of(t)) == n_sites, logical(1))
if (!all(n_ok)) stop("Type assignments do not match the number of observed sites (iterations ",
                     paste(typ$iteration[!n_ok], collapse = ","), ").")
k_ok <- vapply(seq_len(N), function(t) length(tabulate(hard_of(t))) == typ$n_types[t], logical(1))
stopifnot(all(k_ok))

cat(sprintf("Typologies: %d (from %d retained environments) | families: %s\n",
            N, nrow(env_dt), paste(names(table(typ$family)), table(typ$family), sep = " = ", collapse = ", ")))

## 3.2 Site coordinates (spatially constrained family only) ----
if (any(typ$family == "spatial_ward")) {
        path_coords <- find_model_files("data/misc/site_coordinates/")
        if (length(path_coords) != 1L) stop("No site coordinates for ", model_name)
        site_coords <- as.matrix(as.data.frame(readRDS(path_coords))[, c("x.coord", "y.coord")])
        if (nrow(site_coords) != n_sites || anyNA(site_coords)) stop("Site coordinates do not match the model sites.")
        GEO_XY <- geo_coordinates(site_coords)
}

## 3.3 Fuzzy classification of the original environment ----
# Same call as in section 5.9 of the simulation, but on the unmanipulated S.
# Seeded per iteration for reproducibility; cluster order is arbitrary.
cat("Refitting fuzzy c-means on the original environment... ")
fuzzy_orig <- lapply(seq_len(N), function(t) {
        S <- strsplit(typ$selected_vars[t], ";", fixed = TRUE)[[1]]
        if (!all(S %in% env_cols)) stop("Selected variables missing from XData: ",
                                        paste(setdiff(S, env_cols), collapse = ", "))
        X_S <- X_model_all[, S, drop = FALSE]
        set.seed(1000 + typ$iteration[t])
        if (typ$family[t] == "spatial_ward") {
                fuzzy_classification(X_S, typ$n_types[t], FCM_M, FCM_NSTART, GEO_XY, typ$alpha[t])
        } else {
                fuzzy_classification(X_S, typ$n_types[t], FCM_M, FCM_NSTART)
        }
})
memb_of <- function(t) as.matrix(fuzzy_orig[[t]]$memb)
cat("✅\n")


# 4 COMPUTATIONS -----

out <- list()

## 4.1 Distance matrix ----
# One observed community, so one distance matrix for all typologies.
d <- parallelDist(Y, method = "binary", threads = THREADS)

## 4.2 Classification strength ----
cat("Compute Classification Strength... ")
cs <- vapply(seq_len(N), function(t) {
        unname(unlist(summary(meandist(dist = d, grouping = hard_of(t))))["CS"])
}, numeric(1))
out[[length(out) + 1]] <- render_table(x = cs, variable = "Classification Strength", id = typ$ID)
cat("✅\n")

## 4.3 ANOSIM ----
cat("Compute ANOSIM... ")
outR <- vapply(seq_len(N), function(t) {
        unname(anosim(x = d, grouping = hard_of(t), permutations = 0)$statistic)
}, numeric(1))
out[[length(out) + 1]] <- render_table(x = outR, variable = "ANOSIM R", id = typ$ID)
rm(outR)
cat("✅\n")

## 4.4 AUC Zeta ----
cat("Compute AUC Zeta... ")
# The baseline (all sites) is the same for every typology. Exact zeta of an
# order does not depend on which other orders are requested, so it is
# computed once up to zeta_max and subset to each typology's orders, which
# equals the per-community computation in evaluate_typologies.R.
z_base_all <- exact_zeta(Y, orders = seq_len(min(zeta_max, n_sites)))

zeta_rows <- vector("list", length = N)
for (t in seq_len(N)) {
        g   <- hard_of(t)
        sbt <- split(seq_along(g), g)
        sbt <- sbt[lengths(sbt) >= 2]
        if (length(sbt) == 0) next
        
        zeta_rows[[t]] <- rbindlist(lapply(names(sbt), function(ty) {
                idx <- sbt[[ty]]
                ord <- seq_len(min(length(idx), zeta_max))
                z_w <- exact_zeta(Y[idx, , drop = FALSE], orders = ord)
                data.table(
                        ID         = typ$ID[t],
                        type       = ty,
                        n_sites    = length(idx),
                        k          = length(ord),
                        auc_within = calculate_auc(x = ord, y = z_w),
                        auc_base   = calculate_auc(x = ord, y = z_base_all[ord])
                )
        }))
}
zeta_dt <- rbindlist(zeta_rows)
zeta_dt[, auc_rel := auc_within - auc_base]
zeta_dt[, auc_rel_per_order := auc_rel / (k - 1)]
zeta_sum <- zeta_dt[, .(
        auc_rel           = mean(auc_rel, na.rm = TRUE),
        auc_rel_per_order = mean(auc_rel_per_order, na.rm = TRUE)
), by = "ID"]

out[[length(out) + 1]] <- render_table(zeta_sum$auc_rel,           variable = "AucZeta",           id = zeta_sum$ID)
out[[length(out) + 1]] <- render_table(zeta_sum$auc_rel_per_order, variable = "AucZeta per order", id = zeta_sum$ID)
rm(zeta_dt, zeta_sum, zeta_rows)
cat("✅\n")

## 4.5 PERMANOVA ----
cat("Compute PERMANOVA... ")
outR <- vapply(seq_len(N), function(t) {
        g <- factor(hard_of(t))
        adonis2(d ~ g, permutations = 0)$R2[1]
}, numeric(1))
out[[length(out) + 1]] <- render_table(x = outR, variable = "PERMANOVA R2", id = typ$ID)
rm(outR)
cat("✅\n")

## 4.6 Fuzzy PERMANOVA ----
cat("Compute Fuzzy PERMANOVA... ")
outR <- vapply(seq_len(N), function(t) {
        g <- memb_of(t)
        adonis2(d ~ g, permutations = 0)$R2[1]
}, numeric(1))
out[[length(out) + 1]] <- render_table(x = outR, variable = "PERMANOVA Fuzzy R2", id = typ$ID)
rm(outR)
cat("✅\n")

## 4.7 Fuzzy Mantel ----
cat("Compute Fuzzy Mantel... ")
outF <- vapply(seq_len(N), function(t) {
        hel <- dist(decostand(memb_of(t), "hellinger"))
        mantel(d, hel, permutations = 0)$statistic
}, numeric(1))
out[[length(out) + 1]] <- render_table(x = outF, variable = "Fuzzy Mantel", id = typ$ID)
cat("✅\n")


# 5 SAVE RESULTS -----
# n_retained: number of retained environments that share this typology. Each
# contributes N_POSTERIOR_DRAWS simulated communities, so n_retained x 5 is
# the weight that matches the simulated results per typology.
# Column names follow build_meta() where an equivalent exists (ID, scheme_id
# first). gamma = 1: the original environment is the unmanipulated one.
# env_id / rep_id and the community-predictor columns have no counterpart and
# are filled with NA when binding with rbindlist(fill = TRUE).
i.meta <- copy(typ)
i.meta[, `:=`(scheme_id     = i.ss$scheme_id,
              model_name    = model_name,
              source        = "original",
              realization   = "observed",
              gamma         = 1,
              n_sites       = n_sites,
              n_taxa        = n_taxa,
              n_empty_sites = n_empty)]
setcolorder(i.meta, c("ID", "scheme_id"))

out2 <- rbindlist(out, use.names = TRUE)
out2 <- i.meta[out2, on = "ID"]
saveRDS(out2, file = output_file)

cat("\n========================================\n")
cat("Analysis completed successfully!\n")
cat("Results saved to:", output_file, "\n")
cat("========================================\n\n")
################################################################################
# Script Name:        evaluate_typologies_array.R
# Description:        evaluate coherence between simulated typologies and communities
#                     (Array version - processes one iteration per job)
#
# Author:             Jonathan Jupke 
# Date Created:       2025-09-15
# Last Modified:      2025-11-06
#
# R Version:          R 4.5.1
# Required Packages:  data.table, mvabund, parallel, vegan, parallelDist, 
#                     doParallel, zetadiv, foreach, indicspecies, philentropy, 
#                     labdsv, optparse
#
# Usage:              Rscript evaluate_typologies_array.R --iter_id <ID> --input <file> --output <file>
# Notes:              Designed to run as SLURM array job - one iteration per job
################################################################################


# 1 SETUP -----

## 1.1 Parse command line arguments ----
# Parse command line arguments
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

## 1.2 libraries ----
suppressPackageStartupMessages({
        library(data.table)
        library(mvabund)
        library(parallel)
        library(vegan)
        library(parallelDist)
        library(doParallel)
        library(zetadiv)
        library(foreach)
        library(philentropy)
        library(labdsv)
        library(indicspecies)
        library(robCompositions)
})

## 1.3 custom scripts ----
source("../pulse/code/functions/render_table.R")
source("../pulse/code/functions/prop_sample.R")
source("../pulse/code/functions/calculate_auc.R")
source("../pulse/code/functions/group_sites_by_cluster.R")
#source("../pulse/code/functions/cv_d_squared_V2.R")
#source("code/functions/render_table.R")
#source("code/functions/prop_sample.R")
#source("code/functions/calculate_auc.R")
#source("code/functions/group_sites_by_cluster.R")
#source("code/functions/cv_d_squared_V2.R")

## 1.4 Parameters --------
#source("../pulse/code/helper/run_parameters.R")

# Use all available cores for within-job parallelization
cores_to_use <- 8

#simulated_files <- list.files("data/005_simulated_data/", full.names = T)
ss_files        <- list.files("data/misc/spatial_scale/", full.names = T)
cat("length of ss files is:" ,length(ss_files), "\n")
cat(ss_files[1], "\n")

# 2 LOAD DATA -----
cat("Loading input data...")
# Read your input data
cat("✅\n")
i.file <- readRDS(input_file)
cat("Loading spatial scale data...\n")
i.ss.index <- sub("006_simulated_data", "misc/spatial_scale/", x = input_file)
cat(i.ss.index, "\n")
i.ss.index <- which(ss_files == i.ss.index)
cat("i.ss.index: ",i.ss.index, "\n")
if (length(i.ss.index) == 0) {
        cat("Spatial scale file not found. Exiting.\n")
        quit(save = "no", status = 0) 
}
i.ss <- readRDS(ss_files[i.ss.index])
cat("✅\n")
# 3 PROCESS THIS ITERATION -----
cat("Starting analysis for iteration", iter_id, "...\n")

# list to store all the results of this loop
out <- list()

#- Read the simulated data sets created in script 08_simulate_data.R
#- Each of these files is a list with the following elements:
# data = a predicted community data set,
# number_of_variables = the number of variables the typology is based on
# contraction_points  = the factor by which samples are contracted in environmental space
# contraction_centroids= the factor by which cluster centroids are contracted in environmental space
# variable_importance = cumulative importance of variables on which the typology is based
# asw = average silhouette width of environmental clusters
# hard_cluster_assignment = type assignment of each sample in hard clustering
# fuzzy_cluster_assignment =  type assignment of each sample in fuzzy clustering

# 2. Computations
## 2.1 IO and Distance Matrices ----
i.data <- unlist(i.file$data, recursive = F)
# number of artificial environments
i.N1    <- length(i.file$data)
# Number of simulated communities 
i.N2    <- length(i.data)
# How many communities did each iteration contribute
i.N3    <- i.file$n_valid_per_iteration
i.N3    <- i.N3[which(!is.na(i.N3))]

i.file$quality = i.file$contraction_centroids/i.file$contraction_points


# compute distance matrices
i.distance.matrices <-
        rapply(
                i.data,
                parallelDist,
                method = "binary",
                how = "list",
                threads = 4
        )
## 2.1 CS -----
cat("Start CS Analysis...\n")
i.cs0 <- sapply(
        1:length(i.distance.matrices),
        function(x) {
                y <- meandist(
                        dist     = i.distance.matrices[[x]],
                        grouping = rep(i.file$hard_cluster_assignment, each = 5)[[x]]
                )
                c(unlist(summary(y))[c("W", "B", "CS")])
        }
)
i.cs1 <- render_table(i.cs0[1,], "within-type dissimilarity")
i.cs2 <- render_table(i.cs0[2,], "between-type dissimilarity")
i.cs3 <- render_table(i.cs0[3,], "classification strength")
out[[length(out) + 1]] <- i.cs1
out[[length(out) + 1]] <- i.cs2
out[[length(out) + 1]] <- i.cs3
rm(i.cs0)
rm(i.cs1)
rm(i.cs2)
rm(i.cs3)
## 2.2 ANOSIM -----
cat("Start ANOSIM...\n")
# Preallocate objects
i.out.r.min <- numeric(length = i.N1)
i.out.r.men <- numeric(length = i.N1)
i.out.r.max <- numeric(length = i.N1)
i.out.p.min <- numeric(length = i.N1)
i.out.p.men <- numeric(length = i.N1)
i.out.p.max <- numeric(length = i.N1)

i.combinations <-
        lapply(
                1:i.N1,
                function(x) {
                        i.file$hard_cluster_assignment[[x]] |>
                                unique() |>
                                combn(m = 2) |>
                                t()
                }
        )
# Pre-compute cluster indices
cluster_indices <- list()
# Loop through each unique cluster set
for (ci in 1:i.N1) {
        # Get the cluster assignment for this set
        cluster_assignments <- i.file$hard_cluster_assignment[[ci]]
        # Pre-compute all possible pair indices for this cluster assignment
        # A pair index returns the indices of all samples that are in a given
        # Combination of types.
        pair_indices <- list()
        for (pa in 1:nrow(i.combinations[[ci]])) {
                pair_indices[[pa]] <- which(cluster_assignments %in% i.combinations[[ci]][pa, ])
        }
        cluster_indices[[ci]] <- pair_indices
}
rm(ci)
# Setup parllel backend
cl <- makeCluster(cores_to_use)
registerDoParallel(cl)


results <- foreach(evaluation = 1:i.N2, .combine = "rbind", .packages = c('vegan', 'dplyr')) %dopar% {
        # Get cluster set index
        adj.id <- ceiling(evaluation/5)
        combination_set <- i.combinations[[adj.id]]

        pa_results <-
                t(
                        sapply(
                                1:nrow(combination_set),
                                function(pa) {
                                        pa_id       <- cluster_indices[[adj.id]][[pa]]
                                        dist_matrix <- as.matrix(i.distance.matrices[[evaluation]])
                                        pa_dist      <- as.dist(dist_matrix[pa_id, pa_id])

                                        pa_ano <- anosim(
                                                x = pa_dist,
                                                grouping = i.file$hard_cluster_assignment[[adj.id]][pa_id],
                                                permutations = 999
                                        )

                                        return(c(pa_ano$statistic, pa_ano$signif, evaluation))
                                }
                        )
                )
        return(pa_results)
}
stopCluster(cl)
for (evaluation in 1:i.N2) {
        eval_results <- results[which(results[,3] == evaluation), 1:2]

        # Get min, mean, max for this evaluation
        if (is.null(nrow(eval_results))){
                i.out.r.min[evaluation] <-
                        i.out.r.men[evaluation] <-
                        i.out.r.max[evaluation] <-eval_results[1]
                i.out.p.min[evaluation] <-
                        i.out.p.men[evaluation] <-
                        i.out.p.max[evaluation] <- eval_results[1]
        } else {
                i.out.r.min[evaluation] <- min(eval_results[,1])
                i.out.r.men[evaluation] <- mean(eval_results[,1])
                i.out.r.max[evaluation] <- max(eval_results[,1])
                i.out.p.min[evaluation] <- min(eval_results[,2])
                i.out.p.men[evaluation] <- mean(eval_results[,2])
                i.out.p.max[evaluation] <- max(eval_results[,2])
        }
}
# Move into final format
out[[length(out) + 1]] <- render_table(i.out.r.min, "ANOSIM R min")
out[[length(out) + 1]] <- render_table(i.out.r.men, "ANOSIM R mean")
out[[length(out) + 1]] <- render_table(i.out.r.max, "ANOSIM R max")
out[[length(out) + 1]] <- render_table(i.out.p.min, "ANOSIM p min")
out[[length(out) + 1]] <- render_table(i.out.p.men, "ANOSIM p mean")
out[[length(out) + 1]] <- render_table(i.out.p.max, "ANOSIM p max")

rm(cluster_assignments)
rm(cluster_indices)
rm(results)
rm(eval_results)


## 2.3 AUC Zeta ----
cat("Start AUC Zeta...\n")
# The AUC baseline is the inter-type AUCζ that intra-type values can be compared to
baseline <-numeric(i.N2)
# Constants
zeta_orders <- 1:10
zeta_x      <- 1:10
n_reps      <- 10

# Setup parallel backend
# cl <- makeCluster(cores_to_use)
# registerDoParallel(cl)
# results <- foreach(evaluation = 1:(i.N2), .combine = "rbind", .packages = c('zetadiv', 'dplyr')) %dopar% {

for (evaluation in seq_len(i.N2)) {

        adj.id <- ceiling(evaluation/5)
        # Get cluster assignments
        ev_types <- i.file$hard_cluster_assignment[[adj.id]]

        # Pre-calculate the mean sample size (avoid recalculating)
        type_counts <- table(ev_types)
        ev_N <- round(mean(type_counts))

        # Allocate for inner loop results
        ev_baseline2 <- numeric(n_reps)

        # Get the data this evaluation
        eval_data <- i.data[[evaluation]]

        for (evaluation2 in seq_len(n_reps)) {
                # Sample indices
                ev2_id <- prop_sample(x = ev_types, N = ev_N)

                # Sample data and compute AUCζ
                ev2d <- eval_data[ev2_id, ]
                if (sum(ev2d) == 0) {
                        ev_baseline2[evaluation2] <- 0   
                } else {
                        ev_baseline2[evaluation2] <-
                                ev2d |>
                                Zeta.decline.ex(
                                        orders = zeta_orders,
                                        plot = FALSE,
                                        rescale = TRUE
                                ) |>
                                _$zeta.val |>
                                calculate_auc(x = zeta_x)    
                }
                
        }
        baseline[evaluation] <- mean(ev_baseline2)
}

baseline        <- data.frame(baseline, ID = 1:i.N2)
names(baseline) <- c("baseline", "ID")
# pre-allocate object
cluster_indices <- list()
# apply the custom function group_sites_by_cluster to create a list
# for each type that shows which sites are part of that type
for (evaluation in seq_along(i.file$hard_cluster_assignment)) {
        cluster_indices[[evaluation]] <- group_sites_by_cluster(i.file$hard_cluster_assignment[[evaluation]])
}
# Setup parallel backend
cl <- makeCluster(cores_to_use)
registerDoParallel(cl)
results <- foreach(evaluation = 1:i.N2, .combine = "rbind", 
                   .packages = c('zetadiv', 'dplyr'),
                   .export = c("calculate_auc")) %dopar% {

        cluster_set     <- ceiling(evaluation/5)
        combination_set <- cluster_indices[[cluster_set]]
        if (any(sapply(combination_set, length)==1)){
                dropID <- which(sapply(combination_set, length)==1)
                combination_set <- combination_set[-dropID]
        }
        orders = sapply(combination_set, function(x) min(length(x), 10))
        zetaout <- sapply(
                1:length(combination_set),
                #2,
                function(pa) {
                        pa_id <- combination_set[[pa]]
                        if (sum(i.data[[evaluation]][pa_id,]) == 0){
                                zeta = 0
                                zeta
                        } else {
                                zeta <-
                                        Zeta.decline.ex(
                                                data.spec = i.data[[evaluation]][pa_id,],
                                                orders     = 1:orders[pa],
                                                plot       = FALSE,
                                                rescale    = TRUE)$zeta.val |>
                                        calculate_auc(x = 1:orders[pa])

                                zeta
                        }
                }
        )
        matrix(c(zetaout, rep(evaluation, times = length(zetaout))), ncol = 2)
}
stopCluster(cl)
results2 <- as.data.frame(results)
names(results2)[which(names(results2) == "V1")] <- "auczeta_raw"
names(results2)[which(names(results2) == "V2")] <- "ID"
setDT(results2)
setDT(baseline)
results2 <- baseline[results2, on = "ID"]
results2[, auczeta_relative := auczeta_raw - baseline]
results2[, min :=  min(auczeta_relative), by = "ID"]
results2[, mean := mean(auczeta_relative), by = "ID"]
results2[, max :=  max(auczeta_relative), by = "ID"]
results2 <- unique(results2, by = "ID")
out[[length(out) + 1]] <- render_table(results2$min, "AucZeta min")
out[[length(out) + 1]] <- render_table(results2$mean, "AucZeta mean")
out[[length(out) + 1]] <- render_table(results2$max, "AucZeta max")
rm(baseline)
# ================== #
# 2.4 PERMANOVA ---- #
# ================== #
cat("Start PERMANOVA ...\n")
i.r <- i.p  <- i.F  <-
        vector(mode = "numeric", length = i.N2)
for (evaluation in 1:i.N2){
        evaluation2 <- ceiling(evaluation/5)
        types <- i.file$hard_cluster_assignment[[evaluation2]]
        perma <- adonis2(formula = i.distance.matrices[[evaluation]] ~ types, permutations = 999)
        i.r[evaluation]  <- perma$R2[1]
        i.p[evaluation]  <- perma$`Pr(>F)`[1]
}
rm(perma, evaluation)

# Add to output list
out[[length(out) + 1]] <- render_table(i.r, "PERMANOVA R2")
out[[length(out) + 1]] <- render_table(i.p, "PERMANOVA p")

cat("Start Fuzzy PERMANOVA ...\n")
i.r <- i.p  <- i.F  <-
        vector(mode = "numeric", length = i.N2)
for (evaluation in 1:i.N2){
        evaluation2 <- ceiling(evaluation/5)
        types <- i.file$fuzzy_cluster_assignment[[evaluation2]]$memb
        types <- as.matrix(types)
        perma <- adonis2(formula = i.distance.matrices[[evaluation]] ~ types, permutations = 999)
        i.r[evaluation]  <- perma$R2[1]
        i.p[evaluation]  <- perma$`Pr(>F)`[1]
}
rm(perma, evaluation)

# Add to output list
out[[length(out) + 1]] <- render_table(i.r, "PERMANOVA Fuzzy R2")
out[[length(out) + 1]] <- render_table(i.p, "PERMANOVA Fuzzy p")

## 2.5 Indicator Value ----
# Setup parallel backend
cat("Start Indicator Value ...\n")
cl <- makeCluster(cores_to_use)
registerDoParallel(cl)

# Start parallel loop 
results <- foreach(evaluation = 1:i.N1, .combine = "rbind", .packages = c('indicspecies', 'permute')) %dopar% {
        
        out1 <- c()
        out2 <- c()
        
        isa <- lapply(
                i.file$data[[evaluation]], 
                function(x) 
                        multipatt(x, 
                                i.file$hard_cluster_assignment[[evaluation]],
                                func = "indval",
                                control = permute::how(nperm = 999),
                                duleg = T
                        )
        )
        holm_p <- lapply(
                isa, 
                function(x)
                p.adjust(x$sign$p.value, 
                        method = "holm")       
        )                 

        for (ii in 1:length(holm_p)){
                if (any(is.na(holm_p[[ii]]))) holm_p[[ii]] <- holm_p[[ii]][-which(is.na(holm_p[[ii]]))]
                out1[ii] <-  sum(holm_p[[ii]]<=0.05)/length(holm_p[[ii]])
                out2[ii] <- mean(holm_p[[ii]])
        }       

        cbind(out1, out2)
}

out[[length(out) + 1]] <- render_table(results[,1], "isa_number")
out[[length(out) + 1]] <- render_table(results[,2], "isa_avg_p")
rm(results)

## 2.6 ISAMIC ----
cat("Start ISAMIC...\n")
j.isamic <- vector(mode = "numeric", length = i.N2)
for (evaluation in 1:i.N2){
        isa <- isamic(
                comm = i.data[[evaluation]],
                clustering = i.file$hard_cluster_assignment[[ceiling(evaluation/5)]])
        j.isamic[evaluation] <- mean(isa)
}
rm(isa, evaluation)
out[[length(out) + 1]] <- render_table(j.isamic, "isamic")
rm(j.isamic)

## 2.7 Fuzzy Evaluations ----
cat("Start FUZZY... \n")
i.fuzzy.out1 <-
        i.fuzzy.out2 <- numeric(i.N2)
# set family for MVGLM
i.fam <- "binary"
#fuzzy_assignments <- unlist(i.file$fuzzy_cluster_assignment, recursive = FALSE)
cluster_indices     <- ceiling(1:i.N2 / 5)

for (evaluation in 1:i.N2){
        
        # Create matrix for fuzzy type memberships for this iteration
        current_memb <- as.matrix(i.file$fuzzy_cluster_assignment[[cluster_indices[evaluation]]]$memb)
        # Compute Cross validated Squared Divergence
        # i.fuzzy.out1[evaluation] <- cv_d_squared(
        #         species_data = i.data[[evaluation]],
        #         fuzzy_memberships =  current_memb,
        #         family = i.fam
        # )
        # Compute Mantel test
        i.fuzzy.out2[evaluation] <- mantel(
                xdis = i.distance.matrices[[evaluation]],
                ydis = robCompositions::aDist(current_memb),
                permutations = 1)[["statistic"]]

}
# out[[length(out) + 1]] <- render_table(i.fuzzy.out1, "fuzzy_divergence")
out[[length(out) + 1]] <- render_table(i.fuzzy.out2, "fuzzy_mantel")


# # 4 SAVE RESULTS -----
out2 <- rbindlist(out, fill =T )
saveRDS(out2, file = output_file)

cat("\n========================================\n")
cat("Analysis completed successfully!\n")
cat("Results saved to:", output_file, "\n")
cat("========================================\n\n")
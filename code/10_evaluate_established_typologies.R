################################################################################
# Script Name:        evaluate_originals_array.R
# Description:        Evaluate coherence between simulated typologies and communities 
#                     (Array version - processes one iteration per job)
#
# Author:             Jonathan Jupke 
# Date Created:       2025-09-15
# Last Modified:      2025-11-12 (Array adaptation)
#
# R Version:          R 4.5.1
# Required Packages:  data.table, mvabund, parallel, vegan, parallelDist, 
#                     doParallel, zetadiv, foreach, indicspecies, philentropy, 
#                     labdsv, robCompositions
#
# Usage:              Rscript evaluate_originals_array.R --iter_id <ID> --input <file> --output <file>
# Notes:              Designed to run as SLURM array job
################################################################################

# 1 SETUP -------------------------------------------------------------------

## 1.1 Parse command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

# Default values
iter_id <- NULL
input_file <- NULL
output_file <- NULL

# Manual parsing
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

## 1.2 libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
        library(data.table)
        library(mvabund)
        library(parallel)
        library(vegan)
        library(parallelDist)
        library(doParallel)
        library(zetadiv)
        library(foreach)
        library(indicspecies)
        library(philentropy)
        library(labdsv)
        library(robCompositions)
})

## 1.3 custom scripts ----------------------------------------------------------
# Adjust paths if necessary relative to where the script is executed
source("../pulse/code/functions/prop_sample.R")
source("../pulse/code/functions/calculate_auc.R")
source("../pulse/code/functions/group_sites_by_cluster.R")
source("../pulse/code/functions/balance_clusters.R")
source("../pulse/code/functions/cv_d_squared_V2.R")

### =================== ###
# ---- 1.4 Parameters -----
### =================== ###
# Set the number of cores to use 
# Check if we are running inside a SLURM job

cores_to_use <- 2

# slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK")
# 
# if (slurm_cpus != "") {
#         # If on SLURM, use the specific number of CPUs assigned to this job
#         #cores_to_use <- as.numeric(slurm_cpus)
#         #cores_to_use <- as.numeric(slurm_cpus)
#         cores_to_use <- 5
# } else {
#         # If running locally/interactively, fall back to hardware detection
#         cores_to_use <- parallel::detectCores() - 2
# }
# 
# # Safety check: Ensure we have at least 1 core
# if (is.na(cores_to_use) || cores_to_use < 1) cores_to_use <- 1

cat("Using", cores_to_use, "cores for parallel processing.\n")

cores_to_use <- detectCores() - 2
# Lists for matching auxiliary files
scheme_types    <- list.files("data/misc/scheme_types/", full.names = T)
ss_files        <- list.files("data/misc/spatial_scale/", full.names = T)
set.seed(1)

# 2 LOAD DATA ---------------------------------------------------------
cat("Loading input data...\n")

i.file <- try(readRDS(input_file))

# Error Handling for Data Loading
if (inherits(i.file, "try-error")) {
        cat("ERROR: Couldn't read model file.\n")
        quit(save = "no", status = 1)
}
if (ncol(i.file$Y) == 1) {
        cat("WARNING: Only one taxon present. Aborting this iteration.\n")
        # Save empty result or NULL to indicate skip
        saveRDS(NULL, output_file)
        quit(save = "no", status = 0)
}

# 3 MATCH AUXILIARY FILES --------------------------------------------
# Identify corresponding spatial scale file
# This logic assumes the input_file path mirrors the structure expected (data/001_unfitted_hmsc_models/...)
# If input_file is just a filename, you may need to adjust the pattern matching.
i.ss.index <- sub("001_unfitted_hmsc_models", "misc/spatial_scale", x = input_file)
# If input path is relative/different, try matching by basename
if (!any(ss_files == i.ss.index)) {
        cat("Exact path match failed, trying basename match...\n")
        target_base <- basename(i.ss.index)
        i.ss.index <- ss_files[basename(ss_files) == target_base]
} else {
        i.ss.index <- ss_files[ss_files == i.ss.index]
}

if (length(i.ss.index) == 0) {
        stop("ERROR: could not find spatial scale file.")
}
i.ss <- readRDS(i.ss.index)

# Identify corresponding scheme types file
i.typ.index <- sub("001_unfitted_hmsc_models", "misc/scheme_types", x = input_file)
if (!any(scheme_types == i.typ.index)) {
        target_base <- basename(i.typ.index)
        i.typ.index <- scheme_types[basename(scheme_types) == target_base]
} else {
        i.typ.index <- scheme_types[scheme_types == i.typ.index]
}

if (length(i.typ.index) == 0) {
        stop("ERROR: could not find scheme types file.")
}
i.typ <- readRDS(i.typ.index)

## Consistency Check
if (nrow(i.typ) != nrow(i.file$Y)){
        cat("WARNING: Mismatch in rows between Typology and Y. Aborting.\n")
        saveRDS(NULL, output_file)
        quit(save = "no", status = 0)
}

# 4 PROCESSING -------------------------------------------------------
cat("Processing Model... \n")

i.out <- list()

# Extract biotic communities 
i.Y <- i.file$Y 
i.D <- vegdist(i.Y, method = "jaccard")
# Extract environmental variables 
i.X <- i.file$XData
# remove MEMs
if (any(grepl("MEM", names(i.X)))) i.X <- i.X[, -which(grepl("MEM", names(i.X)))]

# 4.1 K-means clustering --------------------------------------------
cat("Running K-means...\n")
i.kmeans.asw <-  i.smallest.cluster <- i.n.small.clusters <- c()
for (numberOfClusters in 2:10) {
        nc.clust <- kmeans(
                x = i.X,
                centers = numberOfClusters,
                nstart = 100,
                iter.max = 100
        )
        i.smallest.cluster[numberOfClusters - 1] <- min(nc.clust$size)
        i.n.small.clusters[numberOfClusters - 1] <- sum(nc.clust$size < 10)
        if (any(nc.clust$size < 10)){
                rm(list = ls()[grepl("^nc\\.", x = ls())])
                next()
        }
        
        nc.asw <- cluster::silhouette(nc.clust$cluster, dist(i.X))
        i.kmeans.asw[length(i.kmeans.asw) + 1] <- mean(nc.asw[, 3])
        names(i.kmeans.asw)[length(i.kmeans.asw)] <- paste(numberOfClusters)
        rm(list = ls()[grepl("^nc\\.", x = ls())])
}

# if there is no solution where the smallest cluster has at least 10 samples, balance manually
if (is.null(i.kmeans.asw)) {
        names(i.smallest.cluster) <- 2:10
        i.smallest.cluster <- i.smallest.cluster[i.n.small.clusters == 1]
        
        if(length(i.smallest.cluster) == 0) {
                # Fallback if extremely poor fit, pick simplest
                i.nc <- 2
        } else {
                i.nc <- as.numeric(names(i.smallest.cluster[max(which(
                        i.smallest.cluster == max(i.smallest.cluster)
                ))]))
        }
        
        i.clusters <- kmeans(
                x = i.X,
                centers = i.nc,
                nstart = 100,
                iter.max = 100
        )
        i.clusters <- i.clusters$cluster
        i.clusters <- balance_clusters(
                data = i.X,
                clusters = i.clusters,
                min_size = 10
        )
        
} else {
        i.nc <- as.numeric(names(which.max(i.kmeans.asw)))
        i.clusters <- kmeans(
                x = i.X,
                centers = i.nc,
                nstart = 100,
                iter.max = 100
        )
        i.clusters <- i.clusters$cluster
}

## Create Fuzzy Clusters
fuzzyClusters <- vegclust::vegclust(i.X, mobileCenters = i.nc, method = "FCM", m = 1.5)

# 4.2 Prepare Typology List (i.D.list) --------------------------------
cat("Preparing Typology Lists...\n")
i.D.list <- list()
for (ii in 1:ncol(i.typ)){
        
        i.D.list[[ii]] <- list(distance = NA, types = NA, Y = NA)
        ii.typ  <- names(i.typ)[ii]
        ii.typ  <- i.typ[, ..ii.typ]
        ii.typ  <- unlist(ii.typ)
        ii.typ  <- as.character(ii.typ)
        ii.drop <- which(is.na(ii.typ))
        ii.tab  <- table(ii.typ)
        
        # only one type? No comparisons possible
        if (length(ii.tab) == 1) next()
        
        if (any(ii.tab < 5)){
                ii.small.id <- names(which(ii.tab < 5))
                
                if (length(ii.drop) == 0){
                        ii.drop <- which(ii.typ %in% ii.small.id)
                } else {
                        ii.drop <- append(ii.drop, which(ii.typ %in% ii.small.id))        
                }
                
                # is there more than one type left after dropping 
                ii.type.dropped <- ii.typ[-ii.drop]
                ii.tab <- table(ii.type.dropped)
                if (length(ii.tab) == 1) next()
                
                ii.D <- as.matrix(i.D)
                ii.D <- ii.D[-ii.drop, -ii.drop]
                ii.D <- as.dist(ii.D)
                i.D.list[[ii]]$distance <- ii.D
                i.D.list[[ii]]$types <- ii.type.dropped
                i.D.list[[ii]]$Y <- i.Y[-ii.drop, ]
                
        } else {
                if (length(ii.drop) > 0){
                        ii.typ <- ii.typ[-ii.drop]
                        ii.D <- as.matrix(i.D)
                        ii.D <- ii.D[-ii.drop, -ii.drop]
                        ii.D <- as.dist(ii.D)   
                        i.D.list[[ii]]$Y <- i.Y[-ii.drop, ]
                        i.D.list[[ii]]$distance <- ii.D
                } else {
                        i.D.list[[ii]]$Y <- i.Y
                        i.D.list[[ii]]$distance <- i.D
                }
                i.D.list[[ii]]$types <- unlist(ii.typ)
        }
        rm(list = ls()[grepl(x = ls(), pattern = "^ii\\.")])
};rm(ii);gc()

# ADD kmeans to list 
i.D.list[[length(i.D.list) + 1]] <- list(distance = NA, types = NA, Y = NA)
i.D.list[[length(i.D.list)]]$distance <- i.D
i.D.list[[length(i.D.list)]]$types     <- i.clusters
i.D.list[[length(i.D.list)]]$Y         <- i.Y
names(i.D.list) <- c(names(i.typ), "kmeans")

# drop empty list elements
i.id <- sapply(i.D.list, function(x) !all(is.na(x$distance)))
i.D.list <- i.D.list[i.id]

# 5 EVALUATIONS -------------------------------------------------------------

## 5.1 Classification Strength ----------------------------------------------
cat("Evaluating: Classification Strength\n")
i.cs <- sapply(i.D.list, function(x) {
        y <- meandist(dist = x$distance, grouping = x$types)
        y <- unlist(summary(y))["CS"]
        y <- round(y, 3)
        y
})
i.ntypes <- sapply(i.D.list, function(x) { length(unique(x$types)) })
i.types.list <- lapply(i.D.list, function(x) { unique(x$types) })

i.out[[length(i.out)+1]] <- data.table(
        scheme = i.ss$scheme_id, 
        value = i.cs, 
        typology = names(i.D.list),
        n_types = i.ntypes,
        types = i.types.list,
        metric = "cs")

## 5.2 ANOSIM ---------------------------------------------------------------
cat("Evaluating: ANOSIM\n")
i.N = length(i.D.list)
i.out.r.min <- i.out.r.men <- i.out.r.max <- numeric(length = i.N)

i.combinations <- lapply(i.D.list, function(x) {
        x$types |> unique() |> combn(m = 2) |> t()
})  

cluster_indices <- list()
for (ci in 1:length(i.D.list)) {
        cluster_assignments <- i.D.list[[ci]]$types
        pair_indices <- list()
        for (pa in 1:nrow(i.combinations[[ci]])) {
                pair_indices[[pa]] <- which(cluster_assignments %in% i.combinations[[ci]][pa, ])
        }
        cluster_indices[[ci]] <- pair_indices
}
rm(ci)

cl <- makeCluster(cores_to_use)
registerDoParallel(cl)

results <- foreach(evaluation = 1:i.N, .combine = "rbind", .packages = c('vegan', 'dplyr')) %dopar% {
        combination_set <- i.combinations[[evaluation]]
        pa_results <- t(sapply(1:nrow(combination_set), function(pa) {
                pa_id <- cluster_indices[[evaluation]][[pa]]
                dist_matrix <- as.matrix(i.D.list[[evaluation]]$distance)
                pa_dist <- as.dist(dist_matrix[pa_id, pa_id])
                
                pa_ano <- anosim(
                        x = pa_dist,
                        grouping = i.D.list[[evaluation]]$types[pa_id],
                        permutations = 5
                )
                return(c(pa_ano$statistic, evaluation))
        }))
        return(pa_results)
}
stopCluster(cl)
rm(evaluation)

for (evaluation in 1:i.N) {
        eval_results <- results[which(results[,2] == evaluation), 1:2]
        if (is.null(nrow(eval_results))){
                i.out.r.min[evaluation] <- i.out.r.men[evaluation] <- i.out.r.max[evaluation] <- eval_results[1]
        } else {
                i.out.r.min[evaluation] <- min(eval_results[,1])
                i.out.r.men[evaluation] <- mean(eval_results[,1])
                i.out.r.max[evaluation] <- max(eval_results[,1])
        }
}
rm(eval_results, cluster_assignments, cluster_indices, results)

i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.out.r.min, typology = names(i.D.list), n_types = i.ntypes, types = i.types.list, metric = "ANOSIM R min")
i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.out.r.men, typology = names(i.D.list), n_types = i.ntypes, types = i.types.list, metric = "ANOSIM R mean")
i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.out.r.max, typology = names(i.D.list), n_types = i.ntypes, types = i.types.list, metric = "ANOSIM R max")

## 5.3 AUC Zeta ------------------------------------------------------------
cat("Evaluating: AUC Zeta\n")
i.baseline <- vector(mode = "numeric", length = i.N)

for (evaluation in 1:i.N){
        ev.types           <- i.D.list[[evaluation]]$types
        ev.baseline2       <- vector(mode = "numeric", length = 10)
        for (evaluation2 in 1:10){
                ev2.N  <- ev.types |> table() |> mean() |> round()
                ev2.id     <- prop_sample(x = ev.types, N = ev2.N, string = T)
                ev2.baseline.data <- i.D.list[[evaluation]]$Y[ev2.id, ]
                ev.baseline2[evaluation2] <- Zeta.decline.ex(
                        data.spec  = ev2.baseline.data,
                        orders     = 1:min(10, ev2.N),
                        plot       = FALSE,
                        rescale    = TRUE)$zeta.val |> calculate_auc(x = 1:min(10, ev2.N))    
        }
        i.baseline[evaluation] <- mean(ev.baseline2)
}
rm(evaluation)

i.baseline        <- data.frame(i.baseline, 1:i.N)
names(i.baseline) <- c("baseline", "ID")

i.cluster_indices <- list()
for (ii in 1:i.N) {
        i.cluster_indices[[ii]] <- group_sites_by_cluster(i.D.list[[ii]]$types)
}
rm(ii)

i.zeta <- list()
for (ii in 1:i.N){
        ii.out <- list()
        ii.combination_set <- i.cluster_indices[[ii]]  
        if (any(sapply(ii.combination_set, length) == 1)){
                dropID <- which(sapply(ii.combination_set, length) == 1)
                ii.combination_set <- ii.combination_set[-dropID]
        }
        ii.orders <- sapply(ii.combination_set, function(x) min(length(x), 10))
        
        for (iii in 1:length(ii.combination_set)){
                pa.id <- ii.combination_set[[iii]]
                pa.zeta  <- Zeta.decline.ex(
                        data.spec = i.D.list[[ii]]$Y[pa.id,],
                        orders     = 1:ii.orders[iii],
                        plot       = FALSE,
                        rescale    = TRUE)
                pa.zeta <- pa.zeta$zeta.val
                pa.zeta <- calculate_auc(pa.zeta, x = 1:ii.orders[iii])
                pa.zeta <- data.table(raw_zeta = pa.zeta, ID = ii)
                ii.out[[length(ii.out) + 1]] <- pa.zeta
        }
        i.zeta[[length(i.zeta) + 1]] <- rbindlist(ii.out)
}
rm(ii)
i.zeta <- rbindlist(i.zeta)
setDT(i.baseline)
i.zeta2 <- i.baseline[i.zeta, on = "ID"]
i.zeta2[, auczeta_relative := raw_zeta - baseline]
i.zeta2[, min  :=  min(auczeta_relative), by = "ID"]
i.zeta2[, mean := mean(auczeta_relative), by = "ID"]
i.zeta2[, max  :=  max(auczeta_relative), by = "ID"]
i.zeta2 <- unique(i.zeta2, by = "ID")

i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.zeta2$min,  typology = names(i.D.list), n_types = i.ntypes, types = i.types.list, metric = "AUCζ min" )
i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.zeta2$mean, typology = names(i.D.list), n_types = i.ntypes, types = i.types.list, metric = "AUCζ mean")
i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.zeta2$max,  typology = names(i.D.list), n_types = i.ntypes, types = i.types.list, metric = "AUCζ max" )

## 5.4 PERMANOVA ----------------------------------------------------------
cat("Evaluating: PERMANOVA\n")
i.r <- i.F <- vector(mode = "numeric", length = i.N)
for (ii in 1:i.N){
        ii.types <- i.D.list[[ii]]$types	
        ii.perma <- adonis2(formula = i.D.list[[ii]]$distance ~ ii.types, permutations = 1)
        i.r[ii]  <- ii.perma$R2[1]
        i.F[ii]  <- ii.perma$F[1]
}
rm(ii)
i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.r, typology = names(i.D.list), n_types = i.ntypes, types = i.types.list, metric = "PERMANOVA R2")
i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.F, typology = names(i.D.list), n_types = i.ntypes, types = i.types.list, metric = "PERMANOVA F")

## PERMANOVA for fuzzy types 
i.fuzzyNOVA <- adonis2(formula = i.D ~ as.matrix(fuzzyClusters$memb), permutations = 1)
i.fnR <- i.fuzzyNOVA$R2[1]
i.fnF <- i.fuzzyNOVA$F[1]
i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.fnR, typology = "fuzzy C means", n_types = i.ntypes["kmeans"], types = i.types.list["kmeans"], metric = "PERMANOVA R2")
i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.fnF, typology = "fuzzy C means", n_types = i.ntypes["kmeans"], types = i.types.list["kmeans"], metric = "PERMANOVA F")

## 5.5 Indicator Value -----------------------------------------------------
cat("Evaluating: IndVal\n")
cl <- makeCluster(cores_to_use)
registerDoParallel(cl)
results <- foreach(evaluation = 1:i.N, .combine = "rbind", .packages = c('indicspecies', 'permute')) %dopar% {
        isa <- multipatt(i.D.list[[evaluation]]$Y, 
                         i.D.list[[evaluation]]$types,
                         func = "indval",
                         control = permute::how(nperm = 999),
                         duleg = T
        )
        holm_p <- p.adjust(isa$sign$p.value, method = "holm")
        if (any(is.na(holm_p))) holm_p <- holm_p[-which(is.na(holm_p))]
        out1 <-  sum(holm_p<=0.05)/length(holm_p)
        out2 <- mean(holm_p)
        cbind(out1, out2)
}
stopCluster(cl)
i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = results[,1], typology = names(i.D.list), n_types = i.ntypes, types = i.types.list, metric = "isa_number")
i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = results[,2], typology = names(i.D.list), n_types = i.ntypes, types = i.types.list, metric = "isa_avg_p")

## 5.6 ISAMIC --------------------------------------------------------------
cat("Evaluating: ISAMIC\n")
j.isamic <- vector(mode = "numeric", length = i.N)
for (evaluation in 1:i.N){
        isa <- isamic(
                comm = i.D.list[[evaluation]]$Y, 
                clustering = i.D.list[[evaluation]]$types)
        j.isamic[evaluation] <- mean(isa)
}
rm(isa, evaluation)
i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = j.isamic, typology = names(i.D.list), n_types = i.ntypes, types = i.types.list, metric = "isamic")
rm(j.isamic)

## 5.7 Fuzzy tests ----------------------------------------------------------
cat("Evaluating: Fuzzy Tests\n")
fuzzy_assignments <- as.matrix(fuzzyClusters$memb)
i.fuzzy.out1 <- cv_d_squared(
        species_data = i.Y,
        fuzzy_memberships =  fuzzy_assignments,
        family = "binary"
)
i.fuzzy.out2 <- mantel(
        xdis = i.D,
        ydis = robCompositions::aDist(fuzzy_assignments),
        permutations = 1
)[["statistic"]]

i.out[[length(i.out) + 1]] <- data.table(scheme = i.ss$scheme_id, value = i.fuzzy.out1, typology = "fuzzy C means", n_types = i.ntypes["kmeans"], types = i.types.list["kmeans"],metric = "fuzzy_divergence")
i.out[[length(i.out) + 1]] <- data.table(scheme = i.ss$scheme_id, value = i.fuzzy.out2, typology = "fuzzy C means", n_types = i.ntypes["kmeans"], types = i.types.list["kmeans"],metric = "fuzzy_mantel")


# 6 WRAP UP AND SAVE -------------------------------------------------------
cat("Analyses complete. Saving data ... \n")
i.out <- rbindlist(i.out)
saveRDS(i.out, output_file)

cat("\n========================================\n")
cat("Analysis completed successfully!\n")
cat("Results saved to:", output_file, "\n")
cat("========================================\n\n")
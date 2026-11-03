################################################################################*
# Script Name:        evaluate_typologies.R
# Description:        evaluate coherence between simulated typologies and communities
#                     (Array version - processes one iteration per job)
# Usage:              Rscript evaluate_typologies_array.R --iter_id <ID> --input <file> --output <file>
# Notes:              Designed to run as SLURM array job - one iteration per job
#                     One array equals one scheme (e.g., "fish_0001") so ty-
#                     pically corresponds too 100 simulated communities. 
################################################################################*


# 1 SETUP -----

## 1.1 Parse command line arguments ----
# Parse command line arguments
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
cat("Processing Iteration:", iter_id, "\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n")
cat("========================================\n\n")

## 1.2 libraries ----
suppressPackageStartupMessages({
        library(data.table)
        library(vegan)
        library(parallelDist)
        library(robCompositions)
})


## 1.3 custom scripts ----
source("../pulse/code/functions/render_table.R")  # build_meta() and render_table()
source("../pulse/code/functions/calculate_auc.R")
source("../pulse/code/functions/exact_zeta.R")

## 1.4 Parameters --------
zeta_max <- 10

# 2 LOAD DATA -----
cat("Loading input data...")

# Read your input data

i.file <- readRDS(input_file)
cat("✅\n")

cat("Loading spatial scale data...\n")
ss_path <- file.path("data/misc/spatial_scale/", basename(input_file))
if (!file.exists(ss_path)) {
        stop("Spatial scale file not found: ", ss_path) 
}
i.ss <- readRDS(ss_path)
cat("✅\n")



# 3 Indexing -----
cat("Starting analysis for iteration", iter_id, "...\n")

i.data <- unlist(i.file$data, recursive = F)
i.N1   <- length(i.file$data) # number of modified environments
i.N2   <- length(i.data) # Number of simulated communities 

# How many communities did each iteration contribute
env_id <- rep(seq_len(i.N1), lengths(i.file$data))
rep_id <- sequence(lengths(i.file$data))
stopifnot(length(env_id) == i.N2)


# Check clusterings
stopifnot(
        length(i.file$hard_cluster_assignment) == i.N1 & 
        length(i.file$fuzzy_cluster_assignment) == i.N1
)

n_rows_ok <- vapply(
  seq_len(i.N2), 
  function(e){
    nrow(i.data[[e]]) == length(i.file$hard_cluster_assignment[[env_id[[e]]]])
  },
  logical(1)
  )

if (!all(n_rows_ok)){
    stop("Site count / type assignment vectors do not match. Issues in communities: ", paste(which(!n_ok_rows), collapse = ","))
}

# one meta data row per community
i.meta <- build_meta(file = i.file, ss = i.ss, env_id = env_id, rep_id = rep_id)

# helper: type labels / memberships for community e
hard_of <- function(e) i.file$hard_cluster_assignment[[env_id[e]]]
memb_of <- function(e) as.matrix(i.file$fuzzy_cluster_assignment[[env_id[e]]]$memb)
 
out <- list()

# 4 COMPUTATIONS -----
 
## 4.1 Distance matrices ----

# make sure they are all matrices otherwise rapply fails silently. 
stopifnot(all(sapply(i.data, function(x) "matrix" %in% class(x))))

# compute distance matrices
i.distance.matrices <-
        rapply(
                i.data,
                parallelDist,
                method = "binary",
                how = "list",
                threads = 4
        )

## 4.2 Classification strength ----
cat("Compute Classification Strength ...\n")
cs <- sapply(
        seq_len(i.N2),
        function(e) {
            y <- meandist(dist = i.distance.matrices[[e]], grouping = hard_of(e))
            unlist(summary(y))[c("CS")]
        }
)

out[[length(out) + 1]] <- render_table(x = cs, variable ="Classification Strength")

cat("✅\n")

## 4.3 ANOSIM ----

cat("Start ANOSIM...\n")
outR <- vapply(
        seq_len(i.N2),
        function(e) {
                unname(anosim(x = i.distance.matrices[[e]],
                              grouping = hard_of(e),
                              permutations = 0)$statistic)
        },
        numeric(1)
)
out[[length(out) + 1]] <- render_table(outR, "ANOSIM R")
rm(outR)

cat("✅\n")

## 4.3 AUC Zeta ----

cat("Compute AUC Zeta ...\n")


zeta_rows <- vector("list", length = i.N2)

for (e in seq_len(i.N2)){

  pa <- as.matrix(i.data[[e]])
  g  <- hard_of(e)

  sbt <- split(seq_along(g), g)
  sbt <- sbt[lengths(sbt) >= 2]

  if (length(sbt) == 0) next

  k_max <- min(zeta_max, max(lengths(sbt)))

  # Compute once 
  z_base <- exact_zeta(pa, orders = seq_len(k_max))

  zeta_rows[[e]] <- rbindlist(lapply(
    names(sbt), 
    function(t){
      idx <- sbt[[t]]
      ord <- seq_len(min(length(idx), zeta_max))
      z_w <- exact_zeta(pa[idx, ,drop = FALSE], orders = ord)
      data.table(
        ID = e, 
        type = t, 
        n_sites = length(idx),
        k = length(ord),
        auc_within = calculate_auc(x = ord, y = z_w),
        auc_base   = calculate_auc(x = ord, y = z_base[ord])
      )
      }
    )
    )
}

zeta_dt <- rbindlist(zeta_rows)
zeta_dt[, auc_rel := auc_within - auc_base]
zeta_dt[, auc_rel_per_oder := auc_rel / (k-1)]
zeta_sum <- zeta_dt[, .(
  auc_rel = mean(auc_rel, na.rm = TRUE),
  auc_rel_per_order = mean(auc_rel_per_oder, na.rm = TRUE)
  ), by = "ID"]

out[[length(out)+1]] <- render_table(zeta_sum$auc_rel, variable = "AucZeta", id = zeta_sum$ID)
out[[length(out)+1]] <- render_table(zeta_sum$auc_rel_per_order, variable = "AucZeta per order", id = zeta_sum$ID)

rm(zeta_dt, zeta_sum, zeta_rows)

cat("✅\n")

## 4.5 PERMANOVA ----

cat("Compute PERMANOVA ...\n")

outR <- vapply(
  seq_len(i.N2), 
  function(e){
    d <- as.matrix(i.distance.matrices[[e]])
    g <- factor(hard_of(e))
    adonis2(d ~ g, permutation = 0)$R2[1]
  },
  numeric(1)
)

# Add to output list
out[[length(out) + 1]] <- render_table(outR, variable =  "PERMANOVA R2")
rm(outR)


cat("✅\n")

## 4.6 Fuzzy PERMANOVA ----

cat("Compute Fuzzy PERMANOVA ...\n")

outR <- vapply(
  seq_len(i.N2), 
  function(e){
    d <- as.matrix(i.distance.matrices[[e]])
    g <- memb_of(e)
    adonis2(d ~ g, permutation = 0)$R2[1]
    },
  numeric(1)
)

# Add to output list

out[[length(out) + 1]] <- render_table(outR, variable =  "PERMANOVA Fuzzy R2")
rm(outR)
cat("✅\n")


## 2.7 Fuzzy Mantel ----
cat("Compute Fuzzy Mantel ... \n")
outF <- vapply( 
  seq_len(i.N2),
  function(e){
    d_comm <- i.distance.matrices[[e]]
    memb   <- memb_of(e)
    hel <- dist(decostand(memb, "hellinger"))
    out <- mantel(d_comm, hel, permutations = 0)$statistic
    },
  numeric(1)
)

out[[length(out) + 1]] <- render_table(x = outF, variable = "Fuzzy Mantel")
cat("✅\n")

# 5 SAVE RESULTS -----
out2 <- rbindlist(out, use.names = TRUE)
out2 <- i.meta[out2, on = "ID"]          # attach metadata once, by ID
saveRDS(out2, file = output_file)
 
cat("\n========================================\n")
cat("Analysis completed successfully!\n")
cat("Results saved to:", output_file, "\n")
cat("========================================\n\n")

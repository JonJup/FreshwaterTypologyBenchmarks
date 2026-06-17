################################################################################*
# Description: Builds and saves unfitted + initialized HMSC models for a single
#              taxon group and scheme row, designed for HPC array job execution.
#              Reads biological data and scheme definitions, performs data
#              filtering, matrix preparation, MEM selection, and model setup,
#              then serializes both the unfitted model and the HPC-initialized
#              MCMC state to disk.
# Notes:       Intended to be called via Slurm array jobs:
#              Rscript 05_build_hmsc_hpc_models.R <TAXON> <ROW_INDEX>
#              TAXON: one of "diatom", "fish", "macroinvertebrate", "macrophyte"
#              ROW_INDEX: integer row of the scheme table to process
################################################################################*

cat("\n========================================\n")
cat("Starting Script Execution (Base R Args)\n")
cat("========================================\n\n")

# 1. Load Packages ------------------------------------------------
suppressPackageStartupMessages({
        library(data.table)   
        library(lubridate)    
        library(sf)           
        library(Hmsc)         
        library(adespatial)   
        library(jsonify)      
        library(mice)         
})

# External functions — halt early if not found rather than failing mid-run
if (file.exists("../R/determine_spatial_scale.R")) {
        source("../R/determine_spatial_scale.R")
} else {
        stop("Custom function 'determine_spatial_scale.R' not found.")
}
if (file.exists("../R/remove_collinearity.R")) {
        source("../R/remove_collinearity.R")
} else {
        stop("Custom function 'remove_collinearity.R' not found.")
}

# 2. Parse Command Line Arguments -------------------------------------------
# Retrieve arguments passed after the script name
args <- commandArgs(trailingOnly = TRUE)

# We expect exactly 2 arguments: [1] Taxon, [2] Row_Index
if (length(args) < 2) {
        stop("Error: Missing arguments. Usage: Rscript script.R <TAXON> <ROW_INDEX>", call. = FALSE)
}

taxon   <- args[1]
row_idx <- as.integer(args[2])

# Validation
if (is.na(row_idx)) {
        stop("Error: Row index (2nd argument) must be an integer.", call. = FALSE)
}

cat("Defining Model For:", taxon, "\n")
cat("Processing Row Index:", row_idx, "\n")

# 3. Helper Functions -----------------------------------------------------
process_scheme <- function(o, b.scheme, b.bio, bio.names, taxon) {

        # -- 3.1 Setup ----
        o.scheme        <- b.scheme[o, ]             # Extract the single scheme row for this array job
        o.scheme.number <- o.scheme$scheme_id        # Unique scheme identifier used in all output filenames
        set.seed(o)                                  # Reproducibility: seed tied to row index

        cat(paste("Processing scheme", o, "of", nrow(b.scheme), "-", bio.names, "\n"))

        # -- Data Selection --
        # Filter biological data to the dataset specified by this scheme
        o.data <- b.bio[data.set == o.scheme$data.set]

        # Apply temporal filter: by month+year if dates are present, by year only otherwise
        focal_months <- as.numeric(o.scheme$focal_months[[1]])
        if (!all(is.na(o.data$eventDate)) && length(focal_months) > 0 && !all(is.na(focal_months))) {
                o.data <- o.data[eventYear == o.scheme$eventYear & month(eventDate) %in% focal_months]
        } else if (!all(is.na(o.data$eventYear))) {
                o.data <- o.data[eventYear == o.scheme$eventYear]
        }

        # Subsample sites if the scheme requests fewer than the full dataset
        if (o.scheme$sample_type != "full") {
                available_ids <- unique(o.data$eventID)
                if (length(available_ids) > o.scheme$samples) {
                        o.sample.ids <- sample(available_ids, size = o.scheme$samples, replace = FALSE)
                        o.data <- o.data[eventID %in% o.sample.ids]
                }
        }

        # ======================================================= *
        # -- 3.2 Taxon Filtering ----
        # ======================================================= *
        # 1. Remove rare taxa: must appear in >= 5 sites AND >= 5% of sites
        unique_sites <- uniqueN(o.data$eventID)
        site_counts  <- o.data[, .(n_sites = uniqueN(eventID)), by = working.taxon]
        valid_taxa   <- site_counts[n_sites >= max(5, 0.05 * unique_sites)]$working.taxon
        o.data       <- o.data[working.taxon %in% valid_taxa]

        # 2. Cap taxa count to sample size to prevent p > n in HMSC
        current_taxa_count <- uniqueN(o.data$working.taxon)
        if (current_taxa_count > o.scheme$samples) {
                # Retain the most widespread (highest site frequency) taxa
                top_taxa <- site_counts[working.taxon %in% valid_taxa][order(-n_sites)][1:o.scheme$samples]$working.taxon
                o.data   <- o.data[working.taxon %in% top_taxa]
        }

        # ======================================================= *
        # -- 3.3 Extract and Save Taxa Count ----
        # ======================================================= *
        n_taxa_final <- uniqueN(o.data$working.taxon)  # Final column count of the species matrix

        taxa_summary <- data.table(
                scheme_id   = paste0(bio.names, "_", o.scheme.number),
                taxon_group = taxon,
                data_set    = o.scheme$data.set,
                n_taxa      = n_taxa_final
        )

        dir.create(paste0(bio.names, "_folder/data/misc/taxa_counts/"), recursive = TRUE, showWarnings = FALSE)
        saveRDS(taxa_summary, paste0("data/misc/taxa_counts/", o.scheme.number, ".rds"))
        # ======================================================= *

        # ======================================================= *
        # -- 3.4 Spatial Scale & Metadata Saving ----
        # ======================================================= *
        # Determine spatial grain of the sampling scheme (e.g. reach / catchment)
        o.sf <- determine_spatial_scale(o.data, o.scheme)
        dir.create("data/misc/spatial_scale/", recursive = TRUE, showWarnings = FALSE)
        saveRDS(o.sf, paste0("data/misc/spatial_scale", o.scheme.number, ".rds"))

        # Record the proportion of records at each taxonomic resolution
        o.tx <- data.table(
                scheme_id    = paste0(bio.names, "_", o.scheme.number),
                taxon        = o.scheme$taxon,
                data.set     = o.scheme$data.set,
                year         = o.scheme$eventYear,
                samples      = o.scheme$samples,
                species_rank = mean(!is.na(o.data$species)),
                genus_rank   = mean(is.na(o.data$species) & !is.na(o.data$genus)),
                family_rank  = mean(is.na(o.data$species) & is.na(o.data$genus) & !is.na(o.data$family)),
                higher_rank  = mean(is.na(o.data$species) & is.na(o.data$genus) & is.na(o.data$family))
        )
        dir.create("data/misc/taxonomic_resolution/", recursive = TRUE, showWarnings = FALSE)
        saveRDS(o.tx, paste0("data/misc/taxonomic_resolution/", o.scheme.number, ".rds"))
        # ======================================================= *

        # ======================================================= *
        # 3.5 Matrix Preparation ----
        # ======================================================= *

        # Pivot to wide format: rows = sites, columns = taxa, values = PA
        o.data2 <- dcast(o.data, eventID ~ working.taxon, value.var = "PA", fun.aggregate = sum, fill = 0)

        # Clamp any values > 1 to 1: aggregation of duplicate records can
        # produce counts > 1, which would break the probit likelihood
        if (any(as.matrix(as.data.frame(o.data2))[, -1] > 1)) {
                for (col in setdiff(names(o.data2), "eventID")) {
                        rows_to_change <- which(o.data2[[col]] > 1)  # Find cells with value > 1 (not just == 2)
                        if (length(rows_to_change) > 0) {
                                set(o.data2, i = rows_to_change, j = col, value = 1L)
                        }
                }
        }

        # Split off the site ID vector and convert the species block to a plain matrix
        site_ids <- o.data2$eventID
        
        # Use column name exclusion instead to drop eventID safely.
        o.data3 <- as.matrix(o.data2[, .SD, .SDcols = setdiff(names(o.data2), "eventID")])

        # -- Environmental Variables --
        # Extract one row per site (unique() collapses repeated taxon rows)
        o.env.samples <- unique(o.data[, .(eventID, x.coord, y.coord, slope, lake_area, soil_oc, soil_pH,
                                           inundation_depth, flooded_area, spi_mean, area_calcareous,
                                           area_siliceous, area_sediment, glacial_area, bioclim05, bioclim06,
                                           bioclim13, bioclim14, Rfactor_max, Rfactor_avg, Rfactor_min,
                                           roughness, elevation, mean_snow_equivalent, mean_discharge,
                                           max_discharge, min_discharge, groundwater_table,
                                           upstream_catchment_area, saturated_soil_water_content, VBM_mean)],
                                by = "eventID")

        # Align env rows to the site order of the species matrix
        o.env.samples <- o.env.samples[match(site_ids, eventID)]
        o.env.samples[, eventID := NULL]  # Drop ID column before numerical processing

        # Drop any constant columns (zero variance → uninformative for modelling)
        o.env.samples <- o.env.samples[, .SD, .SDcols = which(sapply(o.env.samples, uniqueN) > 1)]

        # Pull out coordinates as a separate table before scaling / VIF reduction
        o.coords <- o.env.samples[, .(x.coord, y.coord)]
        o.env.samples[, c("x.coord", "y.coord") := NULL]

        # ======================================================= *
        # -- 3.6 Collinearity & Imputation ----
        # ======================================================= *
        # Impute missing values with predictive mean matching before VIF screening
        if (any(is.na(o.env.samples))) {
                capture.output(o.env.imp <- mice(o.env.samples, method = "pmm", m = 1, maxit = 20, printFlag = FALSE))
                o.env.samples <- complete(o.env.imp)  # Extract the single completed dataset
        }

        # Remove collinear predictors using a VIF threshold of 10
        # NOTE: verify that remove_collinearity.R exports a function named remove_collinearity_vif()
        o.env.samples <- suppressWarnings(remove_collinearity_vif(o.env.samples, threshold = 10))

        # Save the unscaled environment for later back-transformation / reporting
        dir.create("data/misc/unscaled_environments/", recursive = TRUE, showWarnings = FALSE)
        saveRDS(o.env.samples, paste0("data/misc/unscaled_environments/", o.scheme.number, ".rds"))
        
        # Z-standardise all predictors so HMSC priors are on a comparable scale
        o.env.samples <- as.data.table(scale(o.env.samples))

        # ======================================================= *
        # -- 3.7 MEMs (Spatial Eigenvectors) ----
        # # ======================================================= *
        # Compute distance-based Moran's Eigenvector Maps capturing spatial structure
        o.x <- adespatial::dbmem(o.coords, MEM.autocor = "non-null", store.listw = FALSE)

        # Count how often each MEM is significantly correlated with GLM residuals across taxa
        mem_counts <- numeric(ncol(o.x))
        names(mem_counts) <- colnames(o.x)

        for (j in 1:ncol(o.data3)) {
                y <- (o.data3[, j] > 0) * 1L  # Binarise column for robust GLM-based MEM screening

                # Skip degenerate taxa: all-absent, all-present, or very rare (< 5% prevalence)
                if (var(y) == 0 || mean(y) < 0.05) next

                curr_df <- cbind(data.table(y = y), o.env.samples)
                mod <- try(glm(y ~ ., data = curr_df, family = "binomial"), silent = TRUE)
                if (inherits(mod, "try-error")) next

                resids <- residuals(mod)

                # Test each MEM for correlation with env-model residuals (unexplained spatial signal)
                p_vals <- apply(o.x, 2, function(m) cor.test(m, resids)$p.value)
                p_adj  <- p.adjust(p_vals, method = "holm")  # Holm correction for multiple comparisons

                sig_idx <- which(p_adj < 0.005)
                if (length(sig_idx) > 0) {
                        mem_counts[sig_idx] <- mem_counts[sig_idx] + 1  # Tally votes for this MEM
                }
        }

        # Select the top 5 most consistently significant MEMs across all taxa
        top_mems <- names(sort(mem_counts[mem_counts > 0], decreasing = TRUE)[1:5])
        top_mems <- top_mems[!is.na(top_mems)]  # Drop NAs if fewer than 5 MEMs were significant

        # Append selected MEMs to the final predictor set (or use env only if none significant)
        if (length(top_mems) > 0) {
                o.env.final <- cbind(o.env.samples, o.x[, top_mems, drop = FALSE])
        } else {
                o.env.final <- o.env.samples
        }
        # ======================================================= *
        # -- 3.8 HMSC Model Definition ----
        # ======================================================= *
        # 
        # Each site is treated as an independent sample-level random effect
        studyDesign <- data.frame(sample = as.factor(1:nrow(o.data3)))
        rL       <- Hmsc::HmscRandomLevel(units = studyDesign$sample)
        XFormula <- as.formula(paste("~", paste(names(o.env.final), collapse = "+")))

        # Inner helper: construct and save unfitted model, then produce the HPC-engine init object
        fit_and_save_hmsc <- function(Y_matrix, distr_name) {

                mod_obj <- Hmsc(
                        Y           = Y_matrix,
                        XData       = as.data.frame(o.env.final),
                        XFormula    = XFormula,
                        studyDesign = studyDesign,
                        ranLevels   = list("sample" = rL),
                        distr       = distr_name
                )

                # Persist the unfitted model so the HPC fitting step can be re-run independently
                dir.create("data/001_unfitted_hmsc_models/", recursive = TRUE, showWarnings = FALSE)
                saveRDS(mod_obj, paste0("data/001_unfitted_hmsc_models/", o.scheme.number, ".rds"))

                # Run sampleMcmc with engine = "HPC" to generate initialisation state only
                # (actual MCMC sampling is offloaded to the HPC step; this creates the init object)
                init_obj <- sampleMcmc(
                        mod_obj,
                        samples   = 2150,     # Post-thinning samples to retain per chain
                        thin      = 150,      # Keep every 150th sample to reduce autocorrelation
                        transient = 10000,    # Burn-in iterations discarded before sampling
                        nChains   = 2,        # Number of independent MCMC chains
                        verbose   = 1000,     # Print progress every 1000 iterations
                        engine    = "HPC"     # Returns init state for submission to HPC queue
                )

                # Serialise the initialised model as JSON for the HPC fitting pipeline
                json_out <- to_json(init_obj)
                dir.create("data/002_initialized_hmsc_models/", recursive = TRUE, showWarnings = FALSE)
                saveRDS(json_out, paste0("data/002_initialized_hmsc_models/", o.scheme.number, ".rds"))
        }

        # Fit with probit distribution for binary presence/absence data
        fit_and_save_hmsc(o.data3, "probit")

        return(NULL)
}

# 4. Execution -------------------------------------------------------

# Discover available biological data and scheme definition files
bio_file   <- list.files("data/biota/", full.names = TRUE, pattern = "02_")
scheme_file <- list.files("data/biota/", full.names = TRUE, pattern = "03_")

cat("Loading Bio Data:",    bio_file,    "\n")
cat("Loading Scheme Data:", scheme_file, "\n")

b.bio    <- readRDS(bio_file)
b.scheme <- readRDS(scheme_file)
bio.names <- taxon

# Diatom-specific: strip "+" characters appended to some taxon names
if (taxon == "diatom") {
        b.bio[, working.taxon := gsub(x = working.taxon, pattern = "\\+", replacement = "")]
}

# Guard against an out-of-bounds Slurm array index
if (row_idx > nrow(b.scheme)) {
        stop(paste("Row index", row_idx, "exceeds number of schemes", nrow(b.scheme)))
}

# Process the single scheme row assigned to this array job
process_scheme(row_idx, b.scheme, b.bio, bio.names, taxon)

cat("\nDone.\n")

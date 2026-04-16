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

# External functions
if(file.exists("../pulse/code/functions/determine_spatial_scale.R")) {
        source("../pulse/code/functions/determine_spatial_scale.R")
} else {
        stop("Custom function 'determine_spatial_scale.R' not found.")
}
if(file.exists("../pulse/code/functions/remove_collinearity.R")) {
        source("../pulse/code/functions/remove_collinearity.R")
} else {
        stop("Custom function 'remove_collinearity.R' not found.")
}

# 2. Parse Command Line Arguments -------------------------------
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
        
        # -- Setup --
        o.scheme <- b.scheme[o, ]
        o.scheme.number <- o.scheme$scheme_id
        # Reproducibility
        set.seed(o)
        
        cat(paste("Processing scheme", o, "of", nrow(b.scheme), "-", bio.names, "\n"))
        
        # -- Data Selection --
        # Filter by dataset, year, and focal months
        o.data <- b.bio[data.set == o.scheme$data.set]
        
        # If focal_months is not NA/Empty, filter by it
        # If month and year are NA just use all the data
        focal_months <- as.numeric(o.scheme$focal_months[[1]])
        if (!all(is.na(o.data$eventDate)) && length(focal_months) > 0 && !all(is.na(focal_months))) {
                o.data <- o.data[eventYear == o.scheme$eventYear & month(eventDate) %in% focal_months]
        } else if (!all(is.na(o.data$eventYear))){
                o.data <- o.data[eventYear == o.scheme$eventYear]
        } 

        # Subsampling
        if (o.scheme$sample_type != "full") {
                available_ids <- unique(o.data$eventID)
                if (length(available_ids) > o.scheme$samples) {
                        o.sample.ids <- sample(available_ids, size = o.scheme$samples, replace = FALSE)
                        o.data <- o.data[eventID %in% o.sample.ids]
                }
        }
        
        # =======================================================
        # -- Taxon Filtering --
        # =======================================================
        # 1. Remove rare taxa (< 5/5%  sites)
        unique_sites <- uniqueN(o.data$eventID)
        site_counts <- o.data[, .(n_sites = uniqueN(eventID)), by = working.taxon]
        valid_taxa <- site_counts[n_sites >= max(5, 0.05 * unique_sites)]$working.taxon
        o.data <- o.data[working.taxon %in% valid_taxa]
        
        # 2. Cap taxa count to sample size (prevent p > n)
        current_taxa_count <- uniqueN(o.data$working.taxon)
        if (current_taxa_count > o.scheme$samples) {
                # Keep the most abundant/frequent taxa
                top_taxa <- site_counts[working.taxon %in% valid_taxa][order(-n_sites)][1:o.scheme$samples]$working.taxon
                o.data <- o.data[working.taxon %in% top_taxa]
        }
        
        # =======================================================
        # -- Extract and Save Taxa Count --
        # =======================================================
        
        # Calculate final number of taxa (columns in the matrix)
        n_taxa_final <- uniqueN(o.data$working.taxon)
        
        # Create a lightweight summary table
        taxa_summary <- data.table(
                scheme_id   = paste0(bio.names, "_", o.scheme.number),
                taxon_group = taxon,
                data_set    = o.scheme$data.set,
                n_taxa      = n_taxa_final
        )
        
        dir.create("data/misc/taxa_counts/", recursive = TRUE, showWarnings = FALSE)
        saveRDS(taxa_summary, paste0("data/misc/taxa_counts/", o.scheme.number, ".rds"))
        # =======================================================
        
        # =======================================================
        # -- Spatial Scale & Metadata Saving --
        # =======================================================
        o.sf <- determine_spatial_scale(o.data, o.scheme)
        dir.create("data/misc/spatial_scale/", recursive = TRUE, showWarnings = FALSE)
        saveRDS(o.sf, paste0("data/misc/spatial_scale/", o.scheme.number, ".rds"))
        
        # Save Taxonomic Resolution
        o.tx <- data.table(
                scheme_id = paste0(bio.names, "_", o.scheme.number),
                taxon = o.scheme$taxon,
                data.set = o.scheme$data.set,
                year = o.scheme$eventYear,
                samples = o.scheme$samples,
                species_rank = mean(!is.na(o.data$species)),
                genus_rank = mean(is.na(o.data$species) & !is.na(o.data$genus)),
                family_rank = mean(is.na(o.data$species) & is.na(o.data$genus) & !is.na(o.data$family)),
                higher_rank = mean(is.na(o.data$species) & is.na(o.data$genus) & is.na(o.data$family))
        )
        dir.create("data/misc/taxonomic_resolution/", recursive = TRUE, showWarnings = FALSE)
        saveRDS(o.tx, paste0("data/misc/taxonomic_resolution/", o.scheme.number, ".rds"))
        
        # Save Types
        o.types <- unique(o.data[, .(eventID, FEOW, GLORiC, BGR, BRT, EnZ, HER, IFE)], by = "eventID")
        dir.create("data/misc/scheme_types/", recursive = TRUE, showWarnings = FALSE)
        saveRDS(o.types, paste0("data/misc/scheme_types/", o.scheme.number, ".rds"))
        
        # =======================================================
        
        # =======================================================
        # 3. Matrix Preparation --
        # =======================================================
        
        # Reshape data to Wide format
        o.data2 <- dcast(o.data, eventID ~ working.taxon, value.var = "PA", fun.aggregate = sum, fill = 0)
        
        # Remove any >1 numbers which might occur due to aggregation.
        # >1 numbers would lead to errors in the probit model fitting 
        if (any(as.matrix(as.data.frame(o.data2))[,-1] > 1)){
                for (col in setdiff(names(o.data2), "eventID")) {
                        # Find the row indices where the value is 2
                        rows_to_change <- which(o.data2[[col]] > 1)
                        
                        # Update those rows in-place
                        if (length(rows_to_change) > 0) {
                                set(o.data2, i = rows_to_change, j = col, value = 1)
                        }
                }
        }

        # Split ID and Matrix
        site_ids <- o.data2$eventID
        o.data3 <- as.matrix(o.data2[, -1]) # -1 removes eventID assuming it's first
        
        # -- Environmental Variables --
        o.env.samples <- unique(o.data[, .(eventID, x.coord, y.coord, slope, lake_area, soil_oc, soil_pH, 
                                           inundation_depth, flooded_area, spi_mean, area_calcareous, 
                                           area_siliceous, area_sediment, glacial_area, bioclim05, bioclim06, 
                                           bioclim13, bioclim14, Rfactor_max, Rfactor_avg, Rfactor_min, 
                                           roughness, elevation, mean_snow_equivalent, mean_discharge, 
                                           max_discharge, min_discharge, groundwater_table, 
                                           upstream_catchment_area, saturated_soil_water_content, VBM_mean)], 
                                by = "eventID")
        
        # Ensure order matches biological data
        o.env.samples <- o.env.samples[match(site_ids, eventID)]
        o.env.samples[, eventID := NULL] # Remove ID for correlation checks
        
        # Remove constant columns
        o.env.samples <- o.env.samples[, .SD, .SDcols = which(sapply(o.env.samples, uniqueN) > 1)]
        
        # Separate coordinates
        o.coords <- o.env.samples[, .(x.coord, y.coord)]
        o.env.samples[, c("x.coord", "y.coord") := NULL]
        
        # =======================================================
        # -- 3.1 Collinearity & Imputation --
        # =======================================================
        if (any(is.na(o.env.samples))) {
                capture.output(o.env.imp <- mice(o.env.samples, method = "pmm", m = 1, maxit = 20, printFlag = FALSE))
                o.env.samples <- complete(o.env.imp)
        }
        
        o.env.samples <- suppressWarnings(remove_collinearity_vif(o.env.samples, threshold = 10)) 
        
        dir.create("data/misc/unscaled_environments/", recursive = TRUE, showWarnings = FALSE)
        saveRDS(o.env.samples, paste0("data/misc/unscaled_environments/", o.scheme.number, ".rds"))
        
        # -- 3.3 Scaling --
        o.env.samples <- as.data.table(scale(o.env.samples))
        
        # -- 3.4 MEMs (Spatial Eigenvectors) --
        # For MEM selection, we use a simple presence/absence GLM to avoid crashes
        o.x <- adespatial::dbmem(o.coords, MEM.autocor = "non-null", store.listw = FALSE)
        
        mem_counts <- numeric(ncol(o.x))
        names(mem_counts) <- colnames(o.x)
        
        for (j in 1:ncol(o.data3)) {
                # Force binary for MEM selection to be robust
                y <- (o.data3[, j] > 0) * 1 
                
                # skip if all present or all absent or less than 5% present. 
                if (var(y) == 0 || (mean(y > 0) < 0.05)) next
                
                curr_df <- cbind(data.table(y = y), o.env.samples)
                mod <- try(glm(y ~ ., data = curr_df, family = "binomial"), silent = TRUE)
                if (inherits(mod, "try-error")) next
                
                resids <- residuals(mod)
                p_vals <- apply(o.x, 2, function(m) cor.test(m, resids)$p.value)
                p_adj <- p.adjust(p_vals, method = "holm")
                
                sig_idx <- which(p_adj < 0.005)
                if(length(sig_idx) > 0) {
                        mem_counts[sig_idx] <- mem_counts[sig_idx] + 1
                }
        }
        
        top_mems <- names(sort(mem_counts[mem_counts > 0], decreasing = TRUE)[1:5])
        top_mems <- top_mems[!is.na(top_mems)]
        
        if (length(top_mems) > 0) {
                o.env.final <- cbind(o.env.samples, o.x[, top_mems, drop=FALSE])
        } else {
                o.env.final <- o.env.samples
        }
        
        # -- HMSC Model Definition --
        studyDesign <- data.frame(sample = as.factor(1:nrow(o.data3)))
        rL <- HmscRandomLevel(units = studyDesign$sample)
        XFormula <- as.formula(paste("~", paste(names(o.env.final), collapse = "+")))
        
        # Helper to fit and save a model
        fit_and_save_hmsc <- function(Y_matrix, distr_name) {
                
                mod_obj <- Hmsc(
                        Y = Y_matrix,
                        XData = as.data.frame(o.env.final),
                        XFormula = XFormula,
                        studyDesign = studyDesign,
                        ranLevels = list("sample" = rL),
                        distr = distr_name
                )
                
                # Save Unfitted
                dir.create("data/001_unfitted_hmsc_models/", recursive = TRUE, showWarnings = FALSE)
                saveRDS(mod_obj, paste0("data/001_unfitted_hmsc_models/", o.scheme.number, ".rds"))
                
                # Initialize
                init_obj <- sampleMcmc(
                        mod_obj,
                        samples = 2150,
                        thin = 150,
                        transient = 10000,
                        nChains = 2,
                        verbose = 1000,
                        engine = "HPC"
                )
                
                # Save Initialized JSON
                json_out <- to_json(init_obj)
                dir.create("data/002_initialized_hmsc_models/", recursive = TRUE, showWarnings = FALSE)
                saveRDS(json_out, paste0("data/002_initialized_hmsc_models/", o.scheme.number, ".rds"))
        }
        fit_and_save_hmsc(o.data3, "probit")
        return(NULL)
}
# ...

# 4. Execution -------------------------------------------------------

# Load Data
files_all   <- list.files("data/000_biota", pattern = "02_", full.names = TRUE)
schemes_all <- list.files("data/000_biota", pattern = "03_", full.names = TRUE)

# Filter files by the passed 'taxon' variable
bio_file <- grep(taxon, files_all, value = TRUE)
scheme_file <- grep(taxon, schemes_all, value = TRUE)

# Fallback if specific file not found
if(length(bio_file) == 0) bio_file <- files_all[1]
if(length(scheme_file) == 0) scheme_file <- schemes_all[1]

cat("Loading Bio Data:", bio_file, "\n")
cat("Loading Scheme Data:", scheme_file, "\n")

b.bio <- readRDS(bio_file)
b.scheme <- readRDS(scheme_file)
bio.names <- taxon 

# Pre-processing specific to diatoms
if (taxon == "diatoms") {
        b.bio[, working.taxon := gsub(x = working.taxon, pattern = "\\+", replacement = "")]
}

# Validation for row index
if (row_idx > nrow(b.scheme)) {
        stop(paste("Row index", row_idx, "exceeds number of schemes", nrow(b.scheme)))
}

# Run processing for the SINGLE row passed by Slurm
process_scheme(row_idx, b.scheme, b.bio, bio.names, taxon)

cat("\nDone.\n")
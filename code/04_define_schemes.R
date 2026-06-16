################################################################################
 # Script Name: define_schemes.R
 # Description: Each original data set is decomposed into different schemes.
#               Different schemes are created for each year.
#               Within each year only the three consecutive month with the most samples are used.
#               Within these three month further subsets are created based on the number of samples.
#
#               SAMPLING STRATEGY:
#               Uses geometric progession to sample more densely
#               at low N where metrics stabilize rapidly, and more sparsely at high N
#               where changes plateau. Multiplier of 1.5 generates sequence starting at 100:
#               100   150   225   338   506   759  1139  1708  2563  3844  5765
#
#               If abundance data is available, different schemes are run
#               for abundance and for presence/absence.
#               
# Notes: Modified sampling strategy from 100-step to log-scale based on
#        red team analysis identifying non-linear metric stabilization patterns.
################################################################################

# 1. setup ----
## 1.1 wd, packages, scripts ----
setwd(rstudioapi::getActiveProject())
library(data.table)
source("R/find_max_consequtive_sum.R")

## 1.2 load data ----
# List all processed biota files (prefix 02_)
files     <- list.files("data/biota", full.names = T, pattern = "02_")
# Strip path prefix and suffix to get clean taxon names
bio.names <- sapply(files,     function (x) sub("data/biota/02_", "", x))
bio.names <- sapply(bio.names, function (x) sub("_w_environment.rds", "", x))

bio.list       <- lapply(files, readRDS)       # Load all biota datasets into a list
n_bio_datasets <- length(bio.list)             # Number of taxonomic groups

# Lookup table: catchment ID → Environmental Zone name
id.to.enz <- readRDS("data/eu_hydro_dem_w_enz.rds")

# 2. Creating the schemes ----
## 2.1 defining variables before the loop ----
# Number of unique data sets within each taxonomic group
n.data.sets <- lapply(bio.list, function(x) uniqueN(x$data.set))
# Sorted unique data set identifiers per taxonomic group
ud          <- lapply(bio.list, function(x) sort(unique(x$data.set)))

#### 2.1.1 geometric sampling scheme ----
# Rationale: Ecological metrics typically stabilize non-linearly.
# Log-spacing samples more densely where change is rapid (100-300 range)
# and more efficiently where change plateaus (>500 range).
# This provides better characterization of minimum viable N and
# more accurate RF model fitting in critical transition zones.
min_samples <- 100   # Minimum number of samples required to include a year/dataset
multiplier  <- 1.5   # Geometric growth factor between consecutive scheme sizes
max_iter    <- 10    # Number of steps in the geometric sequence

counter <- vector("list", n_bio_datasets)
counter <- lapply(counter, function(x) 0L)  # Initialise all counters to 0

# Generate geometric sequence of candidate sample sizes: min_samples * multiplier^(0..max_iter)
log_sequence <- min_samples * (multiplier ^ (0:max_iter))
# Round to integers and drop any duplicates introduced by rounding
log_sequence <- unique(round(log_sequence))

# Initialise result collector (filled inside the loop)
result_list <- vector(mode = "list")

## 2.2 run Loop ----
for (b in 1:n_bio_datasets) {
        print(b)
        for (d in 1:n.data.sets[[b]]) {
                print(d)
                ds.result.list <- list()  # Collect year-level results for this dataset

                #- Select the focal data set within taxonomic group b
                ds.data.set <- ud[[b]][d]
                ds.bio      <- bio.list[[b]][data.set == ds.data.set]

                #- Coerce eventYear from list to atomic vector if needed

                if (is.list(ds.bio$eventYear)) {
                        ds.bio[, eventYear := unlist(eventYear)]
                }

                # Count unique samples per year (V1 = sample count)
                ds.n.samples <- ds.bio[, uniqueN(eventID), by = "eventYear"]

                #- Skip entire dataset if no year meets the minimum sample threshold
                if (all(ds.n.samples$V1 < min_samples)) {
                        counter[[b]] <- counter[[b]] + nrow(ds.n.samples)  # Track skipped year-slots
                        print(paste("data set", ds.data.set, "contains less than", min_samples, "samples in each year"))
                        rm(list = ls()[grepl(pattern = "^ds\\.", x = ls())])
                        next()
                }

                #- Count and remove years below the sample threshold
                if (any(ds.n.samples$V1 < min_samples)) {
                        count.sub    <- ds.n.samples[V1 < min_samples]
                        counter[[b]] <- counter[[b]] + nrow(count.sub)  # Track sub-threshold years
                }
                ds.n.samples <- ds.n.samples[V1 >= min_samples]  # Keep only years with enough samples

                #- Remove NA year rows (keep NA only if ALL years are NA, i.e. no date info at all)
                if (!all(is.na(ds.n.samples$eventYear)) & any(is.na(ds.n.samples$eventYear))) {
                        ds.n.samples <- ds.n.samples[-which(is.na(ds.n.samples$eventYear))]
                }

                # START LOOP i OVER years in data set d
                for (i in 1:nrow(ds.n.samples)) {

                        if (all(is.na(ds.bio$eventYear))) {
                                #- No year information at all: use entire dataset as one block
                                i.data <- copy(ds.bio)
                                i.focal.months <- NA  

                        } else if (all(is.na(ds.bio$eventDate))) {
                                #- Year known but no date detail: subset by year, skip seasonal filter
                                i.data <- ds.bio[eventYear == ds.n.samples$eventYear[i]]
                                i.focal.months <- NA 

                        } else {
                                #- Full date information available: subset by year then apply seasonal filter
                                i.data <- ds.bio[eventYear == ds.n.samples$eventYear[i]]

                                #--- Seasonal filter: identify the three consecutive months with the
                                #--- most samples to reduce phenological noise in community composition.
                                i.data[, month := month(eventDate)]  # Extract month from date

                                #- Build a frequency table of samples per month (deduplicated by eventID)
                                i.month_table <- unique(i.data, by = "eventID")
                                i.month_table <- table(i.month_table$month)

                                if (length(i.month_table) > 2) {
                                        #- Three or more months present: find the best consecutive triple
                                        i.max_month    <- find_max_consecutive_sum(i.month_table)
                                        i.focal.months <- names(i.max_month$values)
                                        i.data         <- i.data[month %in% i.focal.months]

                                } else {
                                        #- Fewer than three months: handle edge cases manually
                                        i.monthdiff <- diff(as.numeric(names(i.month_table)))

                                        if (length(i.monthdiff) == 0) {
                                                #- Only one month present
                                                if (any(i.month_table > min_samples)) {
                                                        #- Enough samples: keep the single largest month
                                                        i.focal.months <- names(which.max(i.month_table))
                                                        i.data         <- i.data[month == i.focal.months]
                                                } else {
                                                        #- Too few samples even in the single month: skip
                                                        counter[[b]] <- counter[[b]] + 1
                                                        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
                                                        next()
                                                }

                                        } else if (i.monthdiff > 2) {
                                                #- Two months present but non-consecutive (gap > 2 months)
                                                if (any(i.month_table > min_samples)) {
                                                        #- Keep whichever of the two months has more samples
                                                        i.focal.months <- names(which.max(i.month_table))
                                                        i.data         <- i.data[month == i.focal.months]
                                                } else {
                                                        #- Neither month has enough samples: skip
                                                        counter[[b]] <- counter[[b]] + 1
                                                        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
                                                        next()
                                                }

                                        } else {
                                                #- Two consecutive months: keep both
                                                i.focal.months <- names(i.month_table)
                                                i.data         <- i.data[month %in% i.focal.months]
                                        }
                                }
                        }

                        #- Final sample count after seasonal subsetting
                        i.samples <- uniqueN(i.data$eventID)
                        if (i.samples < min_samples) {
                                #- Seasonal filter reduced samples below threshold: skip this year
                                counter[[b]] <- counter[[b]] + 1
                                rm(list = ls()[grepl("^i\\.", x = ls())])
                                next()
                        }

                        #- Build the sample-size scheme for this year
                        # Keep only geometric-sequence steps ≤ actual sample count
                        i.sample_scheme <- log_sequence[log_sequence <= i.samples]

                        # Always append the true full-dataset size as the last scheme step
                        if (max(i.sample_scheme) < i.samples) {
                                i.sample_scheme <- c(i.sample_scheme, i.samples)
                        }

                        # Drop the second-to-last step if it is very close to the full size
                        # (< 25% difference), to avoid near-duplicate schemes
                        if (length(i.sample_scheme) > 1) {
                                i.fullDiff    <- i.samples - i.sample_scheme[length(i.sample_scheme) - 1]
                                i.fullDiffRel <- i.fullDiff / i.samples
                                if (i.fullDiffRel < 0.25) {
                                        i.sample_scheme <- i.sample_scheme[-(length(i.sample_scheme) - 1)]
                                }
                        }

                        # Safety fallback: ensure at least one scheme step exists
                        if (length(i.sample_scheme) == 0) {
                                i.sample_scheme <- i.samples
                        }

                        # Cap at 6 scheme steps to limit downstream computation;
                        # select endpoints plus evenly-spaced intermediates
                        if (length(i.sample_scheme) > 6) {
                                i.keep_indices  <- unique(c(
                                        1,                                                               # Always keep smallest
                                        round(seq(2, length(i.sample_scheme) - 1, length.out = 4)),     # 4 evenly-spaced intermediates
                                        length(i.sample_scheme)                                          # Always keep full sample
                                ))
                                i.sample_scheme <- i.sample_scheme[i.keep_indices]
                        }

                        #- Identify which Environmental Zones (EnZ) are represented in this subset
                        i.enz <- id.to.enz[ID %in% i.data$ID]
                        i.enz <- sort(unique(i.enz$EnZ_name))

                        #- Assemble output row for each scheme step
                        i.out <- data.table(
                                taxon       = bio.names[b],
                                data.set    = ud[[b]][d],
                                eventYear   = unique(i.data$eventYear),
                                catchments  = list(i.enz),
                                samples     = i.sample_scheme,
                                # Last entry is "full"; all preceding are subsamples
                                sample_type = c(rep("sub", length(i.sample_scheme) - 1), "full")
                        )

                        if (all(is.na(ds.bio$eventYear)) || all(is.na(ds.bio$eventDate))) {
                                i.out[, focal_months := list(list(NA))]
                        } else {
                              
                                i.out[, focal_months := list(list(i.focal.months))]
                        }

                        ds.result.list[[length(ds.result.list) + 1]] <- i.out
                        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])

                } # END loop i over years in data.set d

                # Combine year-level results for this dataset and append to global list
                ds.results.list <- rbindlist(ds.result.list)
                result_list[[length(result_list) + 1]] <- ds.results.list
                rm(list = ls()[grepl(pattern = "^ds\\.", x = ls())])

        } # END loop d over data sets in taxonomic group b
} # END loop b over taxonomic groups

# 3. Post-processing ----
# Combine all results into one data.table (fill = TRUE handles any missing columns)
result.data        <- rbindlist(result_list, fill = TRUE)
result.data$taxon  <- factor(result.data$taxon)

# Assign a unique scheme_id per combination of dataset × year × sample size × sample type,
# scoped within each taxon, then format as "<taxon>_<4-digit-number>"
result.data <- split(result.data, f = result.data$taxon)
result.data <- lapply(result.data, function(x) {
        x[, scheme_id := .GRP, by = c("data.set", "eventYear", "samples", "sample_type")]
})
result.data <- rbindlist(result.data)
result.data[, scheme_id := sprintf("%s_%04d", taxon, scheme_id)]

# Split by taxon for per-taxon saving
split.data <- split(result.data, by = "taxon")

# # data output ---------------------------------------------------------------
# Save one .rds file per taxonomic group, matching the 03_ naming convention
lapply(seq_len(n_bio_datasets), function(x)
        saveRDS(
                split.data[[x]],
                paste0("data/biota/03_", bio.names[x], "_scheme.rds")
        )
)

# Print per-taxon counts of skipped year-slots (below min_samples threshold)
counter
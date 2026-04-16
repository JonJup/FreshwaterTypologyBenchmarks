################################################################################
# Script Name:        define_schemes.R
# Description:        Each original data set is decomposed into different schemes. 
#                     Different schemes are created for each year.
#                     Within each year only the three consecutive month with the most samples are used.
#                     Within these three month further subsets are created based on the number of samples.
#                     
#                     SAMPLING STRATEGY:
#                     Uses geometric progession to sample more densely 
#                     at low N where metrics stabilize rapidly, and more sparsely at high N 
#                     where changes plateau. Multiplier of 1.7 generates sequence:
#                     50    85   144   246   418   710  1207  2052  3488  5929 10080

#                     If abundance data is available, different schemes are run 
#                     for abundance and for presence/absence.
# 
# Author:             Jonathan Jupke
# Date Created:       2025-09-18
# Last Modified:      2025-10-14
#
# R Version:          R 4.5.1
# Required Packages:  data.table
#
# Notes:              Modified sampling strategy from 100-step to log-scale based on
#                     red team analysis identifying non-linear metric stabilization patterns.
################################################################################

# 1. setup ----
## 1.1 wd, packages, scripts ----
setwd(rstudioapi::getActiveProject())
library(data.table)
source("R/find_max_consequtive_sum.R")

## 1.2 load data ----
files     <- list.files("data/biota", full.names = T, pattern = "02_")
bio.names <- sapply(files,     function (x) sub("data/biota/02_", "", x)) 
bio.names <- sapply(bio.names, function (x) sub("_w_environment.rds", "", x)) 

bio.list       <- lapply(files, readRDS)
n_bio_datasets <- length(bio.list)

id.to.enz <- readRDS("data/eu_hydro_dem_w_enz.rds")

# 2. Creating the schemes ----
## 2.1 defining variables before the loop ----
n.data.sets      <- lapply(bio.list, function(x) uniqueN(x$data.set))
ud               <- lapply(bio.list, function(x) sort(unique(x$data.set)))



#### 2.1.1 geometric sampling scheme ----
# Rationale: Ecological metrics typically stabilize non-linearly.
# Log-spacing samples more densely where change is rapid (50-200 range)
# and more efficiently where change plateaus (>500 range).
# This provides better characterization of minimum viable N and 
# more accurate RF model fitting in critical transition zones.
# Generate geometric sequence: 50 × 1.7^(0, 1, 2, 3...)
min_samples <- 100
multiplier <- 1.5
max_iter <- 10  # Sufficient to reach very large N
counter = list(0,0,0,0)
# Generate candidate sample sizes
log_sequence <- min_samples * (multiplier ^ (0:(max_iter)))
# Round to whole numbers and remove duplicates after rounding
log_sequence <- unique(round(log_sequence))

# list to store results of loop 
result_list <- vector(mode = "list")

## 2.2 run Loop ----
for (b in 1:n_bio_datasets){
        print(b)
        for (d in 1:n.data.sets[[b]]) {
                if (b == 1 & d == 1) result_list <- vector(mode = "list")
                print(d)
                ds.result.list <- list()
                #- select data set
                ds.data.set <- ud[[b]][d]
                ds.bio      <- bio.list[[b]][data.set == ds.data.set]
                #- how many years?
                if (class(ds.bio$eventYear) == "list") {
                        ds.bio[, eventYear := unlist(eventYear)]
                }
                ds.n.year <- uniqueN(ds.bio$eventYear)
                #- how many samples per year?
                ds.n.samples <- ds.bio[, uniqueN(eventID), by = "eventYear"]
                
                if (all(ds.n.samples$V1 < min_samples)) {
                        counter[[b]] = counter[[b]] + nrow(ds.n.samples)
                        print(paste("data set", ds.data.set, "contains less than", min_samples, "samples in each year"))
                        rm(list = ls()[grepl(pattern = "^ds\\.", x = ls())])
                        next()
                }
                #- select years with more than 100 samples
                if (any(ds.n.samples$V1 < min_samples)){
                        count.sub <- ds.n.samples[V1 < min_samples]
                        counter[[b]] = counter[[b]] + nrow(count.sub)
                }
                ds.n.samples <- ds.n.samples[V1 >= min_samples]
                #- remove NA row, if one of the entries in the year table is NA
                if (!all(is.na(ds.n.samples$eventYear)) & any(is.na(ds.n.samples$eventYear))){
                        ds.n.samples <- ds.n.samples[-which(is.na(ds.n.samples$eventYear))]
                }
                
                # START LOOP i OVER years in data set
                for (i in 1:nrow(ds.n.samples)) {
                        
                        if (all(is.na(ds.bio$eventYear))) {
                                i.data    <- copy(ds.bio)
                        } else if (all(is.na(ds.bio$eventDate))){
                                i.data <- ds.bio[eventYear == ds.n.samples$eventYear[i]]
                        } else {
                                #- Create subset of focal year.
                                i.data <- ds.bio[eventYear == ds.n.samples$eventYear[i]]
                                #--- Check seasons. To prevent strong seasonal changes from influencing
                                #--- the community composition, we identify the three consecutive month
                                #--- with the most samples.
                                #- create month variable
                                i.data[, month := month(eventDate)]
                                #- count samples per month
                                i.month_table <- unique(i.data, by = "eventID")
                                i.month_table <- i.month_table$month
                                i.month_table <- table(i.month_table)
                                #- Which three consecutive months have the most samples?
                                #- This function is loaded in the beginning of the script
                                if (length(i.month_table) > 2) {
                                        i.max_month    <- find_max_consecutive_sum(i.month_table)
                                        i.focal.months <- names(i.max_month$values)
                                        
                                        #- create subset of i.data only containing the focal months
                                        i.data <- i.data[month %in% i.focal.months]
                                        
                                } else {
                                        #- Are the two months consecutive?
                                        i.monthdiff <- diff(as.numeric(names(
                                                i.month_table
                                        )))
                                        #- No or only one month
                                        if (length(i.monthdiff) == 0){
                                                # Does one of the month have more than 100 samples?
                                                if (any(i.month_table > min_samples)) {
                                                        i.focal.months <- names(which.max(
                                                                i.month_table
                                                        ))
                                                        i.data <- i.data[month == i.focal.months]
                                                } else {
                                                        counter[[b]] <- counter[[b]] + 1
                                                        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
                                                        next()
                                                }
                                        } else if (i.monthdiff > 2) {
                                                # Does one of the month have more than 100 samples?
                                                if (any(i.month_table > min_samples)) {
                                                        i.focal.months <- names(which.max(
                                                                i.month_table
                                                        ))
                                                        i.data <- i.data[month == i.focal.months]
                                                } else {
                                                        counter[[b]] <- counter[[b]] + 1
                                                        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
                                                        next()
                                                }
                                        } else {
                                                #- yes
                                                i.focal.months <- names(i.month_table)
                                                i.data <- i.data[month %in% i.focal.months]
                                        }
                                }
                        }
                        
                        #- how many samples are left?
                        i.samples <- uniqueN(i.data$eventID)
                        if (i.samples < min_samples) {
                                counter[[b]] <- counter[[b]] + 1
                                rm(list = ls()[grepl("^i\\.", x = ls())])
                                next()
                        }
                        

                        
                        # Keep only values ≤ actual sample size
                        i.sample_scheme <- log_sequence[log_sequence <= i.samples]
                        
                        # Always include the full sample size if not already present
                        if (max(i.sample_scheme) < i.samples) {
                                i.sample_scheme <- c(i.sample_scheme, i.samples)
                        }
                        # no need to keep predetermined sample if it is relatively similar in size compared to the full data set
                        i.fullDiff <- i.samples - i.sample_scheme[length(i.sample_scheme)-1]
                        i.fullDiffRel <- i.fullDiff /  i.samples
                        if (length(i.fullDiffRel) != 0){
                                if (i.fullDiffRel < 0.25){
                                        i.sample_scheme <- i.sample_scheme[- (length(i.sample_scheme)-1)]
                                }
                        }


                        # If somehow we have no valid samples (shouldn't happen given i.samples >= 50)
                        if (length(i.sample_scheme) == 0) {
                                i.sample_scheme <- i.samples
                        }
                        
                        # Limit to maximum 6 schemes to prevent excessive computation
                        # (typically won't be reached unless i.samples is very large)
                        if (length(i.sample_scheme) > 6) {
                                # Keep the smallest, largest, and evenly spaced intermediate values
                                i.keep_indices <- unique(c(
                                        1,  # Always keep smallest
                                        round(seq(2, length(i.sample_scheme) - 1, length.out = 4)),
                                        length(i.sample_scheme)  # Always keep largest (full sample)
                                ))
                                i.sample_scheme <- i.sample_scheme[unique(i.keep_indices)]
                        }
                        
                        # check what environmental zones are included
                        i.enz <- id.to.enz[ID %in% i.data$ID]
                        i.enz <- unique(i.enz$EnZ_name)
                        i.enz <- sort(i.enz)
                        
                        #i.sample_scheme
                        i.out <- data.table(
                                taxon    = bio.names[b],
                                data.set = ud[[b]][d],
                                eventYear     = unique(i.data$eventYear),
                                catchments = list(i.enz),
                                samples  = i.sample_scheme,
                                sample_type = c(rep("sub", length(i.sample_scheme)-1), "full")
                        )
                        if (all(is.na(ds.bio$eventYear))| all(is.na(ds.bio$eventDate))) {
                                i.out[, focal_months := list(list(NA))]
                        } else {
                                i.out[, focal_months := list(list(i.focal.months))]
                        }
                        ds.result.list[[length(ds.result.list) + 1]] <- i.out
                        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
                } # END of loop i over years in data.set d
                
                # combine results of all years. 
                ds.results.list <- rbindlist(ds.result.list)
                
                result_list[[length(result_list) + 1]] <- ds.results.list
                
                rm(list = ls()[grepl(pattern = "^ds\\.", x = ls())])
                
        } # END d of loop over data sets in taxonomic group b
} # END of loop b over taxonomic groups
        
result.data <- rbindlist(result_list, fill = T)
result.data$taxon <- factor(result.data$taxon) 
result.data <- split(result.data, f = result.data$taxon)
result.data <- lapply(result.data, function(x) x[, scheme_id := .GRP, by =  c("data.set", "eventYear", "samples", "sample_type")])
result.data <- rbindlist(result.data)
result.data[, scheme_id := sprintf("%s_%04d", taxon, scheme_id)]
split.data <- split(result.data, by = "taxon")

test <- result.data[sample_type == "full"]
split2 <- split(test, by = "taxon")


# # data output -------------------------------------------------------------
lapply(1:n_bio_datasets,
        function(x)
                saveRDS(
                        split.data[[x]],
                        paste0("data/biota/03_",bio.names[x],"_scheme.rds")
                )
)

counter


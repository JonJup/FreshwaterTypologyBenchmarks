# ==============================================================================
# setup
# ==============================================================================

library(data.table)
library(bundle)
library(workflows)
library(ranger)
library(tidymodels)

setwd("E://arbeit/Projekte/02_ongoing/PULSE/wp1/data/")

# ==============================================================================
# prep loop
# ==============================================================================

taxa <- c("diatoms", "fish", "invertebrates", "macrophytes")
results_importance <- list()
results_metrics    <- list()
results <- list()
# ==============================================================================
# loop
# ==============================================================================

for (t in taxa){
        # Print status using sprintf
        cat(sprintf("Starting with %s\n", t))
        
        # Extract the first letter of the taxa and capitalize it (e.g., "fish" -> "F")
        taxa_letter <- toupper(substr(t, 1, 1))
        
        # Dynamically build the directory path and set it
        dir_path <- sprintf("E://arbeit/Projekte/02_ongoing/PULSE/wp1/data/pulse%s", taxa_letter)
        setwd(dir_path)

        all_files <- list.files("data/008_qrf/", pattern = "\\.rds$", full.names = TRUE)
        met_files <- all_files[grepl("_metrics", all_files)]
        imp_files <- all_files[grepl("_importance", all_files)]
        mod_files <- setdiff(all_files, met_files)
        mod_files <- setdiff(mod_files, imp_files)
        
        imp <- lapply(imp_files, readRDS)
        imp <- lapply(imp, setDT)
        for (k in seq_along(imp_files)){
                imp[[k]][, metric := gsub("data/008_qrf/", "", imp_files[k])]
                imp[[k]][, metric := gsub("_importance.rds", "", metric)]
        }
        imp = rbindlist(imp)
        imp$taxon = t
        results_importance[[length(results_importance) + 1]] <- imp
        
        
        met <- lapply(met_files, readRDS)
        met <- lapply(met, setDT)
        for (k in seq_along(met_files)){
                met[[k]][, metric := gsub("data/008_qrf/", "", met_files[k])]
                met[[k]][, metric := gsub("_metrics.rds", "", metric)]
        }
        met = rbindlist(met)
        met$taxon = t
        results_metrics[[length(results_metrics) + 1]] <-met
        
# ==============================================================================
# save to file
# ==============================================================================

        cat("Changed WD \n")
        cat("\rLoaded evaluations \n")
        fileList <- list.files("data/007_evaluations/", full.names = TRUE)
        evl <- lapply(fileList, readRDS)
        evl <- rbindlist(evl)
        cat("\rLoading variation partitioning data...\n")
        vp <- list.files("data/005_variation_partitioning/", full.names = TRUE)
        vp <- lapply(vp, readRDS)
        vp2 <- list()
        for (i in seq_along(vp)) {
                x <- vp[[i]]$VP
                y <- x[, mean(scaled_values), by = c("driver")]
                y <- dcast(y, formula = . ~ driver, value.var = "V1")
                y$. <- unique(x$scheme)
                names(y)[1] <- "scheme_id"
                vp2[[i]] <- y
        }
        vp2 <- rbindlist(vp2)

        cat("\rLoading spatial scale data...\n")
        ss <- list.files("data/misc/spatial_scale/", full.names = TRUE)
        ss <- lapply(ss, readRDS)
        ss <- rbindlist(ss, fill = T)

        cat("\rLoading taxonomic resolution data...\n")
        tx <- list.files("data/misc/taxonomic_resolution/", full.names = TRUE)
        tx <- lapply(tx, readRDS)
        tx <- rbindlist(tx)

        if (length(gregexpr("_", tx$scheme_id[1])[[1]]) == 2) {
                tx$scheme_id <- sub("^([A-Za-z]+_)\\1", "\\1", tx$scheme_id)
        }


        cat("\rLoading number of taxa...\n")
        tn <- list.files("data/misc/taxa_counts/", full.names = TRUE)
        tn <- lapply(tn, readRDS)
        tn <- rbindlist(tn)

        if (length(gregexpr("_", tn$scheme_id[1])[[1]]) == 2) {
                tn[, scheme_id := sub("^([A-Za-z]+_)\\1", "\\1", scheme_id)]
        }
        cat("\rMerging datasets...\n")
        evl2 <- evl[ss, on = "scheme_id"]
        evl2 <- evl2[vp2, on = "scheme_id"]
        evl2 <- evl2[tx, on = "scheme_id"]
        evl2 <- evl2[tn, on = "scheme_id"]
        evl2 <- evl2[scheme_id %in% evl$scheme_id]

        # Remove duplicate columns
        if (any(grepl("^i\\.", names(evl2)))) {
                drop.id <- names(evl2)[grep("^i\\.", names(evl2))]
                evl2[, (drop.id) := NULL]
        }

        for (f in mod_files) {

                metric_name <- tools::file_path_sans_ext(basename(f))

                if (!metric_name %in%
                    c("isa_number",
                      "isa_avg_p",
                      "fuzzy_divergence",
                      "fuzzy_mantel")
                    ) {metric_name <- gsub("_", "\\ ", metric_name)}



                # Filter to specific metric
                cat("Filtering to metric:", metric_name, "\n")
                j.d <- evl2[metric == metric_name, ]

                # Check if we have enough data
                if (nrow(j.d) < 10) {
                        stop(paste("Not enough data for metric:", metric_name, "(only", nrow(j.d), "rows)"))
                }



                #===============================================================================
                # 4. Define predictors and model ----
                #===============================================================================

                predictor_vec <- c(
                        "n_types", "variables", "env_asw", "fuzzy_npe",
                        "samples", "min_distance", "mean_distance",
                        "median_distance", "max_distance", "median_latitude", "median_longitude",
                        "max_latitude", "max_longitude", "min_latitude", "min_longitude",
                        "bio", "env", "space", "stochastic", "species_rank",
                        "genus_rank", "family_rank", "higher_rank", "n_taxa"
                )

                if (!metric_name %in% c("fuzzy_divergence", "fuzzy_mantel")){
                        predictor_vec <- predictor_vec[!predictor_vec %in% c("fuzzy_npe")]
                }

                #===============================================================================
                # 5. Split data and create folds ----
                #===============================================================================

                set.seed(1)
                j.split <- initial_split(data = j.d)
                j.train <- training(j.split)
                j.test  <- testing(j.split)
                j.folds <- vfold_cv(j.train, v = 5)

                cat("\rTraining set:", nrow(j.train), "rows\n")
                cat("\rTest set:", nrow(j.test), "rows\n")

                #===============================================================================
                # 6. Create recipe and preprocess ----
                #===============================================================================
                cat("\rCreating recipe...\n")

                rf_recipe <-
                        recipe(value ~ ., data = j.train) %>%
                        update_role(everything(), new_role = "unused") %>% # Set all to unused first
                        update_role(all_of(predictor_vec), new_role = "predictor") %>%
                        update_role(value, new_role = "outcome") %>%
                        step_zv(all_predictors()) %>%
                        step_naomit(all_predictors(), value)


                # Prep the recipe
                prepped_recipe <- prep(rf_recipe, training = j.train)
                test_processed <- bake(prepped_recipe, new_data = j.test)

                # Unbundle
                j.bundle <- readRDS(f)
                j.wf <- unbundle(j.bundle)
                j.ranger <- extract_fit_engine(j.wf)

                # OOB metrics from ranger

                # You'd also need to reload the test data and predict quantiles:
                pred_q <- predict(j.ranger, data = test_processed,
                                  type = "quantiles",
                                  quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95))

                observed <- test_processed$value

                coverage_90 <- mean(observed >= pred_q$predictions[,1] &
                                            observed <= pred_q$predictions[,5])
                coverage_50 <- mean(observed >= pred_q$predictions[,2] &
                                            observed <= pred_q$predictions[,4])
                cat("\rSaving to list")
                results[[length(results)+1]] <- data.table(
                        taxon = t,
                        metric = metric_name,
                        # oob_rmse = oob_rmse,
                        # oob_r2 = oob_r2,
                        cov_90 = coverage_90,
                        cov_50 = coverage_50
                )
        }   # END OF LOOP OVER METRICS
} # END OF LOOP OVER TAXA

metrics    <- rbindlist(results_metrics)
importance <- rbindlist(results_importance)
setwd(rstudioapi::getActiveProject())
saveRDS(metrics   , "data/results/qrf_metrics.rds")
saveRDS(importance, "data/results/qrf_variable_importance.rds")
results_dt <- rbindlist(results)
saveRDS(results_dt, "data/results/qrf_interval_coverage.rds")

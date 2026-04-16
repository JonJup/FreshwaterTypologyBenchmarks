#' Calculate environmental optima for species
#' @param comm Community matrix (sites x species)
#' @param env_vars Environmental variables (sites x variables)
#' @param focal_species Character vector of species names to check (default: 10 most abundant)
calculate_species_optima <- function(comm, env_vars, focal_species = NULL) {
        
        # Select focal species (most abundant if not specified)
        if(is.null(focal_species)) {
                abundances <- colSums(comm)
                focal_species <- names(sort(abundances, decreasing = TRUE)[1:min(10, ncol(comm))])
        }
        
        optima_list <- list()
        
        for(sp in focal_species) {
                # Weighted mean of environmental values (weighted by abundance)
                sp_abundance <- comm[, sp]
                
                if(sum(sp_abundance) == 0) next  # Skip if species not present
                
                # Calculate optimum for each environmental variable
                optima <- sapply(env_vars, 
                                 function(env_col) {
                                         weighted.mean(env_col,
                                                       w = sp_abundance
                                                       )
                                         }
                                 )
                
                # Calculate tolerance (weighted SD)
                tolerance <- sapply(env_vars, function(env_col) {
                        sqrt(weighted.mean((env_col - weighted.mean(env_col, w = sp_abundance))^2, 
                                           w = sp_abundance))
                })
                
                optima_list[[sp]] <- data.frame(
                        species = sp,
                        variable = names(env_vars),
                        optimum = optima,
                        tolerance = tolerance
                )
        }
        
        do.call(rbind, optima_list)
}

#' Posterior predictive check for species-environment relationships
#' @param obs_comm Observed community matrix
#' @param sim_comm_list List of simulated community matrices
#' @param env_vars Environmental variables
#' @param focal_species Species to check
ppc_species_optima <- function(obs_comm, sim_comm_list, env_vars, focal_species = NULL) {
        
        # Observed optima
        obs_optima <- calculate_species_optima(obs_comm, env_vars, focal_species)
        
        # Simulated optima
        sim_optima_list <- lapply(sim_comm_list, function(sim) {
                calculate_species_optima(sim, env_vars, obs_optima$species)
        })
        
        # For each species-variable combination, calculate deviation
        sim_results <- apply(obs_optima, 1, function(row) {
                sp  <- row[["species"]]
                var <- row[["variable"]]
                opt <- as.numeric(row[["optimum"]])
                
                sim_vals <- sapply(sim_optima_list, function(sim_opt) {
                        sim_opt$optimum[sim_opt$species == sp & sim_opt$variable == var]
                })
                
                sim_mean <- mean(sim_vals)
                sim_sd   <- sd(sim_vals)
                z_score  <- (opt - sim_mean) / sim_sd
                flag     <- abs(z_score) > 2
                
                list(
                        sim_optima_values = sim_vals,
                        sim_mean = sim_mean,
                        sim_sd = sim_sd,
                        z_score = z_score,
                        flag = flag
                )
        })
        results <- obs_optima
        results$sim_optima_values <- lapply(sim_results, `[[`, "sim_optima_values")
        results$sim_mean          <- sapply(sim_results, `[[`, "sim_mean")
        results$sim_sd            <- sapply(sim_results, `[[`, "sim_sd")
        results$z_score           <- sapply(sim_results, `[[`, "z_score")
        results$flag              <- sapply(sim_results, `[[`, "flag")
        
        
        # Summary: what proportion of species-environment relationships are problematic?
        summary_stats <- list(
                n_relationships = nrow(results),
                n_flagged = sum(results$flag),
                prop_flagged = mean(results$flag),
                max_deviation = max(abs(results$z_score)),
                overall_flag = mean(results$flag) > 0.20  # Flag if >20% deviate
        )
        
        return(list(
                details = results,
                summary = summary_stats
        ))
}

# Usage:
# env_data <- data.frame(temperature = ..., pH = ..., elevation = ...)
# result <- ppc_species_optima(observed_community, list_of_simulations, env_data)
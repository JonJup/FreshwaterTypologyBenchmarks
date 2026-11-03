#' Calculate C-score for co-occurrence analysis
#' @param comm Community matrix (sites x species), presence-absence
calculate_cscore <- function(comm) {
        # Convert to presence-absence if needed
        comm_pa <- ifelse(comm > 0, 1, 0)
        
        # Remove species with <2 occurrences (can't have checkerboard patterns)
        species_occurrences <- colSums(comm_pa)
        comm_pa <- comm_pa[, species_occurrences >= 2]
        
        # Calculate C-score 
        n_species <- ncol(comm_pa)
        n_pairs <- 0
        c_score_sum <- 0
        
        for(i in 1:(n_species-1)) {
                for(j in (i+1):n_species) {
                        # For each species pair, count checkerboard units
                        sp1_present <- comm_pa[, i] == 1
                        sp2_present <- comm_pa[, j] == 1
                        
                        # Sites where sp1 present but sp2 absent
                        sp1_only <- sum(sp1_present & !sp2_present)
                        # Sites where sp2 present but sp1 absent  
                        sp2_only <- sum(!sp1_present & sp2_present)
                        # C-score for this pair (number of checkerboard units)
                        c_score_sum <- c_score_sum + (sp1_only * sp2_only)
                        n_pairs <- n_pairs + 1
                }
        }
        
        # Mean C-score across all pairs
        return(c_score_sum / n_pairs)
}


#' Posterior predictive check for co-occurrence patterns
#' @param obs_comm Observed community matrix
#' @param sim_comm_list List of simulated community matrices
ppc_cooccurrence <- function(obs_comm, sim_comm_list) {
        
        # Observed C-score
        obs_cscore <- calculate_cscore(obs_comm)
        
        # Simulated C-scores
        sim_cscores <- sapply(sim_comm_list, calculate_cscore)
        
        # Standardized effect size
        z_score <- (obs_cscore - mean(sim_cscores)) / sd(sim_cscores)
        
        # # Bayesian posterior probability
        # p_extreme <- mean(abs(sim_cscores - mean(sim_cscores)) >= 
        #                           abs(obs_cscore - mean(sim_cscores)))
        
        # Results
        results <- list(
                obs_cscore = obs_cscore,
                sim_mean = mean(sim_cscores),
                sim_sd = sd(sim_cscores),
                sim_quantiles = quantile(sim_cscores, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
                z_score = z_score,
                # p_extreme = p_extreme,
                flag = abs(z_score) > 2  # Flag if outside 95% envelope
        )
        
        return(results)
}

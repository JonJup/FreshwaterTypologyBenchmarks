#' Calculate species prevalence (occupancy)
#' @param comm Community matrix (sites x species)
calculate_prevalence <- function(comm) {
        # Convert to presence-absence
        comm_pa <- ifelse(comm > 0, 1, 0)
        
        # Proportion of sites each species occupies
        prevalence <- colSums(comm_pa) / nrow(comm_pa)
        
        return(prevalence)
}

#' Posterior predictive check for prevalence distribution
ppc_prevalence <- function(obs_comm, sim_comm_list) {
        
        # Observed prevalence
        obs_prev <- calculate_prevalence(obs_comm)
        
        # Simulated prevalence distributions
        sim_prev_list <- lapply(sim_comm_list, calculate_prevalence)
        
        # Key summary statistics for the distribution
        summarize_prevalence <- function(prev) {
                c(
                        n_rare = sum(prev < 0.05),           # Very rare species
                        n_uncommon = sum(prev >= 0.05 & prev < 0.25),
                        n_common = sum(prev >= 0.25 & prev < 0.75),
                        n_ubiquitous = sum(prev >= 0.75),    # Very common species
                        mean_prev = mean(prev),
                        median_prev = median(prev),
                        sd_prev = sd(prev),
                        skewness = moments::skewness(prev)   # Shape of distribution
                )
        }
        
        obs_summary <- summarize_prevalence(obs_prev)
        sim_summaries <- t(sapply(sim_prev_list, summarize_prevalence))
        
        # standard deviations are computed beforehand to catch any variables without deviations
        # this would lead to NA z-scores 
        sds <- apply(sim_summaries, 2, sd)
        if (any(sds == 0)) sds[which(sds == 0)] <- 0.000001
        # Calculate z-scores for each summary statistic
        z_scores <- (obs_summary - colMeans(sim_summaries)) / sds
        
        # Flag if key statistics are problematic
        flags <- abs(z_scores) > 2
        
        # Detailed comparison
        results <- data.frame(
                statistic = names(obs_summary),
                observed = obs_summary,
                sim_mean = colMeans(sim_summaries),
                sim_sd = apply(sim_summaries, 2, sd),
                z_score = z_scores,
                flag = flags
        )
        
        # # Create prevalence histogram comparison data for plotting
        # prevalence_bins <- seq(0, 1, by = 0.05)
        # 
        # obs_hist <- hist(obs_prev, breaks = prevalence_bins, plot = FALSE)
        # sim_hists <- lapply(sim_prev_list, function(sp) {
        #         hist(sp, breaks = prevalence_bins, plot = FALSE)$counts
        # })
        # sim_hist_matrix <- do.call(rbind, sim_hists)
        # 
        # histogram_comparison <- data.frame(
        #         bin_midpoint = obs_hist$mids,
        #         obs_count = obs_hist$counts,
        #         sim_mean = colMeans(sim_hist_matrix),
        #         sim_lower = apply(sim_hist_matrix, 2, quantile, probs = 0.025),
        #         sim_upper = apply(sim_hist_matrix, 2, quantile, probs = 0.975)
        # )
        
        # Overall flag: are rare species or ubiquitous species counts way off?
        summary_flag <- sum(flags[c("n_rare", "n_ubiquitous")]) >= 1
        
        return(list(
                summary_stats = results,
                
                overall_flag = sum(results$flag) > 2
        ))
}

# Plotting function
# plot_prevalence_check <- function(ppc_result) {
#         library(ggplot2)
#         
#         ggplot(ppc_result$histogram, aes(x = bin_midpoint)) +
#                 geom_ribbon(aes(ymin = sim_lower, ymax = sim_upper), 
#                             fill = "grey80", alpha = 0.5) +
#                 geom_line(aes(y = sim_mean), color = "grey40", linewidth = 1) +
#                 geom_line(aes(y = obs_count), color = "red", linewidth = 1.5) +
#                 geom_point(aes(y = obs_count), color = "red", size = 2) +
#                 labs(x = "Species Prevalence (proportion of sites)",
#                      y = "Number of Species",
#                      title = "Species Prevalence Distribution: Observed vs Simulated",
#                      subtitle = "Red = Observed, Grey = 95% envelope of simulations") +
#                 theme_minimal()
# }

# Usage:
# result <- ppc_prevalence(observed_community, simulations)
# print(result$summary_stats)
# plot_prevalence_check(result)
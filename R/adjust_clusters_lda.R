################################################################################
# Function:           adjust_clusters_lda.R
# Description:        LDA-based alternative to adjust_clusters(). Manipulates
#                     cluster separation and within-cluster dispersion along the
#                     discriminant axes — the directions that actually distinguish
#                     the clusters — rather than isotropically in the original
#                     variable space.
#
# Why this is better: Isotropic scaling wastes manipulation budget on dimensions
#                     where clusters already don't overlap (or where they can't
#                     be separated). LDA identifies the directions that carry
#                     the cluster signal, so every unit of perturbation directly
#                     affects cluster distinction. This means:
#                       - More ASW change per unit of environmental perturbation
#                       - Fewer combinations rejected by plausibility filters
#                       - Separation and dispersion factors map cleanly onto
#                         actual between/within-cluster geometry
#
# Author:             Jonathan Jupke
# Date Created:       2026-03-16
#
# Dependencies:       MASS (for lda)
#
# Usage:
#   result <- adjust_clusters_lda(
#     centroids           = centroid_matrix,
#     observations        = obs_matrix,
#     cluster_assignments = cluster_vec,
#     separation_factor   = 1.2,
#     dispersion_factor   = 0.8,
#     constraints         = list(lower = 0)   # optional
#   )
################################################################################

library(MASS)

# ==============================================================================
# MAIN FUNCTION
# ==============================================================================

#' Adjust cluster separation and dispersion along LDA discriminant axes
#'
#' @param centroids           Matrix (k × p). Cluster centroids in original space.
#' @param observations        Matrix (n × p). Observations in original space.
#' @param cluster_assignments Integer vector (n). Cluster ID per observation.
##' @param quality_factor  Numeric. Single control for cluster coherence.
#'                        > 1 = better typology (clusters pulled apart, tightened)
#'                        < 1 = worse typology  (clusters pulled together, loosened)
#'                        1.0 = no change
#'                        Internally maps to:
#'                          separation_factor = quality_factor
#'                          dispersion_factor = 1 / quality_factor
#' @param constraints         Optional list with named elements:
#'                            - lower: scalar or named vector of lower bounds
#'                            - upper: scalar or named vector of upper bounds
#'                            Applied after back-projection.
#'
#'
#' @return List with:
#'   - observations:  Adjusted observation matrix (n × p)
#'   - centroids:     Adjusted centroid matrix (k × p)
#'   - discriminant_axes: The LDA scaling matrix (p × d), useful for diagnostics
#'   - variance_explained: Proportion of between-class variance on each axis
#'   - effective_separation:  Effective separation factor achieved after constraints
#'   - effective_dispersion:  Effective dispersion factor achieved after constraints

adjust_clusters_lda <- function(centroids,
                                observations,
                                cluster_assignments,
                                quality_factor      = 1.0,   
                                constraints         = NULL) {
        
        # ------------------------------------------------------------------
        # 0. Input validation
        # ------------------------------------------------------------------
        observations <- as.matrix(observations)
        centroids    <- as.matrix(centroids)
        n <- nrow(observations)
        p <- ncol(observations)
        k <- nrow(centroids)
        d <- min(k - 1, p)  # discriminant dimensionality
        
        stopifnot(
                ncol(centroids) == p,
                length(cluster_assignments) == n,
                d >= 1
        )
        
        cluster_ids <- sort(unique(cluster_assignments))
        
        # ------------------------------------------------------------------
        # 1. Fit LDA via MASS::lda()
        # ------------------------------------------------------------------
        # MASS::lda tol parameter rejects variables with variance < tol^2.
        # This handles collinearity. For near-singular within-class scatter,
        # lda() is already robust internally.
        lda_fit <- MASS::lda(
                x        = observations,
                grouping = factor(cluster_assignments)
        )
        
        # scaling matrix: p × d (discriminant axes in original variable space)
        scaling <- lda_fit$scaling[, 1:d, drop = FALSE]
        
        # Eigenvalues from the svd slot: proportion of between-class variance
        eigenvalues <- lda_fit$svd^2
        eigenvalues <- eigenvalues[1:d]
        variance_explained <- eigenvalues / max(sum(eigenvalues), 1e-10)
        
        grand_mean <- colMeans(observations)
        
        # ------------------------------------------------------------------
        # 2. Project into discriminant space
        # ------------------------------------------------------------------
        # Z = (X - grand_mean) %*% scaling  →  n × d
        obs_centered   <- sweep(observations, 2, grand_mean)
        Z_obs          <- obs_centered %*% scaling
        
        cent_centered  <- sweep(centroids, 2, grand_mean)
        Z_cent         <- cent_centered %*% scaling
        
        Z_grand_mean   <- colMeans(Z_obs)
        
        # ------------------------------------------------------------------
        # 3. Manipulate in discriminant space
        # ------------------------------------------------------------------
        
        # Derive the two internal factors from the single quality parameter
        separation_factor <- quality_factor
        dispersion_factor <- 1
        
        # 3a. Separation: scale centroid positions relative to grand mean
        Z_cent_new <- sweep(Z_cent, 2, Z_grand_mean)          # center
        Z_cent_new <- Z_cent_new * separation_factor           # scale
        Z_cent_new <- sweep(Z_cent_new, 2, -Z_grand_mean)     # uncenter
        
        # 3b. Dispersion: scale observations relative to their (new) centroid
        Z_obs_new <- Z_obs
        for (cl in cluster_ids) {
                idx <- which(cluster_assignments == cl)
                deviations <- sweep(Z_obs[idx, , drop = FALSE], 2, Z_cent[cl, ])
                Z_obs_new[idx, ] <- sweep(deviations * dispersion_factor, 2, -Z_cent_new[cl, ])
        }
        
        # ------------------------------------------------------------------
        # 4. Back-project to original variable space
        # ------------------------------------------------------------------
        # We only modify the discriminant subspace component.
        # The null-space component (everything LDA can't see) is preserved.
        #
        # delta_Z = Z_new - Z_orig  (change in discriminant coordinates)
        # delta_X = delta_Z %*% ginv(scaling)  (change in original space)
        # X_new = X + delta_X
        #
        # This preserves all non-discriminant structure exactly.
        
        scaling_pinv <- tryCatch(
                # Use pseudoinverse: scaling^+ = (scaling^T scaling)^{-1} scaling^T
                solve(crossprod(scaling), t(scaling)),
                error = function(e) {
                        # Fallback to MASS::ginv if crossprod is singular
                        MASS::ginv(scaling)
                }
        )
        
        delta_Z_obs  <- Z_obs_new - Z_obs
        delta_X_obs  <- delta_Z_obs %*% scaling_pinv
        obs_new      <- observations + delta_X_obs
        
        delta_Z_cent <- Z_cent_new - Z_cent
        delta_X_cent <- delta_Z_cent %*% scaling_pinv
        cent_new     <- centroids + delta_X_cent
        
        # ------------------------------------------------------------------
        # 5. Apply constraints (optional)
        # ------------------------------------------------------------------
        if (!is.null(constraints)) {
                obs_new  <- apply_constraints(obs_new, constraints, colnames(observations))
                cent_new <- apply_constraints(cent_new, constraints, colnames(observations))
        }
        
        # ------------------------------------------------------------------
        # 6. Compute effective factors (in discriminant space)
        # ------------------------------------------------------------------
        # Re-project constrained results back to discriminant space to see
        # what factors were actually achieved after constraint clipping.
        
        Z_obs_final  <- sweep(obs_new, 2, grand_mean) %*% scaling
        Z_cent_final <- sweep(cent_new, 2, grand_mean) %*% scaling
        
        eff_sep  <- compute_effective_separation(Z_cent, Z_cent_final, Z_grand_mean)
        eff_disp <- compute_effective_dispersion(Z_obs, Z_obs_final, Z_cent, Z_cent_final,
                                                 cluster_assignments, cluster_ids)
        
        # ------------------------------------------------------------------
        # 7. Return
        # ------------------------------------------------------------------
        
        # Preserve column names
        colnames(obs_new)  <- colnames(observations)
        colnames(cent_new) <- colnames(centroids)
        
        return(list(
                observations         = obs_new,
                centroids            = cent_new,
                discriminant_axes    = scaling,
                eigenvalues          = eigenvalues,
                variance_explained   = variance_explained,
                effective_separation = eff_sep,
                effective_dispersion = eff_disp,
                effective_quality    = eff_sep / eff_disp   
        ))
}


# ==============================================================================
# HELPER: Apply variable constraints
# ==============================================================================

apply_constraints <- function(mat, constraints, var_names = NULL) {
        
        if (!is.null(constraints$lower)) {
                if (length(constraints$lower) == 1) {
                        # Scalar: apply to all columns
                        mat[mat < constraints$lower] <- constraints$lower
                } else if (!is.null(var_names) && !is.null(names(constraints$lower))) {
                        # Named vector: apply per variable
                        for (v in names(constraints$lower)) {
                                if (v %in% var_names) {
                                        col_idx <- which(var_names == v)
                                        mat[, col_idx] <- pmax(mat[, col_idx], constraints$lower[v])
                                }
                        }
                }
        }
        
        if (!is.null(constraints$upper)) {
                if (length(constraints$upper) == 1) {
                        mat[mat > constraints$upper] <- constraints$upper
                } else if (!is.null(var_names) && !is.null(names(constraints$upper))) {
                        for (v in names(constraints$upper)) {
                                if (v %in% var_names) {
                                        col_idx <- which(var_names == v)
                                        mat[, col_idx] <- pmin(mat[, col_idx], constraints$upper[v])
                                }
                        }
                }
        }
        
        return(mat)
}


# ==============================================================================
# HELPER: Effective separation factor in discriminant space
# ==============================================================================
# Ratio of centroid-to-grand-mean distances: after / before

compute_effective_separation <- function(Z_cent_orig, Z_cent_final, Z_grand_mean) {
        
        dist_orig  <- sqrt(rowSums(sweep(Z_cent_orig, 2, Z_grand_mean)^2))
        dist_final <- sqrt(rowSums(sweep(Z_cent_final, 2, Z_grand_mean)^2))
        
        valid <- dist_orig > 1e-10
        if (sum(valid) == 0) return(1.0)
        
        # Weighted mean by original distance (larger clusters contribute more)
        mean(dist_final[valid] / dist_orig[valid])
}


# ==============================================================================
# HELPER: Effective dispersion factor in discriminant space
# ==============================================================================
# Ratio of point-to-centroid distances: after / before

compute_effective_dispersion <- function(Z_obs_orig, Z_obs_final,
                                         Z_cent_orig, Z_cent_final,
                                         cluster_assignments, cluster_ids) {
        
        ratios <- numeric(0)
        
        for (cl in cluster_ids) {
                idx <- which(cluster_assignments == cl)
                
                # Original distances to original centroid
                dev_orig <- sweep(Z_obs_orig[idx, , drop = FALSE], 2, Z_cent_orig[cl, ])
                d_orig   <- sqrt(rowSums(dev_orig^2))
                
                # Final distances to final centroid
                dev_final <- sweep(Z_obs_final[idx, , drop = FALSE], 2, Z_cent_final[cl, ])
                d_final   <- sqrt(rowSums(dev_final^2))
                
                valid <- d_orig > 1e-10
                if (sum(valid) > 0) {
                        ratios <- c(ratios, d_final[valid] / d_orig[valid])
                }
        }
        
        if (length(ratios) == 0) return(1.0)
        mean(ratios)
}


# ==============================================================================
# DIAGNOSTIC: Visualize what LDA sees
# ==============================================================================
# Call this to understand the discriminant structure before/after manipulation.
# Returns a data.table suitable for ggplot.

discriminant_diagnostic <- function(observations,
                                    cluster_assignments,
                                    result,
                                    max_axes = 3) {
        
        scaling    <- result$discriminant_axes
        grand_mean <- colMeans(observations)
        d          <- min(ncol(scaling), max_axes)
        
        # Project original and adjusted into discriminant space
        Z_orig <- sweep(observations, 2, grand_mean) %*% scaling[, 1:d, drop = FALSE]
        Z_new  <- sweep(result$observations, 2, grand_mean) %*% scaling[, 1:d, drop = FALSE]
        
        # Build data.table for plotting
        dt_orig <- data.table::data.table(
                as.data.frame(Z_orig),
                cluster = factor(cluster_assignments),
                state   = "Before"
        )
        dt_new <- data.table::data.table(
                as.data.frame(Z_new),
                cluster = factor(cluster_assignments),
                state   = "After"
        )
        
        # Standardize column names
        ax_names <- paste0("LD", 1:d)
        names(dt_orig)[1:d] <- ax_names
        names(dt_new)[1:d]  <- ax_names
        
        dt <- rbind(dt_orig, dt_new)
        dt[, state := factor(state, levels = c("Before", "After"))]
        
        # Attach variance explained
        attr(dt, "variance_explained") <- result$variance_explained[1:d]
        attr(dt, "axis_names") <- ax_names
        
        return(dt)
}
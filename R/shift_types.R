################################################################################
# File:         shift_types.R
# Description:  Unit-free manipulation of how distinct the types of a typology
#               are, plus the descriptors recorded for every manipulated
#               environment.
# Dependencies: MASS (only for ld1_dominance(), called as MASS::lda)
#
################################################################################


# ==============================================================================
# type_deviations(): centroid deviations from the grand mean
# ==============================================================================
# X  : numeric matrix (n sites x p variables)
# cl : type label per site (length n)
#
# Returns a k x p matrix with one row per type (row names = type labels)
# holding c_t - mu. mu is colMeans(X), i.e. the size-weighted mean of the
# centroids. Using the size-weighted mean is what guarantees that the grand
# mean is unchanged by shift_types(), because sum_t n_t (c_t - mu) = 0.

type_deviations <- function(X, cl) {
        X   <- as.matrix(X)
        cl  <- as.character(cl)
        tot <- rowsum(X, cl)                           # k x p type sums, rows named by type
        n_t <- as.vector(table(cl)[rownames(tot)])     # type sizes in the same row order
        sweep(tot / n_t, 2, colMeans(X))               # centroids minus grand mean
}


# ==============================================================================
# shift_types(): the manipulation
# ==============================================================================
# X     : numeric matrix (n x p) of the variables to manipulate (the subset S)
# cl    : type label per site
# gamma : non-negative scalar; multiplier of the between-type component
#
# Returns X(gamma) with the dimnames of X.

shift_types <- function(X, cl, gamma) {
        X <- as.matrix(X)
        stopifnot(is.numeric(gamma), length(gamma) == 1L, gamma >= 0,
                  nrow(X) == length(cl))
        dev <- type_deviations(X, cl)
        # Row t(i) of dev is the shift direction of site i; the same vector
        # for all sites of a type, hence a rigid translation of the type.
        out <- X + (gamma - 1) * dev[as.character(cl), , drop = FALSE]
        # Arithmetic would otherwise copy the (duplicated) type labels from
        # dev into the row names of the result.
        dimnames(out) <- dimnames(X)
        out
}


# ==============================================================================
# clip_to_bounds(): enforce physical bounds after shifting
# ==============================================================================
# X     : numeric matrix with column names
# lower : named numeric vector of lower bounds (variables without a bound are
#         simply not listed); NULL for none
# upper : named numeric vector of upper bounds; NULL for none
#
# Only columns of X that appear in lower / upper are touched. Returns a list
#   X         : the clipped matrix
#   n_clipped : number of values that were moved onto a bound
#
# Clipping breaks the exact properties listed in the header for the clipped
# values (the affected sites no longer move rigidly with their type). Record
# n_clipped and achieved_gamma() to quantify this.

clip_to_bounds <- function(X, lower = NULL, upper = NULL) {
        X <- as.matrix(X)
        n_clipped <- 0L
        for (v in intersect(names(lower), colnames(X))) {
                below         <- X[, v] < lower[[v]]
                n_clipped     <- n_clipped + sum(below)
                X[below, v]   <- lower[[v]]
        }
        for (v in intersect(names(upper), colnames(X))) {
                above         <- X[, v] > upper[[v]]
                n_clipped     <- n_clipped + sum(above)
                X[above, v]   <- upper[[v]]
        }
        list(X = X, n_clipped = n_clipped)
}


# ==============================================================================
# between_type_ss(): between-type sum of squares
# ==============================================================================
# sum_t n_t * ||c_t - mu||^2, summed over all columns of X. Pass X in
# standardised (model) units so that every variable carries the same weight;
# in raw units the sum would be dominated by the variables with the largest
# numbers.

between_type_ss <- function(X, cl) {
        dev <- type_deviations(X, cl)
        n_t <- as.vector(table(as.character(cl))[rownames(dev)])
        sum(n_t * rowSums(dev^2))
}


# ==============================================================================
# achieved_gamma(): gamma actually realised after clipping
# ==============================================================================
# X_orig, X_new : the variables in S before and after manipulation, both in
#                 standardised (model) units
# cl            : type label per site
#
# Without clipping, the between-type sum of squares scales by gamma^2 in every
# variable, so this returns the requested gamma exactly (in any units). With
# clipping it returns the square root of the realised between-type SS ratio,
# i.e. the gamma that an unclipped manipulation with the same overall
# between-type change would have had.

achieved_gamma <- function(X_orig, X_new, cl) {
        ss_orig <- between_type_ss(X_orig, cl)
        if (ss_orig < 1e-12) return(NA_real_)
        sqrt(between_type_ss(X_new, cl) / ss_orig)
}


# ==============================================================================
# ld1_dominance(): dimensionality descriptor of a typology
# ==============================================================================
# X  : variables in S (any units; LDA is invariant to linear rescaling)
# cl : type label per site
#
# With r = min(k - 1, |S|) discriminant axes, the share of between-type variance
# on LD1, s1, lies in [1/r, 1]. The normalised dominance
#
#     D = (s1 - 1/r) / (1 - 1/r)
#
# is 0 when the between-type variance is spread evenly over all axes and 1 when
# it lies entirely on LD1 (a typology along a single gradient). With a single
# axis (k = 2) the typology is one-dimensional by definition and D = 1.
#
# shift_types() does not change D, so D is a property of the clustering and is
# recorded once per iteration. It is a descriptor, not a filter.
#
# Returns c(ld1_share, ld1_dominance, n_ld_axes); all NA if the LDA fails (e.g.
# a variable constant within all types).

ld1_dominance <- function(X, cl) {
        out <- c(ld1_share = NA_real_, ld1_dominance = NA_real_, n_ld_axes = NA_real_)
        fit <- tryCatch(suppressWarnings(MASS::lda(x = as.matrix(X), grouping = factor(cl))),
                        error = function(e) NULL)
        if (is.null(fit)) return(out)
        ev <- fit$svd^2                     # canonical F statistics, one per LD axis
        r  <- length(ev)
        s1 <- ev[1] / sum(ev)
        out[] <- c(s1, if (r > 1) (s1 - 1 / r) / (1 - 1 / r) else 1, r)
        out
}

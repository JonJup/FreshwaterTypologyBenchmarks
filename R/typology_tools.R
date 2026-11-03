################################################################################
# File:         typology_tools.R
# Description:  (1) Construction of typologies from three classification
#                   families: spatially constrained Ward clustering and crossed
#                   threshold typologies (k-means is called directly in the
#                   simulation script), and their fuzzy counterparts.
#               (2) Algorithm-independent descriptors of any typology. The same
#                   functions are used for simulated and for real typologies
#                   (e.g. broad river types, HER, FEOW), so both are described
#                   identically.
# Author:       Jonathan Jupke
# Dependencies: stats (base R), cluster (silhouette), vegclust (fuzzy c-means)
################################################################################


# ==============================================================================
# 1. TYPOLOGY CONSTRUCTION
# ==============================================================================

# ------------------------------------------------------------------------------
# geo_dist(): geographic distances between sites
# ------------------------------------------------------------------------------
# coords : n x 2 matrix (x, y). If all values fall within [-180, 180] x
#          [-90, 90], they are treated as longitude / latitude and great-circle
#          distances (km) are returned; otherwise Euclidean distances in the
#          units of the (projected) coordinates.
# Returns a dist object.

geo_dist <- function(coords) {
        coords <- as.matrix(coords)
        is_lonlat <- all(abs(coords[, 1]) <= 180) && all(abs(coords[, 2]) <= 90)
        if (!is_lonlat) return(dist(coords))
        lon  <- coords[, 1] * pi / 180
        lat  <- coords[, 2] * pi / 180
        a    <- sin(outer(lat, lat, "-") / 2)^2 +
                outer(cos(lat), cos(lat)) * sin(outer(lon, lon, "-") / 2)^2
        as.dist(2 * 6371 * asin(pmin(sqrt(a), 1)))            # haversine, km (pmin keeps dims of its first argument)
}


# ------------------------------------------------------------------------------
# geo_coordinates(): planar coordinates for the spatially constrained methods
# ------------------------------------------------------------------------------
# Projected coordinates are returned unchanged. Longitude/latitude are replaced
# by a two-dimensional classical MDS embedding of the great-circle distances
# (max. relative distance error < 1% over a European extent), so that spatial
# Ward and its fuzzy counterpart can both work with Euclidean coordinates and
# receive identical geographic information.

geo_coordinates <- function(coords) {
        coords <- as.matrix(coords)
        is_lonlat <- all(abs(coords[, 1]) <= 180) && all(abs(coords[, 2]) <= 90)
        if (!is_lonlat) return(coords)
        cmdscale(geo_dist(coords), k = 2)
}


# ------------------------------------------------------------------------------
# spatial_feature_space(): environment and geography in one Euclidean space
# ------------------------------------------------------------------------------
# Z = [ sqrt(1 - alpha) * X / max D0 ,  sqrt(alpha) * XY / max D1 ]
#
# with D0 = dist(X) and D1 = dist(XY). Euclidean distances in Z equal the
# combined dissimilarity of ClustGeo,
#
#     D_alpha^2 = (1 - alpha) * (D0 / max D0)^2 + alpha * (D1 / max D1)^2,
#
# exactly (verified numerically). Any method that works in Z therefore uses the
# same environmental and geographic information, weighted by the same alpha.

spatial_feature_space <- function(X, geo_xy, alpha) {
        stopifnot(alpha >= 0, alpha <= 1)
        X      <- as.matrix(X)
        geo_xy <- as.matrix(geo_xy)
        cbind(sqrt(1 - alpha) * X      / max(dist(X)),
              sqrt(alpha)     * geo_xy / max(dist(geo_xy)))
}


# ------------------------------------------------------------------------------
# spatial_ward(): spatially constrained hierarchical clustering
# ------------------------------------------------------------------------------
# Ward clustering of D_alpha (see spatial_feature_space()).
#
# X      : environmental variables (here: standardised S)
# geo_xy : planar site coordinates from geo_coordinates()
# alpha  : in [0, 1]; 0 = purely environmental (Ward), 1 = purely spatial
# k      : number of types
#
# Reproduces ClustGeo::hclustgeo(dist(X), dist(geo_xy), alpha) (Chavent et al.
# 2018) exactly: hclustgeo() runs hclust(method = "ward.D") on
# (1 - alpha) D0^2 / (2n) + alpha D1^2 / (2n) after the same normalisation; the
# constant 1 / (2n) does not change the tree, and "ward.D2" on distances equals
# "ward.D" on their squares. Verified on 200 random data sets.

spatial_ward <- function(X, geo_xy, alpha, k) {
        cutree(hclust(dist(spatial_feature_space(X, geo_xy, alpha)), method = "ward.D2"), k = k)
}


# ------------------------------------------------------------------------------
# merge_small_types(): enforce a minimum type size and a maximum number of types
# ------------------------------------------------------------------------------
# Repeatedly merges the smallest type into the type with the nearest centroid
# (Euclidean, in the columns of X) until every type has >= min_size sites and
# there are <= k_max types. Mimics how rare combinations are pooled in
# threshold typologies. Returns labels relabelled to 1..k.

merge_small_types <- function(X, cl, min_size, k_max) {
        X  <- as.matrix(X)
        cl <- as.integer(factor(cl))
        repeat {
                sizes <- tabulate(cl)
                if (length(sizes) <= 1L) break
                if (min(sizes) >= min_size && length(sizes) <= k_max) break
                s <- which.min(sizes)                                   # smallest type
                C <- rowsum(X, cl) / sizes                              # centroids, rows 1..k
                d <- sqrt(colSums((t(C) - C[s, ])^2))                   # centroid distances to s
                d[s] <- Inf
                cl[cl == s] <- which.min(d)
                cl <- as.integer(factor(cl))
        }
        cl
}


# ------------------------------------------------------------------------------
# crossed_thresholds(): a priori (WFD System A-like) threshold typology
# ------------------------------------------------------------------------------
# Every variable in X is split into classes independently by one-dimensional
# k-means (the thresholds lie midway between adjacent class centres); the
# number of classes per variable is drawn from `classes`. Types are all
# combinations of classes (crossed design, the same thresholds everywhere).
# Cells with fewer than min_size sites are merged with their nearest cell, and
# further merging keeps the number of types <= k_max.
#
# One-dimensional k-means is invariant to linear rescaling, so the thresholds
# do not depend on the units of X.
#
# Returns a list:
#   cl        : type labels 1..k
#   n_classes : number of classes per variable (1 = variable could not be split)
#   n_cells   : number of occupied cells before merging

crossed_thresholds <- function(X, classes = 2:3, min_size = 3L, k_max = 10L) {
        X <- as.matrix(X)
        n_classes <- setNames(rep(1L, ncol(X)), colnames(X))
        labs      <- matrix(1L, nrow(X), ncol(X))
        for (j in seq_len(ncol(X))) {
                x   <- X[, j]
                c_j <- min(classes[sample.int(length(classes), 1L)], length(unique(x)))
                if (c_j < 2L) next                                     # constant: no threshold
                km  <- tryCatch(kmeans(x, centers = c_j, nstart = 10), error = function(e) NULL)
                if (is.null(km)) next
                labs[, j]    <- rank(km$centers)[km$cluster]           # classes ordered along x
                n_classes[j] <- c_j
        }
        cells <- as.integer(factor(apply(labs, 1, paste, collapse = ":")))
        list(cl        = merge_small_types(X, cells, min_size, k_max),
             n_classes = n_classes,
             n_cells   = max(cells))
}


# ------------------------------------------------------------------------------
# fuzzy_classification(): fuzzy counterpart of a typology
# ------------------------------------------------------------------------------
# Fuzzy c-means (Bezdek 1981) with k clusters and random starts. It receives the
# same information as the hard classification of the same family, and nothing
# derived from the hard result:
#   kmeans, crossed_thresholds : the variables in S
#   spatial_ward               : S and the site coordinates, combined with the
#                                same alpha (spatial_feature_space())
# The fuzzy clusters are therefore in arbitrary order and are not matched to the
# hard types.
#
# X      : variables in S (standardised)
# k      : number of clusters (= number of types of the hard typology)
# m      : fuzziness exponent
# nstart : random starts; the run with the lowest objective is kept
# geo_xy, alpha : only for the spatially constrained family

fuzzy_classification <- function(X, k, m = 1.5, nstart = 10L, geo_xy = NULL, alpha = NULL) {
        if (!is.null(geo_xy)) X <- spatial_feature_space(X, geo_xy, alpha)
        vegclust::vegclust(as.matrix(X), mobileCenters = k, method = "FCM", m = m, nstart = nstart)
}


# ==============================================================================
# 2. ALGORITHM-INDEPENDENT DESCRIPTORS
# ==============================================================================
# All descriptors are functions of the labels and a fixed reference space X
# (in the simulations: all environmental variables of the model, standardised
# with the model's scaling). They do not depend on how the typology was made.

# ------------------------------------------------------------------------------
# type_evenness(): Pielou evenness of type sizes, H / ln(k), in [0, 1]
# ------------------------------------------------------------------------------
type_evenness <- function(cl) {
        p <- tabulate(as.integer(factor(cl))) / length(cl)
        if (length(p) < 2L) return(NA_real_)
        -sum(p * log(p)) / log(length(p))
}


# ------------------------------------------------------------------------------
# knn_index() / knn_kappa(): cross-validated classifiability
# ------------------------------------------------------------------------------
# Leave-one-out k-nearest-neighbour classification of the types from the
# reference space: each site is assigned the majority type of its k_nn nearest
# other sites (ties go to the nearer neighbours). Classifiability is Cohen's
# kappa between assigned and true types, which corrects for unequal type sizes.
# kNN makes no assumption about the shape of the types, so threshold and
# spatial boundaries are handled as well as compact clusters.
#
# knn_index() depends only on the distances, so for a fixed environment it can
# be computed once and reused for any number of labellings.
# k_nn should not exceed the minimum type size, otherwise the smallest types
# can never be recovered.

knn_index <- function(Dm, k_nn = 3L) {
        Dm <- as.matrix(Dm)
        diag(Dm) <- Inf                                                 # leave one out
        idx <- apply(Dm, 1, function(d) order(d)[seq_len(k_nn)])
        matrix(idx, nrow = nrow(Dm), ncol = k_nn, byrow = TRUE)         # row i: neighbours of site i
}

knn_kappa <- function(nn, cl) {
        cl  <- as.integer(factor(cl))
        K   <- max(cl)
        n   <- length(cl)
        lab <- matrix(cl[nn], nrow = n)                                 # neighbour labels
        w   <- 1 + rev(seq_len(ncol(nn))) * 1e-6                        # nearer neighbours win ties
        votes <- matrix(vapply(seq_len(K), function(t) as.vector((lab == t) %*% w), numeric(n)),
                        nrow = n)
        pred <- max.col(votes, ties.method = "first")
        p_o  <- mean(pred == cl)
        p_e  <- sum(tabulate(cl, K) * tabulate(pred, K)) / n^2
        if (p_e >= 1) return(NA_real_)
        (p_o - p_e) / (1 - p_e)
}


# ------------------------------------------------------------------------------
# typology_descriptors(): all descriptors of one typology in one reference space
# ------------------------------------------------------------------------------
# X    : n x p matrix, the reference space
# cl   : type labels
# Dm   : optional n x n distance matrix of X (Euclidean); computed if NULL
# nn   : optional knn_index() of Dm; computed if NULL
# k_nn : neighbours for the classifiability
#
# Returns a named numeric vector:
#   n_types       number of types
#   type_evenness Pielou evenness of type sizes
#   env_r2        PERMANOVA R^2 (Anderson 2001): share of the variation in X
#                 explained by the types, SS_between / SS_total. With Euclidean
#                 distances this equals adonis2()'s R^2 exactly.
#   env_ch        Calinski-Harabasz index, (SS_B / (k - 1)) / (SS_W / (n - k)),
#                 i.e. env_r2 corrected for k. Equals the PERMANOVA pseudo-F.
#   disp_mean     mean distance of sites to their type centroid (within-type
#                 dispersion)
#   disp_F        PERMDISP F statistic (Anderson 2006): one-way ANOVA F of the
#                 distances to centroid among types, i.e. heterogeneity of
#                 dispersion. Equals anova(betadisper(type = "centroid")) with
#                 Euclidean distances. Large values flag that env_r2 / env_ch
#                 are partly driven by dispersion differences.
#   knn_kappa     cross-validated classifiability (see above)
#   asw           average silhouette width, for comparison with the literature
#                 only: it assumes compact, convex types
#
# Permutation p-values are deliberately not computed: they scale with n and
# would dominate the run time without adding a descriptor.

typology_descriptors <- function(X, cl, Dm = NULL, nn = NULL, k_nn = 3L) {
        X     <- as.matrix(X)
        cl    <- as.integer(factor(cl))
        n     <- nrow(X)
        k     <- max(cl)
        sizes <- tabulate(cl, k)

        # Sum-of-squares decomposition (PERMANOVA R^2, Calinski-Harabasz)
        C     <- rowsum(X, cl) / sizes                                   # type centroids
        dev_w <- X - C[cl, , drop = FALSE]                               # within-type deviations
        ss_t  <- sum(sweep(X, 2, colMeans(X))^2)
        ss_w  <- sum(dev_w^2)
        ss_b  <- ss_t - ss_w

        # Distances to own centroid and their between-type ANOVA (PERMDISP)
        z      <- sqrt(rowSums(dev_w^2))
        z_type <- as.vector(rowsum(z, cl)) / sizes
        f_disp <- (sum(sizes * (z_type - mean(z))^2) / (k - 1)) /
                  (sum((z - z_type[cl])^2) / (n - k))

        # Distance-based descriptors
        if (is.null(Dm)) Dm <- as.matrix(dist(X))
        if (is.null(nn)) nn <- knn_index(Dm, k_nn)

        c(n_types       = k,
          type_evenness = type_evenness(cl),
          env_r2        = ss_b / ss_t,
          env_ch        = (ss_b / (k - 1)) / (ss_w / (n - k)),
          disp_mean     = mean(z),
          disp_F        = f_disp,
          knn_kappa     = knn_kappa(nn, cl),
          asw           = mean(cluster::silhouette(cl, dmatrix = as.matrix(Dm))[, 3]))
}
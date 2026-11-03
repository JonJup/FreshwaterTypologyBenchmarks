################################################################################
# Script Name:        exact_zeta.R
# Description:        Exact (expected) zeta diversity for orders 1..k from species
#                     occupancies:
#                         zeta_i = sum_s choose(o_s, i) / choose(N, i)
#                     o_s = number of sites occupied by species s, N = number of sites.
#                     This is the same computation as zetadiv::Zeta.order.ex()
#                     (occupancy-based formula, computed on the log scale via lchoose),
#                     without the variance matrix and the decline-model fits of
#                     Zeta.decline.ex().
################################################################################

exact_zeta <- function(pa, orders, rescale = TRUE) {
        pa  <- as.matrix(pa) > 0
        N   <- nrow(pa)
        if (max(orders) > N) stop("Zeta order ", max(orders), " exceeds number of sites (", N, ").")
        occ <- colSums(pa)
        occ <- occ[occ > 0]
        # no species at all -> zeta undefined
        if (length(occ) == 0) return(rep(NA_real_, length(orders)))
        z <- vapply(orders, function(i) sum(exp(lchoose(occ, i) - lchoose(N, i))), numeric(1))
        if (rescale) z <- z / z[1]
        z
}
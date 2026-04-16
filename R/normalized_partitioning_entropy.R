#' Calculate the normalized partitioning entropy
#'
#' @param U Membership matrix (n x c) where n = number of data points, c = number of clusters; or object of class vegclust from which membership matrix can be extracted
#' @return y normalized partitioning entropy (lower is better)


NPE <- function(U){
        if (class(U) == "vegclust") U <- U$memb
        n <- nrow(U)
        y <- -sum(U * log(U + 1e-10)) / n
        y <- y / log(ncol(U))
        y
}
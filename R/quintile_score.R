# This funciton assigns a score between 1 and 5 for each metric based on the quintile it is in.
quintile_score <- function(x) {
        y <- findInterval(x, quantile(x, probs = seq(0, 1, 0.2))) 
        y[which(y == 6)] <- 5
        return(y)
}
return(quintile_score)
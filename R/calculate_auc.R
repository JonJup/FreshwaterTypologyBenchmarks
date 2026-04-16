calculate_auc <- function(x, y) {
        # Sort x and y together (in case they're not already ordered)
        ord <- order(x)
        x <- x[ord]
        y <- y[ord]
        
        # Calculate width of each trapezoid
        dx <- diff(x)
        
        # Calculate mean height of each trapezoid
        mean_height <- (y[-1] + y[-length(y)]) / 2
        
        # Calculate area as sum of trapezoid areas
        area <- sum(dx * mean_height)
        
        return(area)
}
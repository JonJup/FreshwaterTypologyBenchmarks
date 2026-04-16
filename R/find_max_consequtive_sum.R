find_max_consecutive_sum <- function(x, window_size = 3) {
        # Input validation
        if (!is.numeric(x)) stop("Input must be numeric")
        if (length(x) < window_size) stop("Input length must be >= window_size")
        if (window_size < 1) stop("Window size must be positive")
        
        # Create all possible windows
        n_windows <- length(x) - window_size + 1
        window_sums <- numeric(n_windows)
        window_indices <- list()
        
        # Calculate sum for each window
        for (i in 1:n_windows) {
                current_window <- x[i:(i + window_size - 1)]
                window_sums[i] <- sum(current_window)
                window_indices[[i]] <- i:(i + window_size - 1)
        }
        
        # Find the maximum sum and its location
        max_index <- which.max(window_sums)
        
        # Return results as a list
        return(list(
                max_sum = window_sums[max_index],
                values = x[window_indices[[max_index]]],
                positions = window_indices[[max_index]],
                all_sums = window_sums
        ))
}
return(find_max_consecutive_sum)
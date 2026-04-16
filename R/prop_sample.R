prop_sample <- function(x, N, string = F) {
        # Calculate original proportions
        props <- table(x) / length(x)
        
        # Calculate how many of each number we need
        # Using round to ensure we get whole numbers
        n_each <- round(N * props)
        
        # Make sure we get exactly N samples (might be off by 1-2 due to rounding)
        while(sum(n_each) != N) {
                if(sum(n_each) < N) {
                        # Add 1 to the value with the biggest difference from target
                        target_props <- props * N
                        diff <- target_props - n_each
                        n_each[which.max(diff)] <- n_each[which.max(diff)] + 1
                } else {
                        # Subtract 1 from the value with the smallest difference from target
                        target_props <- props * N
                        diff <- n_each - target_props
                        n_each[which.max(diff)] <- n_each[which.max(diff)] - 1
                }
        }
        # Sample from each group
        result <- numeric(0)
        for(num in names(n_each)) {
                if (!string){
                positions <- which(x == as.numeric(num))
                result <- c(result, sample(positions, n_each[num]))
                } else {
                        positions <- which(x == num)
                        result <- c(result, sample(positions, n_each[num]))     
                }
        }
        
        # Return sorted positions
        sort(result)
}
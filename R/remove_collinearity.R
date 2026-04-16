remove_collinearity_vif <- function(env_df, threshold = 10) {
        # Ensure data.table is available (it is in your container)
        requireNamespace("data.table", quietly = TRUE)
        
        # 1. Prepare data for VIF calculation
        # We convert to a standard data.frame to ensure compatibility with 'lm'
        # and handle potential 'sf' objects by dropping geometry for the calculation only.
        df_calc <- as.data.frame(env_df)
        
        # Identify numeric columns (VIF only applies to numeric predictors)
        # This automatically excludes ID columns, character strings, or sf geometry columns
        numeric_cols <- sapply(df_calc, is.numeric)
        
        # We will perform stepwise selection on these variables
        active_vars <- names(df_calc)[numeric_cols]
        
        # 2. Iterative VIF Stepwise Elimination
        repeat {
                # If we have fewer than 2 variables, we can't calculate VIF (requires at least one predictor)
                if (length(active_vars) < 2) break
                
                # Placeholder for VIF values
                vifs <- numeric(length(active_vars))
                names(vifs) <- active_vars
                
                # Calculate VIF for each active variable
                # VIF_i = 1 / (1 - R_squared_i)
                for (var in active_vars) {
                        # Formula: current_var ~ all_other_active_vars
                        form <- as.formula(paste(var, "~ ."))
                        
                        # Subset data to current active variables
                        sub_dat <- df_calc[, active_vars, drop = FALSE]
                        
                        # Run linear model
                        # We use tryCatch to handle perfect collinearity or singular matrices gracefully
                        tryCatch({
                                m <- lm(form, data = sub_dat)
                                r2 <- summary(m)$r.squared
                                
                                # Handle perfect correlation (R2=1 implies infinite VIF)
                                if (r2 == 1) {
                                        vifs[var] <- Inf 
                                } else {
                                        vifs[var] <- 1 / (1 - r2)
                                }
                        }, error = function(e) {
                                # If the model fails (e.g., singular matrix), assign Inf to force removal
                                vifs[var] <- Inf
                        })
                }
                
                # 3. Check Threshold and Remove
                max_vif <- max(vifs, na.rm = TRUE)
                
                if (max_vif > threshold) {
                        # Identify variable with the highest VIF
                        exclude_var <- names(vifs)[which.max(vifs)]
                        
                        # Remove it from the active set
                        active_vars <- setdiff(active_vars, exclude_var)
                } else {
                        # If max VIF is below threshold, we are done
                        break
                }
        }
        
        # 4. Final list of variables to keep
        # This includes the survivors of the VIF process + any non-numeric columns (like IDs or Geometry)
        vars_to_keep <- c(active_vars, names(df_calc)[!numeric_cols])
        
        # 5. Return subset (Handling data.table vs data.frame)
        if (data.table::is.data.table(env_df)) {
                return(env_df[, vars_to_keep, with = FALSE])
        } else {
                return(env_df[, vars_to_keep, drop = FALSE])
        }
}
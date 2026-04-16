group_sites_by_cluster <- function(scenario) {
  # Find all unique cluster types
  unique_types <- sort(unique(scenario))
  
  # Create a list to store indices for each type
  result <- list()
  
  # For each unique type, find all indices
  for (i in seq_along(unique_types)) {
    type <- unique_types[i]
    indices <- which(scenario == type)
    result[[i]] <- indices
  }
  
  return(result)
}

balance_clusters <- function(data, clusters, min_size) {
        # Count members in each cluster
        cluster_counts <- table(clusters)
        small_clusters <- which(cluster_counts < min_size)
        
        # Identify clusters that need more members
        for (small_cluster in small_clusters) {
                needed <- min_size - cluster_counts[small_cluster]
                
                # Find cluster centers
                centers <- aggregate(data, by = list(clusters), FUN = mean)
                
                # For each small cluster, find closest points from large clusters
                large_clusters <- which(cluster_counts > min_size + needed)
                if (length(large_clusters) == 0) next
                
                # Points in large clusters
                large_indices <- which(clusters %in% large_clusters)
                
                # Calculate distances to the small cluster center
                small_center <- centers[small_cluster, -1]
                distances <- apply(data[large_indices, ], 1, function(x) 
                        sqrt(sum((x - small_center)^2)))
                
                # Find the closest points and reassign them
                closest <- large_indices[order(distances)][1:needed]
                clusters[closest] <- small_cluster
        }
        
        return(clusters)
}
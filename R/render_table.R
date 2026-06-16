################################################################################*
# Description:        Custom function that takes output of evaluation functions 
#                     and create harmonized table with relevant metadata. 
################################################################################*

if(!require(data.table)) library(data.table)

render_table <- function(x, variable, rep=5){
        
        rep =  length(x)/length(i.file$number_of_clusters)
        
        x2  <- data.table(
                scheme_id = rep(i.ss$scheme_id, each = length(x)),
                value     = x, 
                n_types   = rep(i.file$number_of_clusters, each = rep),
                metric                = variable, 
                variables             = rep(i.file$number_of_variables, each = rep),
                quality               = rep(i.file$quality, each = rep),
                compactness           = rep(i.file$contraction_points, each = rep),
                separation            = rep(i.file$contraction_centroids,each = rep),
                env_asw               = rep(i.file$asw,each = rep),
                importance            = rep(i.file$variable_importance,each = rep),
                fuzzy_npe             = rep(i.file$npe,each = rep)
        )
        return(x2)
}
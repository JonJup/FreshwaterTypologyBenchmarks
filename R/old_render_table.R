################################################################################
# Script Name:        render_table.R
# Description:        Custom function. Called from evaluate_typologies. 
#                    takes output of evaluation functions and create harmonized table with relevant metadata. 
#
# Author:             Jonathan Jupke
# Date Created:       2025-09-15
# Last Modified:      2025-09-15
#
# R Version:          R 4.5.1
# Required Packages:  data.table
#
# Notes:              Any notes or assumptions
################################################################################

if(!require(data.table)) library(data.table)

render_table <- function(x, variable,rep=5){
        x2  <- data.table(
                scheme_id             = rep(i.ss$scheme_id, each = rep),
                value                 = x, 
                #model                 = model_names[i],
                n_types               = rep(i.file$number_of_clusters, each = rep),
                metric                = variable, 
                variables             = rep(i.file$number_of_variables, each = 5),
                compactness           = rep(i.file$contraction_points, each = 5),
                separation            = rep(i.file$contraction_centroids,each = 5),
                env_asw               = rep(i.file$asw,each = 5),
                importance            = rep(i.file$variable_importance,each = 5),
                fuzzy_npe             = rep(i.file$npe,each = 5)
        )
        return(x2)
}

## alternative version; copied from setup script on the 14.08.25
# render_table <- function(x, variable){
#         x2  <- data.table(
#                 value                 = x, 
#                 simulation            = rep(1:max.ss, times = neval/max.ss),
#                 q                     = rep(1:max.q, each  = neval/max.q),
#                 cluster_adjustment    = rep(1:within.q, times = max.q),
#                 metric                = variable, 
#                 variables             = rep(j.n.variables, each = neval/max.q),
#                 contraction_points    = j.contraction.points,
#                 contraction_centroids = j.contraction.centroids,
#                 env_asw               = j.env.asw,
#                 importance            = rep(j.importance, each = neval/max.q)
#         )
#         return(x2)
# }
# return(render_table)
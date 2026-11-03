################################################################################
# Script Name:        render_table.R
# Description:        Custom functions. Called from evaluate_typologies_array.R.
#                     build_meta():   one metadata row per simulated community, keyed by ID
#                     render_table(): long-format metric table keyed by the same ID
#
# Author:             Jonathan Jupke
# Date Created:       2025-09-15
# Last Modified:      2026-09-30
#
# R Version:          R 4.5.1
# Required Packages:  data.table
#
# Notes:              ID = position of the community in unlist(file$data, recursive = FALSE).
#                     build_meta() is the ONLY place that relies on row order.
################################################################################

library(data.table)

build_meta <- function(file, ss, env_id, rep_id) {
        meta <- file$environments[file$community_predictors, on = "production_index"]
        if (nrow(meta) != length(env_id)) {
                stop("build_meta: ", nrow(meta), " metadata rows but ", length(env_id),
                     " communities. Row order cannot be trusted.")
        }
        meta[, `:=`(ID = seq_len(.N), env_id = env_id, rep_id = rep_id, scheme_id = ss$scheme_id)]
        setcolorder(meta, c("ID", "scheme_id", "env_id", "rep_id"))
        meta[]
}

render_table <- function(x, variable, id = seq_along(x)) {
        if (length(x) != length(id)) {
                stop("render_table(", variable, "): length(x) = ", length(x),
                     " but length(id) = ", length(id))
        }
        data.table(ID = id, metric = variable, value = as.numeric(x))
}
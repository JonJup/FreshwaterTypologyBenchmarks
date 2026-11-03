library(data.table)

run_D <- "//turso-fs.cs.helsinki.fi/turso/home/jonjupke/pulseD/data/008_qrf/"   # where the *_quantile_importance.rds files are
filesD   <- list.files(run_D, pattern = "_quantile_importance\\.rds$", full.names = TRUE)
run_F <- "//turso-fs.cs.helsinki.fi/turso/home/jonjupke/pulseF/data/008_qrf/"   # where the *_quantile_importance.rds files are
filesF   <- list.files(run_F, pattern = "_quantile_importance\\.rds$", full.names = TRUE)
run_I <- "//turso-fs.cs.helsinki.fi/turso/home/jonjupke/pulseI/data/008_qrf/"   # where the *_quantile_importance.rds files are
filesI   <- list.files(run_I, pattern = "_quantile_importance\\.rds$", full.names = TRUE)
run_M <- "//turso-fs.cs.helsinki.fi/turso/home/jonjupke/pulseM/data/008_qrf/"   # where the *_quantile_importance.rds files are
filesM   <- list.files(run_M, pattern = "_quantile_importance\\.rds$", full.names = TRUE)

files <- c(filesD, filesF, filesI, filesM)


bands <- list(lower  = c(0.05, 0.10),
              centre = c(0.40, 0.45, 0.50, 0.55, 0.60),
              upper  = c(0.90, 0.95))

# Relative importance on losses pooled within a band of tau levels
band_rel <- function(L, taus) {
        cols <- paste0("pinball_", taus)          # same naming as group_losses() in 11_qrf_fit.R
        stopifnot(all(cols %in% colnames(L$G_base)))
        base <- sum(L$G_base[, cols])
        vapply(L$G_perm, function(Gp) sum(Gp[, cols]) / base - 1, numeric(1))
}

tail_imp <- rbindlist(lapply(files, function(f) {
        x  <- readRDS(f)
        id <- sub("_quantile_importance\\.rds$", "", basename(f))   # "<taxon>__<metric>"
        rbindlist(lapply(names(bands), function(b) {
                rel <- band_rel(x$losses_cv, bands[[b]])
                data.table(run = id, band = b, variable = names(rel), rel = rel,
                           share = 100 * pmax(rel, 0) / sum(pmax(rel, 0)))
        }))
}))
tail_imp[, c("taxon", "metric") := tstrsplit(run, "__", fixed = TRUE)]

key_vars <- c("vp_env_pure", "env_r2_after", "vp_space_pure", "n_vars")

# 1. Median share across the 24 models, per band
med <- dcast(tail_imp[variable %in% key_vars, .(share = median(share)), by = .(variable, band)],
             variable ~ band, value.var = "share")
print(med[, .(variable, lower = round(lower, 1), centre = round(centre, 1), upper = round(upper, 1))])

# 2. In how many models is the upper tail above / below the centre?
w <- dcast(tail_imp[variable %in% key_vars], run + variable ~ band, value.var = "share")
print(w[, .(upper_gt_centre = sum(upper > centre), lower_gt_centre = sum(lower > centre),
            n = .N), by = variable])

# 3. The logic check: which predictor is largest in each band, per model?
top <- tail_imp[, .SD[which.max(share)], by = .(run, band)]
print(dcast(top[, .N, by = .(band, variable)], variable ~ band, value.var = "N", fill = 0))
print(top[band == "lower" & variable != "vp_env_pure", .(run, variable, share = round(share, 1))])

# 4. Per-metric medians across taxa (the coloured lines of Figure S28), same aggregation
per_tau <- rbindlist(lapply(files, function(f) {
        x <- readRDS(f); id <- sub("_quantile_importance\\.rds$", "", basename(f))
        d <- x$results[measure == "pinball", .(variable, tau, rel)]
        d[, share := 100 * pmax(rel, 0) / sum(pmax(rel, 0)), by = tau][, run := id][]
}))
per_tau[, c("taxon", "metric") := tstrsplit(run, "__", fixed = TRUE)]
fig_s28 <- per_tau[variable %in% key_vars,
                   .(share = median(share)), by = .(metric, variable, tau)]
fig_s28_all <- per_tau[variable %in% key_vars, .(share = median(share)), by = .(variable, tau)]
fwrite(tail_imp, "data/supplement/S7_importance_by_band.csv")
fwrite(per_tau, "data/supplement/S7_importance_by_tau.csv")



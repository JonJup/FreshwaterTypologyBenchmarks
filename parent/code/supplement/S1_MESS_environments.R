## mess_coverage.R
## 
## 
## Environmental coverage of European river catchments (EU-Hydro) by the sampling sites.
##   1. MESS (Elith et al. 2010): univariate extrapolation; MoD = most dissimilar variable.
##   2. ExDet NT2 (Mesgaran et al. 2014): novel combinations within univariate ranges.
##   3. Stressor exclusion: coverage by pressure-filtered site sets vs. random removal
##      of the same number of sites.
##
## Sites and catchments must carry the same natural variables, computed identically
## (e.g. each site takes the values of the catchment it lies in). Keep pressure
## variables out of nat_vars, otherwise the exclusion comparison is circular.


# setup -------------------------------------------------------------------
library(data.table)
library(sf)
library(mapview)

set.seed(42)

## ---- settings -----------------------------------------------------------------
areas        <- "E://Arbeit/data/river_network/eu_hydro_dem/processedData/areas/"
path_eu      <- "E://Arbeit/data/river_network/eu_hydro_dem/processedData/all_catchments_w_environment.rds"  # one row per EU-Hydro catchment
path_sites   <- "data/new_environment.rds"          # one row per site x organism group
out_dir      <- "data/environmental_coverage"
pressures    <- "E://Arbeit/data/my_stressor_db/processedData/PCA_stressors.rds"

nat_vars     <- c("slope", "lake_area", "soil_oc", "soil_pH", 
                  "inundation_depth", "flooded_area", "spi_mean", "area_calcareous", 
                  "area_siliceous", "area_sediment", "glacial_area", "bioclim05", "bioclim06", 
                  "bioclim13", "bioclim14", "Rfactor_max", "Rfactor_avg", "Rfactor_min", 
                  "roughness", "elevation", "mean_snow_equivalent", "mean_discharge", 
                  "max_discharge", "min_discharge", "groundwater_table", 
                  "upstream_catchment_area", "saturated_soil_water_content", "VBM_mean") # already transformed (log area etc.)

id_col       <- "ID"
site_col     <- "siteID"
group_col    <- "taxon"
pressure_col <- "dist_ref" 

keep_q       <- c(1, 0.75, 0.5, 0.25, 0.1)     # keep sites with pressure <= this quantile
n_null       <- 199
n_cores      <- 1


## custom functions --------------------------------------------------------
### ---- MESS --------------------------------------------------------------
## Similarity of values p to reference sample ref for one variable.
## Negative = outside the reference range.
mess_var <- function(p, ref) {
        ref <- sort(ref)
        n   <- length(ref)
        mn  <- ref[1L]; mx <- ref[n]; rng <- mx - mn
        if (rng == 0) return(ifelse(p == mn, 0, -100))
        f  <- 100 * findInterval(p, ref, left.open = TRUE) / n  # % of ref values < p
        s  <- 2 * pmin(f, 100 - f)                               # inside the range
        i0 <- which(f == 0)
        s[i0] <- 100 * (p[i0] - mn) / rng # at or below min
        i1 <- which(f == 100)
        s[i1] <- 100 * (mx - p[i1]) / rng # above max
        s
}

mess <- function(new, ref, vars) {
        S <- vapply(vars, function(v) mess_var(new[[v]], ref[[v]]), numeric(nrow(new)))
        S <- matrix(S, ncol = length(vars), dimnames = list(NULL, vars))
        j <- max.col(-S, ties.method = "first")
        list(mess = S[cbind(seq_len(nrow(S)), j)], mod = vars[j])
}

### ---- ExDet NT2 -------------------------------------------------------
## Mahalanobis distance to the reference centroid, relative to the largest distance
## within the reference sample. > 1 = combination outside the sampled space.
nt2 <- function(new, ref, vars) {
        X     <- as.matrix(ref[, vars, with = FALSE])
        mu    <- colMeans(X)
        Sinv  <- tryCatch(solve(cov(X)), error = function(e) MASS::ginv(cov(X)))
        d_ref <- sqrt(mahalanobis(X, mu, Sinv, inverted = TRUE))
        d_new <- sqrt(mahalanobis(as.matrix(new[, vars, with = FALSE]), mu, Sinv, inverted = TRUE))
        d_new / quantile(d_ref, .95)
}

### ---- classification and shares -------------------------------------------
classify <- function(eu, ref, vars) {
        m  <- mess(eu, ref, vars)
        n2 <- nt2(eu, ref, vars)
        data.table(mess = m$mess, mod = m$mod, nt2 = n2,
                   class = fifelse(m$mess < 0, "univariate",
                                   fifelse(n2 > 1, "combinatorial", "covered")))
}

shares <- function(cl,w) {
        c(within_ranges = sum(w[cl$mess >= 0]) / sum(w),
          covered       = sum(w[cl$class == "covered"]) / sum(w))
}

null_summary <- function(null, obs, prefix) {
        x <- null[, prefix]
        setNames(list(mean(x), quantile(x, 0.025, names = FALSE),
                      quantile(x, 0.975, names = FALSE),
                      (sum(x <= obs[[prefix]]) + 1) / (length(x) + 1)),
                 paste0("null_", prefix, c("_mean", "_lo", "_hi", "_p")))
}




# Load Data ---------------------------------------------------------------
eu    <- setDT(readRDS(path_eu))
sites <- readRDS(path_sites)
for (i in seq_along(sites)) sites[[i]]$taxon <- names(sites)[i]
sites <- rbindlist(sites)
stessors <- readRDS(pressures)
## add a site ID 
sites[, siteID := .GRP, by = c("x.coord", "y.coord")]

# add area to EU as weight for later 
areas2 <- list.files(areas, full.names = T)
areas3 <- lapply(areas2, readRDS)
areas3 <- rbindlist(areas3)
areas3$area <- units::drop_units(areas3$area)

eu <- areas3[eu, on = "ID"]


n0 <- nrow(eu)
eu <- eu[complete.cases(eu[, ..nat_vars])]
sites <- sites[complete.cases(sites[, ..nat_vars])]
message(sprintf("%d of %d catchments with complete natural variables", nrow(eu), n0))

sites2 <- merge(sites, stessors, by = "ID",)


# ## inspect missing ones 
# sites3 <- sites$ID[which(!sites$ID %in% sites2$ID)]
# sites3 <- sites[ID %in% sites3]
# site3 <- st_as_sf(sites3, coords = c("x.coord", "y.coord"), crs = 3035)
# mapview(site3)
# sites4 <- copy(sites3)


## one absolute criterion for all groups: quantiles of pressure across unique sites
thr <- quantile(unique(sites2, by = site_col)[[pressure_col]], keep_q,
                na.rm = TRUE, names = FALSE)

groups <- c(sort(unique(as.character(sites[[group_col]]))))
min_n  <- length(nat_vars) + 2
res    <- list()
w <- eu$area
## ---- main loop ----------------------------------------------------------------
for (g in groups) {
      
        full <- sites2[get(group_col) == g]
        cl_full   <- classify(eu, full, nat_vars)
        cl_strict <- NULL
        ## per-catchment output, ready to join to the EU-Hydro geometry for maps
        pc <- data.table(id = eu[[id_col]], w = w,
                         mess = cl_full$mess,
                         mod = cl_full$mod,
                         nt2 = cl_full$nt2, 
                         class_full = cl_full$class)
        if (!is.null(cl_strict)) {
                pc[, class_strict := cl_strict$class]
                pc[, lost := class_full == "covered" & class_strict != "covered"]
        }
        setnames(pc, "id", id_col)
        saveRDS(pc, file.path(out_dir, sprintf("coverage_per_catchment_%s.rds", g)))
        
        ## which variable limits coverage under the full data (share of total weight)
        mod_tab <- pc[mess < 0, .(share_of_total = sum(w) / sum(pc$w)), by = mod][
                order(-share_of_total)]
        fwrite(mod_tab, file.path(out_dir, sprintf("most_dissimilar_variable_%s.csv", g)))
}
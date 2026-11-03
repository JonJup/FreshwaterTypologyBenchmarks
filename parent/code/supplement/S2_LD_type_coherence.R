## Least-disturbed type coherence test (Appendix S2.2.3): French river fish, 2016
##
## Builds least-disturbed (LD) subsets of increasing stringency from the French
## fish data of the application example. For every LD subset it also draws
## N_RANDOM random subsets from all sites, with the same number of sites per type.
## Every subset runs through pulse_predictors(), so each gets its own predictors
## and observed coherence. The batch file goes into the app's batch tab, which
## returns each subset's percentile within its own benchmark. Part 2 places the
## LD percentile of each level and metric within the null distribution of its
## matched random subsets.
##
## Part 1 writes (OUT_DIR)
##   ld_levels.csv         sites and types per level, before and after dropping small types
##   ld_design.csv         one row per subset: level, kind, replicate, sizes, pressure
##   ld_subsets.rds        site membership of every subset
##   ld_batch.csv          batch file for the app (subset x metric)
## Part 2 (runs once APP_RESULTS exists; reads only the files above)
##   ld_comparison.csv     LD percentile vs. null distribution, per level x metric
##   fig_ld_coherence.png

# setup -------------------------------------------------------------------
library(sf)
library(data.table)
library(parallel)
library(recipes)

# settings ----------------------------------------------------------------
PULSE_SCRIPT   <- "SHINY APP/prepare_for_app.R"   
OUT_DIR        <- "data/supplement/S2_ld_case_study/"

TYPOLOGY_FILE  <- "E://Arbeit/data/typology_systems/broad_river_types/rawData/m_river_fec_broad_type.shp"
FISH_FILE      <- "data/000_biota/02_fish_w_environment.rds"
PRESSURE_FILE  <- "E://Arbeit/data/my_stressor_db/processedData/PCA_stressors.rds"

RUN_PREDICTORS <- TRUE   # FALSE skips the expensive pulse_predictors() step

LEVELS         <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.50)   # LD: dist_ref <= this quantile
THRESHOLD_BASE <- "catchments"   # "catchments": quantiles of dist_ref over French catchments
# "sites":      quantiles over the sampled sites
N_RANDOM       <- 199L
MIN_TYPE_SIZE  <- 3L     # types with fewer LD sites are dropped from that level (training typologies had >= 3)
MIN_SITES      <- 100L   # smallest subset in the BGM training data; levels with fewer LD sites are skipped
N_CORES        <- max(1L, detectCores() - 1L)
CHUNK          <- 10L    # subsets per pulse_predictors() call
SEED           <- 20261022L

TYPE_COL       <- "m_btype12"
TYPOLOGY_VARS  <- c("elevation", "upstream_catchment_area",
                    "area_calcareous", "area_siliceous", "area_sediment")

# Optional second LD criterion of S2.2.1 (no single pressure above its limit):
# a function of the pressures table that returns TRUE for admissible catchments.
LD_EXTRA <- NULL

# IDs of the French catchments in `pressures` (used if THRESHOLD_BASE = "catchments").
# EDIT if `pressures` has no country column, e.g. by intersecting the catchments
# with a France outline.
french_catchments <- function(p) {
        cc <- intersect(c("country", "Country", "CNTR_CODE", "cntr_code", "country_code", "ISO2"), names(p))
        if (!length(cc)) stop("`pressures` has no country column: edit french_catchments() ",
                              "to return the IDs of the French catchments.")
        p$ID[p[[cc[1]]] %in% c("FR", "FRA", "France")]
}

APP_RESULTS <- file.path(OUT_DIR, "pulse_batch_all_results_2026-10-22.csv")  # EDIT: where you save the app's batch results
PCT_COL     <- "observed_percentile"                              # EDIT: their column with the benchmark percentile

# Evaluate only the definitions in PULSE_SCRIPT (assignments and library calls),
# so that its RUN block is not executed.
source_defs <- function(file, envir = globalenv()) {
        for (e in parse(file)) {
                head <- if (is.call(e)) as.character(e[[1]])[1] else ""
                if (head %in% c("<-", "=", "library", "require", "suppressPackageStartupMessages")) eval(e, envir)
        }
        invisible(NULL)
}
source_defs(PULSE_SCRIPT)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================*
# 1. SUBSETS AND BATCH FILE ----
# ==============================================================================*

# load data ---------------------------------------------------------------
broadRiverTypes <- st_read(TYPOLOGY_FILE, quiet = TRUE)
fish_w_env      <- readRDS(FISH_FILE)
pressures       <- as.data.table(drop_geometry(readRDS(PRESSURE_FILE)))
pressures[, ID := as.character(ID)]
if (anyDuplicated(pressures$ID)) stop("`pressures` has several rows per ID.")

# prepare data (as in the application example) ----------------------------
setDT(fish_w_env)
fish_w_env[, month := month(eventDate)]
fishSub <- fish_w_env[data.set == "fish_france_naiades" & eventYear == 2016 & month %in% 7:9]

community <- fishSub[, .(site = as.character(siteID), taxon = working.taxon,
                         value = as.numeric(organismQuantity > 0), date = eventDate)]

env_cols <- names(fishSub)[match("slope", names(fishSub)):match("VBM_max", names(fishSub))]
env      <- unique(fishSub, by = "siteID")[, c(list(site = as.character(siteID)), .SD), .SDcols = env_cols]

sites <- unique(fishSub, by = "siteID")[, .(site = as.character(siteID), x = x.coord, y = y.coord, EPSG)]
crs_sites <- unique(sites$EPSG)
if (length(crs_sites) != 1L) stop("Site coordinates come in several CRS.")
if (grepl("^[0-9]+$", crs_sites)) crs_sites <- as.integer(crs_sites)
sites[, EPSG := NULL]

site_catch <- unique(fishSub[, .(site = as.character(siteID), ID = as.character(ID))])
if (anyDuplicated(site_catch$site)) stop("Some sites are linked to several catchment IDs.")

# site pool: pressure, type, environment ----------------------------------
pts   <- read_sites(sites, crs = crs_sites)
types <- assign_types(broadRiverTypes, sites$site, pts, type_col = TYPE_COL)

pool <- merge(sites[, .(site)], site_catch, by = "site", all.x = TRUE)
pool[pressures, dist_ref := i.dist_ref, on = "ID"]
pool[, type := unname(types[site])]
adm_ids <- if (is.null(LD_EXTRA)) pressures$ID else pressures$ID[which(LD_EXTRA(pressures))]
pool[, admissible := ID %in% adm_ids]

say("Sites: %d | without pressure data: %d | without type: %d | without environment: %d",
    nrow(pool), sum(is.na(pool$dist_ref)), sum(is.na(pool$type)), sum(!pool$site %in% env$site))
# LD and random subsets come from the same population: sites with pressure, type and environment
pool <- pool[!is.na(dist_ref) & !is.na(type) & site %in% env$site & site %in% community$site]
say("Pool: %d sites in %d types", nrow(pool), uniqueN(pool$type))

community <- community[site %in% pool$site]
env       <- env[site %in% pool$site]
sites     <- sites[site %in% pool$site]

# thresholds and levels ---------------------------------------------------
level_names <- sprintf("q%02d", round(100 * LEVELS))

ld_thresholds <- function(base) {
        ref <- switch(base,
                      catchments = pressures[ID %in% french_catchments(pressures), dist_ref],
                      sites      = pool$dist_ref,
                      stop("THRESHOLD_BASE must be 'catchments' or 'sites'."))
        setNames(quantile(ref, LEVELS, na.rm = TRUE, names = FALSE), level_names)
}

level_table <- function(thr) {
        rbindlist(lapply(names(thr), function(lv) {
                ld  <- pool[dist_ref <= thr[[lv]] & admissible]
                n_t <- ld[, .N, by = type]
                big <- n_t[N >= MIN_TYPE_SIZE]
                data.table(level = lv, q = LEVELS[match(lv, names(thr))], threshold = thr[[lv]],
                           n_ld = nrow(ld), n_types_ld = nrow(n_t),
                           n_sites = sum(big$N), n_types = nrow(big),
                           share_of_pool = sum(big$N) / nrow(pool),
                           run = sum(big$N) >= MIN_SITES && nrow(big) >= 2L)
        }))
}

# Both bases are printed so that the choice can be checked against the subset sizes
thr_all <- list(sites      = ld_thresholds("sites"),
                catchments = tryCatch(ld_thresholds("catchments"),
                                      error = function(e) { say("  %s", conditionMessage(e)); NULL }))
for (b in names(thr_all)) if (!is.null(thr_all[[b]])) {
        say("\nLD levels, quantiles over %s:", b)
        print(level_table(thr_all[[b]]))
}

thr <- thr_all[["sites"]]
if (is.null(thr)) stop("No thresholds for THRESHOLD_BASE = '", THRESHOLD_BASE, "'.")
levels_dt <- level_table(thr)[, base := THRESHOLD_BASE]
fwrite(levels_dt, file.path(OUT_DIR, "ld_levels.csv"))
if (any(!levels_dt$run)) {
        say("\nSkipped (fewer than %d sites or 2 types after dropping types with < %d LD sites): %s",
            MIN_SITES, MIN_TYPE_SIZE, paste(levels_dt[run == FALSE, level], collapse = ", "))
}
if (!any(levels_dt$run)) stop("No level has enough least-disturbed sites.")

# LD subsets and matched random subsets ------------------------------------
by_type <- split(pool$site, pool$type)
draw <- function(n_t) rbindlist(lapply(seq_len(nrow(n_t)), function(i)
        data.table(site = sample(by_type[[n_t$type[i]]], n_t$N[i]), type = n_t$type[i])))

subsets <- list()
for (lv in levels_dt[run == TRUE, level]) {
        ld  <- pool[dist_ref <= thr[[lv]] & admissible]
        n_t <- ld[, .N, by = type][N >= MIN_TYPE_SIZE]
        subsets[[paste0(lv, "_ld")]] <- ld[type %in% n_t$type, .(site, type)]
        set.seed(SEED + match(lv, level_names))            # one stream per level
        for (r in seq_len(N_RANDOM)) subsets[[sprintf("%s_r%03d", lv, r)]] <- draw(n_t)
}

dref   <- setNames(pool$dist_ref, pool$site)
design <- rbindlist(lapply(names(subsets), function(id) {
        s  <- subsets[[id]]
        lv <- sub("_.*", "", id)
        is_ld <- grepl("_ld$", id)
        data.table(subset = id, level = lv, q = levels_dt[level == lv, q],
                   kind = if (is_ld) "ld" else "random",
                   replicate = if (is_ld) 0L else as.integer(sub(".*_r", "", id)),
                   n_sites = nrow(s), n_types = uniqueN(s$type),
                   share_ld = mean(s$site %in% subsets[[paste0(lv, "_ld")]]$site),
                   median_dist_ref = median(dref[s$site]))
}))
# fwrite(design, file.path(OUT_DIR, "ld_design.csv"))
# saveRDS(subsets, file.path(OUT_DIR, "ld_subsets.rds"))
say("\n%d subsets: %d levels x (1 LD + %d random)", length(subsets), sum(levels_dt$run), N_RANDOM)
print(design[, .(n_sites = n_sites[1], n_types = n_types[1],
                 dist_ref_ld = median_dist_ref[kind == "ld"],
                 dist_ref_random = median(median_dist_ref[kind == "random"]),
                 share_ld_in_random = mean(share_ld[kind == "random"])), by = level])

# predictors and observed coherence for every subset ------------------------
run_chunk <- function(typ) {
        tryCatch(
                pulse_predictors(community = community, env = env, typology = typ, sites = sites,
                                 typology_vars = TYPOLOGY_VARS,
                                 site_col = "site", taxon_col = "taxon", value_col = "value",
                                 date_col = "date", x_col = "x", y_col = "y", crs = crs_sites,
                                 batch_output = tempfile(fileext = ".csv"))$batch,
                error = function(e) data.table(typology = names(typ), error = conditionMessage(e)))
}

if (RUN_PREDICTORS) {
        typ_list <- lapply(subsets, as.data.frame)      # site-type tables; pulse_predictors() keeps only these sites
        chunks   <- split(typ_list, ceiling(seq_along(typ_list) / CHUNK))
        
        t1 <- system.time(invisible(capture.output(run_chunk(typ_list[1]))))[["elapsed"]]
        say("\nOne subset takes %.0f s: about %.0f min for %d subsets on %d cores",
            t1, t1 * length(typ_list) / N_CORES / 60, length(typ_list), N_CORES)
        
        if (N_CORES > 1L) {
                cl <- makeCluster(N_CORES)
                clusterCall(cl, source_defs, normalizePath(PULSE_SCRIPT))
                clusterExport(cl, c("community", "env", "sites", "crs_sites", "TYPOLOGY_VARS"))
                res <- parLapplyLB(cl, chunks, run_chunk)
                stopCluster(cl)
        } else {
                res <- lapply(chunks, function(ch) { invisible(capture.output(out <- run_chunk(ch))); out })
        }
        
        failed <- vapply(res, function(x) "error" %in% names(x), TRUE)
        if (any(failed)) {
                errs <- rbindlist(res[failed])
                fwrite(errs, file.path(OUT_DIR, "ld_errors.csv"))
                say("%d subsets failed (ld_errors.csv), first error: %s", nrow(errs), errs$error[1])
        }
        batch <- rbindlist(res[!failed])
        batch <- batch[order(match(typology, design$subset), match(metric, METRICS))]
        fwrite(batch, file.path(OUT_DIR, "ld_batch.csv"), na = "")
        say("Written: %s (%d subsets x %d metrics)", file.path(OUT_DIR, "ld_batch.csv"),
            uniqueN(batch$typology), length(METRICS))
        # Load ld_batch.csv in the app's batch tab with fish selected, and save the
        # results as APP_RESULTS.
}


APP_RESULTS <- "data/supplement/SM2_pulse_batch_all_results_LD.csv"

# ==============================================================================*
# 2. COMPARISON (after the app) ----
# ==============================================================================*
#
# For every level and metric: the LD subset's benchmark percentile, the median
# and 95% interval of the matched random subsets, and a two-sided Monte Carlo
# p-value. The levels are nested, so read them as one gradient rather than as
# independent tests.

if (file.exists(APP_RESULTS)) {
        design <- fread(file.path(OUT_DIR, "ld_design.csv"))
        base   <- fread(file.path(OUT_DIR, "ld_levels.csv"))$base[1]
        app    <- fread(APP_RESULTS)
        check_cols(app, c("typology", "metric", PCT_COL), "app results")
        d <- merge(app[, .(subset = typology, metric, value = as.numeric(get(PCT_COL)))], design, by = "subset")
        
        cmp <- d[, {
                v0 <- value[kind == "ld"]
                vr <- value[kind == "random" & !is.na(value)]
                ok <- length(v0) == 1L && !is.na(v0) && length(vr) > 0L
                .(n_sites = n_sites[kind == "ld"][1], n_types = n_types[kind == "ld"][1],
                  dist_ref_ld = median_dist_ref[kind == "ld"][1],
                  dist_ref_random = median(median_dist_ref[kind == "random"]),
                  ld = if (ok) v0 else NA_real_,
                  null_median = median(vr),
                  null_lo = quantile(vr, 0.025, names = FALSE),
                  null_hi = quantile(vr, 0.975, names = FALSE),
                  diff = if (ok) v0 - median(vr) else NA_real_,
                  p = if (ok) min(1, 2 * min(1 + sum(vr <= v0), 1 + sum(vr >= v0)) / (length(vr) + 1)) else NA_real_,
                  n_random = length(vr))
        }, by = .(level, q, metric)]
        setorder(cmp, metric, q)
        fwrite(cmp, file.path(OUT_DIR, "ld_comparison.csv"))
        print(cmp[, .(metric, level, n_sites, ld, null_median, null_lo, null_hi, p)])
        
        if (has_pkg("ggplot2")) {
                library(ggplot2)
                lab <- function(q) factor(sprintf("%g%%", 100 * q), levels = sprintf("%g%%", 100 * sort(unique(design$q))))
                cmp[, level_lab := lab(q)]
                d[, level_lab := lab(q)]
                col_ld   <- "Least-disturbed subset"
                col_null <- "Matched random subsets (median, 95% interval)"
                p_fig <- ggplot() +
                        geom_jitter(data = d[kind == "random"], aes(level_lab, value),
                                    width = 0.15, height = 0, colour = "grey80", size = 0.5) +
                        geom_pointrange(data = cmp, aes(level_lab, null_median, ymin = null_lo, ymax = null_hi,
                                                        colour = col_null)) +
                        geom_point(data = cmp, aes(level_lab, ld, colour = col_ld), size = 2.5) +
                        scale_colour_manual(values = setNames(c("#D55E00", "grey30"), c(col_ld, col_null)), name = NULL) +
                        facet_wrap(~ metric, nrow = 1) +
                        labs(x = sprintf("Least-disturbed threshold (quantile of dist_ref over %s)",
                                         if (base == "sites") "sampled sites" else "French catchments"),
                             y = "Benchmark percentile") +
                        theme_bw(base_size = 10) + theme(legend.position = "bottom")
                ggsave(file.path(OUT_DIR, "S2_fig_ld_coherence.png"), p_fig, width = 9, height = 3.8, dpi = 300)
        }
} else {
        say("\nPart 2 skipped: %s not found.", APP_RESULTS)
}
# fwrite(d, "data/ld_test/results.csv")

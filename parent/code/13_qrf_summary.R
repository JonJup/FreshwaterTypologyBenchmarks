### - QRF summary across taxonomic groups and metrics - ###
#
# Collects the outputs of qrf_fit.R (one *_run_info.rds per taxon x metric run,
# found recursively under input_dir) and writes summary tables to tab_dir,
# figures to fig_dir and everything to <output_dir>/summary_tables.rds.
# Only completed runs are used; of duplicate taxon x metric runs the most
# recently finished one is kept. Every output file is read once.
#
#   4  predictive performance: test set vs. grouped CV, calibration, PIT, pinball
#      skill by quantile level, out-of-sample performance per data set,
#      pooled vs. per-data-set weighting; figures res_2, res_3, res_4
#   5  tuning (bagging unit) and how well the test data sets span the training data sets
#   6  importance: CRPS shares, tails vs. centre, agreement of grouped CV and test set
#   7  effects: ALE summaries and additivity
#   8  data sets: concentration, LODO influence relative to the replicate-fit
#      noise floor, equal-weight sensitivity, stability of the tracked ALE classes
#   9  typology families (LOFO): extrapolation cost relative to the size-matched
#      control, calibration, PIT, extrapolation beyond the training range
#   10 overview: one row per taxon x metric


library(data.table)
library(ggplot2)


#===============================================================================*
# 1. Setup ----
#===============================================================================*

input_base   <- ""
taxa_suffix  <- c("D", "F", "I", "M")
input_ending <- "/data/008_qrf/"
input_dir    <- paste0(input_base, taxa_suffix, input_ending)
output_dir   <- "~/projects/pulse/01_wp1/MEE_REVISION_1"

tab_dir <- file.path(output_dir, "data/009_qrf_results/")
fig_dir <- file.path(output_dir, "figures")
for (d in c(tab_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## --- Settings ----
taxon_order  <- NULL
metric_order <- NULL
family_order <- NULL

pit_bins         <- 10
tail_tau         <- 0.1   # tails: tau <= tail_tau or tau >= 1 - tail_tau
center_tau       <- 0.1   # centre: |tau - 0.5| <= center_tau
# a training data set is influential if removing it shifts the test
# median predictions by more than influence_factor x the shift of a
# replicate refit (the noise floor rests on a single replicate)
influence_factor <- 2

# Predictor classes for aggregated importance; unlisted predictors -> "Other"
predictor_classes <- c(
  samples         = "Sampling effort",
  n_types         = "Typology structure",       n_vars          = "Typology structure",
  type_evenness   = "Typology structure",
  env_r2_after    = "Environmental separation", env_ch_after    = "Environmental separation",
  asw_after       = "Environmental separation",
  min_distance    = "Spatial extent",           median_distance = "Spatial extent",
  species_rank    = "Taxonomic resolution",     genus_rank      = "Taxonomic resolution",
  family_rank     = "Taxonomic resolution",     higher_rank     = "Taxonomic resolution",
  n_taxa_used     = "Taxon richness",
  vp_env_pure     = "Variation partitioning",   vp_shared       = "Variation partitioning",
  vp_space_pure   = "Variation partitioning",
  mean_doy        = "Sampling season",          circ_sd_days    = "Sampling season"
)
extra_classes <- c(value = "Outcome", any = "Any predictor")   # split report, LOFO range table

eval_levels <- c("grouped CV", "test set")
lofo_models <- c("LOFO", "control", "full")

## --- Helpers ----
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
has_rows <- function(dt) !is.null(dt) && NROW(dt) > 0
as_dt    <- function(x) if (is.data.table(x)) copy(x) else as.data.table(x)
num      <- function(x) if (is.null(x) || !length(x)) NA_real_ else as.numeric(x[1])
chr      <- function(x) if (is.null(x) || !length(x)) NA_character_ else as.character(x[1])
recode_eval <- function(x) {                          # set labels of _intervals.rds -> evaluation
  u <- unique(as.character(x))
  r <- fifelse(grepl("^test", u, ignore.case = TRUE), "test set",
               fifelse(grepl("cv", u, ignore.case = TRUE), "grouped CV", u))
  r[match(as.character(x), u)]
}

save_fig <- function(plot, name, width, height) {
  f <- file.path(fig_dir, name)
  ggsave(f, plot, width = width, height = height, limitsize = FALSE)
  cat("  saved:", f, "\n")
}
save_table <- function(dt, name) {
  if (!has_rows(dt)) return(invisible())
  f <- file.path(tab_dir, name)
  fwrite(dt, f)
  cat("  saved:", f, "\n")
}

# Histogram of PIT values per group (density 1 = calibrated) and PIT moments.
# Dispersion = 12 var(PIT): 1 = calibrated, > 1 = predictive distributions too
# narrow (U-shaped PIT), < 1 = too wide (hump). Mean != 0.5 = bias.
bin_pit <- function(d, grp) {
  br <- seq(0, 1, length.out = pit_bins + 1)
  d  <- d[, c(grp, "pit"), with = FALSE]
  d[, bin := findInterval(pit, br, rightmost.closed = TRUE, all.inside = TRUE)]
  h  <- d[, .(count = .N), by = c(grp, "bin")]
  full <- unique(d[, grp, with = FALSE])[, .(bin = seq_len(pit_bins)), by = grp]
  h  <- h[full, on = c(grp, "bin")][is.na(count), count := 0L]
  h[, `:=`(lo = br[bin], hi = br[bin + 1], mid = (br[bin] + br[bin + 1]) / 2)]
  h[, density := count / sum(count) / (hi - lo), by = grp]
  h[]
}
pit_moments <- function(d, grp) {
  d[, .(n = .N, pit_mean = mean(pit), pit_dispersion = 12 * var(pit)), by = grp]
}


#===============================================================================*
# 2. Discover runs ----
#===============================================================================*
info_files <- list.files(input_dir, pattern = "_run_info\\.rds$", recursive = TRUE, full.names = TRUE)
if (!length(info_files)) stop("No *_run_info.rds files found under ", input_dir)

runs <- rbindlist(lapply(info_files, function(f) {
  ri  <- readRDS(f)
  fin <- ri$finished
  data.table(
    prefix           = sub("_run_info\\.rds$", "", f),
    folder           = basename(dirname(f)),
    run_id           = chr(ri$run_id),
    taxon_recorded   = chr(ri$taxon),
    metric           = chr(ri$metric),
    status           = chr(ri$status),
    finished         = if (inherits(fin, "POSIXct")) fin else .POSIXct(NA_real_),
    n_rows           = num(ri$n_rows),
    n_rows_full      = num(ri$n_rows_full),
    n_datasets       = num(ri$n_datasets),
    n_train          = num(ri$n_train),
    n_datasets_train = num(ri$n_datasets_train),
    n_test           = num(ri$n_test),
    n_datasets_test  = num(ri$n_datasets_test),
    n_schemes        = num(ri$n_schemes),
    n_schemes_full   = num(ri$n_schemes_full),
    scheme_cap       = num(ri$scheme_cap),
    max_rows         = num(ri$max_rows),
    main_level       = num(ri$settings$main_level %||% 0.9),
    lodo             = as.logical(num(ri$settings$lodo)),
    lofo             = as.logical(num(ri$settings$lofo)),
    n_predictors     = length(ri$predictors),
    families         = paste(ri$families, collapse = ";"),
    mtry             = num(ri$best$mtry),
    min_n            = num(ri$best$min_n),
    bagging          = chr(ri$best$bagging),
    ranger           = chr(ri$versions$ranger)
  )
}), fill = TRUE)

runs[, taxon := stringr::str_to_title(taxon_recorded)]
runs[, run_label := paste(taxon, "|", metric)]

incomplete <- runs[is.na(status) | status != "complete"]
if (nrow(incomplete)) {
  warning(nrow(incomplete), " incomplete run(s) skipped: ",
          paste(incomplete$run_label, collapse = ", "), call. = FALSE)
}
runs <- runs[status == "complete"]
setorder(runs, taxon, metric, -finished)
dups <- runs[duplicated(runs, by = c("taxon", "metric"))]
if (nrow(dups)) {
  warning("Duplicate taxon x metric runs; keeping the most recently finished. Dropped: ",
          paste(dups$prefix, collapse = ", "), call. = FALSE)
  runs <- runs[!duplicated(runs, by = c("taxon", "metric"))]
}
if (!nrow(runs)) stop("No completed runs found.")
if (any(is.finite(runs$max_rows))) {
  warning("Runs fitted with --max-rows (testing subsample): ",
          paste(runs[is.finite(max_rows), run_label], collapse = ", "), call. = FALSE)
}

taxon_levels  <- taxon_order  %||% sort(unique(runs$taxon))
metric_levels <- metric_order %||% sort(unique(runs$metric))
if (uniqueN(runs$main_level) > 1) {
  warning("Runs differ in main_level; headline coverage is not comparable across runs.", call. = FALSE)
}

cat("\nCompleted runs:", nrow(runs), "\n")
print(runs[, .(n_metrics = .N, recorded_taxon = paste(unique(taxon_recorded), collapse = ", ")),
           by = .(folder, taxon)])
print(dcast(runs, taxon ~ metric, value.var = "n_rows", fun.aggregate = sum, fill = NA))

## --- Inventory of output files ----
components <- c(split = "_split.rds", tuning = "_tuning.rds", metrics = "_metrics.rds",
                intervals = "_intervals.rds", importance = "_quantile_importance.rds",
                ale = "_ale.rds", datasets = "_datasets.rds", lofo = "_lofo.rds")
inventory <- runs[, c(.(taxon = taxon, metric = metric, lodo = lodo, lofo_requested = lofo),
                      lapply(components, function(s) file.exists(paste0(prefix, s))))]
missing_core <- inventory[!(split & tuning & metrics & intervals & importance & ale & datasets)]
if (nrow(missing_core)) {
  warning("Completed runs with missing output files:\n",
          paste(capture.output(print(missing_core)), collapse = "\n"), call. = FALSE)
}
no_lofo <- inventory[lofo_requested %in% TRUE & !lofo]
if (nrow(no_lofo)) {
  cat("  LOFO requested but no _lofo.rds (fewer than two families or no evaluable family):",
      paste(no_lofo$taxon, no_lofo$metric, sep = " | ", collapse = ", "), "\n")
}


#===============================================================================*
# 3. Read run outputs ----
#===============================================================================*

# read_parts(): reads <prefix><suffix> of every run once and applies each
# extractor (function(x, r), r = the run's row of `runs`). Results are bound
# across runs, with taxon and metric added.

read_parts <- function(suffix, extractors) {
  out <- lapply(extractors, function(e) list())
  for (i in seq_len(nrow(runs))) {
    f <- paste0(runs$prefix[i], suffix)
    if (!file.exists(f)) next
    x <- readRDS(f)
    for (nm in names(extractors)) {
      d <- tryCatch(extractors[[nm]](x, runs[i]), error = function(e) {
        warning(sprintf("%s | %s: %s", runs$run_label[i], nm, conditionMessage(e)), call. = FALSE)
        NULL
      })
      if (!has_rows(d)) next
      d <- as_dt(d)
      d[, `:=`(taxon = runs$taxon[i], metric = runs$metric[i])]
      setcolorder(d, c("taxon", "metric"))
      out[[nm]] <- c(out[[nm]], list(d))
    }
    rm(x)
  }
  lapply(out, function(l) if (length(l)) rbindlist(l, fill = TRUE) else NULL)
}

cat("\nReading run outputs...\n")
S <- c(
  read_parts("_metrics.rds", list(
    point = function(x, r) as.data.table(setNames(as.list(x$.estimate), x$.metric))
  )),
  read_parts("_intervals.rds", list(
    headline = function(x, r) {
      lev <- r$main_level
      ti  <- as.data.table(x$test_intervals)[abs(nominal - lev) < 1e-6]
      ci  <- as.data.table(x$cv_intervals)
      nf  <- unique(ci[, .(fold, n)])
      cm  <- ci[abs(nominal - lev) < 1e-6, .(coverage = weighted.mean(coverage, n),
                                             width = weighted.mean(mean_width, n),
                                             interval_score = weighted.mean(interval_score, n))]
      cv_cal <- ci[, .(coverage = weighted.mean(coverage, n)), by = nominal]
      cp  <- nf[as.data.table(x$cv_pinball), on = "fold"][, .(pb = weighted.mean(pinball, n),
                                                              pr = weighted.mean(pinball_ref, n)), by = tau]
      tp  <- as.data.table(x$test_pinball)
      data.table(
        coverage_test       = ti$coverage,
        coverage_error_test = ti$coverage_error,
        width_test          = ti$mean_width,
        rel_width_test      = ti$mean_width / sd(x$test_pred$value),   # interval width / outcome SD
        interval_score_test = ti$interval_score,
        mace_test           = x$test_mace,
        crps_test           = 2 * mean(tp$pinball),
        skill_test          = 1 - mean(tp$pinball) / mean(tp$pinball_ref),
        coverage_cv         = cm$coverage,
        coverage_error_cv   = cm$coverage - lev,
        width_cv            = cm$width,
        interval_score_cv   = cm$interval_score,
        mace_cv             = mean(abs(cv_cal$coverage - cv_cal$nominal)),
        crps_cv             = 2 * mean(cp$pb),
        skill_cv            = 1 - mean(cp$pb) / mean(cp$pr),
        skill_cv_fold_sd    = sd(x$cv_fold_table$skill),
        n_folds             = nrow(nf))
    },
    # Calibration curves; grouped CV pooled over folds (each row counts once)
    calibration = function(x, r) {
      te <- as.data.table(x$test_intervals)[, .(nominal, coverage, mean_width, interval_score, n,
                                                coverage_fold_min = NA_real_, coverage_fold_max = NA_real_)]
      cv <- as.data.table(x$cv_intervals)[, .(coverage = weighted.mean(coverage, n),
                                              mean_width = weighted.mean(mean_width, n),
                                              interval_score = weighted.mean(interval_score, n),
                                              n = sum(n), coverage_fold_min = min(coverage),
                                              coverage_fold_max = max(coverage)), by = nominal]
      rbind(te[, evaluation := "test set"], cv[, evaluation := "grouped CV"], use.names = TRUE
      )[, coverage_error := coverage - nominal][]
    },
    pinball = function(x, r) {
      nf <- unique(as.data.table(x$cv_intervals)[, .(fold, n)])
      cv <- nf[as.data.table(x$cv_pinball), on = "fold"][, .(pinball = weighted.mean(pinball, n),
                                                             pinball_ref = weighted.mean(pinball_ref, n)),
                                                         by = tau][, skill := 1 - pinball / pinball_ref]
      te <- as.data.table(x$test_pinball)[, .(tau, pinball, pinball_ref, skill)]
      rbind(te[, evaluation := "test set"], cv[, evaluation := "grouped CV"])
    },
    pit_hist    = function(x, r) bin_pit(as.data.table(x$pit)[, .(evaluation = recode_eval(set), pit)], "evaluation"),
    pit_stats   = function(x, r) pit_moments(as.data.table(x$pit)[, .(evaluation = recode_eval(set), pit)], "evaluation"),
    # Fallback for the per-data-set table if _datasets.rds is missing
    by_dataset  = function(x, r) rbind(as_dt(x$test_by_dataset)[, source := "test set"],
                                       as_dt(x$cv_by_dataset)[, fold := NULL][, source := "grouped CV"],
                                       fill = TRUE)
  )),
  read_parts("_tuning.rds", list(
    tuning = function(x, r) {
      s <- as_dt(x$summary)
      b <- x$best
      s[, bagging := as.character(bagging)]
      s[, selected := mtry == b$mtry & min_n == b$min_n & bagging == as.character(b$bagging)]
      s[, `:=`(crps_rel = crps / min(crps) - 1, se_rel = crps_se / min(crps))]
      s[]
    }
  )),
  read_parts("_split.rds", list(
    split_coverage = function(x, r) setnames(as_dt(x$coverage), "feature", "variable"),
    split_groups   = function(x, r) rbind(data.table(data.set = x$test_groups, role = "test"),
                                          data.table(data.set = x$protected_groups,
                                                     role = "protected (training)"))
  )),
  read_parts("_quantile_importance.rds", list(
    importance = function(x, r) rbind(as_dt(x$results)[, evaluation := "grouped CV"],
                                      as_dt(x$results_test)[, evaluation := "test set"], fill = TRUE)
  )),
  read_parts("_ale.rds", list(
    ale_summary    = function(x, r) x$summary,
    ale_curves     = function(x, r) x$curves,
    ale_additivity = function(x, r) x$additivity
  )),
  read_parts("_datasets.rds", list(
    contribution   = function(x, r) x$contribution,
    concentration  = function(x, r) x$concentration,
    programmes     = function(x, r) x$programmes,
    ds_performance = function(x, r) x$performance_by_dataset,
    weighting      = function(x, r) x$weighting,
    influence      = function(x, r) x$influence,
    noise_floor    = function(x, r) x$noise_floor,
    equal_weight   = function(x, r) x$equal_weight,
    ale_stability  = function(x, r) x$lodo_ale_shifts,
    size_rho       = function(x, r) if (length(x$size_influence_rho)) as.data.table(as.list(x$size_influence_rho)),
    tracked        = function(x, r) data.table(variable = x$tracked_predictors,
                                               rank = seq_along(x$tracked_predictors))
  )),
  read_parts("_lofo.rds", list(
    lofo_summary  = function(x, r) x$summary,
    # Pooled contrast of qrf_fit.R plus the range of the CRPS change over CV folds
    lofo_contrast = function(x, r) {
      fc <- x$rows[model %in% c("LOFO", "control"), .(crps = mean(crps)),
                   by = .(evaluation, fold, family, model)]
      fc <- dcast(fc, evaluation + fold + family ~ model, value.var = "crps")[, change := LOFO / control - 1]
      fr <- fc[, .(n_folds = .N, crps_change_fold_min = if (.N > 1) min(change) else NA_real_,
                   crps_change_fold_max = if (.N > 1) max(change) else NA_real_), by = .(evaluation, family)]
      fr[as_dt(x$contrast), on = .(evaluation, family)]
    },
    lofo_calibration = function(x, r) x$calibration,
    lofo_pit_hist    = function(x, r) bin_pit(x$rows, c("evaluation", "family", "model")),
    lofo_pit_stats   = function(x, r) pit_moments(x$rows, c("evaluation", "family", "model")),
    lofo_range       = function(x, r) x$outside_range[, .(share = weighted.mean(share, n_eval)),
                                                      by = .(evaluation, family, model, variable)],
    lofo_by_dataset  = function(x, r) x$by_dataset,
    lofo_families    = function(x, r) x$families
  ))
)

## --- Factor levels ----
family_levels <- family_order %||% sort(unique(unlist(lapply(S, function(d)
  if (has_rows(d) && "family" %in% names(d)) unique(as.character(d$family))))))
class_levels <- unique(c(unname(predictor_classes), "Other", unname(extra_classes)))

set_levels <- function(dt) {
  if (!has_rows(dt)) return(dt)
  dt[, taxon  := factor(taxon,  levels = taxon_levels)]
  dt[, metric := factor(metric, levels = metric_levels)]
  if ("evaluation" %in% names(dt) && all(dt$evaluation %in% eval_levels)) {
    dt[, evaluation := factor(evaluation, levels = eval_levels)]
  }
  if ("family" %in% names(dt)) dt[, family := factor(family, levels = family_levels)]
  if ("model" %in% names(dt) && all(dt$model %in% lofo_models)) {
    dt[, model := factor(model, levels = lofo_models)]
  }
  if ("variable" %in% names(dt)) {
    cls <- c(predictor_classes, extra_classes)
    v   <- as.character(dt$variable)
    dt[, pclass := factor(fifelse(v %in% names(cls), cls[v], "Other"), levels = class_levels)]
  }
  dt[]
}
S <- lapply(S, set_levels)
runs_f <- set_levels(copy(runs))

has_part <- function(nm) has_rows(S[[nm]])
n_tax <- length(taxon_levels)

# Vertical dodge: taxa on the y axis, one offset row per level of `by`
dodge_y <- function(dt, by = "metric", width = 0.6) {
  f   <- dt[[by]]
  if (!is.factor(f)) f <- factor(f)
  nl  <- nlevels(f)
  off <- if (nl > 1) (as.integer(f) - (nl + 1) / 2) * width / (nl - 1) else 0
  dt[, y_pos := (n_tax + 1 - as.integer(taxon)) - off]
  dt[]
}
scale_y_taxa <- function() {
  scale_y_continuous(breaks = rev(seq_along(taxon_levels)), labels = taxon_levels,
                     expand = expansion(add = 0.5))
}
theme_summary <- function(...) theme_bw() + theme(legend.position = "bottom", ...)


#===============================================================================*
# 4. Predictive performance ----
#===============================================================================*
cat("\nPredictive performance...\n")

## --- One row per model ----
perf <- runs_f[, .(taxon, metric, n_rows, n_datasets_train, n_datasets_test,
                   bagging, mtry, min_n, main_level)]
for (nm in c("point", "headline")) {
  if (has_part(nm)) perf <- merge(perf, S[[nm]], by = c("taxon", "metric"), all.x = TRUE)
}
if (has_part("pit_stats")) {
  ps <- dcast(S$pit_stats, taxon + metric ~ evaluation, value.var = c("pit_mean", "pit_dispersion"))
  setnames(ps, gsub(" ", "_", names(ps)))
  setnames(ps, sub("_grouped_CV$", "_cv", sub("_test_set$", "_test", names(ps))))
  perf <- merge(perf, ps, by = c("taxon", "metric"), all.x = TRUE)
}
if (all(c("skill_test", "skill_cv") %in% names(perf))) perf[, skill_gap := skill_test - skill_cv]

# Out-of-sample performance of every data set: test sets from the main model,
# training sets from their LODO refit (or grouped CV if LODO was skipped)
ds_perf <- if (has_part("ds_performance")) S$ds_performance else S$by_dataset
if (has_rows(ds_perf)) {
  ds_perf <- runs_f[, .(taxon, metric, main_level)][ds_perf, on = .(taxon, metric)]
  ds_sum  <- ds_perf[, .(skill_ds_median   = median(skill),
                         share_ds_neg_skill = mean(skill < 0),
                         coverage_ds_median = median(coverage),
                         coverage_ds_min    = min(coverage),
                         share_ds_undercov  = mean(coverage < main_level - 0.1)),
                     by = .(taxon, metric)]
  perf <- merge(perf, ds_sum, by = c("taxon", "metric"), all.x = TRUE)
}
if (has_part("weighting")) {
  pw <- S$weighting[evaluation == "all data sets, out-of-sample" & weighting == "per data set",
                    .(taxon, metric, skill_per_ds = skill, coverage_per_ds = coverage)]
  perf <- merge(perf, pw, by = c("taxon", "metric"), all.x = TRUE)
}
setcolorder(perf, intersect(c("taxon", "metric", "n_rows", "n_datasets_train", "n_datasets_test",
                              "rmse", "rsq", "mae", "crps_test", "skill_test", "crps_cv", "skill_cv",
                              "skill_cv_fold_sd", "skill_gap", "coverage_test", "coverage_cv",
                              "coverage_error_test", "coverage_error_cv", "mace_test", "mace_cv",
                              "width_test", "rel_width_test", "interval_score_test",
                              "pit_mean_test", "pit_dispersion_test", "pit_mean_cv",
                              "pit_dispersion_cv", "skill_per_ds", "coverage_per_ds",
                              "skill_ds_median", "share_ds_neg_skill", "coverage_ds_median",
                              "coverage_ds_min", "share_ds_undercov"), names(perf)))
save_table(perf, "performance.csv")
save_table(S$calibration, "calibration.csv")
save_table(S$pinball, "pinball_by_tau.csv")
save_table(ds_perf, "performance_by_dataset.csv")
save_table(S$weighting, "weighting_pooled_vs_per_dataset.csv")
print(perf[, intersect(c("taxon", "metric", "rsq", "skill_test", "skill_cv", "coverage_test",
                         "coverage_cv", "mace_test", "pit_dispersion_test", "skill_per_ds"),
                       names(perf)), with = FALSE], digits = 3)

## --- Figure settings ----
metric_levels_fig <- c("PERMANOVA R2", "Classification Strength", "AucZeta",
                       "ANOSIM R", "Fuzzy Mantel", "PERMANOVA Fuzzy R2")
metric_colours <- c(`PERMANOVA R2`            = "#1b9e77",
                    `Classification Strength` = "#d95f02",
                    AucZeta                   = "#7570b3",
                    `ANOSIM R`                = "#e7298a",
                    `Fuzzy Mantel`            = "#66a61e",
                    `PERMANOVA Fuzzy R2`      = "#e6ab02")

## --- Overview figure: R2 (test set) and CRPS skill (test set vs. grouped CV) ----
pm <- c(rsq = "R2", skill = "CRPS skill")
pl <- rbind(
  perf[, .(taxon, metric, evaluation = "test set",   skill = skill_test, rsq)],
  perf[, .(taxon, metric, evaluation = "Grouped-CV", skill = skill_cv)],
  fill = TRUE
)
pl <- melt(pl, id.vars = c("taxon", "metric", "evaluation"), variable.name = "measure", na.rm = TRUE)
pl[, `:=`(measure = factor(measure, levels = names(pm), labels = pm),
          metric  = factor(metric,  levels = metric_levels_fig))]

ref  <- data.table(measure = factor(pm, levels = pm), x = 0)
tags <- data.frame(measure = pm, label = c("b", "a"))

p <- ggplot(dodge_y(pl), aes(x = value, y = y_pos, colour = metric, shape = evaluation)) +
  geom_vline(data = ref, aes(xintercept = x), linetype = 2, colour = "grey50") +
  geom_point(size = 4, stroke = 0.8) +
  scale_color_manual(values = metric_colours) +
  facet_wrap(~ measure, scales = "free_x", nrow = 1) +
  scale_y_taxa() +
  geom_text(data = tags, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.5, vjust = 1.5, fontface = "bold", size = 5, inherit.aes = FALSE) +
  labs(x = NULL, y = NULL, colour = "Metric") +
  theme_summary()
save_fig(p, "main_text/Results/res_1.png", width = 8, height = 4)

## --- Pinball skill by quantile level (grouped CV) ----
pb <- S$pinball[evaluation == "grouped CV"][, metric := factor(metric, levels = metric_levels_fig)][]

p <- ggplot(pb, aes(tau, skill, colour = metric, group = interaction(taxon, metric, evaluation))) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_line(linewidth = 0.5, alpha = 0.85) +
  scale_colour_manual(values = metric_colours, name = NULL) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), expand = expansion(add = 0.01)) +
  labs(x = expression(tau), y = "Pinball skill") +
  facet_wrap(~ taxon) +
  theme_bw() +
  theme(legend.position = "bottom", panel.spacing = unit(1, "lines"))
save_fig(p, "main_text/Results/res_2.png", width = 5, height = 4)

## --- Calibration: coverage error by nominal level ----
cal <- copy(S$calibration)
cal[, metric := factor(metric, levels = metric_levels_fig)]   # same order and colours as the pinball figure

p <- ggplot(cal, aes(nominal, coverage - nominal, colour = metric,
                     group = interaction(metric, evaluation))) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_line(aes(linetype = evaluation), linewidth = 0.5) +
  geom_point(size = 0.9) +
  scale_colour_manual(values = metric_colours, name = NULL) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  facet_wrap(~ taxon) +
  labs(x = "Nominal coverage of the central interval",
       y = "Coverage error (empirical − nominal)") +
  guides(linetype = guide_legend(ncol = 1, title = NULL)) +
  theme_bw() +
  theme(legend.position = "bottom", panel.spacing = unit(1, "lines"))
save_fig(p, "main_text/Results/res_3.png", width = 7, height = 7)


#===============================================================================*
# 5. Tuning and test-set representativeness ----
#===============================================================================*
cat("\nTuning and split...\n")

tu <- copy(S$tuning)[, bagging := factor(bagging, levels = intersect(c("row", "group"),
                                                                     unique(bagging)))]
best_bag <- tu[, .SD[which.min(crps)], by = .(taxon, metric, bagging)]
bagging_cmp <- dcast(best_bag, taxon + metric ~ bagging,
                     value.var = c("crps", "crps_se", "mtry", "min_n", "coverage"))
if (all(c("crps_group", "crps_row") %in% names(bagging_cmp))) {
  bagging_cmp[, `:=`(crps_group_vs_row = crps_group / crps_row - 1,
                     diff_in_se = (crps_group - crps_row) / sqrt(crps_se_group^2 + crps_se_row^2))]
}
bagging_cmp <- runs_f[, .(taxon, metric, selected = bagging)][bagging_cmp, on = .(taxon, metric)]

save_table(S$tuning, "tuning_candidates.csv")
save_table(bagging_cmp, "tuning_bagging_comparison.csv")
save_table(S$split_coverage, "split_coverage.csv")
save_table(S$split_groups, "split_groups.csv")


#===============================================================================*
# 6. Predictor importance ----
#     Primary: grouped CV (held-out training data sets); confirmatory: test set.
#     `relevant` is the relevance rule of qrf_fit.R.
#===============================================================================*
cat("\nImportance...\n")
imp <- S$importance
share_of <- function(d, by) d[, share := 100 * pmax(rel, 0) / sum(pmax(rel, 0)), by = by][!is.finite(share), share := 0][]
imp_crps      <- share_of(imp[evaluation == "grouped CV" & measure == "crps"], c("taxon", "metric"))
imp_crps_test <- share_of(imp[evaluation == "test set"   & measure == "crps"], c("taxon", "metric"))

imp_pb   <- share_of(imp[evaluation == "grouped CV" & measure == "pinball"], c("taxon", "metric", "tau"))
imp_tail <- imp_pb[, .(tail   = mean(share[tau <= tail_tau + 1e-9 | tau >= 1 - tail_tau - 1e-9]),
                       centre = mean(share[abs(tau - 0.5) <= center_tau + 1e-9])),
                   by = .(taxon, metric, variable, pclass)][, contrast := tail - centre][]

imp_agree <- merge(
  imp_crps[, .(taxon, metric, variable, pclass,
               rel_cv = rel, share_cv = share, relevant_cv = relevant)],
  imp_crps_test[, .(taxon, metric, variable,
                    rel_test = rel, share_test = share, relevant_test = relevant)],
  by = c("taxon", "metric", "variable")
)
imp_rho <- imp_agree[, .(
  rho = suppressWarnings(cor(rel_cv, rel_test, method = "spearman")),
  n_relevant_cv = sum(relevant_cv), n_relevant_test = sum(relevant_test),
  n_relevant_both = sum(relevant_cv & relevant_test)
), by = .(taxon, metric)]

save_table(imp_crps, "importance_crps_cv.csv")
save_table(imp_crps_test, "importance_crps_test.csv")
save_table(imp_tail, "importance_tail_vs_centre.csv")
save_table(imp_rho, "importance_cv_vs_test.csv")


#===============================================================================*
# 7. Effects (ALE) ----
#===============================================================================*
cat("\nEffects...\n")
save_table(S$ale_summary, "ale_summary.csv")
save_table(S$ale_additivity, "ale_additivity.csv")

# ALE output used for additivity in the overview: the quantile output closest
# to the median (q0.5, q50, q_0.5, ...); NULL if none can be parsed
q_outs  <- grep("^q", unique(as.character(S$ale_summary$output)), value = TRUE)
q_tau   <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", q_outs)))
q_tau   <- ifelse(q_tau > 1, q_tau / 100, q_tau)
ale_loc <- if (any(is.finite(q_tau))) q_outs[which.min(abs(q_tau - 0.5))] else NULL


#===============================================================================*
# 8. Data-set contribution and influence ----
#===============================================================================*
cat("\nData sets...\n")
save_table(S$contribution, "dataset_contribution.csv")
save_table(S$concentration, "dataset_concentration.csv")
save_table(S$programmes, "dataset_programmes.csv")

## --- LODO influence relative to the replicate-fit noise floor ----
nf <- S$noise_floor[, .(taxon, metric,
                        floor_shift = pred_shift,
                        floor_crps  = test_crps_change,
                        floor_ale   = ale_max_shift_rel)]

# Equal data-set weights
eq <- nf[S$equal_weight[model == "equal data-set weights"], on = .(taxon, metric)]
eq[, shift_ratio := pred_shift / pmax(floor_shift, .Machine$double.eps)]

# Leave one data set out
infl <- nf[S$influence, on = .(taxon, metric)]
infl[, shift_ratio := pred_shift / pmax(floor_shift, .Machine$double.eps)]
infl[, influential := shift_ratio > influence_factor]
infl_sum <- infl[, .(n_train_ds          = .N,
                     n_influential       = sum(influential),
                     max_shift_ratio     = max(shift_ratio),
                     most_influential    = data.set[which.max(shift_ratio)],
                     share_most_infl     = share_train[which.max(shift_ratio)],
                     max_crps_change     = max(test_crps_change),
                     n_ale_class_changes = sum(ale_class_changes)),
                 by = .(taxon, metric)]
infl_sum <- merge(infl_sum,
                  S$size_rho[, .(taxon, metric, rho_size_shift = pred_shift, rho_size_crps = test_crps)],
                  by = c("taxon", "metric"), all.x = TRUE)
infl_sum <- merge(infl_sum,
                  eq[, .(taxon, metric,
                         eq_shift_ratio       = shift_ratio,
                         eq_crps_change       = test_crps_change,
                         eq_crps_change_dw    = test_crps_change_dw,
                         eq_ale_class_changes = ale_class_changes)],
                  by = c("taxon", "metric"), all = TRUE)

save_table(infl, "dataset_influence.csv")
save_table(infl_sum, "dataset_influence_summary.csv")

## --- Stability of the tracked ALE classes under LODO ----
ale_stab <- S$ale_stability[, .(n_refits = .N, share_class_changed = mean(changed),
                                max_shift_rel = max(shift_rel), median_shift_rel = median(shift_rel),
                                class_full = class_full[1]),
                            by = .(taxon, metric, variable, pclass)]
save_table(ale_stab, "dataset_ale_stability.csv")


#===============================================================================*
# 9. Transfer between typology families: leave one family out ----
#     Extrapolation cost of family f = CRPS of the model trained without f
#     (LOFO) relative to a size-matched control that saw all families, on
#     identical held-out rows. Primary: grouped CV; confirmatory: test set.
#===============================================================================*
lc <- S$lofo_contrast
save_table(lc, "lofo_contrast.csv")
save_table(S$lofo_summary, "lofo_summary.csv")
save_table(S$lofo_calibration, "lofo_calibration.csv")
save_table(S$lofo_by_dataset, "lofo_by_dataset.csv")
save_table(S$lofo_pit_stats, "lofo_pit.csv")
save_table(S$lofo_families, "lofo_families.csv")

## --- Which predictors leave the training range without the family? ----
lr <- S$lofo_range[model %in% c("LOFO", "control")]
lofo_range_diff <- dcast(lr, taxon + metric + evaluation + family + variable + pclass ~ model,
                         value.var = "share")[, diff := LOFO - control][]
save_table(lofo_range_diff, "lofo_outside_range.csv")


#===============================================================================*
# 10. Overview: one row per taxon x metric ----
#===============================================================================*
cat("\nOverview...\n")
ov <- perf[, intersect(c("taxon", "metric", "n_rows", "n_datasets_train", "n_datasets_test",
                         "bagging", "mtry", "min_n", "rsq", "skill_test", "skill_cv", "skill_gap",
                         "coverage_test", "coverage_cv", "mace_test", "mace_cv", "rel_width_test",
                         "pit_dispersion_test", "skill_per_ds", "share_ds_neg_skill"), names(perf)),
           with = FALSE]
add_ov <- function(ov, d) if (has_rows(d)) merge(ov, d, by = c("taxon", "metric"), all.x = TRUE) else ov

if (has_part("concentration")) {
  ov <- add_ov(ov, S$concentration[data_sets == "all", .(taxon, metric, n_eff_datasets = n_eff_rows,
                                                         largest_ds_share = largest_share)])
}
if (has_rows(bagging_cmp) && "crps_group_vs_row" %in% names(bagging_cmp)) {
  ov <- add_ov(ov, bagging_cmp[, .(taxon, metric, crps_group_vs_row)])
}
if (has_rows(imp_crps)) {
  ov <- add_ov(ov, imp_crps[order(-share), .(n_relevant = sum(relevant), top_predictor = as.character(variable[1]),
                                             top_share = share[1], top3_share = sum(head(share, 3))),
                            by = .(taxon, metric)])
  ov <- add_ov(ov, imp_rho[, .(taxon, metric, importance_rho_cv_test = rho)])
}
if (has_part("ale_additivity") && !is.null(ale_loc)) {
  ov <- add_ov(ov, S$ale_additivity[output == ale_loc, .(taxon, metric, additivity_median = r2_main_effects)])
}
if (has_rows(infl_sum)) {
  ov <- add_ov(ov, infl_sum[, intersect(c("taxon", "metric", "n_influential", "max_shift_ratio",
                                          "most_influential", "eq_shift_ratio", "eq_crps_change"),
                                        names(infl_sum)), with = FALSE])
}
if (has_rows(lc)) {
  for (e in eval_levels) {
    w <- lc[evaluation == e]
    if (!nrow(w)) next
    w <- dcast(w, taxon + metric ~ family, value.var = c("crps_change", "coverage_change"))
    tag <- if (e == "grouped CV") "cv" else "test"
    setnames(w, setdiff(names(w), c("taxon", "metric")),
             paste0("lofo_", tag, "_", setdiff(names(w), c("taxon", "metric"))))
    ov <- add_ov(ov, w)
  }
}
setorder(ov, metric, taxon)
save_table(ov, "overview.csv")
print(ov[, lapply(.SD, function(z) if (is.double(z)) signif(z, 3) else z)], nrows = 200)

saveRDS(list(runs = runs, inventory = inventory, overview = ov, performance = perf,
             performance_by_dataset = ds_perf, bagging = bagging_cmp,
             importance_crps = imp_crps, importance_tail = imp_tail,
             importance_agreement = imp_agree, importance_rho = imp_rho,
             influence = infl, influence_summary = infl_sum, ale_stability = ale_stab,
             lofo_outside_range = lofo_range_diff, parts = S,
             settings = list(influence_factor = influence_factor, pit_bins = pit_bins,
                             tail_tau = tail_tau, center_tau = center_tau)),
        file.path(output_dir, "summary_tables.rds"))

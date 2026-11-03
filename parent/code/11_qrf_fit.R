### - Quantile Random Forest: fitting and evaluation - ###
#
# Usage:
#   Rscript qrf_fit.R --metric "<metric name>" --taxon "<taxon group>" --output <dir>
#                     [--scheme-cap K] [--max-rows N] [--no-lodo] [--no-lofo] [--allow-dropped-subsets]
#
#   --taxon           taxonomic group to model (filters on column `taxon_group`);
#                     if omitted, all taxa are pooled and the run is labelled "all"
#   --scheme-cap K    at most K subsets (schemes) per data set, evenly spaced along
#                     subset size (`samples`); default 20, 0 = no cap
#   --max-rows N      random subsample of whole schemes up to N rows (testing only)
#   --no-lodo         skip the leave-one-data-set-out refits in section 13
#   --no-lofo         skip the leave-one-family-out refits in section 14
#   --allow-dropped-subsets  continue if whole subsets have missing predictors
#                     (default: stop; dropped rows are always written to <run>_dropped.rds)
#   Requires R/prepare_qrf_inputs.R 
#
# Pipeline: load + merge -> designed grouped split -> grouped-CV tuning on CRPS
# (mtry, min_n, bagging unit) -> ONE final fit (saved immediately) -> evaluation:
#    grouped CV (out-of-fold, from the tuning fits of the selected candidate) and
#    test set: coverage, width, interval score, PIT, pinball/CRPS;
#    quantile-specific permutation importance (test set), accumulated local effects;
#    data-set influence (LODO); transfer between typology families (LOFO).
# All fits use ranger directly (quantile regression forest) with all threads.

#===============================================================================*
# 1. Setup ----
#===============================================================================*
suppressPackageStartupMessages({
  library(data.table)
  library(tidymodels)   # recipes, rsample, dials, dplyr, tibble, purrr, ggplot2
  library(ranger)
})

## --- Command-line arguments ----
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value for ", flag)
  args[i + 1]
}
metric_name     <- get_arg("--metric")
taxon_name      <- get_arg("--taxon")
output_dir      <- get_arg("--output")
max_rows        <- as.numeric(get_arg("--max-rows", Inf))
scheme_cap      <- as.numeric(get_arg("--scheme-cap", 20))
if (is.na(scheme_cap) || scheme_cap <= 0) scheme_cap <- Inf
do_lodo         <- !("--no-lodo" %in% args)
do_lofo         <- !("--no-lofo" %in% args)
allow_dropped   <- "--allow-dropped-subsets" %in% args

if (is.null(metric_name) || is.null(output_dir)) {
  stop("Missing --metric or --output argument.")
}

## --- Output paths ----
# Files are named <taxon>__<metric><suffix>; run_id is set in section 4 once the
# taxon label is known. The summary script finds runs via their _run_info.rds.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
clean_name <- function(x) gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", x))

out_path <- function(suffix) file.path(output_dir, paste0(run_id, suffix))
save_out <- function(obj, suffix) {
  saveRDS(obj, out_path(suffix))
  cat("  saved:", out_path(suffix), "\n")
}
save_fig <- function(plot, suffix, width, height) {
  ggsave(out_path(suffix), plot, width = width, height = height)
  cat("  saved:", out_path(suffix), "\n")
}
log_step <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

## --- Cores: respect the SLURM allocation, not the whole node ----
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores)) n_cores <- max(parallel::detectCores(logical = FALSE), 2L, na.rm = TRUE)

## --- Analysis settings ----
n_trees        <- 3000
n_trees_tune   <- 1500
tune_grid_size <- 10                    # (mtry, min_n) combinations
tune_bagging   <- c("row", "group")    # bootstrap unit: observations or whole data sets
cv_folds       <- 5

# Designed train/test split (fixed rule, identical for every model)
split_test_prop <- 0.25                  # share of DATA SETS held out for testing
split_n_pc      <- 3                     # leading PCs of data-set predictor means
split_protect   <- c(paste0("PC", seq_len(split_n_pc)), "outcome")  # extremes stay in training

# Predictive evaluation
nominal_levels <- c(seq(0.1, 0.9, by = 0.1), 0.95)   # central interval levels
main_level     <- 0.9                                # headline interval / Winkler level
pit_step       <- 0.005
eval_probs     <- round(seq(pit_step, 1 - pit_step, by = pit_step), 6)

# Permutation importance (test set)
imp_probs    <- round(seq(0.01, 0.99, by = 0.01), 6) # uniform grid -> CRPS
imp_taus     <- round(seq(0.05, 0.95, by = 0.05), 6) # heatmap levels
imp_n_perm   <- 3
imp_n_boot   <- 2000
imp_max_rows_per_ds <- 3000   # rows per data set used for permutation importance

# Relevance rule for importance (used for heatmap marks and ALE "negligible" class):
# loss increase >= imp_min_rel AND permuting hurt in >= imp_min_consistency of
# held-out data sets. Bootstrap p-values are only reported with enough groups.
imp_min_rel         <- 0.01
imp_min_consistency <- 0.75
imp_min_groups_p    <- 8

# Data-set contribution and influence (section 13)
lodo_top_k       <- 3      # most important predictors whose ALE stability is tracked
ds_upstream_file <- NULL   # optional CSV: data.set + upstream counts (e.g. n_subsets, n_hmsc, n_communities)
programme_of     <- function(ds) vapply(strsplit(ds, "_"), function(z) z[min(2, length(z))], "")

# ALE (training data)
ale_taus           <- c(0.1, 0.5, 0.9)               # width = last - first
ale_bins           <- 20
ale_max_rows       <- 10000
ale_n_boot         <- 1000
ale_mono_threshold <- 0.8

a_main <- 1 - main_level
stopifnot(
  all(round(c((1 - nominal_levels) / 2, 1 - (1 - nominal_levels) / 2), 6) %in% eval_probs),
  all(imp_taus %in% imp_probs),
  all(round(c(a_main / 2, 1 - a_main / 2), 6) %in% imp_probs),
  length(ale_taus) >= 2, !is.unsorted(ale_taus)
)

cat("\n========================================\n")
cat("Metric:          ", metric_name, "\n")
cat("Taxon:           ", if (is.null(taxon_name)) "(all pooled)" else taxon_name, "\n")
cat("Output directory:", output_dir, "\n")
cat("Cores:           ", n_cores, "\n")
cat("Scheme cap:      ", scheme_cap, "subsets per data set\n")
cat("Max rows:        ", max_rows, "\n")
cat("Bagging units:   ", paste(tune_bagging, collapse = ", "), "\n")
cat("LODO refits:     ", do_lodo, "\n")
cat("LOFO refits:     ", do_lofo, "\n")
cat("========================================\n\n")


#===============================================================================*
# 2. Helper functions ----
#===============================================================================*

## --- Cap the number of subsets (schemes) per data set ----
# A subset S is identified by data.set + scheme_id; its size is `samples`.
# Data sets with more than k subsets keep k of them at evenly spaced positions
# along their subsets sorted by size (smallest and largest always kept); ties in
# size are broken at random (fixed seed). All rows of a kept subset are kept.
cap_schemes <- function(d, k, seed = 1) {
  s <- unique(d[, .(data.set, scheme_id, samples)])
  if (anyDuplicated(s[, .(data.set, scheme_id)])) {
    stop("`samples` is not constant within data.set x scheme_id.")
  }
  if (!is.finite(k)) return(d)
  set.seed(seed)
  s <- s[sample(.N)][order(data.set, samples)]
  keep <- s[, .SD[if (.N <= k) seq_len(.N) else round(seq(1, .N, length.out = k))],
            by = data.set]
  d[keep[, .(data.set, scheme_id)], on = .(data.set, scheme_id), nomatch = NULL]
}
scheme_summary <- function(d) {
  d[, .(n_rows = .N, n_schemes = uniqueN(scheme_id),
        samples_min = min(samples), samples_max = max(samples)), by = data.set]
}

## --- Designed grouped split ----
# Each data set is profiled by its mean predictors (reduced to the leading PCs)
# and its mean outcome. Data sets at the extremes of the `protect` axes stay in
# training (QRF cannot extrapolate). The remaining data sets are clustered into
# n_test clusters; the most central data set of each cluster goes to the test set.
dataset_profile <- function(data, group, features) {
  prof <- as.data.table(data)[, c(lapply(.SD, mean, na.rm = TRUE), .(n_rows = .N)),
                              by = group, .SDcols = features]
  setnames(prof, group, "group")
  setorder(prof, group)          # deterministic order, independent of row order
  prof[]
}

designed_group_split <- function(data, group, predictors, outcome, test_prop = 0.25,
                                 n_pc = 3, protect = c("PC1", "outcome"),
                                 nstart = 100, seed = 1) {
  prof     <- dataset_profile(data, group, c(predictors, outcome))
  n_groups <- nrow(prof)
  n_test   <- max(2L, as.integer(round(test_prop * n_groups)))
  if (n_groups - n_test < 3) stop("Too few data sets (", n_groups, ") for a train/test split.")
  
  P <- as.matrix(prof[, ..predictors])
  P <- P[, apply(P, 2, function(z) sd(z) > 0), drop = FALSE]
  k <- min(n_pc, ncol(P), n_groups - 1)
  pcs <- prcomp(P, scale. = TRUE)$x[, seq_len(k), drop = FALSE]
  Z <- scale(cbind(pcs, outcome = prof[[outcome]]))
  Z[!is.finite(Z)] <- 0
  
  extremes <- function(axes) unique(unlist(lapply(intersect(axes, colnames(Z)), function(a)
    c(which.min(Z[, a]), which.max(Z[, a])))))
  # Protect extremes of all axes; drop trailing axes only if too few candidates remain
  protected <- integer(0)
  for (m in rev(seq_along(protect))) {
    pr <- extremes(protect[seq_len(m)])
    if (n_groups - length(pr) >= n_test + 1) { protected <- pr; break }
  }
  if (length(protected) < length(extremes(protect))) {
    warning("Too few candidate data sets to protect all axes; protection was reduced.")
  }
  candidates <- setdiff(seq_len(n_groups), protected)
  
  set.seed(seed)
  km <- kmeans(Z[candidates, , drop = FALSE], centers = n_test, nstart = nstart)
  test_rows <- vapply(seq_len(n_test), function(j) {
    members <- candidates[km$cluster == j]
    members[which.min(rowSums(sweep(Z[members, , drop = FALSE], 2, km$centers[j, ])^2))]
  }, integer(1))
  
  in_test <- data[[group]] %in% prof$group[test_rows]
  split <- make_splits(list(analysis = which(!in_test), assessment = which(in_test)),
                       data = data)
  attr(split, "test_groups")      <- prof$group[test_rows]
  attr(split, "protected_groups") <- prof$group[protected]
  attr(split, "profile")          <- cbind(prof, as.data.table(Z))
  split
}

# How well do the test data sets span the training data sets (data-set means)?
split_coverage_report <- function(split, group, features) {
  tr <- dataset_profile(training(split), group, features)
  te <- dataset_profile(testing(split),  group, features)
  rbindlist(lapply(features, function(f) {
    rng <- max(tr[[f]]) - min(tr[[f]])
    data.table(feature       = f,
               train_min     = min(tr[[f]]), train_max = max(tr[[f]]),
               test_min      = min(te[[f]]), test_max  = max(te[[f]]),
               range_spanned = if (rng > 0) (max(te[[f]]) - min(te[[f]])) / rng else NA_real_,
               test_in_range = mean(te[[f]] >= min(tr[[f]]) & te[[f]] <= max(tr[[f]])))
  }))
}

## --- Quantile regression forest fit (row or data-set bagging) ----
# Group bagging: each tree is grown on a bootstrap sample of whole data sets
# (in-bag count of a row = number of times its data set was drawn).
# equalize = TRUE: every drawn data set contributes the same number of rows
# (n / number of data sets, sampled with replacement within the data set).
make_group_inbag <- function(groups, n_trees, equalize = FALSE) {
  gi <- as.integer(factor(groups)); ng <- max(gi); n <- length(gi)
  if (!equalize) {
    return(lapply(seq_len(n_trees), function(t)
      tabulate(sample.int(ng, ng, replace = TRUE), ng)[gi]))
  }
  rows_by_g <- split(seq_len(n), gi)
  m <- max(1L, round(n / ng))
  lapply(seq_len(n_trees), function(t) {
    cnt <- tabulate(sample.int(ng, ng, replace = TRUE), ng)
    idx <- unlist(lapply(which(cnt > 0), function(g) {
      r <- rows_by_g[[g]]
      r[sample.int(length(r), cnt[g] * m, replace = TRUE)]
    }))
    tabulate(idx, n)
  })
}
# equalize = TRUE gives every data set equal total influence:
#   row bagging   -> case weights 1 / n_rows(data set)
#   group bagging -> equal rows per drawn data set
fit_qrf <- function(X, y, groups, mtry, min_n, bagging = c("row", "group"), seed = 1,
                    equalize = FALSE, trees = n_trees) {
  bagging <- match.arg(bagging)
  set.seed(seed)
  inbag <- if (bagging == "group") make_group_inbag(groups, trees, equalize) else NULL
  cw    <- if (equalize && bagging == "row") as.numeric(1 / table(groups)[as.character(groups)]) else NULL
  ranger(x = X, y = y, num.trees = trees, mtry = min(mtry, ncol(X)),
         min.node.size = min_n, quantreg = TRUE, inbag = inbag, case.weights = cw,
         num.threads = n_cores, seed = seed, verbose = FALSE)
}

## --- Row- and data-set-level scores ----
# Per-row CRPS (quantile grid), reference CRPS (unconditional quantiles ref_q),
# interval hit, width and residual; column loop keeps memory low for large n.
row_scores <- function(y, q, probs, level, ref_q) {
  crps <- crps_ref <- numeric(length(y))
  for (j in seq_along(probs)) {
    u <- y - q[, j];     crps     <- crps     + u * (probs[j] - (u < 0))
    u <- y - ref_q[j];   crps_ref <- crps_ref + u * (probs[j] - (u < 0))
  }
  a  <- 1 - level
  lo <- get_q(q, probs, a / 2); hi <- get_q(q, probs, 1 - a / 2)
  data.table(crps = 2 * crps / length(probs), crps_ref = 2 * crps_ref / length(probs),
             inside = y >= lo & y <= hi, width = hi - lo,
             resid = y - get_q(q, probs, 0.5))
}
dataset_scores <- function(rs, groups) {
  rs <- copy(rs)[, data.set := groups]
  rs[, .(n = .N, crps = mean(crps), crps_ref = mean(crps_ref),
         coverage = mean(inside), width = mean(width), bias = median(resid)),
     by = data.set][, skill := 1 - crps / crps_ref][]
}
# Pooled (each row counts once) vs. per-data-set (each data set counts once)
weighting_table <- function(d, label) {
  data.table(evaluation = label,
             weighting  = c("pooled (rows)", "per data set"),
             n_datasets = nrow(d),
             coverage   = c(weighted.mean(d$coverage, d$n), mean(d$coverage)),
             crps       = c(weighted.mean(d$crps, d$n), mean(d$crps)),
             skill      = c(1 - weighted.mean(d$crps, d$n) / weighted.mean(d$crps_ref, d$n),
                            1 - mean(d$crps) / mean(d$crps_ref)))
}
n_eff    <- function(p) 1 / sum((p / sum(p))^2)          # effective number of data sets
short_ds <- function(x) sub("^[^_]+_", "", x)             # drop the taxon prefix for labels

## --- Point metrics (same format as yardstick::metrics) ----
point_metrics <- function(obs, pred) {
  tibble(.metric    = c("rmse", "rsq", "mae"),
         .estimator = "standard",
         .estimate  = c(sqrt(mean((obs - pred)^2)), cor(obs, pred)^2, mean(abs(obs - pred))))
}

## --- Evaluate one fitted model on held-out rows (tuning folds) ----
evaluate_heldout <- function(fit, X, y, groups, y_ref, fold) {
  q  <- predict_q(fit, X, eval_probs)
  pt <- predict(fit, data = X, num.threads = n_cores)$predictions
  pb <- pinball_metrics(y, q, eval_probs, y_ref)
  lo <- get_q(q, eval_probs, a_main / 2)
  hi <- get_q(q, eval_probs, 1 - a_main / 2)
  inside <- y >= lo & y <= hi
  list(
    crps          = 2 * mean(pb$pinball),        # uniform grid -> CRPS
    rmse          = sqrt(mean((y - pt)^2)),
    coverage_main = mean(inside),
    width_main    = mean(hi - lo),
    metrics       = interval_metrics(y, q, eval_probs, nominal_levels) %>% mutate(fold = fold),
    pinball       = pb %>% mutate(fold = fold),
    pit           = tibble(fold = fold, data.set = groups, pit = pit_values(y, q, pit_step)),
    by_dataset    = dataset_scores(row_scores(y, q, eval_probs, main_level,
                                              quantile(y_ref, eval_probs, names = FALSE)),
                                   groups)[, fold := fold][]
  )
}

## --- Quantile predictions ----
predict_q <- function(model, X, probs) {
  predict(model, data = X, type = "quantiles", quantiles = probs,
          num.threads = n_cores)$predictions
}
get_q <- function(q, probs, p) q[, which.min(abs(probs - p))]

## --- Interval metrics: coverage, width, interval (Winkler) score ----
interval_metrics <- function(obs, q, probs, levels) {
  purrr::map_dfr(levels, function(lev) {
    a  <- 1 - lev
    lo <- get_q(q, probs, a / 2)
    hi <- get_q(q, probs, 1 - a / 2)
    inside <- obs >= lo & obs <= hi
    tibble(
      nominal        = lev,
      coverage       = mean(inside),
      coverage_error = mean(inside) - lev,          # < 0 = too narrow
      mean_width     = mean(hi - lo),
      interval_score = mean((hi - lo) +
                              (2 / a) * (lo - obs) * (obs < lo) +
                              (2 / a) * (obs - hi) * (obs > hi)),
      n              = length(obs)
    )
  })
}

## --- Randomized PIT from a quantile grid (ties handled by uniform draw) ----
pit_values <- function(obs, q, step) {
  k_lo <- rowSums(q <  obs)
  k_hi <- rowSums(q <= obs)
  runif(length(obs), min = step * k_lo, max = step * (k_hi + 1))
}

## --- Pinball loss, skill vs. unconditional quantiles, CRPS approximation ----
pinball_loss <- function(obs, q, probs) {
  tau  <- matrix(probs, nrow = length(obs), ncol = length(probs), byrow = TRUE)
  diff <- obs - q
  colMeans(diff * (tau - (diff < 0)))
}
pinball_metrics <- function(obs, q, probs, ref_values) {
  q_ref <- matrix(quantile(ref_values, probs, names = FALSE),
                  nrow = length(obs), ncol = length(probs), byrow = TRUE)
  tibble(tau         = probs,
         pinball     = pinball_loss(obs, q, probs),
         pinball_ref = pinball_loss(obs, q_ref, probs)) %>%
    mutate(skill = 1 - pinball / pinball_ref)
}
summarise_pinball <- function(pb) {
  tibble(mean_pinball     = mean(pb$pinball),
         crps_approx      = 2 * mean(pb$pinball),
         mean_pinball_ref = mean(pb$pinball_ref),
         skill            = 1 - mean(pb$pinball) / mean(pb$pinball_ref))
}

## --- Row losses summed per group (for permutation importance) ----
group_losses <- function(y, q, groups, probs, taus, level) {
  a   <- 1 - level
  tau <- matrix(probs, nrow = length(y), ncol = length(probs), byrow = TRUE)
  u   <- y - q
  pl  <- u * (tau - (u < 0))
  crps <- 2 * rowMeans(pl)
  lo  <- q[, match(round(a / 2, 6), probs)]
  hi  <- q[, match(round(1 - a / 2, 6), probs)]
  wk  <- (hi - lo) + (2 / a) * (lo - y) * (y < lo) + (2 / a) * (y - hi) * (y > hi)
  M <- cbind(pl[, match(taus, probs), drop = FALSE], crps, wk)
  colnames(M) <- c(paste0("pinball_", taus), "crps", "winkler")
  rowsum(M, groups, reorder = TRUE)
}

## --- ALE from bin sums / counts; columns = bootstrap replicates (col 1 = point) ----
ale_from_bins <- function(Ssum, Ncnt, fallback) {
  K  <- nrow(Ssum)
  mu <- Ssum / Ncnt
  miss <- !is.finite(mu)                                 # bin empty in a replicate
  if (any(miss)) mu[miss] <- fallback[row(mu)[miss]]
  acc <- rbind(0, matrix(apply(mu, 2, cumsum), nrow = K))
  mid <- (acc[-1, , drop = FALSE] + acc[-(K + 1), , drop = FALSE]) / 2
  Nc  <- Ncnt; Nc[miss] <- 0
  ctr <- colSums(mid * Nc) / colSums(Nc)
  acc <- sweep(acc, 2, ctr)
  mid <- sweep(mid, 2, ctr)
  list(ale       = acc,
       amplitude = apply(acc, 2, max) - apply(acc, 2, min),
       sd_effect = sqrt(colSums(Nc * mid^2) / colSums(Nc)),
       net       = colSums(mu),
       mono      = abs(colSums(mu)) / colSums(abs(mu)))
}

## --- ALE for one variable: curves, summaries, per-row effects ----
ale_variable <- function(v, X, g, W, predict_fun, outputs, bins) {
  x <- X[[v]]
  edges <- unique(quantile(x, probs = seq(0, 1, length.out = bins + 1),
                           names = FALSE, type = 1))
  if (length(edges) < 2) return(NULL)                    # constant predictor
  K   <- length(edges) - 1
  G   <- nrow(W)
  bin <- pmax(findInterval(x, edges, left.open = TRUE, rightmost.closed = TRUE), 1L)
  
  Xlo <- X; Xlo[[v]] <- edges[bin]
  Xhi <- X; Xhi[[v]] <- edges[bin + 1]
  n <- nrow(X)
  P <- predict_fun(rbind(Xlo, Xhi))
  D <- P[(n + 1):(2 * n), , drop = FALSE] - P[1:n, , drop = FALSE]
  
  key  <- (g - 1L) * K + bin
  Ncnt <- matrix(tabulate(key, K * G), K, G) %*% W
  rs   <- rowsum(D, key)
  S    <- matrix(0, K * G, ncol(D))
  S[as.integer(rownames(rs)), ] <- rs
  
  bq <- function(z, p) unname(quantile(z[-1], p, na.rm = TRUE))
  curves <- list(); summ <- list(); centred <- list()
  for (o in seq_along(outputs)) {
    Ssum <- matrix(S[, o], K, G) %*% W
    r    <- ale_from_bins(Ssum, Ncnt, fallback = Ssum[, 1] / Ncnt[, 1])
    centred[[o]] <- r$ale[, 1]
    curves[[o]] <- data.table(
      variable = v, output = outputs[o], x = edges, ale = r$ale[, 1],
      lo = apply(r$ale[, -1, drop = FALSE], 1, quantile, 0.025),
      hi = apply(r$ale[, -1, drop = FALSE], 1, quantile, 0.975))
    summ[[o]] <- data.table(
      variable  = v, output = outputs[o],
      amplitude = r$amplitude[1],
      amp_lo    = bq(r$amplitude, 0.025), amp_hi = bq(r$amplitude, 0.975),
      sd_effect = r$sd_effect[1],
      net       = r$net[1],
      net_lo    = bq(r$net, 0.025), net_hi = bq(r$net, 0.975),
      mono      = r$mono[1],
      mono_lo   = bq(r$mono, 0.025), mono_hi = bq(r$mono, 0.975),
      p_dir     = mean(sign(r$net[-1]) == sign(r$net[1])))
  }
  row_eff <- sapply(centred, function(ce)
    approx(edges, ce, xout = x, rule = 2, ties = "ordered")$y)
  colnames(row_eff) <- outputs
  list(curves = rbindlist(curves), summary = rbindlist(summ), row_eff = row_eff)
}

## --- Leave-one-family-out helpers (section 14) ----
# Size-matched control for a LOFO training set: within each scheme, as many rows
# as the LOFO set keeps (`keep`), drawn at random from all rows of that scheme.
matched_rows <- function(idx, scheme, keep, seed) {
  set.seed(seed)
  d <- data.table(row = idx, s = scheme, k = keep)
  sort(d[, .(row = row[sample.int(.N, sum(k))]), by = s]$row)
}
# Scores of one fit on held-out rows. y_ref: outcome whose unconditional quantiles
# are the CRPS reference (skill); pass the same y_ref to all models being compared.
transfer_scores <- function(fit, X, y, g, y_ref) {
  q  <- predict_q(fit, X, eval_probs)
  rs <- row_scores(y, q, eval_probs, main_level, quantile(y_ref, eval_probs, names = FALSE))
  rs[, `:=`(data.set = g, pit = pit_values(y, q, pit_step))]
  list(rows = rs, intervals = interval_metrics(y, q, eval_probs, nominal_levels))
}
# Share of new rows with any (and with each) predictor outside the range of the fitting rows
outside_range <- function(X_fit, X_new) {
  M <- vapply(names(X_new), function(v) {
    r <- range(X_fit[[v]])
    X_new[[v]] < r[1] | X_new[[v]] > r[2]
  }, logical(nrow(X_new)))
  M <- matrix(M, nrow = nrow(X_new))
  data.table(variable = c("any", names(X_new)), share = c(mean(rowSums(M) > 0), colMeans(M)))
}
add_cols <- function(d, cols) {
  d <- as.data.table(d)
  d[, (names(cols)) := cols]
  setcolorder(d, names(cols))
  d[]
}

## --- Figures ----
plot_pit <- function(pit_df, bins = 10, title = NULL) {
  ggplot(pit_df, aes(x = pit)) +
    geom_histogram(aes(y = after_stat(density)),
                   breaks = seq(0, 1, length.out = bins + 1),
                   fill = "grey70", colour = "white") +
    geom_hline(yintercept = 1, linetype = 2) +
    facet_wrap(~ set) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "PIT", y = "Density", title = title) +
    theme_bw()
}

plot_pinball_skill <- function(test_pb, cv_pb = NULL, title = NULL) {
  p <- ggplot() + geom_hline(yintercept = 0, linetype = 2)
  if (!is.null(cv_pb)) {
    p <- p + geom_line(data = cv_pb, aes(x = tau, y = skill, group = fold),
                       colour = "grey60", linewidth = 0.4)
  }
  p + geom_line(data = test_pb, aes(x = tau, y = skill), linewidth = 0.8) +
    labs(x = expression(tau), y = "Pinball skill vs. unconditional quantiles",
         title = title,
         subtitle = if (is.null(cv_pb)) "Test set" else "grey = CV folds, black = test set") +
    theme_bw()
}

relevance_label <- function() {
  paste0("loss increase >= ", 100 * imp_min_rel, "% and worse in >= ",
         100 * imp_min_consistency, "% of held-out data sets")
}

plot_importance_heatmap <- function(res, value = c("rel", "share"), title = NULL) {
  value <- match.arg(value)
  d <- copy(res[measure == "pinball"])
  d[, variable := factor(variable, levels = res[measure == "crps"][order(rel), variable])]
  d[, sig := relevant]
  if (value == "rel") {
    d[, v := rel * 100]
    lim <- max(abs(d$v), na.rm = TRUE)
    fill_scale <- scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                                       midpoint = 0, limits = c(-lim, lim),
                                       name = "Pinball\nincrease (%)")
  } else {
    d[, v := 100 * pmax(rel, 0) / sum(pmax(rel, 0)), by = tau]
    fill_scale <- scale_fill_gradient(low = "white", high = "#b2182b", limits = c(0, NA),
                                      name = "Share of\nimportance\nat tau (%)")
  }
  ggplot(d, aes(x = tau, y = variable, fill = v)) +
    geom_tile(aes(alpha = sig), colour = "white", linewidth = 0.3) +
    geom_point(data = d[sig == TRUE], size = 0.8, colour = "black", show.legend = FALSE) +
    fill_scale +
    scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 1), guide = "none") +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2), expand = c(0, 0)) +
    labs(x = expression("Quantile level " * tau), y = NULL, title = title,
         subtitle = paste0(res$source[1], ". Dots / full colour: ", relevance_label())) +
    theme_bw() +
    theme(panel.grid = element_blank())
}

plot_importance_distribution <- function(res, level, title = NULL) {
  d <- copy(res[measure %in% c("crps", "winkler")])
  d[, measure := factor(measure, levels = c("crps", "winkler"),
                        labels = c("CRPS", paste0("Winkler (", level * 100, "%)")))]
  d[, variable := factor(variable, levels = res[measure == "crps"][order(rel), variable])]
  d[, sig := relevant]
  ggplot(d, aes(y = variable)) +
    geom_vline(xintercept = 0, colour = "grey50") +
    geom_errorbarh(aes(xmin = ci_lo * 100, xmax = ci_hi * 100),
                   height = 0, colour = "grey55", linewidth = 0.6) +
    geom_point(aes(x = rel * 100, shape = sig), size = 2) +
    geom_text(aes(x = ci_hi * 100, label = paste0(groups_pos, "/", n_groups)),
              hjust = -0.25, size = 2.6, colour = "grey30") +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "Relevant") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    facet_wrap(~ measure, scales = "free_x") +
    labs(x = "Loss increase when permuted (%)", y = NULL, title = title,
         subtitle = paste0("Bars: 95% group-bootstrap CI; labels: held-out groups where permuting hurt.\n",
                           "Relevant: ", relevance_label())) +
    theme_bw()
}

plot_ale <- function(curves, rug, order, outputs, out_levels, out_labels, title, subtitle) {
  d <- curves[output %in% outputs]
  d[, output := factor(output, levels = out_levels, labels = out_labels)]
  d[, variable := factor(variable, levels = order)]
  rug <- copy(rug)[, variable := factor(variable, levels = order)]
  ggplot(d, aes(x = x)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = output), alpha = 0.2) +
    geom_line(aes(y = ale, colour = output), linewidth = 0.6) +
    geom_rug(data = rug, aes(x = x), inherit.aes = FALSE,
             sides = "b", colour = "grey40", length = unit(0.02, "npc")) +
    facet_wrap(~ variable, scales = "free_x") +
    labs(x = NULL, y = "Accumulated local effect (centred, outcome units)",
         colour = NULL, fill = NULL, title = title, subtitle = subtitle) +
    theme_bw() +
    theme(legend.position = "bottom", strip.text = element_text(size = 8))
}

plot_ale_summary <- function(summary, order, outputs, out_levels, out_labels,
                             threshold, title) {
  d <- copy(summary[output %in% outputs])
  d[, output := factor(output, levels = out_levels, labels = out_labels)]
  d[, variable := factor(variable, levels = rev(order))]
  d[, class := factor(class, levels = c("increasing", "decreasing", "non-monotone", "negligible"))]
  ggplot(d, aes(y = variable)) +
    geom_vline(xintercept = 0, colour = "grey50") +
    geom_errorbarh(aes(xmin = amp_lo, xmax = amp_hi), height = 0, colour = "grey55") +
    geom_point(aes(x = amplitude, shape = class, colour = class), size = 2.2) +
    scale_shape_manual(values = c(increasing = 17, decreasing = 25,
                                  `non-monotone` = 4, negligible = 1), drop = FALSE) +
    scale_colour_manual(values = c(increasing = "#b2182b", decreasing = "#2166ac",
                                   `non-monotone` = "black", negligible = "grey60"),
                        drop = FALSE) +
    facet_wrap(~ output, scales = "free_x") +
    labs(x = "ALE amplitude: max - min of the effect curve (outcome units)", y = NULL,
         shape = NULL, colour = NULL, title = title,
         subtitle = paste0("Symbol: direction if monotonicity >= ", threshold,
                           ", else non-monotone.\nBars: 95% group-bootstrap CI")) +
    theme_bw()
}


#===============================================================================*
# 3. Load data ----
#===============================================================================*
# Subsets are keyed by FILE NAME (canonical "<taxon>_<0000>"), not by the scheme_id
# stored in the files, whose format deviates between sources. The previous join on
# the stored variable silently dropped whole subsets.
source("../pulse/code/functions/prepare_qrf_inputs.R")
date_vars <- QRF_DATE_VARS
inp <- load_qrf_inputs(log = log_step)

id_problems <- qrf_id_problems(inp$id_check)       # ids that cannot be resolved by name
if (nrow(id_problems)) {
  print(id_problems)
  stop(nrow(id_problems), " identifier problem(s); run 11a_check_qrf_inputs.R for details.")
}
fd <- inp$id_check[, .(n = sum(format_differs)), by = source]
cat("  stored ids differing from the file name in format only:",
    paste(fd$source, fd$n, sep = " = ", collapse = ", "), "\n")
cat("  evaluated subsets:", nrow(inp$coverage),
    "| without spatial scale:", sum(!inp$coverage$in_spatial),
    "| without taxonomic resolution:", sum(!inp$coverage$in_taxonomy),
    "| without date variables:", sum(!inp$coverage$in_dates), "\n")


log_step("Typology families")
sim_dir            <- "data/006_simulated_data/"
spatial_scale_dir  <- "data/misc/spatial_scale/"


#===============================================================================*
# 4. Prepare data ----
#===============================================================================*

# Merge done in load_qrf_inputs(): left joins on the file-name id, row count checked.
evl2 <- inp$evl2
rm(inp); invisible(gc())

predictor_vec <- qrf_predictor_vec                # defined in prepare_qrf_inputs.R

log_step("Filtering to metric:", metric_name, "| taxon:", taxon_name)
j.d <- evl2[metric == metric_name]
if (!is.null(taxon_name)) {                       # --taxon was documented but not applied
  j.d <- j.d[taxon_key == tolower(taxon_name)]
  if (!nrow(j.d)) {
    stop("No rows for taxon '", taxon_name, "'. Available: ",
         paste(sort(unique(evl2$taxon_key)), collapse = ", "))
  }
}
taxon_label <- if (!is.null(taxon_name)) {
  taxon_name
} else if ("taxon_group" %in% names(j.d) &&
           uniqueN(j.d$taxon_group) == 1) {
  unique(j.d$taxon_group)
} else {
  "all"
}
run_id    <- paste0(clean_name(taxon_label), "__", clean_name(metric_name))
run_label <- paste(taxon_label, "|", metric_name)

model_cols  <- qrf_model_cols
drop_report <- qrf_drop_report(j.d, model_cols)   # subsets with incomplete rows (prepare_qrf_inputs.R)
if (nrow(drop_report)) {
  save_out(drop_report, "_dropped.rds")
  n_whole <- drop_report[whole_subset == TRUE, .N]
  cat("  incomplete rows dropped:", sum(drop_report$n_dropped), "in", nrow(drop_report),
      "subsets, of which", n_whole, "entirely\n")
  if (n_whole > 0) {
    print(drop_report[whole_subset == TRUE, .(subsets = .N), by = .(data.set, na_columns)])
    if (!allow_dropped) {
      stop(n_whole, " subsets would be dropped entirely for ", run_label,
           "; fix the inputs (see 11a_check_qrf_inputs.R) or rerun with --allow-dropped-subsets.")
    }
  }
}
# `family` is kept but not required: rows without it are only excluded from section 14
j.d <- na.omit(j.d[, c(model_cols, "family"), with = FALSE], cols = model_cols)
if (nrow(j.d) < 10) {
  stop("Not enough complete rows for ", run_label, " (", nrow(j.d), ")")
}
ds_before <- scheme_summary(j.d)
j.d <- cap_schemes(j.d, scheme_cap)
if (nrow(j.d) > max_rows) {                       # testing only: whole schemes
  set.seed(1)
  sk <- j.d[, .N, by = .(data.set, scheme_id)][sample(.N)][cumsum(N) <= max_rows]
  j.d <- j.d[sk[, .(data.set, scheme_id)], on = .(data.set, scheme_id)]
}
ds_after <- scheme_summary(j.d)
cap_table <- merge(ds_before, ds_after, by = "data.set", suffixes = c("_full", ""))
setorder(cap_table, -n_schemes_full)
cat("  rows:", nrow(j.d), "(", sum(cap_table$n_rows_full), "before cap ) | data sets:",
    uniqueN(j.d$data.set), "| subsets:", sum(cap_table$n_schemes), "(",
    sum(cap_table$n_schemes_full), "before cap )\n")
if (is.finite(scheme_cap) && any(cap_table$n_schemes < cap_table$n_schemes_full)) {
  cat("  subsets per data set capped at", scheme_cap, "(evenly spaced along `samples`):\n")
  print(cap_table[n_schemes < n_schemes_full,
                  .(data.set, subsets = paste(n_schemes_full, "->", n_schemes),
                    samples_full = paste(samples_min_full, "-", samples_max_full),
                    samples_kept = paste(samples_min, "-", samples_max))])
}
if (do_lofo) {
  cat("  families:", paste(names(table(j.d$family, useNA = "ifany")),
                           table(j.d$family, useNA = "ifany"), sep = " = ", collapse = ", "), "\n")
}


#===============================================================================*
# 5. Designed grouped split and folds ----
#===============================================================================*
log_step("Designed train/test split")
j.split <- designed_group_split(j.d, group = "data.set", predictors = predictor_vec,
                                outcome = "value", test_prop = split_test_prop,
                                n_pc = split_n_pc, protect = split_protect)
j.train <- training(j.split)
j.test  <- testing(j.split)

n_folds <- min(cv_folds, uniqueN(j.train$data.set))
set.seed(1)
j.folds <- group_vfold_cv(j.train, group = "data.set", v = n_folds)

split_report <- split_coverage_report(j.split, "data.set", c(predictor_vec, "value"))
cat("  test data sets:     ", paste(attr(j.split, "test_groups"), collapse = ", "), "\n")
cat("  protected (train):  ", paste(attr(j.split, "protected_groups"), collapse = ", "), "\n")
cat("  median share of training range spanned by test data sets:",
    round(median(split_report$range_spanned, na.rm = TRUE), 2), "\n")
if (any(split_report$test_in_range < 1, na.rm = TRUE)) {
  cat("  test data sets outside the training range on:",
      paste(split_report[test_in_range < 1, feature], collapse = ", "), "\n")
}
save_out(list(test_groups      = attr(j.split, "test_groups"),
              protected_groups = attr(j.split, "protected_groups"),
              profile          = attr(j.split, "profile"),
              coverage         = split_report,
              settings         = list(test_prop = split_test_prop, n_pc = split_n_pc,
                                      protect = split_protect, cv_folds = n_folds)),
         "_split.rds")

stopifnot(length(intersect(j.train$data.set, j.test$data.set)) == 0,
          length(intersect(paste(j.train$data.set, j.train$scheme_id),
                           paste(j.test$data.set, j.test$scheme_id))) == 0)
cat("  training:", nrow(j.train), "rows,", uniqueN(j.train$data.set), "data sets\n")
cat("  test:    ", nrow(j.test), "rows,", uniqueN(j.test$data.set), "data sets\n")

## --- Run record (read by the summary script) ----
run_info <- list(
  run_id   = run_id,
  taxon    = taxon_label,
  metric   = metric_name,
  status   = "started",
  started  = Sys.time(),
  finished = NA,
  n_rows   = nrow(j.d),   n_datasets = uniqueN(j.d$data.set),
  n_train  = nrow(j.train), n_datasets_train = uniqueN(j.train$data.set),
  n_test   = nrow(j.test),  n_datasets_test  = uniqueN(j.test$data.set),
  max_rows = max_rows,
  scheme_cap = scheme_cap,
  n_schemes = sum(cap_table$n_schemes), n_schemes_full = sum(cap_table$n_schemes_full),
  n_rows_full = sum(cap_table$n_rows_full),
  n_rows_dropped_na = sum(drop_report$n_dropped), n_schemes_dropped_na = sum(drop_report$whole_subset),
  test_datasets      = attr(j.split, "test_groups"),
  protected_datasets = attr(j.split, "protected_groups"),
  settings = list(n_trees = n_trees, main_level = main_level,
                  nominal_levels = nominal_levels, imp_taus = imp_taus,
                  ale_taus = ale_taus, ale_bins = ale_bins,
                  tune_bagging = tune_bagging, imp_min_rel = imp_min_rel,
                  imp_min_consistency = imp_min_consistency, lodo = do_lodo,
                  lofo = do_lofo, date_vars = date_vars),
  versions = list(R = R.version.string,
                  ranger = as.character(packageVersion("ranger")))
)
save_out(run_info, "_run_info.rds")


#===============================================================================*
# 6. Recipe and model-space data ----
#===============================================================================*
rf_recipe <- recipe(value ~ ., data = j.train) %>%
  update_role(everything(), new_role = "unused") %>%
  update_role(all_of(predictor_vec), new_role = "predictor") %>%
  update_role(value, new_role = "outcome") %>%
  update_role_requirements("unused", bake = FALSE) %>%   # ID columns not needed for new data
  step_zv(all_predictors())
rec <- prep(rf_recipe, training = j.train)

final_predictors <- summary(rec) %>% filter(role == "predictor") %>% pull(variable)
if (length(final_predictors) == 0) {
  stop("No predictors remain after preprocessing for ", run_label)
}
cat("  predictors (", length(final_predictors), "): ",
    paste(final_predictors, collapse = ", "), "\n", sep = "")

to_model_space <- function(df) as.data.frame(bake(rec, new_data = df, all_predictors()))
X_train <- to_model_space(j.train); y_train <- j.train$value; g_train <- j.train$data.set
X_test  <- to_model_space(j.test);  y_test  <- j.test$value;  g_test  <- j.test$data.set
stopifnot(nrow(X_train) == nrow(j.train), nrow(X_test) == nrow(j.test))


#===============================================================================*
# 7. Tune hyperparameters: grouped CV on the training data sets, CRPS ----
#     Candidates: (mtry, min_n) space-filling grid x bagging unit (row / group).
#     Every fit runs sequentially with all threads. The out-of-fold results of
#     the selected candidate are reused as the grouped-CV evaluation (section 10).
#===============================================================================*
grid_fun <- if (exists("grid_space_filling")) grid_space_filling else grid_latin_hypercube
set.seed(1)
tune_grid  <- as.data.frame(grid_fun(mtry(range = c(1L, length(final_predictors))),
                                     min_n(), size = tune_grid_size))
candidates <- merge(unique(tune_grid[, c("mtry", "min_n")]),
                    data.frame(bagging = tune_bagging), by = NULL)
candidates$cand <- seq_len(nrow(candidates))

fold_idx <- lapply(j.folds$splits, function(s)
  list(an = s$in_id, as = setdiff(seq_len(nrow(j.train)), s$in_id)))
n_fits <- nrow(candidates) * length(fold_idx)
log_step("Tuning:", nrow(candidates), "candidates x", length(fold_idx), "folds =",
         n_fits, "fits (", n_cores, "threads each )")

tune_eval  <- vector("list", nrow(candidates))
tune_folds <- vector("list", n_fits)
i_fit <- 0
for (ci in seq_len(nrow(candidates))) {
  cand <- candidates[ci, ]
  tune_eval[[ci]] <- vector("list", length(fold_idx))
  for (k in seq_along(fold_idx)) {
    i_fit <- i_fit + 1
    an <- fold_idx[[k]]$an
    as <- fold_idx[[k]]$as
    t0 <- Sys.time()
    fit <- fit_qrf(X_train[an, , drop = FALSE], y_train[an], g_train[an],
                   cand$mtry, cand$min_n, cand$bagging, seed = 100 * ci + k, trees = n_trees_tune)
    ev <- evaluate_heldout(fit, X_train[as, , drop = FALSE], y_train[as], g_train[as],
                           y_ref = y_train[an], fold = j.folds$id[k])
    rm(fit); invisible(gc(verbose = FALSE))
    tune_eval[[ci]][[k]] <- ev
    tune_folds[[i_fit]]  <- data.table(cand = ci, fold = j.folds$id[k], crps = ev$crps,
                                       rmse = ev$rmse, coverage = ev$coverage_main,
                                       width = ev$width_main)
    log_step(sprintf("  fit %2d/%d | mtry = %2d, min_n = %2d, %-5s | %s: CRPS = %.4f, RMSE = %.4f, coverage = %.2f (%.0f s)",
                     i_fit, n_fits, cand$mtry, cand$min_n, cand$bagging, j.folds$id[k],
                     ev$crps, ev$rmse, ev$coverage_main,
                     as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}
tune_folds   <- rbindlist(tune_folds)
tune_summary <- tune_folds[, .(crps = mean(crps), crps_se = sd(crps) / sqrt(.N),
                               rmse = mean(rmse), coverage = mean(coverage),
                               width = mean(width)), by = cand]
tune_summary <- merge(as.data.table(candidates), tune_summary, by = "cand")
setorder(tune_summary, crps)
best_ci <- tune_summary$cand[1]
j.best  <- candidates[best_ci, c("mtry", "min_n", "bagging")]

cat("\n  Tuning results (mean over folds, sorted by CRPS):\n")
print(tune_summary[, .(mtry, min_n, bagging, crps = signif(crps, 4), crps_se = signif(crps_se, 2),
                       rmse = signif(rmse, 4), coverage = round(coverage, 3),
                       width = signif(width, 3))])
bag_cmp <- tune_summary[, .SD[1], by = bagging]
if (nrow(bag_cmp) > 1) {
  cat("  Best per bagging unit: ",
      paste0(bag_cmp$bagging, " CRPS = ", signif(bag_cmp$crps, 4), collapse = " | "), "\n")
}
cat("  Selected: mtry =", j.best$mtry, "| min_n =", j.best$min_n,
    "| bagging =", j.best$bagging, "\n")
save_out(list(candidates = candidates, folds = tune_folds, summary = tune_summary,
              best = j.best), "_tuning.rds")


#===============================================================================*
# 8. Final fit (the only full-data fit) -- saved immediately ----
#     Saved as list(model = ranger object, recipe = prepped recipe, best, predictors).
#     The recipe only drops zero-variance predictors, so for new data:
#       m <- readRDS(file)
#       predict(m$model, data = as.data.frame(new_data)[, m$predictors],
#               type = "quantiles", quantiles = c(.05, .5, .95))$predictions
#===============================================================================*
log_step("Final fit on the training set (", j.best$bagging, "bagging )")
eng <- fit_qrf(X_train, y_train, g_train, j.best$mtry, j.best$min_n, j.best$bagging, seed = 1)
save_out(list(model = eng, recipe = rec, best = j.best, predictors = final_predictors), ".rds")

log_step("Point metrics on the test set")
test_point_metrics <- point_metrics(y_test, predict(eng, data = X_test, num.threads = n_cores)$predictions)
print(test_point_metrics)
save_out(test_point_metrics, "_metrics.rds")


#===============================================================================*
# 10. Predictive evaluation ----
#===============================================================================*
log_step("Evaluating predictive distributions on the test set")
q_test <- predict_q(eng, X_test, eval_probs)

test_intervals       <- interval_metrics(y_test, q_test, eval_probs, nominal_levels)
test_mace            <- mean(abs(test_intervals$coverage_error))
test_pinball         <- pinball_metrics(y_test, q_test, eval_probs, y_train)
test_pinball_summary <- summarise_pinball(test_pinball)
set.seed(1)
test_pit <- pit_values(y_test, q_test, pit_step)

test_pred <- tibble(
  data.set  = g_test,
  scheme_id = j.test$scheme_id,
  value     = y_test,
  q_median  = get_q(q_test, eval_probs, 0.5),
  lo_main   = get_q(q_test, eval_probs, a_main / 2),
  hi_main   = get_q(q_test, eval_probs, 1 - a_main / 2),
  pit       = test_pit
)
test_by_group <- test_pred %>%
  group_by(data.set) %>%
  summarise(coverage = mean(value >= lo_main & value <= hi_main),
            n = n(), .groups = "drop")
test_ds <- dataset_scores(row_scores(y_test, q_test, eval_probs, main_level,
                                     quantile(y_train, eval_probs, names = FALSE)), g_test)

cat("\n  TEST SET (", uniqueN(g_test), "data sets ) -- main interval level:\n")
print(as.data.table(test_intervals)[abs(nominal - main_level) < 1e-9])
cat("  mean absolute coverage error:", round(test_mace, 3), "\n")
print(as.data.table(test_pinball_summary))
cat("  by test data set:\n")
print(test_ds[order(coverage), .(data.set, n, coverage = round(coverage, 3), crps = signif(crps, 3),
                                 skill = round(skill, 3), bias = signif(bias, 3))])

## --- Grouped CV: out-of-fold results of the selected tuning candidate ----
log_step("Grouped CV evaluation (selected candidate, from tuning fits)")
cv_res        <- tune_eval[[best_ci]]
cv_intervals  <- purrr::map_dfr(cv_res, "metrics")
cv_pit        <- purrr::map_dfr(cv_res, "pit")
cv_pinball    <- purrr::map_dfr(cv_res, "pinball")
cv_by_dataset <- rbindlist(purrr::map(cv_res, "by_dataset"))
cv_summary    <- cv_intervals %>%
  group_by(nominal) %>%
  summarise(across(c(coverage, coverage_error, mean_width, interval_score),
                   list(mean = mean, sd = sd)), .groups = "drop")
cv_mace <- cv_intervals %>%
  group_by(fold) %>%
  summarise(mace = mean(abs(coverage_error)), .groups = "drop")
cv_fold_table <- as.data.table(cv_intervals)[abs(nominal - main_level) < 1e-9,
                                             .(fold, coverage, mean_width, interval_score)
][as.data.table(cv_mace), on = "fold"
][as.data.table(cv_pinball)[, .(crps = 2 * mean(pinball),
                                skill = 1 - mean(pinball) / mean(pinball_ref)), by = fold],
  on = "fold"]
rm(tune_eval)

cat("\n  GROUPED CV (", uniqueN(g_train), "training data sets,", length(cv_res),
    "folds, out-of-fold) -- main interval level:\n")
print(as.data.table(cv_summary)[abs(nominal - main_level) < 1e-9])
cat("  per fold:\n")
print(cv_fold_table[, .(fold, coverage = round(coverage, 3), mace = round(mace, 3),
                        width = signif(mean_width, 3), crps = signif(crps, 4),
                        skill = round(skill, 3))])
cat("  by training data set (out-of-fold):\n")
print(cv_by_dataset[order(coverage), .(fold, data.set, n, coverage = round(coverage, 3),
                                       crps = signif(crps, 3), skill = round(skill, 3),
                                       bias = signif(bias, 3))])
cat("  note: CV folds were also used for tuning; CV results are slightly optimistic.\n")

## --- Save ----
pit_df <- bind_rows(
  tibble(set = "Test set", data.set = g_test, pit = test_pit),
  if (!is.null(cv_pit)) transmute(cv_pit, set = "Grouped CV (out-of-fold)", data.set, pit)
)
save_out(list(
  probs                = eval_probs,
  test_intervals       = test_intervals,
  test_mace            = test_mace,
  test_by_group        = test_by_group,
  test_pinball         = test_pinball,
  test_pinball_summary = test_pinball_summary,
  test_pred            = test_pred,
  pit                  = pit_df,
  cv_intervals         = cv_intervals,
  cv_summary           = cv_summary,
  cv_mace              = cv_mace,
  cv_pinball           = cv_pinball,
  cv_fold_table        = cv_fold_table,
  cv_by_dataset        = cv_by_dataset,
  test_by_dataset      = test_ds
), "_intervals.rds")
save_fig(plot_pit(pit_df, title = paste("PIT:", run_label)),
         "_pit.png", width = if (is.null(cv_pit)) 4 else 7, height = 3.5)
save_fig(plot_pinball_skill(test_pinball, cv_pinball,
                            title = paste("Pinball skill:", run_label)),
         "_pinball.png", width = 6, height = 4)
rm(q_test)


#===============================================================================*
# 11. Permutation importance ----
#     Primary: grouped CV folds -- the selected model refitted on each fold (same
#     seed as its tuning fit, i.e. the identical model), permutations evaluated on
#     the held-out data sets; every training data set is held out once.
#     Confirmatory: the final model on the test data sets.
#     Rows are capped per data set, so that no data set dominates the losses.
#     Uncertainty: bootstrap over held-out data sets.
#===============================================================================*
cap_rows_per_group <- function(g, max_n, seed = 1) {
  set.seed(seed)
  idx <- split(seq_along(g), g)
  sort(unlist(lapply(idx, function(i)
    if (length(i) > max_n) i[sample.int(length(i), max_n)] else i), use.names = FALSE))
}
imp_losses <- function(fit, X, y, g) {
  group_losses(y, predict_q(fit, X, imp_probs), g, imp_probs, imp_taus, main_level)
}
# Summed losses per held-out data set: baseline and after permuting each variable
perm_group_losses <- function(fit, X, y, g, vars) {
  G_base <- imp_losses(fit, X, y, g)
  G_perm <- lapply(vars, function(v) {
    Reduce(`+`, lapply(seq_len(imp_n_perm), function(b) {
      Xp <- X
      Xp[[v]] <- sample(Xp[[v]])
      imp_losses(fit, Xp, y, g)
    })) / imp_n_perm
  })
  names(G_perm) <- vars
  list(G_base = G_base, G_perm = G_perm, n = as.vector(table(g)[rownames(G_base)]))
}
summarise_importance <- function(L, source) {
  G_base <- L$G_base; n_g <- L$n; ng <- nrow(G_base)
  set.seed(1)
  B_w       <- t(rmultinom(imp_n_boot, size = ng, prob = rep(1, ng)))
  boot_base <- B_w %*% G_base
  res <- rbindlist(lapply(names(L$G_perm), function(v) {
    Gp   <- L$G_perm[[v]]
    relb <- (B_w %*% Gp) / boot_base - 1
    data.table(
      variable   = v,
      mkey       = colnames(G_base),
      rel        = colSums(Gp) / colSums(G_base) - 1,
      rel_dw     = colSums(Gp / n_g) / colSums(G_base / n_g) - 1,   # each data set equal
      diff       = (colSums(Gp) - colSums(G_base)) / sum(n_g),
      ci_lo      = apply(relb, 2, quantile, 0.025),
      ci_hi      = apply(relb, 2, quantile, 0.975),
      p_boot     = (1 + colSums(relb <= 0)) / (1 + imp_n_boot),
      groups_pos = colSums(Gp / G_base - 1 > 0),
      n_groups   = ng
    )
  }))
  res[, measure := sub("_.*", "", mkey)]
  res[, tau := suppressWarnings(as.numeric(sub("^pinball_", "", mkey)))]
  res[measure == "winkler", tau := main_level]
  res[, mkey := NULL]
  res[, p_adj := p.adjust(p_boot, method = "BH"), by = .(measure, tau)]
  res[, relevant := rel >= imp_min_rel & groups_pos / n_groups >= imp_min_consistency]
  if (ng < imp_min_groups_p) res[, c("p_boot", "p_adj") := NA_real_]
  res[, source := source]
  res[]
}
print_importance <- function(res) {
  print(res[measure == "crps"][order(-rel),
                               .(variable, rel = round(rel, 3), rel_dw = round(rel_dw, 3), ci_lo = round(ci_lo, 3),
                                 ci_hi = round(ci_hi, 3), p_adj = round(p_adj, 3),
                                 worse_in = paste0(groups_pos, "/", n_groups), relevant)])
}
imp_vars <- final_predictors

## --- Primary: grouped CV folds ----
log_step("Permutation importance on grouped CV folds (primary;", length(fold_idx), "refits )")
L_cv <- list(G_base = NULL, G_perm = setNames(vector("list", length(imp_vars)), imp_vars), n = NULL)
for (k in seq_along(fold_idx)) {
  t0 <- Sys.time()
  an  <- fold_idx[[k]]$an
  as  <- fold_idx[[k]]$as
  as  <- as[cap_rows_per_group(g_train[as], imp_max_rows_per_ds, seed = k)]
  fit <- fit_qrf(X_train[an, , drop = FALSE], y_train[an], g_train[an],
                 j.best$mtry, j.best$min_n, j.best$bagging, seed = 100 * best_ci + k)
  Lk  <- perm_group_losses(fit, X_train[as, , drop = FALSE], y_train[as], g_train[as], imp_vars)
  rm(fit); invisible(gc(verbose = FALSE))
  L_cv$G_base <- rbind(L_cv$G_base, Lk$G_base)
  for (v in imp_vars) L_cv$G_perm[[v]] <- rbind(L_cv$G_perm[[v]], Lk$G_perm[[v]])
  L_cv$n <- c(L_cv$n, Lk$n)
  log_step(sprintf("  %s: %d held-out data sets, %d rows (%.0f s)", j.folds$id[k],
                   nrow(Lk$G_base), sum(Lk$n),
                   as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
imp_res <- summarise_importance(L_cv, "Grouped CV (held-out training data sets)")

## --- Confirmatory: final model on the test set ----
log_step("Permutation importance on the test set (confirmatory)")
te_rows  <- cap_rows_per_group(g_test, imp_max_rows_per_ds)
L_te     <- perm_group_losses(eng, X_test[te_rows, , drop = FALSE], y_test[te_rows],
                              g_test[te_rows], imp_vars)
imp_test <- summarise_importance(L_te, "Test set")

cat("\n  relevant =", relevance_label(), "\n")
cat("  GROUPED CV (", nrow(L_cv$G_base), "held-out data sets ) -- CRPS importance:\n")
print_importance(imp_res)
cat("  TEST SET (", nrow(L_te$G_base), "data sets; confirmatory ) -- CRPS importance:\n")
if (nrow(L_te$G_base) < imp_min_groups_p) {
  cat("  (p-values not reported: fewer than", imp_min_groups_p, "held-out data sets)\n")
}
print_importance(imp_test)
imp_agree <- merge(imp_res[measure == "crps", .(variable, rel_cv = rel)],
                   imp_test[measure == "crps", .(variable, rel_test = rel)], by = "variable")
cat("  Agreement CV vs. test (Spearman, CRPS importance):",
    round(cor(imp_agree$rel_cv, imp_agree$rel_test, method = "spearman"), 2), "\n")

save_out(list(
  settings = list(probs = imp_probs, taus = imp_taus, level = main_level,
                  n_perm = imp_n_perm, n_boot = imp_n_boot,
                  max_rows_per_ds = imp_max_rows_per_ds,
                  n_groups_cv = nrow(L_cv$G_base), n_groups_test = nrow(L_te$G_base),
                  min_rel = imp_min_rel, min_consistency = imp_min_consistency),
  results      = imp_res,       # primary (grouped CV)
  results_test = imp_test,      # confirmatory (test set)
  losses_cv    = L_cv,
  losses_test  = L_te
), "_quantile_importance.rds")

fig_h <- max(3, 0.25 * length(imp_vars) + 1.5)
save_fig(plot_importance_heatmap(imp_res, "rel",
                                 title = paste("Quantile importance:", run_label)),
         "_imp_heatmap.png", width = 8, height = fig_h)
save_fig(plot_importance_heatmap(imp_res, "share",
                                 title = paste("Quantile importance (share per tau):", run_label)),
         "_imp_heatmap_share.png", width = 8, height = fig_h)
save_fig(plot_importance_distribution(imp_res, main_level,
                                      title = paste("Distribution-level importance (grouped CV):", run_label)),
         "_imp_distribution.png", width = 9, height = fig_h)
save_fig(plot_importance_distribution(imp_test, main_level,
                                      title = paste("Distribution-level importance (test set):", run_label)),
         "_imp_distribution_test.png", width = 9, height = fig_h)
rm(L_cv, L_te)


#===============================================================================*
# 12. Accumulated local effects (training data, no refits) ----
#     Replaces the Pearson-correlation "direction" of effects.
#===============================================================================*
log_step("Accumulated local effects")
set.seed(1)
ale_rows <- seq_len(nrow(X_train))
if (length(ale_rows) > ale_max_rows) ale_rows <- sort(sample(ale_rows, ale_max_rows))
X_ale <- X_train[ale_rows, , drop = FALSE]
y_ale <- y_train[ale_rows]
g_ale <- as.integer(factor(g_train[ale_rows]))
G_ale <- max(g_ale)

ale_out       <- c(paste0("q", ale_taus), "width")
ale_out_label <- c(paste0("tau = ", ale_taus),
                   paste0("width (", ale_taus[length(ale_taus)], " - ", ale_taus[1], ")"))
ale_predict <- function(X) {
  q <- predict_q(eng, X, ale_taus)
  q <- cbind(q, q[, ncol(q)] - q[, 1])
  colnames(q) <- ale_out
  q
}

set.seed(1)
ale_W <- cbind(1, rmultinom(ale_n_boot, size = G_ale, prob = rep(1, G_ale)))

ale_res <- lapply(intersect(final_predictors, names(X_ale)), function(v)
  ale_variable(v, X_ale, g_ale, ale_W, ale_predict, ale_out, ale_bins))
names(ale_res) <- intersect(final_predictors, names(X_ale))
ale_res <- Filter(Negate(is.null), ale_res)

ale_curves  <- rbindlist(lapply(ale_res, `[[`, "curves"))
ale_summary <- rbindlist(lapply(ale_res, `[[`, "summary"))

ale_summary[, class := fifelse(is.na(mono) | mono < ale_mono_threshold, "non-monotone",
                               fifelse(net > 0, "increasing", "decreasing"))]
ale_summary <- imp_res[measure == "crps", .(variable, crps_rel = rel, crps_relevant = relevant)
][ale_summary, on = "variable"]
ale_summary[crps_relevant %in% FALSE, class := "negligible"]   # below the relevance rule
ale_summary[, cor_marginal := sapply(variable, function(v)
  suppressWarnings(cor(X_ale[[v]], y_ale)))]          # for comparison only

base_pred <- ale_predict(X_ale)
sum_eff   <- Reduce(`+`, lapply(ale_res, `[[`, "row_eff"))
ale_additivity <- data.table(
  output = ale_out,
  r2_main_effects = sapply(ale_out, function(o) {
    f <- base_pred[, o]
    1 - var(f - mean(f) - sum_eff[, o]) / var(f)
  }))
print(ale_additivity)
print(ale_summary[output == ale_out[ceiling(length(ale_taus) / 2)]][order(-amplitude),
                                                                    .(variable, amplitude = signif(amplitude, 3), net = signif(net, 3),
                                                                      mono = round(mono, 2), class, cor_marginal = round(cor_marginal, 2))])

ale_order <- intersect(imp_res[measure == "crps"][order(-rel), variable], names(ale_res))
ale_rug   <- rbindlist(lapply(names(ale_res), function(v)
  data.table(variable = v, x = quantile(X_ale[[v]], seq(0.1, 0.9, 0.1), names = FALSE))))

save_out(list(
  settings   = list(taus = ale_taus, bins = ale_bins, n_rows = nrow(X_ale),
                    n_groups = G_ale, n_boot = ale_n_boot,
                    mono_threshold = ale_mono_threshold),
  curves     = ale_curves,
  summary    = ale_summary,
  additivity = ale_additivity,
  order      = ale_order,       # panel order of the ALE figures
  rug        = ale_rug          # predictor deciles of the ALE rows (figure rug)
), "_ale.rds")
n_col <- min(4, length(ale_order))
ale_h <- 2.2 * ceiling(length(ale_order) / n_col) + 1.2

save_fig(plot_ale(ale_curves, ale_rug, ale_order, ale_out[seq_along(ale_taus)],
                  ale_out, ale_out_label,
                  title = paste("ALE of predicted quantiles:", run_label),
                  subtitle = paste0("Parallel curves = location shift; diverging curves = change in spread.\n",
                                    "Ribbons: 95% group-bootstrap CI; rug: deciles of the predictor")),
         "_ale_quantiles.png", width = 2.6 * n_col, height = ale_h)
save_fig(plot_ale(ale_curves, ale_rug, ale_order, "width", ale_out, ale_out_label,
                  title = paste("ALE of prediction-interval width:", run_label),
                  subtitle = "Above zero = wider predictive distribution. Ribbons: 95% group-bootstrap CI"),
         "_ale_width.png", width = 2.6 * n_col, height = ale_h)
save_fig(plot_ale_summary(ale_summary, ale_order,
                          c(ale_out[ceiling(length(ale_taus) / 2)], "width"),
                          ale_out, ale_out_label, ale_mono_threshold,
                          title = paste("Size and shape of effects:", run_label)),
         "_ale_summary.png", width = 8, height = max(3, 0.25 * length(ale_order) + 1.5))


#===============================================================================*
# 13. Data-set contribution and influence ----
#     (1) distribution of observations among data sets; effective number of data sets
#     (2) pooled (row-weighted) vs. per-data-set (each data set equal) performance
#     (3) leave-one-data-set-out (LODO) refits: out-of-sample performance of every
#         training data set, and how much removing it changes the model
#         (test predictions, test CRPS, ALE curves of the key predictors)
#     (4) equal-weight sensitivity fit (every data set equal total influence)
#     A replicate fit of the full model (different seed) gives the noise floor:
#     a data set is only influential if removing it changes more than refitting does.
#===============================================================================*
log_step("Data-set contribution and influence")

## --- (1) Contribution ----
test_groups <- attr(j.split, "test_groups")
contrib <- cap_table[, .(data.set, n_rows_full, n_schemes_full, n_rows, n_schemes)]
contrib[, `:=`(programme = programme_of(data.set),
               set       = fifelse(data.set %in% test_groups, "test", "train"),
               protected = data.set %in% attr(j.split, "protected_groups"),
               share_all = n_rows / sum(n_rows), share_full = n_rows_full / sum(n_rows_full))]
contrib[set == "train", share_train := n_rows / sum(n_rows)]
if (!is.null(ds_upstream_file) && file.exists(ds_upstream_file)) {
  contrib <- merge(contrib, fread(ds_upstream_file), by = "data.set", all.x = TRUE)
}
setorder(contrib, -n_rows)

conc_row <- function(label, n_rows, n_schemes) {
  data.table(data_sets = label, n_datasets = length(n_rows), n_rows = sum(n_rows),
             n_schemes = sum(n_schemes),
             n_eff_rows = n_eff(n_rows), n_eff_schemes = n_eff(n_schemes),
             top3_share = sum(head(sort(n_rows, decreasing = TRUE), 3)) / sum(n_rows),
             largest_share = max(n_rows) / sum(n_rows))
}
concentration <- rbind(
  conc_row("all, before cap", contrib$n_rows_full, contrib$n_schemes_full),
  conc_row("all", contrib$n_rows, contrib$n_schemes),
  conc_row("train", contrib[set == "train", n_rows], contrib[set == "train", n_schemes]))
programme_tab <- contrib[, .(n_datasets = .N, n_rows = sum(n_rows), n_schemes = sum(n_schemes)),
                         by = programme][, share_rows := n_rows / sum(n_rows)][order(-n_rows)]

cat("\n  Rows per data set:\n")
print(contrib[, .(data.set, programme, set, subsets = paste(n_schemes_full, "->", n_schemes),
                  n_rows, share_full = round(share_full, 3), share = round(share_all, 3))])
cat("  Concentration (effective number of data sets = 1 / sum(share^2)):\n")
print(concentration[, .(data_sets, n_datasets, n_schemes, n_rows, n_eff_rows = round(n_eff_rows, 1),
                        top3_share = round(top3_share, 3), largest_share = round(largest_share, 3))])

## --- Reference quantities of the main model on the test set ----
ref_q_test <- quantile(y_train, imp_probs, names = FALSE)
full_q     <- predict_q(eng, X_test, imp_probs)
full_med   <- get_q(full_q, imp_probs, 0.5)
full_rs    <- row_scores(y_test, full_q, imp_probs, main_level, ref_q_test)
full_test  <- dataset_scores(full_rs, g_test)
rm(full_q)

# Predictors tracked for ALE stability: relevant ones first, by CRPS importance
top_vars <- head(imp_res[measure == "crps"][order(-relevant, -rel), variable], lodo_top_k)
med_out  <- ale_out[which.min(abs(ale_taus - 0.5))]
W_point  <- matrix(1, G_ale, 1)
full_ale <- ale_curves[output == med_out & variable %in% top_vars, .(variable, x, ale_full = ale)
][, idx := seq_len(.N), by = variable]
full_cls <- ale_summary[output == med_out & variable %in% top_vars,
                        .(variable, amp_full = amplitude,
                          class_full = fifelse(mono < ale_mono_threshold, "non-monotone",
                                               fifelse(net > 0, "increasing", "decreasing")))]

ale_of <- function(fit) {
  pf <- function(X) { q <- predict_q(fit, X, 0.5); colnames(q) <- med_out; q }
  rbindlist(lapply(top_vars, function(v) {
    r <- ale_variable(v, X_ale, g_ale, W_point, pf, med_out, ale_bins)
    r$curves[, .(variable, x, ale)][, `:=`(idx = seq_len(.N), net = r$summary$net,
                                           mono = r$summary$mono)]
  }))
}

# How different is a refitted model from the main model?
compare_fit <- function(fit) {
  q   <- predict_q(fit, X_test, imp_probs)
  rs  <- row_scores(y_test, q, imp_probs, main_level, ref_q_test)
  dsc <- dataset_scores(rs, g_test)
  ale <- full_ale[ale_of(fit), on = .(variable, idx)]
  shifts <- ale[, .(shift = max(abs(ale - ale_full)), net = net[1], mono = mono[1]), by = variable
  ][full_cls, on = "variable"]
  shifts[, `:=`(shift_rel  = shift / amp_full,
                class_new  = fifelse(mono < ale_mono_threshold, "non-monotone",
                                     fifelse(net > 0, "increasing", "decreasing")))]
  shifts[, changed := class_new != class_full]
  list(pred_shift      = mean(abs(get_q(q, imp_probs, 0.5) - full_med)) / sd(full_med),
       width_change    = mean(rs$width) / mean(full_rs$width) - 1,
       crps_change     = mean(rs$crps) / mean(full_rs$crps) - 1,
       crps_change_dw  = mean(dsc$crps) / mean(full_test$crps) - 1,
       test_scores     = dsc,
       ale_shifts      = shifts,
       ale_curves      = ale[, .(variable, x, ale)])
}

## --- Noise floor: replicate fit of the main model with a different seed ----
log_step("  Replicate fit (noise floor)")
fit_rep   <- fit_qrf(X_train, y_train, g_train, j.best$mtry, j.best$min_n, j.best$bagging, seed = 2)
cmp_floor <- compare_fit(fit_rep)
rm(fit_rep); invisible(gc(verbose = FALSE))

## --- (3) Leave one training data set out ----
train_ds <- contrib[set == "train", data.set]
lodo <- list()
if (do_lodo) {
  log_step("  Leave-one-data-set-out:", length(train_ds), "refits")
  for (i in seq_along(train_ds)) {
    d <- train_ds[i]; keep <- g_train != d
    t0 <- Sys.time()
    fit <- fit_qrf(X_train[keep, , drop = FALSE], y_train[keep], g_train[keep],
                   j.best$mtry, j.best$min_n, j.best$bagging, seed = 1000 + i)
    qd  <- predict_q(fit, X_train[!keep, , drop = FALSE], imp_probs)
    val <- dataset_scores(row_scores(y_train[!keep], qd, imp_probs, main_level,
                                     quantile(y_train[keep], imp_probs, names = FALSE)),
                          g_train[!keep])
    cmp <- compare_fit(fit)
    rm(fit, qd); invisible(gc(verbose = FALSE))
    lodo[[d]] <- list(validation = val, compare = cmp)
    log_step(sprintf(paste0("    %2d/%d %-40s %5.1f%% of training rows | own: coverage = %.2f,",
                            " skill = %6.3f | shift = %.3f SD, test CRPS %+5.1f%% (%.0f s)"),
                     i, length(train_ds), substr(short_ds(d), 1, 40),
                     100 * contrib[data.set == d, share_train], val$coverage, val$skill,
                     cmp$pred_shift, 100 * cmp$crps_change,
                     as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}

## --- (4) Equal data-set weights ----
log_step("  Equal-weight sensitivity fit")
fit_eq <- fit_qrf(X_train, y_train, g_train, j.best$mtry, j.best$min_n, j.best$bagging,
                  seed = 3, equalize = TRUE)
cmp_eq <- compare_fit(fit_eq)
rm(fit_eq); invisible(gc(verbose = FALSE))

## --- Assemble tables ----
# Out-of-sample performance of EVERY data set: test sets from the main model,
# training sets from their LODO model (or from the grouped CV if LODO is skipped)
train_perf <- if (do_lodo) {
  rbindlist(lapply(lodo, `[[`, "validation"))[, source := "LODO"]
} else {
  cv_by_dataset[, !"fold"][, source := "grouped CV"]
}
ds_perf <- rbind(full_test[, source := "test set"], train_perf, fill = TRUE)
ds_perf <- contrib[, .(data.set, programme, set, n_rows, share_all)][ds_perf, on = "data.set"]

weighting <- rbind(weighting_table(full_test, "test set"),
                   weighting_table(ds_perf, "all data sets, out-of-sample"),
                   weighting_table(cmp_eq$test_scores, "test set, equal-weight model"))

influence <- if (do_lodo) {
  rbindlist(lapply(names(lodo), function(d) {
    v <- lodo[[d]]$validation; cm <- lodo[[d]]$compare
    data.table(data.set = d, share_train = contrib[data.set == d, share_train],
               n_rows = contrib[data.set == d, n_rows],
               lodo_coverage = v$coverage, lodo_skill = v$skill, lodo_bias = v$bias,
               pred_shift = cm$pred_shift, width_change = cm$width_change,
               test_crps_change = cm$crps_change, test_crps_change_dw = cm$crps_change_dw,
               ale_max_shift_rel = max(cm$ale_shifts$shift_rel),
               ale_class_changes = sum(cm$ale_shifts$changed))
  }))
} else NULL
noise_floor <- data.table(pred_shift = cmp_floor$pred_shift, width_change = cmp_floor$width_change,
                          test_crps_change = cmp_floor$crps_change,
                          ale_max_shift_rel = max(cmp_floor$ale_shifts$shift_rel))
lodo_ale_shifts <- if (do_lodo) {
  rbindlist(lapply(names(lodo), function(d) lodo[[d]]$compare$ale_shifts[, data.set := d]))
} else NULL
eq_comparison <- data.table(
  model             = c("replicate fit (noise floor)", "equal data-set weights"),
  pred_shift        = c(cmp_floor$pred_shift, cmp_eq$pred_shift),
  width_change      = c(cmp_floor$width_change, cmp_eq$width_change),
  test_crps_change  = c(cmp_floor$crps_change, cmp_eq$crps_change),
  test_crps_change_dw = c(cmp_floor$crps_change_dw, cmp_eq$crps_change_dw),
  ale_max_shift_rel = c(max(cmp_floor$ale_shifts$shift_rel), max(cmp_eq$ale_shifts$shift_rel)),
  ale_class_changes = c(sum(cmp_floor$ale_shifts$changed), sum(cmp_eq$ale_shifts$changed)))
size_influence_rho <- if (!is.null(influence) && nrow(influence) > 2) {
  c(pred_shift = cor(influence$share_train, influence$pred_shift, method = "spearman"),
    test_crps  = cor(influence$share_train, influence$test_crps_change, method = "spearman"))
} else NULL

## --- Log ----
cat("\n  Pooled (each row counts once) vs. per data set (each data set counts once):\n")
print(weighting[, .(evaluation, weighting, n_datasets, coverage = round(coverage, 3),
                    crps = signif(crps, 4), skill = round(skill, 3))])
if (!is.null(influence)) {
  cat("  Influence of each training data set (LODO; shift in SD of test median predictions):\n")
  print(influence[order(-share_train), .(data.set = substr(short_ds(data.set), 1, 40),
                                         share = round(share_train, 3),
                                         own_cov = round(lodo_coverage, 2), own_skill = round(lodo_skill, 3),
                                         shift = round(pred_shift, 3),
                                         test_crps = sprintf("%+.1f%%", 100 * test_crps_change),
                                         ale_shift = round(ale_max_shift_rel, 2),
                                         class_changes = ale_class_changes)])
  cat("  Noise floor (replicate fit): shift =", round(noise_floor$pred_shift, 3),
      "| test CRPS", sprintf("%+.1f%%", 100 * noise_floor$test_crps_change),
      "| ALE shift =", round(noise_floor$ale_max_shift_rel, 2), "\n")
  cat("  Spearman correlation of data-set size with influence: prediction shift =",
      round(size_influence_rho["pred_shift"], 2), "| test CRPS change =",
      round(size_influence_rho["test_crps"], 2), "\n")
}
cat("  Equal-weight model vs. main model:\n")
print(eq_comparison[, .(model, pred_shift = round(pred_shift, 3),
                        test_crps = sprintf("%+.1f%%", 100 * test_crps_change),
                        test_crps_per_ds = sprintf("%+.1f%%", 100 * test_crps_change_dw),
                        ale_shift = round(ale_max_shift_rel, 2), class_changes = ale_class_changes)])

## --- Save ----
save_out(list(contribution = contrib, concentration = concentration, programmes = programme_tab,
              performance_by_dataset = ds_perf, weighting = weighting,
              influence = influence, noise_floor = noise_floor, lodo_ale_shifts = lodo_ale_shifts,
              equal_weight = eq_comparison, size_influence_rho = size_influence_rho,
              tracked_predictors = top_vars,
              ale_curves = rbind(
                full_ale[, .(variable, x, ale = ale_full, model = "main")],
                cmp_eq$ale_curves[, model := "equal weights"],
                if (do_lodo) rbindlist(lapply(names(lodo), function(d)
                  lodo[[d]]$compare$ale_curves[, model := paste("without", d)])))),
         "_datasets.rds")

## --- Figures ----
p <- ggplot(contrib, aes(x = reorder(short_ds(data.set), n_rows_full))) +
  geom_col(aes(y = n_rows_full), fill = "grey88", width = 0.75) +
  geom_col(aes(y = n_rows, fill = set), width = 0.75) +
  geom_text(aes(y = n_rows_full, label = sprintf("%.1f%% (%.1f%%)", 100 * share_all, 100 * share_full)),
            hjust = -0.1, size = 2.5) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c(train = "grey55", test = "#b2182b"), name = NULL) +
  labs(x = NULL, y = "Observations (rows)", title = paste("Observations per data set:", run_label),
       subtitle = sprintf(paste0("Grey: before the subset cap. Effective number of data sets %.1f ",
                                 "(before cap %.1f) of %d.\nLabels: share of rows after (before) the cap"),
                          concentration[data_sets == "all", n_eff_rows],
                          concentration[data_sets == "all, before cap", n_eff_rows],
                          concentration[data_sets == "all", n_datasets])) +
  theme_bw()
save_fig(p, "_dataset_contribution.png", width = 7.5, height = 0.28 * nrow(contrib) + 1.8)

if (!is.null(influence)) {
  il <- melt(influence[, .(data.set, share_train, `Shift of test median predictions (SD)` = pred_shift,
                           `Change in test CRPS (%)` = 100 * test_crps_change)],
             id.vars = c("data.set", "share_train"), variable.name = "measure")
  fl <- data.table(measure = factor(levels(il$measure), levels = levels(il$measure)),
                   floor = c(noise_floor$pred_shift, 100 * noise_floor$test_crps_change))
  p <- ggplot(il, aes(x = 100 * share_train, y = value)) +
    geom_hline(data = fl, aes(yintercept = floor), linetype = 2, colour = "grey50") +
    geom_point(size = 2) +
    geom_text(aes(label = substr(short_ds(data.set), 1, 25)), size = 2.3, vjust = -0.8,
              check_overlap = TRUE) +
    facet_wrap(~ measure, scales = "free_y") +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.25))) +
    labs(x = "Share of training rows (%)", y = NULL,
         title = paste("Influence of individual data sets (leave one out):", run_label),
         subtitle = paste0("Dashed: replicate refit of the full model (noise floor). ",
                           "Spearman (size vs. shift) = ",
                           round(size_influence_rho["pred_shift"], 2))) +
    theme_bw()
  save_fig(p, "_dataset_influence.png", width = 10, height = 4.5)
}

ale_plot_data <- rbind(
  cmp_eq$ale_curves[, .(variable, x, ale, model = "equal data-set weights", ds = "equal")],
  if (do_lodo) rbindlist(lapply(names(lodo), function(d)
    lodo[[d]]$compare$ale_curves[, .(variable, x, ale, model = "leave one data set out", ds = d)])),
  full_ale[, .(variable, x, ale = ale_full, model = "main model", ds = "main")])
ale_plot_data[, model := factor(model, levels = c("leave one data set out", "equal data-set weights", "main model"))]
p <- ggplot(ale_plot_data, aes(x = x, y = ale, group = ds, colour = model, linewidth = model,
                               linetype = model)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line() +
  scale_colour_manual(values = c(`leave one data set out` = "grey70",
                                 `equal data-set weights` = "#b2182b", `main model` = "black")) +
  scale_linewidth_manual(values = c(`leave one data set out` = 0.4,
                                    `equal data-set weights` = 0.7, `main model` = 0.8)) +
  scale_linetype_manual(values = c(`leave one data set out` = 1,
                                   `equal data-set weights` = 2, `main model` = 1)) +
  facet_wrap(~ factor(variable, levels = top_vars), scales = "free_x") +
  labs(x = NULL, y = "ALE on the median prediction", colour = NULL, linewidth = NULL, linetype = NULL,
       title = paste("Stability of the main effects:", run_label),
       subtitle = "Grey: one refit per removed training data set; red dashed: equal data-set weights") +
  theme_bw() + theme(legend.position = "bottom")
save_fig(p, "_dataset_ale_stability.png", width = 3.2 * length(top_vars) + 1, height = 4.2)


#===============================================================================*
# 14. Transfer between typology families: leave one family out ----
#     Can the model benchmark a typology built by a method it has never seen
#     (e.g. an expert typology)? Held-out rows of family f are predicted by
#       LOFO     a model trained without family f, on other data sets;
#       control  a model trained on the same data sets with the same number of
#                rows per scheme, drawn at random from all families. It differs
#                from LOFO only in having seen family f, so the size of the
#                training set does not enter the comparison;
#       full     (test set only) the main model.
#     Extrapolation cost of family f = CRPS of LOFO relative to control.
#     Primary: grouped CV folds (every training data set held out once);
#     confirmatory: the test data sets. Hyperparameters were tuned with all
#     families present, so the LOFO results are slightly optimistic.
#===============================================================================*
f_train  <- j.train$family
f_test   <- j.test$family
s_train  <- j.train$scheme_id
families <- sort(setdiff(unique(c(f_train, f_test)), NA))
run_lofo <- do_lofo && length(families) >= 2
if (do_lofo && !run_lofo) log_step("Leave-one-family-out skipped: fewer than two families")

if (run_lofo) {
  log_step("Leave-one-family-out:", paste(families, collapse = ", "))
  fam_tab <- rbind(
    j.train[, .(set = "train", n_rows = .N, n_schemes = uniqueN(scheme_id),
                n_datasets = uniqueN(data.set)), by = family],
    j.test[,  .(set = "test",  n_rows = .N, n_schemes = uniqueN(scheme_id),
                n_datasets = uniqueN(data.set)), by = family])
  print(dcast(fam_tab, family ~ set, value.var = c("n_rows", "n_datasets")))
  if (anyNA(c(f_train, f_test))) {
    cat("  rows without a family (excluded):", sum(is.na(c(f_train, f_test))), "\n")
  }
  
  # Evaluation designs: one per grouped CV fold, plus the test set
  lofo_splits <- c(
    lapply(seq_along(fold_idx), function(k) {
      as <- fold_idx[[k]]$as
      list(evaluation = "grouped CV", fold = j.folds$id[k], an = fold_idx[[k]]$an,
           X = X_train[as, , drop = FALSE], y = y_train[as], g = g_train[as], f = f_train[as])
    }),
    list(list(evaluation = "test set", fold = "test", an = seq_len(nrow(j.train)),
              X = X_test, y = y_test, g = g_test, f = f_test)))
  
  lofo_out <- list()
  n_runs   <- (2 * length(lofo_splits) + 1) * length(families)
  i_run    <- 0
  for (si in seq_along(lofo_splits)) {
    sp <- lofo_splits[[si]]
    an <- sp$an[!is.na(f_train[sp$an])]
    for (fi in seq_along(families)) {
      fam  <- families[fi]
      ev   <- which(sp$f == fam)
      keep <- f_train[an] != fam
      if (!length(ev) || all(keep) || sum(keep) < 10) {
        log_step("  skipped", sp$fold, "|", fam, ": family absent from evaluation or training rows")
        next
      }
      seed <- 5000 + 10 * si + fi
      runs <- list(LOFO = an[keep], control = matched_rows(an, s_train[an], keep, seed))
      if (sp$evaluation == "test set") runs$full <- seq_len(nrow(j.train))
      for (m in names(runs)) {
        i_run <- i_run + 1
        t0    <- Sys.time()
        tr    <- runs[[m]]
        fit   <- if (m == "full") eng else
          fit_qrf(X_train[tr, , drop = FALSE], y_train[tr], g_train[tr],
                  j.best$mtry, j.best$min_n, j.best$bagging, seed = seed)
        set.seed(seed)
        # Common skill reference: all training rows of the split (as in section 10)
        sc  <- transfer_scores(fit, sp$X[ev, , drop = FALSE], sp$y[ev], sp$g[ev], y_train[sp$an])
        rng <- outside_range(X_train[tr, , drop = FALSE], sp$X[ev, , drop = FALSE])
        rm(fit); invisible(gc(verbose = FALSE))
        tag <- list(evaluation = sp$evaluation, fold = sp$fold, family = fam, model = m)
        lofo_out[[length(lofo_out) + 1]] <- list(
          rows      = add_cols(sc$rows, tag),
          intervals = add_cols(sc$intervals, tag),
          range     = add_cols(rng, c(tag, list(n_train = length(tr), n_eval = length(ev)))))
        log_step(sprintf(paste0("  %2d/%d | %-6s | %-17s | %-7s | %6d training rows |",
                                " CRPS = %.4f, coverage = %.2f, outside range = %.2f (%.0f s)"),
                         i_run, n_runs, sp$fold, fam, m, length(tr), mean(sc$rows$crps),
                         mean(sc$rows$inside), rng[variable == "any", share],
                         as.numeric(difftime(Sys.time(), t0, units = "secs"))))
      }
    }
  }
  rm(lofo_splits)
  run_lofo <- length(lofo_out) > 0
}

if (run_lofo) {
  lofo_rows <- rbindlist(lapply(lofo_out, `[[`, "rows"))
  lofo_iv   <- rbindlist(lapply(lofo_out, `[[`, "intervals"))
  lofo_rng  <- rbindlist(lapply(lofo_out, `[[`, "range"))
  rm(lofo_out)
  keys <- c("evaluation", "family", "model")
  
  ## --- Pooled over folds (each row counts once) ----
  lofo_cal <- lofo_iv[, .(coverage       = weighted.mean(coverage, n),
                          mean_width     = weighted.mean(mean_width, n),
                          interval_score = weighted.mean(interval_score, n)),
                      by = c(keys, "nominal")]
  lofo_ds  <- lofo_rows[, .(n = .N, crps = mean(crps), coverage = mean(inside)),
                        by = c(keys, "data.set")]
  lofo_summary <- lofo_rows[, .(n_rows = .N, n_datasets = uniqueN(data.set), crps = mean(crps),
                                skill = 1 - sum(crps) / sum(crps_ref), coverage = mean(inside),
                                width = mean(width), bias = median(resid)), by = keys
  ][lofo_ds[, .(crps_dw = mean(crps)), by = keys], on = keys
  ][lofo_cal[, .(mace = mean(abs(coverage - nominal))), by = keys], on = keys
  ][lofo_rng[variable == "any", .(outside_range = weighted.mean(share, n_eval)), by = keys],
    on = keys]
  
  ## --- LOFO vs. control (paired: identical evaluation rows) ----
  lofo_contrast <- dcast(lofo_summary[model %in% c("LOFO", "control")], evaluation + family ~ model,
                         value.var = c("n_rows", "crps", "crps_dw", "coverage", "width", "mace",
                                       "skill", "bias", "outside_range"))
  lofo_contrast[, `:=`(crps_change     = crps_LOFO / crps_control - 1,
                       crps_change_dw  = crps_dw_LOFO / crps_dw_control - 1,
                       coverage_change = coverage_LOFO - coverage_control,
                       width_change    = width_LOFO / width_control - 1)]
  lofo_worse <- dcast(lofo_ds[model %in% c("LOFO", "control")],
                      evaluation + family + data.set ~ model, value.var = "crps"
  )[, .(worse_in = sum(LOFO > control), n_ds = .N), by = .(evaluation, family)]
  lofo_contrast <- lofo_worse[lofo_contrast, on = .(evaluation, family)]
  
  cat("\n  LEAVE ONE FAMILY OUT vs. size-matched control (", 100 * main_level, "% intervals):\n")
  print(lofo_contrast[order(evaluation, family),
                      .(evaluation, family, n = n_rows_LOFO,
                        crps_change        = sprintf("%+.1f%%", 100 * crps_change),
                        crps_change_per_ds = sprintf("%+.1f%%", 100 * crps_change_dw),
                        worse_in           = paste0(worse_in, "/", n_ds),
                        coverage           = sprintf("%.3f vs %.3f", coverage_LOFO, coverage_control),
                        mace               = sprintf("%.3f vs %.3f", mace_LOFO, mace_control),
                        width_change       = sprintf("%+.1f%%", 100 * width_change),
                        bias               = sprintf("%.3g vs %.3g", bias_LOFO, bias_control),
                        outside_range      = sprintf("%.2f vs %.2f", outside_range_LOFO,
                                                     outside_range_control))])
  cat("  All models (pooled over folds):\n")
  print(lofo_summary[order(evaluation, family, model),
                     .(evaluation, family, model, n_rows, n_datasets, crps = signif(crps, 4),
                       skill = round(skill, 3), coverage = round(coverage, 3), mace = round(mace, 3),
                       width = signif(width, 3), bias = signif(bias, 3),
                       outside_range = round(outside_range, 2))])
  
  save_out(list(settings      = list(families = families, n_trees = n_trees, best = j.best,
                                     main_level = main_level, nominal_levels = nominal_levels),
                families      = fam_tab,
                summary       = lofo_summary,
                contrast      = lofo_contrast,
                calibration   = lofo_cal,
                by_dataset    = lofo_ds,
                outside_range = lofo_rng,
                rows          = lofo_rows), "_lofo.rds")
  
  ## --- Figures ----
  model_cols_fig <- c(LOFO = "#b2182b", control = "grey40", full = "black")
  d <- copy(lofo_cal)[, model := factor(model, levels = names(model_cols_fig))]
  p <- ggplot(d, aes(x = nominal, y = coverage, colour = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
    geom_line() + geom_point(size = 1.2) +
    facet_grid(evaluation ~ family) +
    scale_colour_manual(values = model_cols_fig, drop = TRUE) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Nominal central interval level", y = "Empirical coverage", colour = NULL,
         title = paste("Calibration by held-out typology family:", run_label),
         subtitle = "LOFO: family not seen in training; control: same rows per scheme, all families") +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(p, "_lofo_calibration.png", width = 2.8 * length(families) + 1, height = 6.5)
  
  d <- lofo_rows[evaluation == "grouped CV"][, model := factor(model, levels = c("LOFO", "control"))]
  p <- ggplot(d, aes(x = pit)) +
    geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, length.out = 11),
                   fill = "grey70", colour = "white") +
    geom_hline(yintercept = 1, linetype = 2) +
    facet_grid(model ~ family) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "PIT", y = "Density",
         title = paste("PIT by held-out typology family (grouped CV):", run_label)) +
    theme_bw()
  save_fig(p, "_lofo_pit.png", width = 2.8 * length(families) + 1, height = 4.5)
}


## --- Mark run as complete ----
run_info$status     <- "complete"
run_info$finished   <- Sys.time()
run_info$predictors <- final_predictors
run_info$best       <- as.list(j.best)
run_info$families   <- if (run_lofo) families else character(0)
save_out(run_info, "_run_info.rds")

log_step("Processing complete:", run_label)
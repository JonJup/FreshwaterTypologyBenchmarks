################################################################################
# compile_simulation_diagnostics.R
#
# Compiles the diagnostics written by finalize() in 09_simulate_data.R into
# tables, tests and figures that document
#   A  which schemes (HMSC models) yield retained environments, and why the
#      others do not
#   B  at which stage each classification family fails, and under which
#      typology settings and manipulated variables
#   C  how plausibility depends on gamma, and how the two checks interact
#   D  whether the plausibility filter selects particular typologies
#   E  how retention would change under other tolerance rules (post hoc)
#   F  optional: production outputs (manipulation check, VP status)
#
# Units of inference. Iterations and environments of one scheme share its data,
# baselines and selection history, so they are not independent. Tests use the
# scheme as unit (rank tests on per-scheme rates, Friedman tests with schemes
# as blocks, Mantel-Haenszel with scheme strata) or GLMs with scheme-clustered
# sandwich SEs. Rates in figures give every scheme the same total weight, so
# schemes that ran many iterations do not dominate; pooled (environment-
# weighted) rates are reported alongside in the tables.
#
# Usage
#   Rscript compile_simulation_diagnostics.R --diagDir <DIAG_DIR of script 09>
#          [--outDir <dir>]     default: <diagDir>/report
#          [--modelDir <dir>]   unfitted models (.rds); schemes without
#                               diagnostics are reported as "no diagnostics"
#          [--prodDir <dir>]    production outputs of script 09 (section F)
#
# Packages: base R, data.table, ggplot2. Memory scales with the number of
# manipulated environments (one row each); only the columns used are kept.
################################################################################

suppressPackageStartupMessages({
        library(data.table)
        library(ggplot2)
})
options(width = 200, datatable.print.class = FALSE)


# ==============================================================================*
# 0. SETTINGS
# ==============================================================================*

## 0.1 Command-line arguments ---------------------------------------------------

OUT_DIR   <- "data/supplement/simulation_diagnostics/"


## 0.2 Settings of 09_simulate_data.R (must match) -------------------------------
N_TARGET_RETAINED <- 100L
GAMMA_BINS        <- c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)
FAMILIES          <- c("kmeans", "spatial_ward", "crossed_thresholds")
FAMILY_PATIENCE   <- 50L
TOL_MULTIPLIER    <- 1.5
TOL_RANGE         <- c(0.1, 0.5)
GEO_VARS          <- c("area_calcareous", "area_sediment", "area_siliceous")
ITER_STATUSES     <- c("FAIL_too_few_types", "FAIL_small_type", "FAIL_isoforest_all",
                       "FAIL_mahalanobis_all", "OK")
ENV_DESCRIPTORS   <- c("env_r2", "env_ch", "disp_mean", "disp_F", "knn_kappa", "asw")

n_bins     <- length(GAMMA_BINS) - 1L
BIN_QUOTA  <- N_TARGET_RETAINED %/% n_bins + as.integer(seq_len(n_bins) <= N_TARGET_RETAINED %% n_bins)
BIN_LABELS <- paste0("[", head(GAMMA_BINS, -1), ",", GAMMA_BINS[-1], ")")

## 0.3 Report settings -------------------------------------------------------------
FIXED_TOL   <- c(0.05, 0.10, 0.20, 0.30, 0.50)  # fixed tolerated anomalous shares (section E)
GAMMA_STEP  <- 0.05   # resolution of the gamma curves
N_QBINS     <- 5L     # quantile bins for continuous settings in binned-rate figures
MIN_CELL    <- 30L    # minimum units per plotted cell in binned-rate figures
MIN_EPV     <- 5      # minimum events per parameter for the scheme-level logistic GLM
MIN_ENV_RHO <- 5L     # minimum retained environments per scheme and family for a Spearman rho

# Scheme-level covariates; LOG_COVARIATES are log-transformed before standardising.
SCHEME_COVARIATES     <- c("n_sites", "n_taxa", "n_enz", "n_candidate_vars", "n_reference_catchments",
                           "vp_environmental", "vp_spatial", "vp_codistribution",
                           "iso_tol_share", "maha_tol_share")
SCHEME_GLM_COVARIATES <- c("n_sites", "n_taxa", "n_enz", "n_candidate_vars",
                           "vp_environmental", "vp_spatial", "iso_tol_share", "maha_tol_share")
LOG_COVARIATES        <- c("n_sites", "n_taxa", "n_reference_catchments")

# Typology settings (predictors of plausibility) and descriptors (section D) of
# an iteration, as stored in the iteration log of script 09.
TYPOLOGY_SETTINGS    <- c("n_vars", "n_types", "ld1_dominance", "type_evenness", "env_r2_orig", "vp_importance")
TYPOLOGY_DESCRIPTORS <- c("n_vars", "n_types", "type_evenness", "ld1_share", "ld1_dominance",
                          "asw_orig_sub", "env_r2_orig", "env_ch_orig", "disp_mean_orig",
                          "disp_F_orig", "knn_kappa_orig", "asw_orig", "vp_importance")
LOG_DESCRIPTORS      <- c("env_ch_orig", "disp_F_orig")

# Okabe-Ito colours
FAMILY_COLOURS  <- c(kmeans = "#0072B2", spatial_ward = "#E69F00", crossed_thresholds = "#009E73")
OUTCOME_COLOURS <- c("complete"                = "#009E73",
                     "partial"                 = "#56B4E9",
                     "no retained environment" = "#D55E00",
                     "skipped: < 4 variables"  = "#E69F00",
                     "skipped: model quality"  = "#999999",
                     "no diagnostics"          = "#000000")
STATUS_COLOURS  <- c(FAIL_too_few_types   = "#CC79A7",
                     FAIL_small_type      = "#E69F00",
                     FAIL_isoforest_all   = "#D55E00",
                     FAIL_mahalanobis_all = "#F0E442",
                     OK                   = "#009E73")

# ==============================================================================*
# 1. OUTPUT AND HELPERS
# ==============================================================================*

## 1.1 Output -------------------------------------------------------------------
for (d in c("tables", "figures")) dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT_DIR, "report.txt"), split = TRUE)
options(error = function() { while (sink.number() > 0L) sink(); quit(save = "no", status = 1L) })

TABLES     <- list()
FIGURES    <- list()
GLM_TABLES <- list()

section <- function(title) cat("\n", strrep("=", 78), "\n", title, "\n", strrep("=", 78), "\n", sep = "")

# Evaluate one analysis block (in the global environment); an error is reported
# and the report continues with the next block.
run_block <- function(label, expr) {
        cat("\n---", label, "---\n")
        tryCatch(expr, error = function(e) cat("!! Block '", label, "' failed: ", conditionMessage(e), "\n", sep = ""))
        invisible(NULL)
}

save_table <- function(x, name, show = TRUE, digits = 3L) {
        x <- as.data.table(x)
        TABLES[[name]] <<- x
        fwrite(x, file.path(OUT_DIR, "tables", paste0(name, ".csv")))
        if (show) {
                y   <- copy(x)
                num <- names(y)[vapply(y, is.double, logical(1))]
                if (length(num)) y[, (num) := lapply(.SD, signif, digits), .SDcols = num]
                cat("\n[", name, "]\n", sep = "")
                print(y, nrows = 200L)
        }
        invisible(x)
}

save_figure <- function(p, name, width = 8, height = 5) {
        FIGURES[[name]] <<- p
        f <- file.path(OUT_DIR, "figures", name)
        tryCatch({
                ggsave(paste0(f, ".pdf"), p, width = width, height = height)
                if (capabilities("cairo")) ggsave(paste0(f, ".png"), p, width = width, height = height,
                                                  dpi = 200, type = "cairo")
        }, error = function(e) cat("!! Figure", name, "failed:", conditionMessage(e), "\n"))
        invisible(p)
}

theme_set(theme_bw(base_size = 9) +
                  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                        strip.background = element_rect(fill = "grey95", colour = NA)))
scale_family_colour <- function() scale_colour_manual(values = FAMILY_COLOURS, name = NULL, drop = FALSE)
scale_family_fill   <- function() scale_fill_manual(values = FAMILY_COLOURS, name = NULL, drop = FALSE)
pct_axis            <- function() scale_y_continuous(labels = function(x) paste0(100 * x, "%"))
bin_edges           <- function() geom_vline(xintercept = GAMMA_BINS, colour = "grey85")

## 1.2 Table helpers ---------------------------------------------------------------

# Median as double: median() of an integer vector is integer for odd lengths and
# double for even lengths, which breaks grouped data.table results.
med <- function(x, na.rm = TRUE) stats::median(as.numeric(x), na.rm = na.rm)

# Append copies in which the grouping columns are set to "all", for every
# combination of cols (marginal totals).
with_totals <- function(d, cols) {
        d <- copy(d)
        for (cc in cols) if (!is.character(d[[cc]])) set(d, j = cc, value = as.character(d[[cc]]))
        out <- list(d)
        for (cc in cols) out <- c(out, lapply(out, function(x) copy(x)[, (cc) := "all"]))
        rbindlist(out)
}

# Bin a continuous setting: exact values if few distinct values, otherwise
# quantile bins represented by their median.
bin_x <- function(x, nb = N_QBINS) {
        if (uniqueN(x[!is.na(x)]) <= 12L) return(as.numeric(x))
        br <- unique(quantile(x, seq(0, 1, length.out = nb + 1L), na.rm = TRUE, names = FALSE))
        g  <- cut(x, br, include.lowest = TRUE)
        ave(as.numeric(x), g, FUN = function(v) med(v, na.rm = TRUE))
}

# Rate of a binary outcome y by group. Each trial is weighted 1 / (number of
# trials of its scheme, column n_total), so every scheme has the same total
# weight. CI from a scheme-clustered ratio estimator; the unweighted (pooled)
# rate is returned for comparison.
scheme_weighted_rate <- function(d, y, by, n = NULL) {
        s <- if (is.null(n)) {
                d[, .(k = sum(.SD[[1]]), n = .N, n_total = .SD[[2]][1]),
                  by = c(by, "model_name"), .SDcols = c(y, "n_total")]
        } else {
                d[, .(k = sum(.SD[[1]]), n = sum(.SD[[2]]), n_total = .SD[[3]][1]),
                  by = c(by, "model_name"), .SDcols = c(y, n, "n_total")]
        }
        s[, {
                K  <- k / n_total
                N  <- n / n_total
                p  <- sum(K) / sum(N)
                G  <- .N
                se <- if (G > 1L) sqrt(G / (G - 1) * sum((K - p * N)^2)) / sum(N) else NA_real_
                .(n_schemes = G, n = sum(n), pooled = sum(k) / sum(n), rate = p,
                  lo = max(0, p - 1.96 * se), hi = min(1, p + 1.96 * se))
        }, by = by]
}

## 1.3 Test helpers ----------------------------------------------------------------

# Standardise covariates in place (z_<name>); constants become 0 and are then
# dropped by fit_glm().
add_z <- function(d, vars, log_vars = character(0)) {
        for (v in vars) {
                x <- as.numeric(d[[v]])
                if (v %in% log_vars) x <- log(x)
                s <- sd(x, na.rm = TRUE)
                set(d, j = paste0("z_", v),
                    value = if (is.finite(s) && s > 0) (x - mean(x, na.rm = TRUE)) / s else 0)
        }
        invisible(d)
}

# Two-group comparison of a scheme-level covariate: medians, rank-biserial
# correlation (positive = larger in the TRUE group) and Wilcoxon rank-sum test.
compare_groups <- function(d, group, vars) {
        out <- rbindlist(lapply(vars, function(v) {
                x  <- as.numeric(d[[v]]); g <- d[[group]]
                ok <- !is.na(x) & !is.na(g); x <- x[ok]; g <- g[ok]
                n1 <- sum(g); n0 <- sum(!g)
                res <- data.table(variable = v, n_true = n1, n_false = n0,
                                  median_true = NA_real_, median_false = NA_real_, r_rank_biserial = NA_real_, p = NA_real_)
                if (n1 >= 3L && n0 >= 3L && uniqueN(x) > 1L) {
                        wt <- wilcox.test(x[g], x[!g], exact = FALSE)
                        res[, `:=`(median_true = med(x[g]), median_false = med(x[!g]),
                                   r_rank_biserial = unname(2 * wt$statistic / (n1 * n0) - 1), p = wt$p.value)]
                }
                res
        }))
        out[, p_holm := p.adjust(p, "holm")]
}

# Binomial / quasibinomial GLM with scheme-clustered sandwich SEs (CR1).
# Predictors without variation (incl. one-level factors) are dropped first.
fit_glm <- function(response, rhs, data, family, label, cluster = "model_name") {
        resp_vars <- all.vars(str2lang(response))
        rhs_vars  <- unique(unlist(lapply(rhs, function(t) all.vars(str2lang(t)))))
        d <- na.omit(data[, unique(c(resp_vars, rhs_vars, cluster)), with = FALSE])
        fac <- names(d)[vapply(d, is.factor, logical(1))]
        if (length(fac)) d[, (fac) := lapply(.SD, droplevels), .SDcols = fac]
        const <- rhs_vars[vapply(rhs_vars, function(v) uniqueN(d[[v]]) < 2L, logical(1))]
        if (length(const)) {
                rhs <- rhs[!vapply(rhs, function(t) any(all.vars(str2lang(t)) %in% const), logical(1))]
                cat(sprintf("  [%s] constant predictors dropped: %s\n", label, paste(const, collapse = ", ")))
        }
        if (!length(rhs)) rhs <- "1"
        fit <- glm(as.formula(paste(response, "~", paste(rhs, collapse = " + "))), family = family, data = d)
        b   <- coef(fit); ok <- !is.na(b)
        X   <- model.matrix(fit)[, ok, drop = FALSE]
        U   <- rowsum(X * (residuals(fit, "working") * fit$weights), d[[cluster]], reorder = FALSE)
        G   <- nrow(U); n <- nrow(X); k <- ncol(X)
        B   <- summary(fit)$cov.unscaled
        V   <- G / (G - 1) * (n - 1) / (n - k) * (B %*% crossprod(U) %*% B)
        se  <- sqrt(diag(V)); est <- b[ok]
        sep <- family$family == "binomial" && any(fit$fitted.values < 1e-8 | fit$fitted.values > 1 - 1e-8)
        if (!fit$converged) cat(sprintf("  [%s] WARNING: GLM did not converge.\n", label))
        if (sep) cat(sprintf("  [%s] WARNING: fitted probabilities of 0 or 1 (separation); extreme ORs are unreliable.\n", label))
        tab <- data.table(model = label, term = names(est), estimate = est,
                          se_cluster = se, se_naive = sqrt(diag(B) * summary(fit)$dispersion),
                          OR = exp(est), OR_lo = exp(est - 1.96 * se), OR_hi = exp(est + 1.96 * se),
                          p = 2 * pnorm(-abs(est / se)), n_obs = n, n_schemes = G)
        GLM_TABLES[[label]] <<- tab
        list(fit = fit, table = tab)
}

# Variance inflation factors of numeric predictors.
vif_table <- function(d, vars) {
        X <- na.omit(as.matrix(d[, vars, with = FALSE]))
        X <- X[, apply(X, 2, sd) > 0, drop = FALSE]
        v <- tryCatch(diag(solve(cor(X))), error = function(e) rep(NA_real_, ncol(X)))
        data.table(variable = colnames(X), vif = unname(v))
}

# Friedman test across groups with schemes as blocks (complete blocks only),
# plus paired Wilcoxon signed-rank tests (Holm) for all or adjacent pairs.
friedman_blocks <- function(s, group, value, pairs = c("all", "adjacent")) {
        pairs <- match.arg(pairs)
        w <- dcast(s, as.formula(paste("model_name ~", group)), value.var = value)
        m <- as.matrix(w[, -1L, with = FALSE])
        m <- m[stats::complete.cases(m), , drop = FALSE]
        if (nrow(m) < 3L || ncol(m) < 2L) return(NULL)
        ft <- friedman.test(m)
        pr <- if (pairs == "all") combn(colnames(m), 2L, simplify = FALSE) else
                lapply(seq_len(ncol(m) - 1L), function(i) colnames(m)[c(i, i + 1L)])
        pw <- rbindlist(lapply(pr, function(p) {
                dl <- m[, p[1]] - m[, p[2]]
                data.table(comparison = paste(p, collapse = " vs "), median_difference = med(dl),
                           share_first_higher = mean(dl > 0),
                           p = if (all(dl == 0)) NA_real_ else
                                   suppressWarnings(wilcox.test(m[, p[1]], m[, p[2]], paired = TRUE, exact = FALSE)$p.value))
        }))
        pw[, p_holm := p.adjust(p, "holm")]
        list(test = data.table(n_blocks = nrow(m), n_groups = ncol(m), chisq = unname(ft$statistic),
                               df = unname(ft$parameter), p = ft$p.value),
             pairwise = pw)
}

# Mantel-Haenszel common odds ratio of y given x across strata, with the
# Robins-Breslow-Greenland CI and the CMH chi-square test (no continuity
# correction). Strata without variation in x or y carry no information.
mh_or <- function(stratum, x, y) {
        t <- data.table(s = stratum, x = x, y = y)[, .(a = as.numeric(sum(x & y)), b = as.numeric(sum(x & !y)),
                                                       c = as.numeric(sum(!x & y)), d = as.numeric(sum(!x & !y))), by = s]
        t[, n := a + b + c + d]
        t <- t[(a + b) > 0 & (c + d) > 0 & (a + c) > 0 & (b + d) > 0]
        empty <- data.table(n_strata = nrow(t), n_units = sum(t$n), OR_MH = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_)
        if (!nrow(t)) return(empty)
        R <- t$a * t$d / t$n; S <- t$b * t$c / t$n
        if (sum(R) == 0 || sum(S) == 0) return(empty)
        P <- (t$a + t$d) / t$n; Q <- (t$b + t$c) / t$n
        or <- sum(R) / sum(S)
        v  <- sum(P * R) / (2 * sum(R)^2) + sum(P * S + Q * R) / (2 * sum(R) * sum(S)) + sum(Q * S) / (2 * sum(S)^2)
        Ea <- (t$a + t$b) * (t$a + t$c) / t$n
        Va <- (t$a + t$b) * (t$c + t$d) * (t$a + t$c) * (t$b + t$d) / (t$n^2 * (t$n - 1))
        data.table(n_strata = nrow(t), n_units = sum(t$n), OR_MH = or,
                   lo = exp(log(or) - 1.96 * sqrt(v)), hi = exp(log(or) + 1.96 * sqrt(v)),
                   p = pchisq(sum(t$a - Ea)^2 / sum(Va), 1, lower.tail = FALSE))
}


# ==============================================================================*
# 2. DATA
# ==============================================================================*
section("2. DATA")

## 2.1 Read diagnostics ---------------------------------------------------------
ENV_COLS <- c("model_name", "taxon", "iteration", "family", "gamma", "gamma_bin", "gamma_achieved",
              "share_clipped", "iso_n_exceed", "maha_n_exceed", "passed_isoforest",
              "passed_mahalanobis", "plausible", "retained")

read_dir <- function(dir, cols = NULL) {
        f <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
        if (!length(f)) stop("No .rds files in ", dir)
        rbindlist(lapply(f, function(x) {
                d <- as.data.table(readRDS(x))
                if (!is.null(cols) && ncol(d)) d <- d[, intersect(cols, names(d)), with = FALSE]
                d
        }), use.names = TRUE, fill = TRUE)
}

DIAG_DIR1 <- paste0("pulse", c("D","F", "I", "M"), "/")
DIAG_DIR2 <- "data/misc/simulation_diagnostics/"

models <- 
        iters <-
        envs <- vector(mode = "list", length = 4)

for (i in 1:4){
        models[[i]] <- read_dir(file.path(paste0(DIAG_DIR1[i],DIAG_DIR2), "models"))
        envs[[i]]   <- read_dir(file.path(paste0(DIAG_DIR1[i],DIAG_DIR2), "environments"), ENV_COLS)
        iters[[i]]  <- read_dir(file.path(paste0(DIAG_DIR1[i],DIAG_DIR2), "iterations"))
        
}
models <- rbindlist(models, use.names = TRUE)
iters  <- rbindlist(iters, use.names = TRUE)
envs   <- rbindlist(envs, use.names = TRUE)


cat(sprintf("Read %d schemes, %d iterations, %d manipulated environments.\n",
            nrow(models), nrow(iters), nrow(envs)))

## 2.2 Schemes without diagnostics ---------------------------------------------------
# if (!is.null(MODEL_DIR)) {
#         expected <- sub("\\.rds$", "", list.files(MODEL_DIR, pattern = "\\.rds$"))
#         missing  <- setdiff(expected, models$model_name)
#         cat(sprintf("Unfitted models: %d | without diagnostics: %d\n", length(expected), length(missing)))
#         if (length(missing)) {
#                 models <- rbind(models, data.table(model_name = missing, taxon = sub("_[^_]*$", "", missing),
#                                                    run_status = "NO_DIAGNOSTICS"), fill = TRUE)
#         }
# }

## 2.3 Scheme outcome ----------------------------------------------------------------
# complete = target reached; partial = >= 1 but < N_TARGET_RETAINED retained.
OUTCOME_OF_STATUS <- c(SKIPPED_model_quality     = "skipped: model quality",
                       SKIPPED_too_few_variables = "skipped: < 4 variables",
                       NO_RETAINED_ENV           = "no retained environment",
                       NO_DIAGNOSTICS            = "no diagnostics")
models[, outcome := unname(OUTCOME_OF_STATUS[run_status])]
models[run_status == "OK", outcome := fifelse(n_env_retained >= N_TARGET_RETAINED, "complete", "partial")]
models[, outcome := factor(outcome, levels = names(OUTCOME_COLOURS))]
models[, entered := run_status %in% c("OK", "NO_RETAINED_ENV")]   # schemes that entered the loop
models[, `:=`(any_retained = entered & n_env_retained > 0,
              complete     = entered & n_env_retained >= N_TARGET_RETAINED)]

## 2.4 Joins and derived columns ------------------------------------------------------
iters[, status := factor(status, levels = ITER_STATUSES)]
iters[, family := factor(family, levels = FAMILIES)]
envs[,  family := factor(family, levels = FAMILIES)]
iters[models, `:=`(n_sites = i.n_sites, n_candidate_vars = i.n_candidate_vars), on = "model_name"]
envs[models, `:=`(n_sites = i.n_sites, iso_tol_share = i.iso_tol_share,
                  maha_tol_share = i.maha_tol_share), on = "model_name"]
iters[, typology_valid := status %in% c("FAIL_isoforest_all", "FAIL_mahalanobis_all", "OK")]
iters[, n_total := .N, by = model_name]
envs[, bin := factor(BIN_LABELS[gamma_bin], levels = BIN_LABELS)]
envs[, gamma_fine := (floor(gamma / GAMMA_STEP) + 0.5) * GAMMA_STEP]

# Iteration x gamma-bin units: an open bin of an iteration with a valid typology.
# Retention happens at this level (one environment per open bin and iteration),
# so it is the unit for the plausibility-stage analyses.
units <- envs[, .(n_env = .N, n_plaus = sum(plausible), retained = any(retained), gamma_mean = mean(gamma)),
              by = .(model_name, taxon, iteration, family, gamma_bin, bin)]
units[, any_plaus := n_plaus > 0L]
units[, n_total := .N, by = model_name]
desc_cols <- intersect(c(TYPOLOGY_DESCRIPTORS, "alpha"), names(iters))
units[iters, (desc_cols) := mget(paste0("i.", desc_cols)), on = .(model_name, iteration)]
units[models, n_sites := i.n_sites, on = "model_name"]
units[, gamma_side := factor(fifelse(gamma_mean < 1, "gamma < 1", "gamma > 1"), levels = c("gamma < 1", "gamma > 1"))]

## 2.5 Integrity checks -----------------------------------------------------------------
run_block("2.5 integrity checks", {
        cnt <- models[entered == TRUE, .(model_name, n_iterations, n_env_generated, n_env_retained)]
        cnt[iters[, .N, by = model_name], n_it := i.N, on = "model_name"]
        cnt[envs[, .(N = .N, R = sum(retained)), by = model_name], `:=`(n_en = i.N, n_re = i.R), on = "model_name"]
        for (v in c("n_it", "n_en", "n_re")) set(cnt, which(is.na(cnt[[v]])), v, 0L)
        reproduced <- envs[, (iso_n_exceed < iso_tol_share * n_sites) & (maha_n_exceed < maha_tol_share * n_sites)]
        integrity <- data.table(
                check = c("schemes with an unknown run_status",
                          "schemes: iteration log length != n_iterations",
                          "schemes: environment log length != n_env_generated",
                          "schemes: retained environments != n_env_retained",
                          "schemes that entered the loop without n_sites or tolerances",
                          "environments: plausibility not reproduced from counts and tolerances",
                          "environments: retained but not plausible",
                          "iteration x bin units with > 1 retained environment",
                          "schemes: candidate variables seen in iterations != n_candidate_vars"),
                n = c(sum(is.na(models$outcome)),
                      cnt[, sum(n_it != n_iterations)],
                      cnt[, sum(n_en != n_env_generated)],
                      cnt[, sum(n_re != n_env_retained)],
                      models[entered == TRUE, sum(is.na(n_sites) | is.na(iso_tol_share) | is.na(maha_tol_share))],
                      sum(is.na(reproduced) | reproduced != envs$plausible),
                      envs[, sum(retained & !plausible)],
                      envs[retained == TRUE, .N, by = .(model_name, iteration, gamma_bin)][, sum(N > 1L)],
                      NA_integer_))
        sv <- strsplit(iters$selected_vars, ";", fixed = TRUE)
        seen <- data.table(model_name = rep(iters$model_name, lengths(sv)), variable = unlist(sv))
        seen <- unique(seen)[, .(n_seen = .N), by = model_name][models, on = "model_name", nomatch = NULL]
        integrity[.N, n := seen[, sum(n_seen != n_candidate_vars)]]
        save_table(integrity, "00_integrity")
        # The last check is informative only: a candidate never drawn is not an error.
})


# ==============================================================================*
# A. SCHEMES: WHICH PASS THE SIMULATION
# ==============================================================================*
section("A. SCHEMES")

run_block("A1 outcome by taxon", {
        a1 <- with_totals(models[, .(taxon, outcome)], "taxon")[, .N, by = .(taxon, outcome)]
        a1[, share := N / sum(N), by = taxon]
        setorder(a1, taxon, outcome)
        save_table(a1, "A1_outcome_by_taxon", show = FALSE)
        print(dcast(a1, taxon ~ outcome, value.var = "N", fill = 0L))
        
        a1b <- with_totals(models[entered == TRUE, .(taxon, stop_reason, n_iterations, n_env_retained, runtime_min)], "taxon")[
                , .(n_schemes = .N, median_iterations = med(n_iterations), median_retained = med(n_env_retained),
                    median_runtime_min = med(runtime_min), max_runtime_min = max(runtime_min)),
                by = .(taxon, stop_reason)][order(taxon, stop_reason)]
        save_table(a1b, "A1b_stop_reasons")
        
        n_lab <- models[, .N, by = taxon]
        p <- ggplot(a1[taxon != "all"], aes(taxon, share, fill = outcome)) +
                geom_col(width = 0.7) +
                geom_text(data = n_lab, aes(taxon, 1.04, label = paste0("n = ", N)), inherit.aes = FALSE, size = 2.8) +
                scale_fill_manual(values = OUTCOME_COLOURS, drop = FALSE, name = NULL) + pct_axis() +
                labs(x = NULL, y = "Share of schemes", title = "Simulation outcome per scheme")
        save_figure(p, "A1_outcome_by_taxon", 7, 4.5)
        
        p <- ggplot(models[entered == TRUE], aes(n_sites, runtime_min, colour = stop_reason)) +
                geom_point(size = 1, alpha = 0.7) + scale_x_log10() + scale_y_log10() +
                facet_wrap(~ taxon, nrow = 1) + labs(x = "Sites", y = "Runtime (min)", colour = NULL,
                                                     title = "Runtime and stop reason")
        save_figure(p, "A1b_runtime", 9, 3.5)
})

run_block("A2 scheme covariates of success", {
        ent  <- models[entered == TRUE]
        covs <- intersect(SCHEME_COVARIATES, names(ent))
        # Two hurdles: (1) any retained environment; (2) given (1), the full target.
        a2 <- rbind(data.table(hurdle = "1: >= 1 retained (TRUE) vs none", compare_groups(ent, "any_retained", covs)),
                    data.table(hurdle = "2: target reached (TRUE) vs partial", compare_groups(ent[any_retained == TRUE], "complete", covs)))
        save_table(a2, "A2_scheme_covariates_wilcoxon")
        
        gcovs <- intersect(SCHEME_GLM_COVARIATES, names(ent))
        add_z(ent, gcovs, LOG_COVARIATES)
        rhs   <- c("taxon", paste0("z_", gcovs))
        n_par <- length(gcovs) + uniqueN(ent$taxon)
        epv   <- min(sum(ent$any_retained), sum(!ent$any_retained)) / n_par
        if (epv >= MIN_EPV) {
                save_table(fit_glm("any_retained", rhs, ent, binomial(), "schemes: >= 1 retained (logit)")$table,
                           "A2b_glm_any_retained")
        } else {
                cat(sprintf("  Logistic GLM skipped: %.1f events per parameter (< %g).\n", epv, MIN_EPV))
        }
        ent[, `:=`(n_ret = pmin(n_env_retained, N_TARGET_RETAINED), n_not = N_TARGET_RETAINED - pmin(n_env_retained, N_TARGET_RETAINED))]
        save_table(fit_glm("cbind(n_ret, n_not)", rhs, ent, quasibinomial(), "schemes: share of target retained (quasi-logit)")$table,
                   "A2c_glm_share_retained")
        save_table(vif_table(ent, paste0("z_", gcovs)), "A2d_vif_scheme_covariates")
        
        wide <- ent[, c("model_name", "outcome", covs), with = FALSE]
        wide[, (covs) := lapply(.SD, as.numeric), .SDcols = covs]
        long <- melt(wide, id.vars = c("model_name", "outcome"), variable.factor = FALSE, value.name = "value")
        long[variable %in% LOG_COVARIATES, `:=`(value = log10(value), variable = paste0("log10(", variable, ")"))]
        long[, outcome := droplevels(outcome)]
        p <- ggplot(long, aes(outcome, value, colour = outcome)) +
                geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.6, alpha = 0.5) +
                facet_wrap(~ variable, scales = "free_y", nrow = 2) +
                scale_colour_manual(values = OUTCOME_COLOURS, guide = "none") +
                labs(x = NULL, y = NULL, title = "Scheme covariates by outcome") +
                theme(axis.text.x = element_text(angle = 30, hjust = 1))
        save_figure(p, "A2_scheme_covariates", 10, 5.5)
})

run_block("A3 gamma-bin coverage", {
        bin_grid <- CJ(model_name = models[entered == TRUE, model_name], gamma_bin = seq_len(n_bins))
        bin_grid[envs[retained == TRUE, .N, by = .(model_name, gamma_bin)], n_ret := i.N, on = .(model_name, gamma_bin)]
        bin_grid[is.na(n_ret), n_ret := 0L]
        bin_grid[, `:=`(quota = BIN_QUOTA[gamma_bin], bin = factor(BIN_LABELS[gamma_bin], levels = BIN_LABELS))]
        bin_grid[models, taxon := i.taxon, on = "model_name"]
        a3 <- with_totals(bin_grid, "taxon")[, .(n_schemes = .N, share_quota_met = mean(n_ret >= quota),
                                                 share_empty = mean(n_ret == 0L), mean_fill = mean(pmin(n_ret / quota, 1))),
                                             by = .(taxon, bin)][order(taxon, bin)]
        save_table(a3, "A3_bin_coverage")
        p <- ggplot(a3, aes(bin, taxon, fill = share_quota_met)) + geom_tile(colour = "white") +
                geom_text(aes(label = sprintf("%.0f%%", 100 * share_quota_met)), size = 2.8) +
                scale_fill_gradient(low = "#F7F7F7", high = "#0072B2", limits = c(0, 1), name = "Quota met") +
                labs(x = "gamma bin", y = NULL, title = "Share of schemes that filled the quota of each gamma bin")
        save_figure(p, "A3_bin_coverage", 7, 3.5)
})

run_block("A4 plausibility baselines", {
        # The tolerance is anchored to the observed sites' anomalous share s:
        # tol = clamp(TOL_MULTIPLIER * s, TOL_RANGE). At the ceiling, s >= 1/3;
        # if s >= 0.5, even the unmanipulated environment fails.
        b <- melt(models[entered == TRUE, .(model_name, taxon, outcome, any_retained, complete, iso_tol_share, maha_tol_share)],
                  measure.vars = c("iso_tol_share", "maha_tol_share"), variable.name = "check", value.name = "tol_share")
        b[, check := fifelse(check == "iso_tol_share", "Isolation Forest", "Mahalanobis")]
        b[, anchor := fifelse(tol_share <= TOL_RANGE[1] + 1e-9, sprintf("floor (observed share < %.3f)", TOL_RANGE[1] / TOL_MULTIPLIER),
                              fifelse(tol_share >= TOL_RANGE[2] - 1e-9, sprintf("ceiling (observed share >= %.3f)", TOL_RANGE[2] / TOL_MULTIPLIER),
                                      "anchored"))]
        a4 <- with_totals(b, "taxon")[, .(n_schemes = .N, share_any_retained = mean(any_retained),
                                          share_complete = mean(complete), median_tol = med(tol_share)),
                                      by = .(check, taxon, anchor)][order(check, taxon, anchor)]
        save_table(a4, "A4_tolerance_anchor")
        a4b <- with_totals(models[entered == TRUE, .(taxon, maha_ridge_used, coords_lonlat)], "taxon")[
                , .(n_schemes = .N, share_ridge = mean(maha_ridge_used, na.rm = TRUE),
                    share_lonlat = mean(coords_lonlat, na.rm = TRUE)), by = taxon]
        save_table(a4b, "A4b_ridge_and_coordinates")
        p <- ggplot(b, aes(tol_share, fill = droplevels(outcome))) + geom_histogram(binwidth = 0.02, boundary = 0) +
                facet_wrap(~ check) + scale_fill_manual(values = OUTCOME_COLOURS, name = NULL) +
                labs(x = "Tolerated share of anomalous sites", y = "Schemes",
                     title = "Anchored tolerances (floor and ceiling at the range limits)")
        save_figure(p, "A4_tolerance_anchor", 8, 4)
})

run_block("A5 near-constant variables", {
        den <- models[run_status %in% c("OK", "NO_RETAINED_ENV", "SKIPPED_too_few_variables"), .(n_schemes = .N), by = taxon]
        nc  <- models[!is.na(near_const_vars) & nzchar(near_const_vars),
                      .(variable = unlist(strsplit(near_const_vars, ";", fixed = TRUE))), by = .(model_name, taxon)]
        if (!nrow(nc)) {
                cat("No near-constant variables.\n")
        } else {
                a5 <- nc[, .N, by = .(taxon, variable)][den, on = "taxon", share_of_schemes := N / i.n_schemes][order(taxon, -N)]
                save_table(a5, "A5_near_constant_variables")
        }
})


# ==============================================================================*
# B. CLASSIFICATION FAMILIES: WHERE AND WHEN THEY FAIL
# ==============================================================================*
section("B. CLASSIFICATION FAMILIES")

run_block("B1 iteration status by family", {
        b1 <- with_totals(iters[, .(taxon, family, status)], c("taxon", "family"))[, .N, by = .(taxon, family, status)]
        b1[, share := N / sum(N), by = .(taxon, family)]
        setorder(b1, taxon, family, status)
        save_table(b1, "B1_iteration_status", show = FALSE)
        b1b <- with_totals(iters[, .(taxon, family, typology_valid, ok = status == "OK")], c("taxon", "family"))[
                , .(n_iterations = .N, p_valid_typology = mean(typology_valid),
                    p_ok_given_valid = sum(ok) / sum(typology_valid)), by = .(taxon, family)][order(taxon, family)]
        save_table(b1b, "B1b_two_stage_rates")
        
        # Within-scheme comparison of families (schemes = blocks). The
        # plausibility-stage rate is standardised over the gamma bins observed
        # for every family of the scheme, because late iterations face only the
        # hardest open bins and families are drawn unevenly over the run.
        s1 <- iters[, .(rate = mean(typology_valid)), by = .(model_name, family)]
        r  <- units[, .(rate = mean(any_plaus)), by = .(model_name, family, gamma_bin)]
        r[, n_fam := uniqueN(family), by = .(model_name, gamma_bin)]
        s2 <- r[n_fam == uniqueN(units$family), .(rate = mean(rate)), by = .(model_name, family)]
        tests <- list("typology valid" = friedman_blocks(s1, "family", "rate"),
                      ">= 1 plausible per open bin (bin-standardised)" = friedman_blocks(s2, "family", "rate"))
        tests <- tests[!vapply(tests, is.null, logical(1))]
        if (length(tests)) {
                save_table(rbindlist(lapply(names(tests), function(nm) data.table(stage = nm, tests[[nm]]$test))), "B1c_family_friedman")
                save_table(rbindlist(lapply(names(tests), function(nm) data.table(stage = nm, tests[[nm]]$pairwise))), "B1d_family_pairwise")
        }
        
        p <- ggplot(b1[taxon != "all" & family != "all"], aes(family, share, fill = status)) + geom_col(width = 0.7) +
                facet_wrap(~ taxon, nrow = 1) + scale_fill_manual(values = STATUS_COLOURS, name = NULL, drop = FALSE) +
                pct_axis() + labs(x = NULL, y = "Share of iterations", title = "Iteration outcome by classification family") +
                theme(axis.text.x = element_text(angle = 30, hjust = 1))
        save_figure(p, "B1_iteration_status", 9, 4.5)
})

run_block("B2 family abandonment", {
        # A family is no longer drawn after FAMILY_PATIENCE consecutive
        # unsuccessful iterations; the streak is counted over that family's iterations.
        setorder(iters, model_name, family, iteration)
        runs <- iters[, {
                fail <- status != "OK"
                rl   <- rle(fail)
                .(n_iter = .N, n_ok = sum(!fail),
                  max_fail_run = if (any(fail)) max(rl$lengths[rl$values]) else 0L)
        }, by = .(model_name, taxon, family)]
        setorder(iters, model_name, iteration)
        runs[, abandoned := max_fail_run >= FAMILY_PATIENCE]
        runs[models, `:=`(n_sites = i.n_sites, n_candidate_vars = i.n_candidate_vars), on = "model_name"]
        b2 <- with_totals(runs, c("taxon", "family"))[, .(n_schemes = .N, share_abandoned = mean(abandoned),
                                                          share_never_ok = mean(n_ok == 0L), median_iterations = med(n_iter)),
                                                      by = .(taxon, family)][order(taxon, family)]
        save_table(b2, "B2_family_abandonment")
        fams <- intersect(FAMILIES, as.character(unique(runs$family)))
        b2b  <- rbindlist(lapply(fams, function(f) data.table(family = f, compare_groups(runs[family == f], "abandoned", c("n_sites", "n_candidate_vars")))))
        save_table(b2b, "B2b_abandonment_covariates")
        
        runs[, x := bin_x(log10(n_sites))]
        fb <- runs[, .(n = .N, k = sum(abandoned)), by = .(family, x)]
        ci <- with(fb, { p <- k / n; z <- 1.96; d <- 1 + z^2 / n; c0 <- (p + z^2 / (2 * n)) / d
        h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d; list(lo = pmax(0, c0 - h), hi = pmin(1, c0 + h)) })
        fb[, `:=`(rate = k / n, lo = ci$lo, hi = ci$hi)]
        p <- ggplot(fb, aes(10^x, rate, colour = family)) + geom_line() + geom_pointrange(aes(ymin = lo, ymax = hi)) +
                scale_x_log10() + scale_family_colour() + pct_axis() +
                labs(x = "Sites (quantile bins, median)", y = "Share of schemes", title = "Families abandoned after the patience limit (Wilson 95% CI)")
        save_figure(p, "B2_family_abandonment", 6, 4)
})

run_block("B3 family share of retained environments", {
        fc <- paste0("n_env_retained_", FAMILIES)
        sh <- melt(models[n_env_retained > 0, c("model_name", "taxon", "n_env_retained", intersect(fc, names(models))), with = FALSE],
                   measure.vars = intersect(fc, names(models)), variable.name = "family", value.name = "n", variable.factor = FALSE)
        sh[, `:=`(family = factor(sub("n_env_retained_", "", family), levels = FAMILIES), share = n / n_env_retained)]
        b3 <- with_totals(sh, "taxon")[, .(n_schemes = .N, mean_share = mean(share), median_share = med(share),
                                           share_schemes_none = mean(n == 0)), by = .(taxon, family)][order(taxon, family)]
        save_table(b3, "B3_family_share_of_retained")
        p <- ggplot(sh, aes(family, share, fill = family)) + geom_boxplot(outlier.size = 0.5) +
                geom_hline(yintercept = 1 / length(FAMILIES), linetype = 2) + facet_wrap(~ taxon, nrow = 1) +
                scale_family_fill() + pct_axis() +
                labs(x = NULL, y = "Share of the scheme's retained environments", title = "Family composition of the benchmark (dashed: equal share)") +
                theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        save_figure(p, "B3_family_share_of_retained", 9, 4)
})

run_block("B4 typology stage", {
        b4 <- iters[, .(n_iterations = .N, p_too_few_types = mean(status == "FAIL_too_few_types"),
                        p_small_type = mean(status == "FAIL_small_type"), median_min_type_size = med(as.numeric(min_type_size))),
                    by = .(family, n_types)][order(family, n_types)]
        save_table(b4, "B4_typology_stage_by_k")
        if (all(c("n_cells", "n_classes") %in% names(iters))) {
                b4b <- iters[family == "crossed_thresholds", .(n_iterations = .N, median_cells = med(as.numeric(n_cells)),
                                                               median_types = med(as.numeric(n_types)),
                                                               p_too_few_types = mean(status == "FAIL_too_few_types"),
                                                               p_small_type = mean(status == "FAIL_small_type")), by = n_vars][order(n_vars)]
                save_table(b4b, "B4b_crossed_thresholds_merging")
        }
        d1 <- iters[, .(model_name, family, typology_valid, n_total, x = bin_x(log10(n_sites)))]
        d2 <- iters[n_types >= 2L, .(model_name, family, typology_valid, n_total, x = as.numeric(pmin(n_types, 10L)))]
        r4 <- rbind(scheme_weighted_rate(d1, "typology_valid", c("family", "x"))[, setting := "log10(sites)"],
                    scheme_weighted_rate(d2, "typology_valid", c("family", "x"))[, setting := "types k (>= 2; 10 = 10+)"])
        p <- ggplot(r4[n >= MIN_CELL], aes(x, rate, colour = family, fill = family)) +
                geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() + geom_point(size = 0.8) +
                facet_wrap(~ setting, scales = "free_x") + scale_family_colour() + scale_family_fill() + pct_axis() +
                labs(x = NULL, y = "P(valid typology)", title = "Typology stage: share of iterations with >= 2 types of >= MIN_TYPE_SIZE sites")
        save_figure(p, "B4_typology_stage", 8, 4)
})

run_block("B5 plausibility stage: typology settings", {
        sets <- intersect(c(TYPOLOGY_SETTINGS, "alpha"), names(units))
        b5r  <- rbindlist(lapply(sets, function(v) {
                keep <- !is.na(units[[v]])
                d    <- units[keep, .(model_name, family, gamma_side, any_plaus, n_total)]
                d[, x := bin_x(units[[v]][keep])]
                scheme_weighted_rate(d, "any_plaus", c("family", "gamma_side", "x"))[, setting := v]
        }))
        save_table(b5r, "B5_binned_rates", show = FALSE)
        p <- ggplot(b5r[n >= MIN_CELL], aes(x, rate, colour = family, fill = family)) +
                geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() + geom_point(size = 0.8) +
                facet_grid(gamma_side ~ setting, scales = "free_x") + scale_family_colour() + scale_family_fill() + pct_axis() +
                labs(x = "Setting (exact value or quantile-bin median)", y = "P(>= 1 plausible environment in an open bin)",
                     title = "Plausibility stage by typology setting (scheme-weighted, clustered 95% CI)")
        save_figure(p, "B5_plausibility_by_setting", 11, 5.5)
        
        # GLM on iteration x bin units. gamma enters piecewise around 1, in
        # units of one bin width (0.25): g_down < 0 below 1, g_up > 0 above 1.
        units[, `:=`(g_down = pmin(gamma_mean - 1, 0) / 0.25, g_up = pmax(gamma_mean - 1, 0) / 0.25)]
        zv <- intersect(TYPOLOGY_SETTINGS, names(units))
        add_z(units, c("n_sites", zv), log_vars = "n_sites")
        rhs <- c("taxon", "family", "g_down", "g_up", "family:g_down", "family:g_up", paste0("z_", c("n_sites", zv)))
        save_table(fit_glm("any_plaus", rhs, units, binomial(), "units: >= 1 plausible in bin")$table, "B5b_glm_plausibility")
        save_table(vif_table(units, paste0("z_", c("n_sites", zv))), "B5c_vif_settings")
        if ("alpha" %in% names(units) && units[family == "spatial_ward", .N] > 0L) {
                uw <- units[family == "spatial_ward"]
                add_z(uw, c("n_sites", zv, "alpha"), log_vars = "n_sites")
                save_table(fit_glm("any_plaus", c("taxon", "g_down", "g_up", paste0("z_", c("n_sites", zv, "alpha"))),
                                   uw, binomial(), "units, spatial Ward: alpha")$table, "B5d_glm_spatial_ward_alpha")
        }
})

run_block("B6 manipulated variables", {
        # Mantel-Haenszel OR of >= 1 plausible environment in an open bin when a
        # variable is in S, stratified by scheme x gamma bin x |S|, so schemes,
        # bins and the number of shifted variables are held constant. The three
        # geology shares are always selected together and enter as one term.
        sv <- strsplit(iters$selected_vars, ";", fixed = TRUE)
        iter_vars <- data.table(model_name = rep(iters$model_name, lengths(sv)),
                                iteration  = rep(iters$iteration, lengths(sv)), variable = unlist(sv))
        iter_vars[variable %in% GEO_VARS, variable := "geology (3 shares)"]
        iter_vars <- unique(iter_vars)
        cand <- unique(iter_vars[, .(model_name, variable)])
        u <- units[, .(model_name, iteration, gamma_bin, n_vars, y = any_plaus, retained)]
        u[, stratum := .GRP, by = .(model_name, gamma_bin, n_vars)]
        b6 <- rbindlist(lapply(sort(unique(cand$variable)), function(v) {
                uv <- u[cand[variable == v, .(model_name)], on = "model_name", nomatch = NULL]
                uv[, x := FALSE]
                uv[iter_vars[variable == v, .(model_name, iteration)], x := TRUE, on = .(model_name, iteration)]
                data.table(variable = v, n_schemes = uniqueN(uv$model_name), share_in_S_proposed = mean(uv$x),
                           share_in_S_retained = if (any(uv$retained)) mean(uv$x[uv$retained]) else NA_real_,
                           mh_or(uv$stratum, uv$x, uv$y))
        }))
        b6[, `:=`(retained_over_proposed = share_in_S_retained / share_in_S_proposed, p_bh = p.adjust(p, "BH"))]
        setorder(b6, OR_MH)
        save_table(b6, "B6_variables_mantel_haenszel")
        pd <- b6[is.finite(OR_MH) & is.finite(lo) & is.finite(hi)]
        pd[, variable := factor(variable, levels = variable)]
        p <- ggplot(pd, aes(OR_MH, variable, colour = p_bh < 0.05)) + geom_vline(xintercept = 1, linetype = 2) +
                geom_pointrange(aes(xmin = lo, xmax = hi)) + scale_x_log10() +
                scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "#D55E00"), name = "BH-adjusted p < 0.05") +
                labs(x = "Mantel-Haenszel OR (>= 1 plausible environment | variable in S)", y = NULL,
                     title = "Manipulated variables and plausibility (strata: scheme x bin x |S|)")
        save_figure(p, "B6_variables_mantel_haenszel", 7, max(3, 0.22 * nrow(pd) + 1.5))
})

run_block("B7 GLM forest plot", {
        g <- rbindlist(GLM_TABLES)
        g <- g[term != "(Intercept)" & !grepl("^taxon", term) & is.finite(OR_lo) & is.finite(OR_hi)]
        save_table(rbindlist(GLM_TABLES), "B7_all_glm_coefficients", show = FALSE)
        p <- ggplot(g, aes(OR, term)) + geom_vline(xintercept = 1, linetype = 2) +
                geom_pointrange(aes(xmin = OR_lo, xmax = OR_hi)) + scale_x_log10() +
                facet_wrap(~ model, scales = "free_y") +
                labs(x = "Odds ratio (scheme-clustered 95% CI; z_ terms per SD, gamma per 0.25)", y = NULL,
                     title = "GLM effects (taxon terms omitted)")
        save_figure(p, "B7_glm_forest", 11, 6)
})


# ==============================================================================*
# C. GAMMA AND THE PLAUSIBILITY CHECKS
# ==============================================================================*
section("C. GAMMA AND PLAUSIBILITY")

run_block("C1 environment funnel", {
        agg <- envs[, .(generated = .N, pass_iso = sum(passed_isoforest), pass_maha = sum(passed_mahalanobis),
                        plausible = sum(plausible), retained = sum(retained)), by = .(model_name, taxon, family)]
        ft  <- with_totals(agg, c("taxon", "family"))[, lapply(.SD, sum), by = .(taxon, family, model_name),
                                                      .SDcols = c("generated", "pass_iso", "pass_maha", "plausible", "retained")]
        c1 <- ft[, .(n_schemes = .N, generated = sum(generated),
                     iso_rejected = 1 - sum(pass_iso) / sum(generated),
                     maha_rejected = 1 - sum(pass_maha) / sum(generated),
                     maha_rejected_after_iso = 1 - sum(plausible) / sum(pass_iso),
                     plausible_share = sum(plausible) / sum(generated),
                     retained_of_plausible = sum(retained) / sum(plausible),
                     plausible_share_scheme_mean = mean(plausible / generated),
                     iso_rejected_scheme_mean = mean(1 - pass_iso / generated)),
                 by = .(taxon, family)][order(taxon, family)]
        save_table(c1, "C1_environment_funnel")
        cat("Note: pooled shares weight schemes by their number of environments (i.e. by iterations run);\n",
            "      *_scheme_mean columns weight every scheme equally.\n")
})

run_block("C2 plausibility along gamma", {
        envs[, n_total := .N, by = model_name]
        cg <- envs[, .(n = .N, plausible = sum(plausible), pass_both = sum(plausible),
                       fail_iso_only = sum(!passed_isoforest & passed_mahalanobis),
                       fail_maha_only = sum(passed_isoforest & !passed_mahalanobis),
                       fail_both = sum(!passed_isoforest & !passed_mahalanobis), n_total = n_total[1]),
                   by = .(model_name, taxon, family, gamma_fine)]
        envs[, n_total := NULL]
        c2 <- with_totals(cg, c("taxon", "family"))
        c2 <- c2[, lapply(.SD, sum), by = .(taxon, family, gamma_fine, model_name),
                 .SDcols = c("n", "plausible", "pass_both", "fail_iso_only", "fail_maha_only", "fail_both")][
                         cg[, .(n_total = n_total[1]), by = model_name], on = "model_name"]
        classes <- c("pass_both", "fail_iso_only", "fail_maha_only", "fail_both")
        c2r <- rbindlist(lapply(c("plausible", classes), function(k) {
                scheme_weighted_rate(c2, k, c("taxon", "family", "gamma_fine"), n = "n")[, outcome := k]
        }))
        save_table(c2r[outcome == "plausible"], "C2_plausibility_by_gamma", show = FALSE)
        save_table(c2r[outcome != "plausible"], "C2b_check_classes_by_gamma", show = FALSE)
        
        pd <- c2r[outcome == "plausible" & taxon != "all" & family != "all" & n >= MIN_CELL]
        pd[, family := factor(family, levels = FAMILIES)]
        p <- ggplot(pd, aes(gamma_fine, rate, colour = family, fill = family)) + bin_edges() +
                geom_vline(xintercept = 1, linetype = 2) +
                geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() +
                facet_wrap(~ taxon, nrow = 1) + scale_family_colour() + scale_family_fill() + pct_axis() +
                labs(x = "gamma", y = "Plausible environments", title = "Plausibility along gamma (scheme-weighted, clustered 95% CI)")
        save_figure(p, "C2_plausibility_by_gamma", 10, 4)
        
        pd <- c2r[outcome %in% classes & family == "all" & n >= MIN_CELL]
        pd[, outcome := factor(outcome, levels = rev(classes),
                               labels = rev(c("pass both", "fail Isolation Forest only", "fail Mahalanobis only", "fail both")))]
        p <- ggplot(pd, aes(gamma_fine, rate, fill = outcome)) + geom_area(position = "stack") +
                geom_vline(xintercept = 1, linetype = 2) + facet_wrap(~ taxon, nrow = 1) +
                scale_fill_manual(values = c("pass both" = "#009E73", "fail Isolation Forest only" = "#E69F00",
                                             "fail Mahalanobis only" = "#56B4E9", "fail both" = "#D55E00"), name = NULL) +
                pct_axis() + labs(x = "gamma", y = "Share of environments (scheme-weighted)", title = "Which check rejects along gamma")
        save_figure(p, "C2b_check_classes_by_gamma", 10, 4)
        
        # Within-scheme trend across gamma bins (schemes = blocks).
        sb <- envs[, .(rate = mean(plausible)), by = .(model_name, gamma_bin)]
        sb[, bin := factor(BIN_LABELS[gamma_bin], levels = BIN_LABELS)]
        fr <- friedman_blocks(sb, "bin", "rate", pairs = "adjacent")
        if (!is.null(fr)) {
                save_table(fr$test, "C2c_gamma_friedman")
                save_table(fr$pairwise, "C2d_gamma_adjacent_bins")
        }
        rho <- sb[, if (.N >= 3L && uniqueN(rate) > 1L) .(rho = cor(gamma_bin, rate, method = "spearman")), by = model_name]
        if (nrow(rho)) {
                save_table(data.table(n_schemes = nrow(rho), median_rho = med(rho$rho), share_negative = mean(rho$rho < 0),
                                      p_wilcoxon = if (nrow(rho) >= 5L) wilcox.test(rho$rho, mu = 0, exact = FALSE)$p.value else NA_real_),
                           "C2e_within_scheme_spearman")
        }
})

run_block("C3 Isolation Forest vs Mahalanobis", {
        ag <- envs[, .(both_pass = sum(passed_isoforest & passed_mahalanobis),
                       fail_maha_only = sum(passed_isoforest & !passed_mahalanobis),
                       fail_iso_only = sum(!passed_isoforest & passed_mahalanobis),
                       both_fail = sum(!passed_isoforest & !passed_mahalanobis)), by = .(model_name, taxon)]
        c3 <- with_totals(ag, "taxon")[, {
                a <- as.numeric(sum(both_pass)); b <- as.numeric(sum(fail_maha_only))
                c <- as.numeric(sum(fail_iso_only)); d <- as.numeric(sum(both_fail)); n <- a + b + c + d
                po <- (a + d) / n; pe <- ((a + b) * (a + c) + (c + d) * (b + d)) / n^2
                m_str <- sum(fail_maha_only > fail_iso_only); i_str <- sum(fail_iso_only > fail_maha_only)
                .(n_schemes = .N, n_environments = n, share_both_pass = a / n, share_fail_iso_only = c / n,
                  share_fail_maha_only = b / n, share_both_fail = d / n, kappa = (po - pe) / (1 - pe),
                  schemes_maha_stricter = m_str, schemes_iso_stricter = i_str,
                  p_sign_test = if (m_str + i_str > 0L) binom.test(m_str, m_str + i_str)$p.value else NA_real_)
        }, by = taxon]
        save_table(c3, "C3_check_agreement")
        cat("Note: 'stricter' compares a scheme's single-check rejections; the sign test uses schemes as units.\n")
})

run_block("C4 clipping and achieved gamma", {
        r <- envs[retained == TRUE, .(model_name, taxon, family, bin, gamma_bin, gamma, gamma_achieved, share_clipped)]
        r[, outside_bin := gamma_achieved < GAMMA_BINS[gamma_bin] | gamma_achieved >= GAMMA_BINS[gamma_bin + 1L]]
        c4 <- with_totals(r, c("taxon", "family", "bin"))[, .(n_retained = .N, share_clipped_any = mean(share_clipped > 0),
                                                              median_share_clipped = med(share_clipped),
                                                              p95_share_clipped = quantile(share_clipped, 0.95, names = FALSE),
                                                              median_abs_gamma_deviation = med(abs(gamma_achieved - gamma), na.rm = TRUE),
                                                              share_outside_bin = mean(outside_bin, na.rm = TRUE)),
                                                          by = .(taxon, family, bin)]
        c4[, bin := factor(bin, levels = c(BIN_LABELS, "all"))]
        setorder(c4, taxon, family, bin)
        save_table(c4[taxon == "all" | family == "all"], "C4_clipping_retained")
        TABLES$C4_clipping_retained_full <- c4
        fwrite(c4, file.path(OUT_DIR, "tables", "C4_clipping_retained_full.csv"))
        c4b <- envs[, .(n = .N, plausible = mean(plausible)), by = .(clipped = share_clipped > 0, bin)][order(bin, clipped)]
        save_table(c4b, "C4b_plausibility_by_clipping")
        p <- ggplot(r, aes(gamma, gamma_achieved)) + geom_bin2d(bins = 60) +
                geom_abline(slope = 1, intercept = 0, linetype = 2) +
                geom_vline(xintercept = GAMMA_BINS, colour = "grey70") + geom_hline(yintercept = GAMMA_BINS, colour = "grey70") +
                facet_wrap(~ family) + scale_fill_gradient(low = "grey85", high = "#0072B2", trans = "log10", name = "Environments") +
                labs(x = "Nominal gamma", y = "Achieved gamma", title = "Retained environments: achieved vs nominal gamma (lines: bin edges)")
        save_figure(p, "C4_achieved_gamma", 9, 3.8)
})

run_block("C5 distance to the tolerance", {
        # Ratio of a manipulated environment's anomalous share to the tolerated
        # share; < 1 passes. Shows whether rejections are marginal or gross.
        qs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
        qf <- function(r) { q <- quantile(r, qs, names = FALSE, na.rm = TRUE); list(q10 = q[1], q25 = q[2], q50 = q[3], q75 = q[4], q90 = q[5]) }
        c5 <- rbind(envs[, qf(iso_n_exceed / (iso_tol_share * n_sites)), by = .(taxon, gamma_fine)][, check := "Isolation Forest"],
                    envs[, qf(maha_n_exceed / (maha_tol_share * n_sites)), by = .(taxon, gamma_fine)][, check := "Mahalanobis"])
        setorder(c5, check, taxon, gamma_fine)
        save_table(c5, "C5_ratio_to_tolerance", show = FALSE)
        p <- ggplot(c5, aes(gamma_fine)) + bin_edges() +
                geom_ribbon(aes(ymin = q10, ymax = q90), fill = "#0072B2", alpha = 0.15) +
                geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#0072B2", alpha = 0.3) +
                geom_line(aes(y = q50), colour = "#0072B2") + geom_hline(yintercept = 1, linetype = 2) +
                facet_grid(check ~ taxon) + scale_y_sqrt() +
                labs(x = "gamma", y = "Anomalous share / tolerated share (sqrt scale)",
                     title = "Distance to the tolerance (median, 50% and 80% ranges; < 1 passes)")
        save_figure(p, "C5_ratio_to_tolerance", 10, 5)
})


# ==============================================================================*
# D. SELECTION BY THE PLAUSIBILITY FILTER
# ==============================================================================*
section("D. SELECTION BY THE PLAUSIBILITY FILTER")

run_block("D1 retained vs proposed typologies", {
        # Proposed = every open bin of an iteration with a valid typology;
        # retained = the bins in which an environment was retained. Descriptors
        # are z-scored within scheme (over its valid iterations), so the SMD
        # compares typologies within schemes rather than schemes with each other.
        dv <- intersect(TYPOLOGY_DESCRIPTORS, names(iters))
        it <- iters[typology_valid == TRUE, c("model_name", "iteration", dv), with = FALSE]
        for (v in dv) set(it, j = v, value = as.numeric(it[[v]]))
        for (v in intersect(LOG_DESCRIPTORS, dv)) set(it, j = v, value = log1p(pmax(it[[v]], 0)))
        it[, (dv) := lapply(.SD, function(x) {
                s <- sd(x, na.rm = TRUE)
                if (is.na(s) || s == 0) rep(NA_real_, length(x)) else (x - mean(x, na.rm = TRUE)) / s
        }), by = model_name, .SDcols = dv]
        zu <- units[, .(model_name, iteration, family, bin, retained)][it, on = .(model_name, iteration), nomatch = NULL]
        lz <- melt(zu, id.vars = c("model_name", "iteration", "family", "bin", "retained"), measure.vars = dv,
                   variable.name = "descriptor", variable.factor = FALSE)[!is.na(value)]
        smd_j <- quote(.(n_proposed = .N, n_retained = sum(retained),
                         smd = if (sum(retained) >= 2L) (mean(value[retained]) - mean(value)) / sd(value) else NA_real_))
        d1 <- rbind(lz[, eval(smd_j), by = .(family, bin, descriptor)],
                    lz[, eval(smd_j), by = .(family, descriptor)][, bin := "all"])
        d1[, bin := factor(bin, levels = c(BIN_LABELS, "all"))]
        setorder(d1, family, descriptor, bin)
        save_table(d1[bin == "all"], "D1_smd_retained_vs_proposed")
        TABLES$D1_smd_by_bin <- d1
        fwrite(d1, file.path(OUT_DIR, "tables", "D1_smd_by_bin.csv"))
        L <- max(0.2, max(abs(d1$smd), na.rm = TRUE))
        p <- ggplot(d1, aes(bin, descriptor, fill = smd)) + geom_tile(colour = "white") +
                geom_text(aes(label = sprintf("%.2f", smd)), size = 2.3) + facet_wrap(~ family, nrow = 1) +
                scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, limits = c(-L, L), name = "SMD") +
                labs(x = "gamma bin", y = NULL,
                     title = "Retained vs proposed typologies: standardised mean difference (within-scheme z; |SMD| > 0.1 notable)") +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        save_figure(p, "D1_smd_retained_vs_proposed", 11, 5)
        
        key <- intersect(c("env_r2_orig", "asw_orig", "ld1_dominance", "type_evenness", "n_types"), dv)
        raw <- iters[typology_valid == TRUE, c("model_name", "iteration", "family", key), with = FALSE]
        rr  <- envs[retained == TRUE, .(model_name, iteration)][raw, on = .(model_name, iteration), nomatch = NULL]
        dd  <- rbind(melt(raw, id.vars = c("model_name", "iteration", "family"), measure.vars = key, variable.factor = FALSE)[, set := "proposed (valid typologies)"],
                     melt(rr,  id.vars = c("model_name", "iteration", "family"), measure.vars = key, variable.factor = FALSE)[, set := "retained environments"])
        dd[, value := as.numeric(value)]
        p <- ggplot(dd[!is.na(value)], aes(value, colour = set)) + geom_density() +
                facet_grid(family ~ variable, scales = "free") +
                scale_colour_manual(values = c("grey40", "#D55E00"), name = NULL) +
                labs(x = NULL, y = "Density", title = "Typology descriptors: proposed vs retained (raw scale)")
        save_figure(p, "D2_descriptor_densities", 11, 6)
})


# ==============================================================================*
# E. POST-HOC TOLERANCE SENSITIVITY
# ==============================================================================*
# The stored exceedance counts allow re-classifying every generated environment
# under another tolerance. Retention is re-derived as in script 09 (one
# environment per open bin and iteration, capped at the bin quota), but only
# within the iterations that were run: a re-run would stop earlier under a
# permissive rule and continue longer under a strict one, and generated no
# environments for bins already closed. The "as run" scenario reproduces the
# observed retention exactly (check printed below).
section("E. POST-HOC TOLERANCE SENSITIVITY")

run_block("E tolerance scenarios", {
        scen <- list("as run (anchored, both)"         = quote(iso_n_exceed < iso_tol_share * n_sites & maha_n_exceed < maha_tol_share * n_sites),
                     "anchored, Isolation Forest only" = quote(iso_n_exceed < iso_tol_share * n_sites),
                     "anchored, Mahalanobis only"      = quote(maha_n_exceed < maha_tol_share * n_sites))
        for (t in FIXED_TOL) {
                scen[[sprintf("fixed %.2f, both", t)]] <- bquote(iso_n_exceed < .(t) * n_sites & maha_n_exceed < .(t) * n_sites)
        }
        grid0 <- CJ(model_name = models[entered == TRUE, model_name], gamma_bin = seq_len(n_bins))
        grid0[models, taxon := i.taxon, on = "model_name"]
        grid0[, quota := BIN_QUOTA[gamma_bin]]
        observed <- envs[retained == TRUE, .(n_obs = .N), by = .(model_name, gamma_bin)]
        
        res <- lapply(names(scen), function(s) {
                j1 <- bquote(list(n = .N, k = sum(.(scen[[s]]))))
                j2 <- bquote(list(ok = any(.(scen[[s]]))))
                pl <- envs[, eval(j1), by = .(model_name, taxon, gamma_bin)]
                it <- envs[, eval(j2), by = .(model_name, iteration, gamma_bin)][, .(n_ok = sum(ok)), by = .(model_name, gamma_bin)]
                g  <- copy(grid0)[it, n_ok := i.n_ok, on = .(model_name, gamma_bin)]
                g[is.na(n_ok), n_ok := 0L]
                g[, n_ret := pmin(n_ok, quota)]
                if (s == names(scen)[1]) {
                        chk <- g[observed, on = .(model_name, gamma_bin)][, sum(n_ret != n_obs)]
                        cat(sprintf("  As-run check: %d scheme x bin cells differ from the observed retention (expected 0).\n", chk))
                }
                by_bin <- pl[, .(plausible_pooled = sum(k) / sum(n), plausible_scheme_mean = mean(k / n)), by = gamma_bin][
                        g[, .(share_quota_met = mean(n_ret >= quota), mean_fill = mean(n_ret / quota)), by = gamma_bin], on = "gamma_bin"]
                by_scheme <- g[, .(n_ret = sum(n_ret), complete = all(n_ret >= quota)), by = model_name]
                list(bin = by_bin[, `:=`(scenario = s, bin = factor(BIN_LABELS[gamma_bin], levels = BIN_LABELS))],
                     scheme = data.table(scenario = s, share_any_retained = mean(by_scheme$n_ret > 0),
                                         share_complete = mean(by_scheme$complete),
                                         mean_share_of_target = mean(by_scheme$n_ret) / N_TARGET_RETAINED))
        })
        e1 <- rbindlist(lapply(res, `[[`, "bin"))
        e2 <- rbindlist(lapply(res, `[[`, "scheme"))
        e1[, scenario := factor(scenario, levels = names(scen))]
        e2[, scenario := factor(scenario, levels = names(scen))]
        save_table(e1[, .(scenario, bin, plausible_pooled, plausible_scheme_mean, share_quota_met, mean_fill)], "E1_sensitivity_by_bin")
        save_table(e2, "E2_sensitivity_by_scheme")
        
        sc_col <- setNames(c("black", "grey45", "grey70", colorRampPalette(c("#BDD7EE", "#08306B"))(length(FIXED_TOL))), names(scen))
        p <- ggplot(e1, aes(bin, plausible_scheme_mean, colour = scenario, group = scenario)) + geom_line() + geom_point(size = 1) +
                scale_colour_manual(values = sc_col, name = NULL) + pct_axis() +
                labs(x = "gamma bin", y = "Plausible environments (scheme mean)", title = "Plausibility under alternative tolerance rules") +
                guides(colour = guide_legend(ncol = 3))
        save_figure(p, "E1_sensitivity_plausibility", 8, 5)
        p <- ggplot(e1, aes(bin, scenario, fill = share_quota_met)) + geom_tile(colour = "white") +
                geom_text(aes(label = sprintf("%.0f%%", 100 * share_quota_met)), size = 2.6) +
                scale_fill_gradient(low = "#F7F7F7", high = "#0072B2", limits = c(0, 1), name = "Quota met") +
                labs(x = "gamma bin", y = NULL, title = "Schemes filling each bin's quota within the iterations run")
        save_figure(p, "E2_sensitivity_quota", 8, 4)
})


# ==============================================================================*
# F. PRODUCTION OUTPUTS (optional, --prodDir)
# ==============================================================================*
if (!is.null(PROD_DIR)) {
        section("F. PRODUCTION OUTPUTS")
        
        run_block("F0 read production outputs", {
                # Output files are matched to schemes by model_name delimited by
                # non-alphanumeric characters (as in script 09); the longest match wins.
                mnames <- models$model_name
                esc    <- function(s) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", s)
                match_model <- function(f) {
                        b <- sub("\\.rds$", "", basename(f))
                        if (b %in% mnames) return(b)
                        hit <- mnames[vapply(mnames, function(m) grepl(paste0("(^|[^[:alnum:]])", esc(m), "([^[:alnum:]]|$)"), b), logical(1))]
                        if (length(hit)) hit[which.max(nchar(hit))] else NA_character_
                }
                pf   <- list.files(PROD_DIR, pattern = "\\.rds$", full.names = TRUE)
                prod <- lapply(pf, function(f) {
                        mn <- match_model(f)
                        if (is.na(mn)) return(NULL)
                        o <- readRDS(f)
                        if (!is.list(o) || is.null(o$environments)) return(NULL)
                        res <- list(env = data.table(model_name = mn, as.data.table(o$environments)),
                                    cp  = data.table(model_name = mn, as.data.table(o$community_predictors)))
                        rm(o); invisible(gc(FALSE))
                        res
                })
                prod_env <- rbindlist(lapply(prod, `[[`, "env"), fill = TRUE)
                prod_cp  <- rbindlist(lapply(prod, `[[`, "cp"),  fill = TRUE)
                rm(prod)
                prod_env[models, taxon := i.taxon, on = "model_name"]
                prod_cp[prod_env, `:=`(taxon = i.taxon, family = i.family, gamma_bin = i.gamma_bin,
                                       gamma_achieved = i.gamma_achieved), on = .(model_name, production_index)]
                cat(sprintf("Production files: %d | matched schemes: %d | environments: %d | realisations: %d\n",
                            length(pf), uniqueN(prod_env$model_name), nrow(prod_env), nrow(prod_cp)))
        })
        
        run_block("F1 manipulation check", {
                # Change of each descriptor (after - original typology) against the
                # achieved gamma; a working manipulation gives a monotone relation.
                ed <- ENV_DESCRIPTORS[paste0(ENV_DESCRIPTORS, "_after") %in% names(prod_env) &
                                              paste0(ENV_DESCRIPTORS, "_orig")  %in% names(prod_env)]
                for (v in ed) set(prod_env, j = paste0("delta_", v), value = prod_env[[paste0(v, "_after")]] - prod_env[[paste0(v, "_orig")]])
                rho <- rbindlist(lapply(ed, function(v) {
                        prod_env[, if (.N >= MIN_ENV_RHO) .(rho = suppressWarnings(cor(gamma_achieved, .SD[[1]], method = "spearman", use = "complete.obs"))),
                                 by = .(model_name, family), .SDcols = paste0("delta_", v)][, descriptor := v]
                }))
                f1 <- rho[!is.na(rho), .(n_schemes = .N, median_rho = med(rho), share_positive = mean(rho > 0),
                                         p_wilcoxon = if (.N >= 5L) wilcox.test(rho, mu = 0, exact = FALSE)$p.value else NA_real_),
                          by = .(descriptor, family)][order(descriptor, family)]
                save_table(f1, "F1_manipulation_check_rho")
                prod_env[, gamma_fine := (floor(gamma_achieved / GAMMA_STEP) + 0.5) * GAMMA_STEP]
                lm_ <- melt(prod_env[, c("model_name", "family", "gamma_fine", paste0("delta_", ed)), with = FALSE],
                            id.vars = c("model_name", "family", "gamma_fine"), variable.factor = FALSE)
                lm_ <- lm_[, .(value = mean(value, na.rm = TRUE)), by = .(model_name, family, gamma_fine, variable)][
                        , .(med = med(value, na.rm = TRUE), q25 = quantile(value, 0.25, na.rm = TRUE, names = FALSE),
                            q75 = quantile(value, 0.75, na.rm = TRUE, names = FALSE), n_schemes = .N), by = .(family, gamma_fine, variable)]
                lm_[, family := factor(family, levels = FAMILIES)]
                p <- ggplot(lm_[n_schemes >= 3L], aes(gamma_fine, med, colour = family, fill = family)) + bin_edges() +
                        geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 1, linetype = 2) +
                        geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, colour = NA) + geom_line() +
                        facet_wrap(~ variable, scales = "free_y") + scale_family_colour() + scale_family_fill() +
                        labs(x = "Achieved gamma", y = "Descriptor change (after - original), median over schemes and IQR",
                             title = "Manipulation check")
                save_figure(p, "F1_manipulation_check", 10, 6)
        })
        
        run_block("F2 variation partitioning of realisations", {
                prod_cp[, vp_class := sub("^error:.*", "error", vp_status)]
                f2 <- with_totals(prod_cp[, .(taxon, family, vp_class)], c("taxon", "family"))[, .N, by = .(taxon, family, vp_class)]
                f2[, share := N / sum(N), by = .(taxon, family)]
                setorder(f2, taxon, family, vp_class)
                save_table(f2, "F2_vp_status")
                f2b <- prod_cp[, .(n_realisations = .N, median_dbrda_r2 = med(dbrda_r2, na.rm = TRUE),
                                   share_mem_forced = mean(mem_forced, na.rm = TRUE), median_n_mem = med(as.numeric(n_mem_selected))),
                               by = .(family, gamma_bin)][order(family, gamma_bin)]
                save_table(f2b, "F2b_dbrda_by_gamma")
                p <- ggplot(f2[taxon != "all" & family != "all"], aes(family, share, fill = vp_class)) + geom_col(width = 0.7) +
                        facet_wrap(~ taxon, nrow = 1) + pct_axis() +
                        labs(x = NULL, y = "Share of realisations", fill = NULL, title = "Variation-partitioning status of realisations") +
                        theme(axis.text.x = element_text(angle = 30, hjust = 1))
                save_figure(p, "F2_vp_status", 9, 4)
        })
}


# ==============================================================================*
# KEY NUMBERS AND OUTPUT
# ==============================================================================*
section("KEY NUMBERS")

run_block("key numbers", {
        kn  <- character(0)
        add <- function(expr) tryCatch(kn <<- c(kn, expr), error = function(e) NULL)
        pc  <- function(x) sprintf("%.1f%%", 100 * x)
        n_m <- nrow(models); ent <- models[entered == TRUE]
        add(sprintf("Schemes: %d; skipped for model quality %d (%s), for < 4 candidate variables %d (%s); without diagnostics %d.",
                    n_m, sum(models$outcome == "skipped: model quality", na.rm = TRUE), pc(mean(models$outcome == "skipped: model quality", na.rm = TRUE)),
                    sum(models$outcome == "skipped: < 4 variables", na.rm = TRUE), pc(mean(models$outcome == "skipped: < 4 variables", na.rm = TRUE)),
                    sum(models$outcome == "no diagnostics", na.rm = TRUE)))
        add(sprintf("Of %d schemes that entered the simulation, %d (%s) yielded >= 1 retained environment and %d (%s) the full target of %d.",
                    nrow(ent), sum(ent$any_retained), pc(mean(ent$any_retained)), sum(ent$complete), pc(mean(ent$complete)), N_TARGET_RETAINED))
        add({
                b <- TABLES$B1b_two_stage_rates[taxon == "all"]
                sprintf("Iterations: %d; valid typology %s; >= 1 retained environment given a valid typology %s (%s).",
                        b[family == "all", n_iterations], pc(b[family == "all", p_valid_typology]), pc(b[family == "all", p_ok_given_valid]),
                        paste(b[family != "all", sprintf("%s %s", family, pc(p_ok_given_valid))], collapse = ", "))
        })
        add({
                f <- TABLES$C1_environment_funnel[taxon == "all" & family == "all"]
                sprintf("Manipulated environments: %d; rejected by Isolation Forest %s and by Mahalanobis %s (marginal, pooled); plausible %s pooled, %s scheme mean.",
                        f$generated, pc(f$iso_rejected), pc(f$maha_rejected), pc(f$plausible_share), pc(f$plausible_share_scheme_mean))
        })
        add({
                e <- TABLES$E1_sensitivity_by_bin[scenario == "as run (anchored, both)"]
                sprintf("Plausible environments by gamma bin (scheme mean): %s.", paste(sprintf("%s %s", e$bin, pc(e$plausible_scheme_mean)), collapse = ", "))
        })
        add({
                c4 <- TABLES$C4_clipping_retained[taxon == "all" & family == "all" & bin == "all"]
                sprintf("Retained environments with clipping: %s; with achieved gamma outside the nominal bin: %s.",
                        pc(c4$share_clipped_any), pc(c4$share_outside_bin))
        })
        writeLines(kn, file.path(OUT_DIR, "key_numbers.txt"))
        cat(kn, sep = "\n")
})

saveRDS(TABLES, file.path(OUT_DIR, "tables", "all_tables.rds"))
if (length(FIGURES)) {
        pdf(file.path(OUT_DIR, "all_figures.pdf"), width = 11, height = 7)
        for (nm in names(FIGURES)) tryCatch(print(FIGURES[[nm]] + labs(caption = nm)),
                                            error = function(e) cat("!! Page", nm, "failed:", conditionMessage(e), "\n"))
        invisible(dev.off())
}
cat(sprintf("\nReport written to %s (%d tables, %d figures).\n", OUT_DIR, length(TABLES), length(FIGURES)))
sink()
### - Reduced QRF benchmark models: tolerances, memory budget and figures - ###
#
# Runs locally (RStudio, R console or Rscript). Either edit SETTINGS below and source()
# the script, or pass command-line flags, which override SETTINGS:
#   Rscript 12b_summarise_reduction.R --input <dir> --output <dir>
#           [--pod-limit-gib 4] [--cache 1] [--baseline-mb 400] [--headroom 0.8]
#           [--image-cap-gib 5]
#
#   --input          folder with the <stem>_reduction_test.rds files of 12a_test_reduced_qrf.R
#   --pod-limit-gib  memory limit of the app container (Deployment: 12Gi)
#   --cache          models kept in memory by the app (QRF_CACHE)
#   --baseline-mb    memory of the idle app (R, shiny, packages); replace with the
#                    pod's measured memory before any model is loaded
#   --headroom       share of the pod limit the peak may use
#   --image-cap-gib  largest image Rahti's internal registry accepts
#
# Tolerances are predefined in TOL below and applied per input panel (heldout,
# training, offsupport); held-out performance applies to the held-out panel only.
# Memory per model, from the footprint measured in 12a (fresh R process):
#   peak = baseline + cache x resident model + transient copy while predicting
# and the largest tree count each model allows under the pod limit (k_max), next to
# the smallest tree count of the ladder that meets all agreement tolerances (k_min).
# Footprints measured on macOS can understate the Linux container (memory compression).
#
# Output (in the output folder):
#   reduction_gates.csv        every model x criterion x panel: value, tolerance, pass
#   reduction_summary.csv      one row per model: verdict, failed criteria, k_min, k_max, memory
#   reduction_ladder.csv       agreement statistics for every evaluated tree count (incl. the
#                              descriptive quantile shift relative to the 90% width, rel_*)
#   reduction_performance.csv  held-out coverage, MACE, CRPS per tree count vs. the full forest
#   reduction_footprint.csv    memory and timing of the reduced files
#   reduction_anchors.csv      reference cases (if any)
#   fig1_precision_ladder.png, fig2_agreement_bland_altman.png, fig3_compliance_scorecard.png

## --- Settings (edit for interactive use; command-line flags override) ----
SETTINGS <- list(
        input         = "SHINY APP/reduction_tests",    # <stem>_reduction_test.rds from 12a
        output        = "reduction_report",
        pod_limit_gib = 12,                   # memory limit of the app container (Deployment: 12Gi)
        cache         = 2,                    # models the app keeps in memory (QRF_CACHE)
        baseline_mb   = 400,                  # idle app memory; replace with the pod's measured value
        headroom      = 0.8,                  # share of the limit the peak may use
        image_cap_gib = 5                     # Rahti internal registry, per image
)

suppressPackageStartupMessages({
        library(data.table)
        library(ggplot2)
        library(dplyr)
})

## --- Arguments ----
args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, key) {
        i <- match(flag, args)
        if (is.na(i)) return(SETTINGS[[key]])
        if (i == length(args) || startsWith(args[i + 1], "--")) stop("Missing value for ", flag)
        args[i + 1]
}
in_dir      <- arg("--input", "input")
out_dir     <- arg("--output", "output")
pod_gib     <- as.numeric(arg("--pod-limit-gib", "pod_limit_gib"))
cache       <- as.integer(arg("--cache", "cache"))
baseline_mb <- as.numeric(arg("--baseline-mb", "baseline_mb"))
headroom    <- as.numeric(arg("--headroom", "headroom"))
image_gib   <- as.numeric(arg("--image-cap-gib", "image_cap_gib"))
if (is.null(in_dir) || is.null(out_dir)) stop("Missing input or output folder.")
if (!dir.exists(in_dir)) stop("Folder not found: ", in_dir, " (working directory: ", getwd(), ")")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%`   <- function(x, y) if (is.null(x)) y else x
GIB_MB    <- 2^30 / 1e6                              # MB per GiB
budget_mb <- headroom * pod_gib * GIB_MB - baseline_mb

## --- Predefined tolerances ----
# Probability-scale deviations in percentage points (pp). Expected at k = N/3 (1,000 of
# 3,000 trees) under pure tree subsampling: P95 ~2.1 pp, max ~4.3 pp for ~2,000 rows,
# ~3% category changes. The same as reseeding the full forest; below the models' own
# calibration error (median MACE: test 3.0 pp, grouped CV 2.1 pp).
TOL <- list(
        pct_p95     = 3,      # percentile of the observed value, P95 |reduced - full|
        pct_max     = 6,      #   maximum
        tau_p95     = 3,      # level of each benchmark quantile within the full distribution, P95
        tau_max     = 6,      #   maximum
        cat_changed = 0.05,   # share of rows changing app category
        cat_nonadj  = 0,      # rows changing by more than one category
        z_sd        = c(0.8, 1.25),   # SD of deviation / expected Monte Carlo SD
        cov90       = 1,      # held-out |coverage change| of the 90% interval, pp
        cov_max     = 1.5,    # held-out max |coverage change| over nominal levels, pp
        mace        = 1,      # held-out |MACE change|, pp
        crps_rel    = 0.01,   # held-out relative CRPS increase
        fit_agree   = 0.99,   # rebuilt held-out rows reproducing the stored median
        memory      = 1,      # peak / (headroom x pod limit)
        load_s      = NA      # seconds to read the file: reported, not judged (HPC disk is not Rahti disk)
)
STAT_PANELS <- c("heldout", "training", "offsupport")   # anchors are judged by category only

## --- Read test results ----
files <- list.files(in_dir, pattern = "_reduction_test\\.rds$", full.names = TRUE)
if (!length(files)) stop("No *_reduction_test.rds files in ", in_dir)
res <- lapply(files, readRDS)
names(res) <- vapply(res, `[[`, "", "model")

grab <- function(el) rbindlist(lapply(res, `[[`, el), fill = TRUE, use.names = TRUE)
integrity <- grab("integrity")
ladder    <- grab("ladder")
perf      <- grab("performance")
rows      <- grab("rows")
foot      <- grab("footprint")
fitchk    <- rbindlist(lapply(res, function(r) if (!is.null(r$fit_check))
        data.table(model = r$model, matched = r$fit_check$matched %||% NA_real_,
                   agree = r$fit_check$agree %||% NA_real_, rows_rebuilt = r$fit_check$rows_rebuilt,
                   rows_at_fit = r$fit_check$rows_at_fit)), fill = TRUE)
meta <- rbindlist(lapply(res, function(r) data.table(model = r$model, trees_full = r$trees_full, trees = r$trees)))
meta[, `:=`(taxon = sub("__.*$", "", model), metric = gsub("_", " ", sub("^.*?__", "", model)))]
k_sel <- if (uniqueN(meta$trees) == 1) meta$trees[1] else NA_integer_
cat("Models:", nrow(meta), "| trees in the reduced files:", paste(unique(meta$trees), collapse = ", "), "\n")

## --- Held-out performance relative to the full forest ----
perf_delta <- NULL
if (nrow(perf)) {
        pm <- perf[, .(cov90 = coverage[abs(nominal - 0.9) < 1e-9], mace = mean(abs(coverage - nominal)),
                       crps = crps[1]), by = .(model, trees)]
        full <- meta[, .(model, trees_full)]
        ref  <- pm[full, on = .(model, trees = trees_full)][, .(model, cov90_full = cov90, mace_full = mace, crps_full = crps)]
        cm   <- perf[full, on = .(model, trees = trees_full)][, .(model, nominal, cov_full = coverage)]
        covd <- perf[cm, on = .(model, nominal)][, .(cov_max = 100 * max(abs(coverage - cov_full))), by = .(model, trees)]
        perf_delta <- pm[ref, on = "model"][covd, on = .(model, trees)][
                , .(model, trees, cov90 = 100 * (cov90 - cov90_full), mace = 100 * (mace - mace_full),
                    crps_rel = crps / crps_full - 1, cov_max)]
        perf_delta <- perf_delta[meta[, .(model, trees_full)], on = "model"][trees < trees_full][, trees_full := NULL][]
}

## --- Gates for one model at one tree count ----
gate_rows <- function(m, k) {
        L <- ladder[model == m & trees == k & panel %in% STAT_PANELS]
        out <- list()
        add <- function(criterion, panel, value, pass, tol) {
                out[[length(out) + 1]] <<- data.table(model = m, trees = k, criterion, panel, value, tolerance = tol, pass)
        }
        for (i in seq_len(nrow(L))) {
                s <- L[i]
                add("pct_p95",     s$panel, s$pct_p95,     s$pct_p95     <= TOL$pct_p95,     TOL$pct_p95)
                add("pct_max",     s$panel, s$pct_max,     s$pct_max     <= TOL$pct_max,     TOL$pct_max)
                add("tau_p95",     s$panel, s$tau_p95,     s$tau_p95     <= TOL$tau_p95,     TOL$tau_p95)
                add("tau_max",     s$panel, s$tau_max,     s$tau_max     <= TOL$tau_max,     TOL$tau_max)
                add("cat_changed", s$panel, s$cat_changed, s$cat_changed <= TOL$cat_changed, TOL$cat_changed)
                add("cat_nonadj",  s$panel, s$cat_nonadj,  s$cat_nonadj  <= TOL$cat_nonadj,  TOL$cat_nonadj)
                add("z_sd",        s$panel, s$z_sd,        s$z_sd >= TOL$z_sd[1] && s$z_sd <= TOL$z_sd[2],
                    paste(TOL$z_sd, collapse = "-"))
        }
        if (!is.null(perf_delta)) {
                p <- perf_delta[model == m & trees == k]
                if (nrow(p)) {
                        add("cov90",    "heldout", p$cov90,    abs(p$cov90)    <= TOL$cov90,    TOL$cov90)
                        add("cov_max",  "heldout", p$cov_max,  p$cov_max       <= TOL$cov_max,  TOL$cov_max)
                        add("mace",     "heldout", p$mace,     abs(p$mace)     <= TOL$mace,     TOL$mace)
                        add("crps_rel", "heldout", p$crps_rel, p$crps_rel      <= TOL$crps_rel, TOL$crps_rel)
                }
        }
        rbindlist(out)
}

## --- Gates at the tree count of the reduced files ----
gates <- rbindlist(lapply(meta$model, function(m) gate_rows(m, meta[model == m, trees])))
extra <- list(
        integrity[, .(model, trees, criterion = "integrity", panel = "all", value = max_abs_diff,
                      tolerance = "identical", pass = identical)]
)
if (nrow(fitchk)) {
        extra[[length(extra) + 1]] <- fitchk[meta, on = "model", nomatch = 0][
                , .(model, trees, criterion = "fit_agree", panel = "heldout", value = agree,
                    tolerance = TOL$fit_agree, pass = !is.na(agree) & agree >= TOL$fit_agree)]
}
anchors <- rows[panel == "anchor"]
if (nrow(anchors)) {
        extra[[length(extra) + 1]] <- anchors[, .(trees = meta[model == .BY$model, trees], criterion = "anchor_category",
                                                  panel = "anchor", value = sum(cat_red != cat_full),
                                                  tolerance = 0, pass = all(cat_red == cat_full)), by = model]
}
if (nrow(foot)) {
        foot[, `:=`(peak_mb = baseline_mb + cache * load_mb + transient_mb)]
        foot[, `:=`(mem_ratio = (peak_mb - baseline_mb) / budget_mb,
                    k_max = floor(trees * budget_mb / (cache * load_mb + transient_mb)))]
        extra[[length(extra) + 1]] <- foot[, .(model, trees, criterion = "memory", panel = "all",
                                               value = peak_mb, tolerance = round(headroom * pod_gib * GIB_MB),
                                               pass = mem_ratio <= TOL$memory)]
}
gates <- rbind(gates, rbindlist(extra, use.names = TRUE, fill = TRUE), fill = TRUE)
gates[, tolerance := as.character(tolerance)]

## --- Smallest ladder tree count meeting all agreement tolerances (k_min) ----
kmin <- rbindlist(lapply(meta$model, function(m) {
        ks <- sort(unique(ladder[model == m, trees]))
        ok <- vapply(ks, function(k) { g <- gate_rows(m, k); nrow(g) > 0 && all(g$pass) }, NA)
        # smallest k from which every larger evaluated k also passes
        pass_from <- rev(cumprod(rev(ok)) == 1)
        data.table(model = m, k_min = if (any(pass_from)) ks[which(pass_from)[1]] else NA_integer_)
}))

summary_tab <- gates[, .(pass = all(pass, na.rm = TRUE), n_fail = sum(!pass, na.rm = TRUE),
                         failed = paste(unique(paste0(criterion, ifelse(panel == "all", "", paste0(" (", panel, ")")))[!pass]),
                                        collapse = "; ")), by = model]
summary_tab <- meta[summary_tab, on = "model"][kmin, on = "model"]
if (nrow(foot)) summary_tab <- foot[, .(model, file_mb = round(file_mb), resident_mb = round(load_mb),
                                        transient_mb = round(transient_mb), peak_mb = round(peak_mb),
                                        load_s = round(load_s, 1), k_max)][summary_tab, on = "model"]
setcolorder(summary_tab, c("taxon", "metric", "model", "trees", "pass", "n_fail", "k_min"))
setorder(summary_tab, taxon, metric)

## --- Write tables ----
fwrite(gates[order(model, criterion, panel)], file.path(out_dir, "reduction_gates.csv"))
fwrite(summary_tab, file.path(out_dir, "reduction_summary.csv"))
fwrite(ladder, file.path(out_dir, "reduction_ladder.csv"))
if (!is.null(perf_delta)) fwrite(perf_delta, file.path(out_dir, "reduction_performance.csv"))
if (nrow(foot)) fwrite(foot, file.path(out_dir, "reduction_footprint.csv"))
if (nrow(anchors)) fwrite(anchors[, .(model, label = data.set, observed = y, pct_full = 100 * pct_full,
                                      pct_red = 100 * pct_red, cat_full, cat_red)],
                          file.path(out_dir, "reduction_anchors.csv"))

## --- Console report ----
cat("\n  TOLERANCES at the tree count of the reduced files\n")
print(summary_tab[, .(model, trees, pass, failed = ifelse(nzchar(failed), failed, "-"), k_min,
                      k_max = if ("k_max" %in% names(summary_tab)) k_max else NA,
                      peak_mb = if ("peak_mb" %in% names(summary_tab)) peak_mb else NA)])
worst <- gates[panel != "anchor" & !criterion %in% c("integrity", "memory", "load_s", "fit_agree"),
               .(worst = if (criterion[1] %in% c("cov90", "mace")) max(abs(value)) else
                       if (criterion[1] == "z_sd") value[which.max(abs(value - 1))] else max(value)), by = criterion]
cat("\n  Worst value over models and panels:\n")
print(worst[, .(criterion, worst = signif(worst, 3))])
cat(sprintf("\n  Agreement: smallest passing tree count over all models (k_min) = %s\n",
            if (all(is.na(kmin$k_min))) "none in the ladder" else max(kmin$k_min, na.rm = TRUE)))
if (nrow(foot)) {
        big <- foot[which.max(peak_mb)]
        cat(sprintf("  Memory: pod limit %.1f GiB, cache %d, baseline %.0f MB, headroom %.0f%% -> budget %.0f MB for models\n",
                    pod_gib, cache, baseline_mb, 100 * headroom, budget_mb))
        cat(sprintf("          largest peak %s: %.0f MB = %.0f%% of the limit; largest tree count all models allow (k_max) = %d\n",
                    big$model, big$peak_mb, 100 * big$peak_mb / (pod_gib * GIB_MB), min(foot$k_max)))
        if ("method" %in% names(foot)) cat(sprintf("          footprint measured: %s\n", paste(unique(foot$method), collapse = "; ")))
        tot <- sum(foot$file_mb) / 1e3
        cat(sprintf("  Storage: online set %.1f GB (%s the %.0f GiB cap of Rahti's image registry, before the base image)\n",
                    tot, if (tot * 1e9 > image_gib * 2^30) "EXCEEDS" else "within", image_gib))
}
cat(sprintf("\n  VERDICT: %d of %d models meet all tolerances\n", sum(summary_tab$pass), nrow(summary_tab)))


#===============================================================================*
# Figures ----
#===============================================================================*
OI <- c(vermillion = "#D55E00", blue = "#0072B2", sky = "#56B4E9", orange = "#E69F00", grey = "grey45")
panel_labels <- c(heldout = "Held-out data sets", training = "Training data sets", offsupport = "Off-support inputs")
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95")))
save_png <- function(p, name, w, h) {
        ggsave(file.path(out_dir, name), p, width = w, height = h, dpi = 200)
        cat("  saved:", file.path(out_dir, name), "\n")
}
N_all <- unique(meta$trees_full)

## Figure 1: agreement vs. number of trees ----
# Lines: P95 |percentile deviation| per model. Dashed: expected under pure tree
# subsampling, from the full-forest percentiles. Band: reseeding the full forest
# (range over models). Red: tolerance. Vertical: tree count of the reduced files.
L1 <- ladder[panel %in% STAT_PANELS]
if (nrow(L1)) {
        rk <- rows[panel %in% STAT_PANELS][meta[, .(model, trees_full)], on = "model", nomatch = 0]
        set.seed(1); rk[, zabs := abs(rnorm(.N))]
        kgrid <- unique(round(exp(seq(log(min(L1$trees) * 0.8), log(max(N_all) * 0.98), length.out = 60))))
        expct <- rbindlist(lapply(kgrid, function(k) rk[k < trees_full, .(trees = k,
                                                                          p95 = quantile(zabs * 100 * sqrt(sv * (1 / k - 1 / trees_full)), 0.95, names = FALSE)), by = panel]))
        seedb <- L1[, .(lo = min(seed_p95), hi = max(seed_p95), med = median(seed_p95)), by = panel]
        for (d in list(L1, expct, seedb)) d[, panel := factor(panel, levels = names(panel_labels), labels = panel_labels)]
        p1 <- ggplot(dplyr::filter(L1, panel == "Training data sets"), aes(trees, pct_rms)) +
                # geom_rect(data = seedb, aes(xmin = min(L1$trees) * 0.8, xmax = max(N_all), ymin = lo, ymax = hi), inherit.aes = FALSE,
                          # fill = OI["sky"], alpha = 0.25) +
                # geom_hline(data = seedb, aes(yintercept = med), colour = OI["blue"], linewidth = 0.4) +
                # geom_hline(yintercept = TOL$pct_p95, colour = OI["vermillion"], linetype = 2) +
                { if (!is.na(k_sel)) geom_vline(xintercept = k_sel, colour = "grey30", linetype = 3) } +
                 geom_line(aes(group = model), colour = "grey55", linewidth = 0.35) +
                geom_point(size = 0.8, colour = "grey35") +
                #geom_line(data = expct, aes(trees, p95), linetype = 2, linewidth = 0.6) +
                facet_wrap(~ panel, nrow = 1) +
                scale_x_log10(breaks = sort(unique(c(L1$trees, N_all))), labels = scales::comma) +
                scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
                labs(x = sprintf("Trees kept (of %s)", paste(scales::comma(N_all), collapse = "/")),
                     y = "RMS percentile deviation from the full forest",
                     # title = "How many trees are enough?",
                     # subtitle = paste0("Grey: one line per model. Dashed: expected from tree subsampling.\n",
                                       # "Blue: reseeding the full forest (band: range over models). Orange dashed: tolerance. Dotted: reduced files.")
                                       ) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        save_png(p1, "fig1_precision_ladder.png", w = 9, h = 3.8)
}

L1$taxon <- stringr::str_split_i(L1$model, "__", i = 1)
L1$method <- stringr::str_split_i(L1$model, "__", i = 2)

L1 %>% 
        filter(panel == "Training data sets") %>%
        # group_by(taxon) %>%
        group_by(method) %>%
        summarize(mean = mean(pct_rms))

## Figure 2: row-level agreement at the tree count of the reduced files ----
# Bland-Altman: deviation vs. full-forest percentile, held-out rows (training rows if no
# held-out panel). Curves: +/-1.96 Monte Carlo SD expected from tree subsampling.
# Vertical lines: app category boundaries. Orange: rows that change category.
pan2 <- if (nrow(rows[panel == "heldout"])) "heldout" else "training"
R2 <- rows[panel == pan2][meta[, .(model, taxon, trees, trees_full)], on = "model", nomatch = 0]
if (nrow(R2)) {
        R2[, changed := cat_red != cat_full]
        env <- CJ(taxon = unique(R2$taxon), p = seq(0, 1, by = 0.005))[
                unique(R2[, .(taxon, trees, trees_full)]), on = "taxon", allow.cartesian = TRUE][
                        , sd := 100 * sqrt(p * (1 - p) * (1 / trees - 1 / trees_full))]
        env <- env[, .(sd = max(sd)), by = .(taxon, p)]
        lab <- R2[, .(txt = sprintf("category changed: %.1f%%\noutside curves: %.1f%%", 100 * mean(changed),
                                    100 * mean(abs(d_pct) > 1.96 * sig, na.rm = TRUE))), by = taxon]
        ylim <- max(abs(R2$d_pct)) * 1.05
        p2 <- ggplot(R2, aes(100 * pct_full, d_pct)) +
                geom_vline(xintercept = c(10, 25, 75, 90), colour = "grey70", linetype = 2, linewidth = 0.3) +
                geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
                geom_point(data = R2[changed == FALSE], colour = "grey30", alpha = 0.3, size = 0.7, stroke = 0) +
                geom_point(data = R2[changed == TRUE], colour = OI["orange"], size = 0.9) +
                geom_line(data = env, aes(100 * p, 1.96 * sd), colour = OI["blue"], linewidth = 0.5) +
                geom_line(data = env, aes(100 * p, -1.96 * sd), colour = OI["blue"], linewidth = 0.5) +
                geom_text(data = lab, aes(x = 50, y = -ylim, label = txt), size = 2.7, vjust = 0, lineheight = 0.9) +
                facet_wrap(~ taxon, nrow = 1) +
                scale_x_continuous(breaks = c(0, 10, 25, 50, 75, 90, 100)) +
                coord_cartesian(ylim = c(-ylim, ylim)) +
                labs(x = "Percentile of the observed value, full forest", y = "Reduced - full (pp)",
                     title = sprintf("Agreement per typology (%s rows, %s trees)", sub("heldout", "held-out", pan2),
                                     if (is.na(k_sel)) "reduced" else scales::comma(k_sel)),
                     subtitle = "Blue: \u00b11.96 Monte Carlo SD expected from tree subsampling (no ties). Dashed: app category boundaries. Orange: category changed.")
        save_png(p2, "fig2_agreement_bland_altman.png", w = 2.6 * uniqueN(R2$taxon) + 1, h = 3.8)
}

## Figure 3: compliance scorecard ----
# Worst value over panels per model and criterion, relative to its tolerance.
CRIT <- c(integrity = "Identical\ntrees", pct_p95 = "Percentile\nP95\n(pp)", pct_max = "Percentile\nmax\n(pp)",
          tau_p95 = "Quantile\nlevel P95\n(pp)", tau_max = "Quantile\nlevel max\n(pp)",
          cat_changed = "Category\nchanged", cat_nonadj = "Change\n>1 class", z_sd = "SD\nof z",
          cov90 = "\u0394 cover.\n90%\n(pp)", cov_max = "\u0394 cover.\nmax\n(pp)", mace = "\u0394\nMACE\n(pp)",
          crps_rel = "\u0394\nCRPS", anchor_category = "Anchor\ncate-\ngories", memory = "Peak\nmemory\n(MB)",
          load_s = "Load\n(s)")
G <- gates[criterion %in% names(CRIT) & criterion != "fit_agree"]
num_tol <- function(c) switch(c, z_sd = NA_real_, integrity = NA_real_, anchor_category = 0, cat_nonadj = 0,
                              memory = headroom * pod_gib * GIB_MB, as.numeric(TOL[[c]]))
G[, ratio := {
        c <- criterion[1]
        if (c == "z_sd") ifelse(value < 1, (1 - value) / (1 - TOL$z_sd[1]), (value - 1) / (TOL$z_sd[2] - 1))
        else if (c %in% c("integrity", "anchor_category", "cat_nonadj")) ifelse(pass, 0, 2)
        else if (c == "memory") (value - baseline_mb) / budget_mb
        else abs(value) / num_tol(c)
}, by = criterion]
G3 <- G[, .SD[which.max(ratio)], by = .(model, criterion)][meta, on = "model"]
G3[, cls := factor(fifelse(!pass, "fails", fifelse(ratio <= 0.5, "\u2264 50% of tolerance", "50\u2013100%")),
                   levels = c("\u2264 50% of tolerance", "50\u2013100%", "fails"))]
fmt <- function(c, v) {
        switch(c, integrity = ifelse(v == 0, "yes", "no"), cat_changed = sprintf("%.1f%%", 100 * v),
               crps_rel = sprintf("%+.2f%%", 100 * v),
               cat_nonadj = sprintf("%d", as.integer(v)), anchor_category = sprintf("%d", as.integer(v)),
               memory = sprintf("%.0f", v), z_sd = sprintf("%.2f", v),
               cov90 = sprintf("%+.2f", v), mace = sprintf("%+.2f", v), sprintf("%.2g", v))
}
G3[, label := mapply(fmt, criterion, value)]
G3[, `:=`(criterion = factor(criterion, levels = names(CRIT), labels = CRIT),
          row = factor(paste0(taxon, " | ", metric), levels = rev(unique(paste0(meta[order(taxon, metric)]$taxon, " | ",
                                                                                meta[order(taxon, metric)]$metric)))))]
p3 <- ggplot(G3, aes(criterion, row, fill = cls)) +
        geom_tile(colour = "white", linewidth = 0.6) +
        geom_text(aes(label = label, colour = cls == "fails"), size = 2.5) +
        scale_fill_manual(values = c("#E3ECF5", OI[["sky"]], OI[["vermillion"]]), drop = FALSE, name = NULL) +
        scale_colour_manual(values = c(`FALSE` = "grey10", `TRUE` = "white"), guide = "none") +
        scale_x_discrete(position = "top") +
        labs(x = NULL, y = NULL,
             title = sprintf("Compliance of the reduced models (%s trees)", if (is.na(k_sel)) "various" else scales::comma(k_sel)),
             subtitle = sprintf("Worst value over input panels. Memory: peak with cache %d at a %.0f GiB limit (%.0f%% headroom).",
                                cache, pod_gib, 100 * headroom)) +
        theme(panel.grid = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
              axis.text.x.top = element_text(size = 7.5, lineheight = 0.9))
save_png(p3, "fig3_compliance_scorecard.png", w = 12, h = 2 + 0.28 * nrow(meta))
cat("Done\n")

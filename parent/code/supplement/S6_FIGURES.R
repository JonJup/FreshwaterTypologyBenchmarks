################################################################################
# plot_simulation_outcomes.R
#
# Creates three figures from the diagnostics written by finalize() in
# 09_simulate_data.R, pooled over the projects pulseD, pulseF, pulseM, pulseI:
#   A1_outcome_by_taxon.png          simulation outcome per scheme
#   B1_iteration_status.png          iteration outcome by classification family
#   B3_family_share_of_retained.png  family composition of the retained environments
################################################################################

library(data.table)
library(ggplot2)

# ==============================================================================*
# 0. SETTINGS
# ==============================================================================*

BASE_DIR      <- ""
TAXON_FOLDERS <- c("pulseD", "pulseF", "pulseM", "pulseI")
DIAG_SUBDIR   <- "data/misc/simulation_diagnostics"      # contains models/ and iterations/
MODEL_SUBDIR  <- "data/001_unfitted_hmsc_models"         # NULL: do not report "no diagnostics"
OUT_DIR       <- "figures/supplement/S6/"

N_TARGET_RETAINED <- 100L
FAMILIES          <- c("kmeans", "spatial_ward", "crossed_thresholds")


# C2b
GAMMA_STEP <- 0.05   # resolution of the gamma curves
MIN_CELL   <- 30L    # minimum environments per plotted gamma step
ALL_TAXA   <- "All taxa"

## Display labels ----------------------------------------------------------------
# Taxon labels: values of the taxon column -> display name. Taxa not listed are
# shown in title case with underscores replaced by spaces.
TAXON_LABELS <- c(diatoms            = "Diatoms",
                  fish               = "Fish",
                  macrophytes        = "Macrophytes",
                  invertebrates      = "Invertebrates",
                  macroinvertebrates = "Macroinvertebrates")

OUTCOME_LABELS <- c(complete     = sprintf("Target reached (%d environments)", N_TARGET_RETAINED),
                    partial      = sprintf("Target partly reached (1\u2013%d)", N_TARGET_RETAINED - 1L),
                    none         = "No environment retained",
                    skip_vars    = "Skipped: < 4 candidate variables",
                    skip_quality = "Skipped: insufficient model quality",
                    no_diag      = "No diagnostics")

STATUS_LABELS <- c(FAIL_too_few_types   = "Fewer than two types",
                   FAIL_small_type      = "Type below minimum size",
                   FAIL_isoforest_all   = "All rejected by Isolation Forest",
                   FAIL_mahalanobis_all = "All rejected by Mahalanobis distance",
                   OK                   = "\u2265 1 environment retained")

CHECK_LABELS <- c(pass_both      = "Passed both checks",
                  fail_iso_only  = "Rejected by Isolation Forest only",
                  fail_maha_only = "Rejected by Mahalanobis distance only",
                  fail_both      = "Rejected by both checks")

FAMILY_LABELS <- c(kmeans             = "k-means",
                   spatial_ward       = "Spatially constrained Ward",
                   crossed_thresholds = "Crossed thresholds")

## Colours (Okabe-Ito) --------------------------------------------------------------
OUTCOME_COLOURS <- c(complete = "#009E73", partial = "#56B4E9", none = "#D55E00",
                     skip_vars = "#E69F00", skip_quality = "#999999", no_diag = "#000000")
STATUS_COLOURS  <- c(FAIL_too_few_types = "#CC79A7", FAIL_small_type = "#E69F00",
                     FAIL_isoforest_all = "#D55E00", FAIL_mahalanobis_all = "#F0E442", OK = "#009E73")
FAMILY_COLOURS  <- c(kmeans = "#0072B2", spatial_ward = "#E69F00", crossed_thresholds = "#009E73")
CHECK_COLOURS   <- c(pass_both = "#009E73", fail_iso_only = "#E69F00", fail_maha_only = "#56B4E9", fail_both = "#D55E00")


# ==============================================================================*
# 1. HELPERS
# ==============================================================================*

read_rds_dir <- function(dir, cols) {
        f <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
        if (!length(f)) return(NULL)
        rbindlist(lapply(f, function(x) {
                d <- as.data.table(readRDS(x))
                d[, intersect(cols, names(d)), with = FALSE]
        }), use.names = TRUE, fill = TRUE)
}

nice_taxon <- function(x) {
        lab <- unname(TAXON_LABELS[x])
        ifelse(is.na(lab), tools::toTitleCase(gsub("_", " ", x)), lab)
}

save_png <- function(p, name, width, height) {
        f <- file.path(OUT_DIR, paste0(name, ".png"))
        if (capabilities("cairo")) {
                ggsave(f, p, width = width, height = height, dpi = 300, type = "cairo")
        } else {
                ggsave(f, p, width = width, height = height, dpi = 300)
        }
        cat("Wrote", f, "\n")
}

pct_axis <- function() scale_y_continuous(labels = function(x) paste0(100 * x, "%"))

theme_set(theme_bw(base_size = 9) +
                  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                        strip.background = element_rect(fill = "grey95", colour = NA)))


# ==============================================================================*
# 2. DATA
# ==============================================================================*

MODEL_COLS <- c("model_name", "taxon", "run_status", "n_env_retained", paste0("n_env_retained_", FAMILIES))
ITER_COLS  <- c("model_name", "family", "status")
ENV_COLS   <- c("model_name", "gamma", "passed_isoforest", "passed_mahalanobis")

models <- list()
iters  <- list()
envs   <- list()
for (pr in TAXON_FOLDERS) {
        dd <- file.path(BASE_DIR, pr, DIAG_SUBDIR)
        m  <- read_rds_dir(file.path(dd, "models"), MODEL_COLS)
        if (is.null(m)) {
                warning("No model diagnostics in ", file.path(dd, "models"), "; project skipped.")
                next
        }
        m[, project := pr]
        
        # Schemes with an unfitted model but no diagnostics
        if (!is.null(MODEL_SUBDIR)) {
                md <- file.path(BASE_DIR, pr, MODEL_SUBDIR)
                if (dir.exists(md)) {
                        expected <- sub("\\.rds$", "", list.files(md, pattern = "\\.rds$"))
                        missing  <- setdiff(expected, m$model_name)
                        if (length(missing)) {
                                m <- rbind(m, data.table(model_name = missing, taxon = sub("_[^_]*$", "", missing),
                                                         run_status = "NO_DIAGNOSTICS", project = pr), fill = TRUE)
                        }
                }
        }
        models[[pr]] <- m
        
        it <- read_rds_dir(file.path(dd, "iterations"), ITER_COLS)
        if (!is.null(it)) iters[[pr]] <- it[, project := pr]
        en <- read_rds_dir(file.path(dd, "environments"), ENV_COLS)
        if (!is.null(en)) envs[[pr]] <- en[, project := pr]
        cat(sprintf("%s: %d schemes, %d iterations, %d environments\n", pr, nrow(m),
                    if (is.null(it)) 0L else nrow(it), if (is.null(en)) 0L else nrow(en)))
}
if (!length(models)) stop("No diagnostics found in any project.")
models <- rbindlist(models, use.names = TRUE, fill = TRUE)
iters  <- rbindlist(iters,  use.names = TRUE, fill = TRUE)
envs   <- rbindlist(envs,   use.names = TRUE, fill = TRUE)

# Scheme outcome
models[, outcome := fcase(
        run_status == "OK" & n_env_retained >= N_TARGET_RETAINED, "complete",
        run_status == "OK",                                       "partial",
        run_status == "NO_RETAINED_ENV",                          "none",
        run_status == "SKIPPED_too_few_variables",                "skip_vars",
        run_status == "SKIPPED_model_quality",                    "skip_quality",
        run_status == "NO_DIAGNOSTICS",                           "no_diag")]
if (models[, any(is.na(outcome))]) {
        warning(models[, sum(is.na(outcome))], " schemes with unknown run_status are omitted from A1.")
}

# Taxon display labels; iterations get the taxon of their scheme
models[, taxon_lab := nice_taxon(taxon)]
models[, taxon_lab := factor(taxon_lab, levels = sort(unique(taxon_lab)))]
iters[models, taxon_lab := i.taxon_lab, on = .(project, model_name)]
envs[models,  taxon_lab := i.taxon_lab, on = .(project, model_name)]

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================*
# 3. FIGURES
# ==============================================================================*
# 
## B1 Iteration outcome by classification family -----------------------------------
b1 <- iters[!is.na(taxon_lab), .N, by = .(taxon_lab, family, status)]
b1[, share := N / sum(N), by = .(taxon_lab, family)]
b1[, `:=`(family = factor(family, levels = FAMILIES),
          status = factor(status, levels = names(STATUS_LABELS)))]

p <- ggplot(b1, aes(family, share, fill = status)) +
        geom_col(width = 0.7, show.legend = TRUE) +
        facet_wrap(~ taxon_lab, nrow = 1) +
        scale_fill_manual(values = STATUS_COLOURS, labels = STATUS_LABELS, drop = FALSE, name = NULL) +
        scale_x_discrete(labels = FAMILY_LABELS) +
        pct_axis() +
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
        labs(x = NULL, y = "Share of iterations") +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_png(p, "rejection_reason", 9, 4.5)

## B3 Family share of retained environments -----------------------------------------
fc <- intersect(paste0("n_env_retained_", FAMILIES), names(models))
sh <- melt(models[n_env_retained > 0, c("project", "model_name", "taxon_lab", "n_env_retained", fc), with = FALSE],
           measure.vars = fc, variable.name = "family", value.name = "n", variable.factor = FALSE)
sh[, `:=`(family = factor(sub("n_env_retained_", "", family), levels = FAMILIES),
          share  = n / n_env_retained)]

p <- ggplot(sh[!is.na(share)], aes(family, share, fill = family)) +
        geom_boxplot(outlier.size = 0.5, show.legend = TRUE) +
        geom_hline(yintercept = 1 / length(FAMILIES), linetype = 2) +
        facet_wrap(~ taxon_lab, nrow = 1) +
        scale_fill_manual(values = FAMILY_COLOURS, labels = FAMILY_LABELS, drop = FALSE, name = NULL) +
        pct_axis() +
        labs(x = NULL, y = "Share of the scheme's retained environments") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_png(p, "rejection_by_family", 9, 4)

## C2b Which check rejects along gamma ------------------------------------------------
# Share of environments in each check class per gamma step. Each environment is
# weighted 1 / (environments of its scheme), so every scheme has the same total
# weight; schemes that ran many iterations do not dominate.
envs[, `:=`(n_total    = .N,
            gamma_fine = (floor(gamma / GAMMA_STEP) + 0.5) * GAMMA_STEP), by = .(project, model_name)]
cg <- envs[!is.na(taxon_lab), .(n              = .N,
                                pass_both      = sum(passed_isoforest & passed_mahalanobis),
                                fail_iso_only  = sum(!passed_isoforest & passed_mahalanobis),
                                fail_maha_only = sum(passed_isoforest & !passed_mahalanobis),
                                fail_both      = sum(!passed_isoforest & !passed_mahalanobis),
                                n_total        = n_total[1]),
           by = .(project, model_name, taxon_lab, gamma_fine)]
cg[, taxon_lab := as.character(taxon_lab)]
cg <- rbind(copy(cg)[, taxon_lab := ALL_TAXA], cg)

c2b <- melt(cg, id.vars = c("project", "model_name", "taxon_lab", "gamma_fine", "n", "n_total"),
            measure.vars = names(CHECK_LABELS), variable.name = "outcome", value.name = "k", variable.factor = FALSE)
c2b <- c2b[, .(n = sum(n), rate = sum(k / n_total) / sum(n / n_total)), by = .(taxon_lab, gamma_fine, outcome)][n >= MIN_CELL]
c2b[, `:=`(taxon_lab = factor(taxon_lab, levels = c(ALL_TAXA, levels(models$taxon_lab))),
           outcome   = factor(outcome, levels = rev(names(CHECK_LABELS))))]   # rejections stacked on top

p <- ggplot(c2b, aes(gamma_fine, rate, fill = outcome)) +
        geom_area(position = "stack") +
        geom_vline(xintercept = 1, linetype = 2) +
        facet_wrap(~ taxon_lab, nrow = 1) +
        scale_fill_manual(values = CHECK_COLOURS, labels = CHECK_LABELS, breaks = names(CHECK_LABELS), name = NULL) +
        pct_axis() +
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
        labs(x = expression(gamma), y = "Share of environments (scheme-weighted)")
save_png(p, "rejection_by_gamma", 10, 4)

### - Reduced vs. full QRF benchmark models: prediction agreement - ###
#
# Runs locally (RStudio, R console or Rscript) or on a cluster. Either edit SETTINGS below
# and source() the script, or pass command-line flags, which override SETTINGS:
#   Rscript 12a_test_reduced_qrf.R --full <dir> --reduced <dir> --output <dir>
#           [--model <taxon__metric>] [--ladder 250,500,750,1000,1500,2000]
#           [--heldout --code-dir ../pulse/code | --heldout-rows <dir>] [--anchors <csv>]
#           [--n-heldout 2000] [--n-train 1000] [--n-offsupport 500] [--chunk-trees 0]
#           [--threads n] [--no-footprint] [--footprint-rows 240] [--overwrite] [--seed 1]
#   Rscript 12a_test_reduced_qrf.R --export-heldout --full <dir> --output <dir> --code-dir ../pulse/code
#
#   --full, --reduced  folders with the full (11_fit_qrf.R) and reduced (12_reduce_qrf.R)
#                      model files; --full also holds <stem>_run_info.rds and
#                      <stem>_intervals.rds (held-out panel)
#   --model            one model (file stem). Default: all models found in both folders.
#                      Models with a result in --output are skipped unless --overwrite,
#                      so an interrupted run resumes where it stopped.
#   --heldout          rebuild the test-set rows with prepare_qrf_inputs.R (in --code-dir);
#                      needs the pipeline data: run from the folder where 11_fit_qrf.R runs.
#                      The rows used are saved as <output>/<stem>_heldout_rows.rds.
#   --heldout-rows     folder with <stem>_heldout_rows.rds (from --export-heldout or an
#                      earlier --heldout run): the held-out panel without the pipeline data
#   --export-heldout   only build and save the held-out rows; loads no forest
#   --chunk-trees n    predict the full forest in blocks of n trees, cut as in 12_reduce_qrf.R.
#                      Peak memory falls from the forest + ~2 x forest (ranger's working copy)
#                      to the forest + ~2 x one block. 0 = whole forest at once. With blocks,
#                      the integrity check compares the reduced file with the first k trees
#                      cut from the full file; the check against the uncut forest is the one
#                      in 12_reduce_qrf.R.
#   --threads          ranger threads (default: physical cores)
#   --anchors          CSV of reference cases: taxon, metric, label, observed,
#                      optional expected_category (1-5), and one column per predictor
#   (--index i, or SLURM_ARRAY_TASK_ID, still selects one model on a cluster.)
#
# Pipeline data on the cluster, forests on the local machine:
#   cluster: Rscript 12a_test_reduced_qrf.R --export-heldout --full <full> --output heldout_rows
#   local:   copy heldout_rows/, then SETTINGS$heldout_rows <- "heldout_rows" and source()
#
# Memory: the largest full forest needs its own size plus ~2 x its forest while predicting
# (about 3 x in total). On machines with less than ~16 GB RAM use chunk_trees = 500.
#
# Design. The full forest is evaluated once (per-tree predictions, n x N trees). The
# reduced file is evaluated separately and must reproduce the columns of its kept trees
# exactly (integrity). A forest of the first k trees therefore equals the first k
# columns, so every tree count of the ladder is evaluated without further files or fits.
# Reference for "how much is too much": the deviation between two independent full
# forests (reseeding), estimated from the two disjoint halves of the full forest:
# (half A - half B) / sqrt(2) has the variance of the difference between two N-tree fits.
#
# Input panels (rows per model):
#   heldout     test-set typologies of the designed split (unseen data sets), observed values
#   training    training typologies, spread evenly over data sets, observed values
#   offsupport  each predictor drawn independently from its training values (breaks the
#               joint structure: unusual input combinations); value drawn from training values
#   anchor      user-supplied reference cases (e.g. the case study), if anchors are given
#
# Statistics per panel and tree count k (probability scale in percentage points):
#   percentile of the observed value (mid-rank, as app.R): |reduced - full|
#   level of each benchmark quantile (taus) of the reduced forest within the full forest:
#     distance of tau from the probability mass the full forest places at or below that
#     value (zero if the value is a valid tau-quantile of the full forest; per-tree values
#     are training values, so the predictive distributions are discrete)
#   quantile shift relative to the full forest's 90% interval width (descriptive)
#   app categories (<10, <25, 25-75, >75, >90): share changed, changes by > 1 class
#   z = deviation / Monte Carlo SD expected from tree subsampling, sqrt(v (1/k - 1/N)),
#     v = variance of one tree's mid-rank score (p(1-p) without ties)
#   held-out panel: coverage (all nominal levels), MACE, CRPS
# Footprint (reduced file, fresh R process; packages callr and ps): memory after loading,
# peak while predicting footprint_rows rows (sampled every 10 ms; exact peak where the OS
# reports one: Linux, Windows), load and prediction time. macOS compresses memory, so its
# figures can understate the Linux container on Rahti.
#
# Output: <output>/<stem>_reduction_test.rds

## --- Settings (edit for interactive use; command-line flags override) ----
SETTINGS <- list(
        full           = "SHINY APP/models/",       # full models + <stem>_run_info.rds, <stem>_intervals.rds
        reduced        = "SHINY APP/small_models/",     # reduced models (12_reduce_qrf.R)
        output         = "reduction_tests",   # results
        model          = NULL,                # one file stem, e.g. "fish__ANOSIM_R"; NULL = all
        ladder         = c(250, 500, 750, 1000, 1500, 2000),
        heldout        = FALSE,               # rebuild held-out rows from the pipeline inputs
        heldout_rows   = NULL,                # or: folder with exported <stem>_heldout_rows.rds
        export_heldout = FALSE,               # only export held-out rows, then stop
        code_dir       = "/code",     # folder with functions/prepare_qrf_inputs.R
        anchors        = NULL,                # CSV of reference cases
        n_heldout      = 2000,
        n_train        = 1000,
        n_offsupport   = 500,
        chunk_trees    = 0,                   # e.g. 500 with less than ~16 GB RAM
        threads        = NA,                  # NA = physical cores
        footprint      = TRUE,
        footprint_rows = 240,                 # app: 200 reference draws + one-at-a-time rows
        overwrite      = FALSE,
        seed           = 1
)

suppressPackageStartupMessages({
        library(data.table)
        library(ranger)
})

## --- Arguments ----
args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, key) {
        i <- match(flag, args)
        if (is.na(i)) return(SETTINGS[[key]])
        if (i == length(args) || startsWith(args[i + 1], "--")) stop("Missing value for ", flag)
        args[i + 1]
}
has <- function(flag, key) flag %in% args || isTRUE(SETTINGS[[key]])
num <- function(x) if (is.null(x)) NA else as.numeric(x)

full_dir    <- arg("--full", "full")
red_dir     <- arg("--reduced", "reduced")
out_dir     <- arg("--output", "output")
model_arg   <- arg("--model", "model")
index_arg   <- if (!is.na(i <- match("--index", args))) args[i + 1] else Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)
ladder      <- arg("--ladder", "ladder")
ladder      <- as.integer(if (is.character(ladder)) strsplit(ladder, ",")[[1]] else ladder)
do_heldout  <- has("--heldout", "heldout")
heldout_dir <- arg("--heldout-rows", "heldout_rows")
export_ho   <- has("--export-heldout", "export_heldout")
code_dir    <- arg("--code-dir", "code_dir")
anchor_file <- arg("--anchors", "anchors")
n_heldout   <- as.integer(arg("--n-heldout", "n_heldout"))
n_train     <- as.integer(arg("--n-train", "n_train"))
n_off       <- as.integer(arg("--n-offsupport", "n_offsupport"))
chunk_trees <- as.integer(arg("--chunk-trees", "chunk_trees"))
do_foot     <- !("--no-footprint" %in% args) && isTRUE(SETTINGS$footprint)
foot_rows   <- as.integer(arg("--footprint-rows", "footprint_rows"))
overwrite   <- has("--overwrite", "overwrite")
seed        <- as.integer(arg("--seed", "seed"))

n_threads <- as.integer(num(arg("--threads", "threads")))
if (is.na(n_threads)) n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_threads)) n_threads <- max(parallel::detectCores(logical = FALSE), 1L, na.rm = TRUE)

if (is.null(full_dir) || is.null(out_dir)) stop("Missing full or output folder.")
if (!export_ho && is.null(red_dir)) stop("Missing reduced folder.")
if (!dir.exists(full_dir)) stop("Folder not found: ", full_dir, " (working directory: ", getwd(), ")")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TAUS       <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)   # quantiles offered by app.R
NOMINAL    <- c(seq(0.1, 0.9, by = 0.1), 0.95)               # as 11_fit_qrf.R
EVAL_PROBS <- round(seq(0.005, 0.995, by = 0.005), 6)        # CRPS grid, as 11_fit_qrf.R

log_step   <- function(...) { cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n"); flush.console() }
clean_name <- function(x) gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", x))


## --- Models ----
AUX_SUFFIXES <- c("ale", "datasets", "intervals", "lofo", "metrics",
                  "quantile_importance", "run_info", "split", "tuning")
discover <- function(dir) {
        s <- sub("\\.rds$", "", list.files(dir, pattern = "__.+\\.rds$"))
        s[!grepl(paste0("_(", paste(AUX_SUFFIXES, collapse = "|"), ")$"), s)]
}
select_stems <- function(stems) {
        if (!is.null(model_arg)) return(intersect(stems, model_arg))
        if (!is.na(index_arg)) {
                i <- as.integer(index_arg)
                if (i < 1 || i > length(stems)) stop("--index ", i, " outside 1..", length(stems))
                return(stems[i])
        }
        stems
}
if (export_ho) {
        stems <- select_stems(sub("_run_info\\.rds$", "", list.files(full_dir, pattern = "__.+_run_info\\.rds$")))
        if (!length(stems)) stop("No <stem>_run_info.rds files in ", full_dir)
} else {
        if (!dir.exists(red_dir)) stop("Folder not found: ", red_dir, " (working directory: ", getwd(), ")")
        stems <- select_stems(intersect(discover(red_dir), discover(full_dir)))
        if (!length(stems)) stop("No requested model present in both ", full_dir, " and ", red_dir)
}
cat(sprintf("\n%s %d model(s) | threads %d | chunk_trees %d | held-out: %s\n\n",
            if (export_ho) "Exporting held-out rows for" else "Testing", length(stems), n_threads, chunk_trees,
            if (!is.null(heldout_dir)) heldout_dir else if (do_heldout) "rebuild" else "none"))


## --- Helpers ----
tree_values <- function(model, X) {                 # n x trees, as app.R
        p <- predict(model, data = X, type = "quantiles", what = function(v) v,
                     num.threads = n_threads)$predictions
        p <- as.matrix(p)
        if (nrow(p) != nrow(X)) p <- t(p)
        unname(p)
}
# Forest cut to some of its trees, as 12_reduce_qrf.R (used for block-wise prediction)
PER_TREE <- c("child.nodeIDs", "split.varIDs", "split.values")
subset_forest <- function(r, keep) {
        n0 <- r$num.trees
        per_tree <- names(r$forest)[vapply(r$forest, function(e) is.list(e) && length(e) == n0, NA)]
        if (!setequal(per_tree, PER_TREE)) stop("Unexpected per-tree forest elements: ", paste(per_tree, collapse = ", "))
        for (nm in PER_TREE) r$forest[[nm]] <- r$forest[[nm]][keep]
        n_new <- length(keep)
        storage.mode(n_new) <- typeof(n0)
        r$forest$num.trees <- n_new
        r$num.trees        <- n_new
        r$random.node.values <- r$random.node.values[seq_len(max(lengths(r$forest$split.varIDs))), keep, drop = FALSE]
        r$predictions <- NULL
        r
}
tree_values_full <- function(r, X, chunk) {         # n x N, optionally in blocks of trees
        N <- as.integer(r$num.trees)
        if (is.na(chunk) || chunk <= 0 || chunk >= N) return(tree_values(r, X))
        out <- matrix(NA_real_, nrow(X), N)
        for (b in split(seq_len(N), ceiling(seq_len(N) / chunk))) {
                out[, b] <- tree_values(subset_forest(r, b), X)
                invisible(gc(verbose = FALSE))
        }
        out
}
mid_pct  <- function(tv, y) (rowSums(tv < y) + 0.5 * rowSums(tv == y)) / ncol(tv)   # app.R mid_pct()
sort_rows <- function(tv) {
        s <- t(apply(tv, 1, sort.int, method = "quick"))
        if (ncol(tv) == 1) t(s) else s
}
q7 <- function(S, probs) {                          # type-7 quantiles of sorted rows (= ranger)
        m <- ncol(S)
        vapply(probs, function(p) {
                h <- (m - 1) * p + 1; lo <- floor(h); hi <- min(m, lo + 1)
                S[, lo] + (h - lo) * (S[, hi] - S[, lo])
        }, numeric(nrow(S)))
}
level_dev <- function(tv, q, tau) {                # generalised inverse, robust to ties
        lo <- rowMeans(tv < q); hi <- rowMeans(tv <= q)
        ifelse(tau < lo, lo - tau, ifelse(tau > hi, hi - tau, 0))
}
cat5 <- function(p) {                               # app.R wording, percent scale
        pc <- 100 * p
        ifelse(pc < 10, 1L, ifelse(pc < 25, 2L, ifelse(pc > 90, 5L, ifelse(pc > 75, 4L, 3L))))
}
spread_sample <- function(groups, n, seed) {        # row indices, evenly over groups
        set.seed(seed)
        dt  <- data.table(g = groups, i = seq_along(groups))
        per <- ceiling(n / uniqueN(dt$g))
        idx <- dt[, .(i = i[sample.int(.N, min(.N, per))]), by = g]$i
        sort(if (length(idx) > n) sample(idx, n) else idx)
}
pinball_crps <- function(y, Q, probs) {             # 2 x mean pinball over a uniform grid
        u <- y - Q
        2 * mean(sweep(u, 2, probs, function(a, p) a * (p - (a < 0))))
}
interval_cov <- function(y, S) {
        vapply(NOMINAL, function(l) {
                b <- q7(S, c((1 - l) / 2, 1 - (1 - l) / 2))
                mean(y >= b[, 1] & y <= b[, 2])
        }, 0)
}

## Held-out rows: rebuilt with the pipeline's input functions, restricted to the test
## data sets and to the (capped) subsets recorded in _intervals.rds; a sample of
## n_heldout rows, spread over data sets, is kept with all model columns.
.evl2 <- NULL
build_heldout <- function(stem) {
        ri <- readRDS(file.path(full_dir, paste0(stem, "_run_info.rds")))
        ip <- readRDS(file.path(full_dir, paste0(stem, "_intervals.rds")))
        if (is.null(.evl2)) {
                if (!exists("load_qrf_inputs")) source(file.path(code_dir, "functions", "prepare_qrf_inputs.R"))
                .evl2 <<- load_qrf_inputs(log = log_step)$evl2
        }
        d <- .evl2[metric == ri$metric]
        if (!identical(ri$taxon, "all")) d <- d[taxon_key == tolower(ri$taxon)]
        d <- d[data.set %in% ri$test_datasets]
        d <- na.omit(d[, unique(qrf_model_cols), with = FALSE], cols = qrf_model_cols)
        tp <- as.data.table(ip$test_pred)
        d  <- d[paste(data.set, scheme_id) %chin% unique(paste(tp$data.set, tp$scheme_id))]
        list(rows = d[spread_sample(d$data.set, n_heldout, seed)], rows_rebuilt = nrow(d),
             rows_at_fit = nrow(tp), n_heldout = n_heldout, seed = seed, created = Sys.time())
}
heldout_file <- function(dir, stem) file.path(dir, paste0(stem, "_heldout_rows.rds"))
get_heldout <- function(stem) {
        if (!is.null(heldout_dir)) {
                f <- heldout_file(heldout_dir, stem)
                if (!file.exists(f)) stop("no ", basename(f), " in ", heldout_dir)
                return(readRDS(f))
        }
        ho <- build_heldout(stem)
        saveRDS(ho, heldout_file(out_dir, stem))                       # reusable via heldout_rows
        ho
}

## Footprint of a model file in a fresh R process. The child reports its memory after
## loading and, where the OS provides it, its exact peak; the parent samples the
## child's resident memory every 10 ms while it predicts (all platforms).
footprint <- function(path, n_rows, threads, poll = 0.01) {
        if (!requireNamespace("callr", quietly = TRUE) || !requireNamespace("ps", quietly = TRUE)) {
                log_step("  footprint skipped: install.packages(c(\"callr\", \"ps\"))")
                return(NULL)
        }
        status <- tempfile()
        child <- function(path, n_rows, threads, status) {
                stage <- function(s) writeLines(s, status)
                rss   <- function() ps::ps_memory_info()[["rss"]] / 1e6
                peak  <- function() {                                         # MB; NA on macOS
                        if (file.exists("/proc/self/status")) {
                                s <- readLines("/proc/self/status")
                                return(as.numeric(gsub("\\D", "", s[startsWith(s, "VmHWM:")])) * 1024 / 1e6)
                        }
                        m <- ps::ps_memory_info()
                        if ("peak_wset" %in% names(m)) m[["peak_wset"]] / 1e6 else NA_real_
                }
                suppressPackageStartupMessages(library(ranger))
                invisible(gc()); base <- rss(); stage("loading")
                t0 <- Sys.time(); o <- readRDS(path); t_load <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
                m <- o$model; p <- m$forest$independent.variable.names
                X <- as.data.frame(o$recipe$template)[, p, drop = FALSE]
                set.seed(1); X <- X[sample.int(nrow(X), min(n_rows, nrow(X)), replace = nrow(X) < n_rows), , drop = FALSE]
                rm(o); invisible(gc()); loaded <- rss(); stage("loaded"); Sys.sleep(0.1)
                stage("predicting"); t0 <- Sys.time()
                v <- predict(m, data = X, type = "quantiles", what = function(v) v, num.threads = threads)$predictions
                t_pred <- as.numeric(difftime(Sys.time(), t0, units = "secs")); stage("done")
                list(rss_base_mb = base, rss_load_mb = loaded, peak_exact_mb = peak(), load_s = t_load,
                     predict_s = t_pred, rows = nrow(X), threads = threads, trees = m$num.trees,
                     object_mb = as.numeric(object.size(m)) / 1e6, forest_mb = as.numeric(object.size(m$forest)) / 1e6)
        }
        proc <- callr::r_bg(child, args = list(path = normalizePath(path), n_rows = n_rows, threads = threads,
                                               status = status), stdout = NULL, stderr = NULL)
        h <- ps::ps_handle(proc$get_pid())
        poll_peak <- NA_real_
        while (proc$is_alive()) {
                st <- if (file.exists(status)) tryCatch(suppressWarnings(readLines(status, warn = FALSE)),
                                                        error = function(e) character(0)) else character(0)
                if (length(st) && identical(st[1], "predicting")) {
                        m <- tryCatch(ps::ps_memory_info(h)[["rss"]] / 1e6, error = function(e) NA_real_)
                        poll_peak <- max(poll_peak, m, na.rm = TRUE)
                }
                Sys.sleep(poll)
        }
        f <- tryCatch(proc$get_result(), error = function(e) {
                log_step("  footprint failed:", conditionMessage(e)); NULL })
        if (is.null(f)) return(NULL)
        f <- as.data.table(f)
        f[, `:=`(peak_poll_mb = poll_peak,
                 peak_proc_mb = suppressWarnings(max(peak_exact_mb, poll_peak, na.rm = TRUE)),
                 method = if (is.na(peak_exact_mb)) "sampled every 10 ms" else "exact peak + sampling")]
        f[!is.finite(peak_proc_mb), peak_proc_mb := NA_real_]
        f[, `:=`(load_mb = rss_load_mb - rss_base_mb,                  # model resident after loading
                 transient_mb = pmax(0, peak_proc_mb - rss_load_mb))]   # ranger's working copy while predicting
        f
}

## Export mode: held-out rows only ----
if (export_ho) {
        for (stem in stems) {
                f <- heldout_file(out_dir, stem)
                if (file.exists(f) && !overwrite) { log_step(stem, ": exists, skipped"); next }
                ho <- build_heldout(stem)
                saveRDS(ho, f)
                log_step(sprintf("%s: %d of %d rebuilt rows (%d at fit time) -> %s", stem, nrow(ho$rows),
                                 ho$rows_rebuilt, ho$rows_at_fit, f))
        }
        log_step("Held-out rows exported; copy them and set heldout_rows to that folder.")
        stems <- character(0)                                          # nothing else to do
}

anchors <- if (!is.null(anchor_file)) fread(anchor_file) else NULL


## --- Per model ----
for (stem in stems) {
        log_step("Model:", stem)
        t_model  <- Sys.time()
        res_file <- file.path(out_dir, paste0(stem, "_reduction_test.rds"))
        if (file.exists(res_file) && !overwrite) {
                log_step("  result exists, skipped (overwrite = TRUE to redo)")
                next
        }
        
        ## Full model: panels and per-tree predictions ----
        obj_full <- readRDS(file.path(full_dir, paste0(stem, ".rds")))
        r_full   <- obj_full$model
        preds    <- r_full$forest$independent.variable.names
        N        <- as.integer(r_full$num.trees)
        tmpl     <- as.data.table(obj_full$recipe$template)
        tmpl     <- tmpl[stats::complete.cases(tmpl[, c(preds, "value", "data.set"), with = FALSE])]
        rm(obj_full)
        
        panels <- list()
        fit_check <- NULL
        if (do_heldout || !is.null(heldout_dir)) {
                ho <- tryCatch(get_heldout(stem), error = function(e) {
                        log_step("  held-out rows unavailable:", conditionMessage(e)); NULL })
                if (!is.null(ho)) {
                        d    <- as.data.table(ho$rows)
                        miss <- setdiff(c(preds, "value", "data.set"), names(d))
                        if (length(miss)) {
                                log_step("  held-out rows lack:", paste(miss, collapse = ", "))
                        } else {
                                panels$heldout <- data.table(panel = "heldout", data.set = d$data.set,
                                                             scheme_id = if ("scheme_id" %in% names(d)) as.character(d$scheme_id) else NA_character_,
                                                             y = d$value, d[, preds, with = FALSE])
                                f_int <- file.path(full_dir, paste0(stem, "_intervals.rds"))
                                fit_check <- list(rows_rebuilt = ho$rows_rebuilt, rows_at_fit = ho$rows_at_fit,
                                                  tp = if (file.exists(f_int)) as.data.table(readRDS(f_int)$test_pred) else NULL)
                                log_step(sprintf("  held-out: %d of %d rebuilt rows (%d at fit time), %d data sets",
                                                 nrow(d), ho$rows_rebuilt, ho$rows_at_fit, uniqueN(d$data.set)))
                        }
                }
        }
        idx <- spread_sample(tmpl$data.set, n_train, seed + 1)
        panels$training <- data.table(panel = "training", data.set = tmpl$data.set[idx], scheme_id = NA_character_,
                                      y = tmpl$value[idx], tmpl[idx, preds, with = FALSE])
        set.seed(seed + 2)
        off <- as.data.table(lapply(setNames(preds, preds), function(p) sample(tmpl[[p]], n_off, replace = TRUE)))
        panels$offsupport <- data.table(panel = "offsupport", data.set = NA_character_, scheme_id = NA_character_,
                                        y = sample(tmpl$value, n_off, replace = TRUE), off)
        if (!is.null(anchors)) {
                a <- anchors[paste0(clean_name(taxon), "__", clean_name(metric)) == stem]
                if (nrow(a)) {
                        miss <- setdiff(preds, names(a))
                        if (length(miss)) {
                                log_step("  anchors skipped, missing predictors:", paste(miss, collapse = ", "))
                        } else {
                                panels$anchor <- data.table(panel = "anchor", data.set = a$label, scheme_id = NA_character_,
                                                            y = a$observed, a[, preds, with = FALSE])
                        }
                }
        }
        P <- rbindlist(panels, use.names = TRUE)
        P[, row := .I]
        X <- as.data.frame(P[, preds, with = FALSE])
        y <- P$y
        
        log_step(sprintf("  full forest: %d trees, %d rows (%s)%s", N, nrow(P),
                         paste(names(table(P$panel)), table(P$panel), sep = " = ", collapse = ", "),
                         if (chunk_trees > 0 && chunk_trees < N) sprintf(", blocks of %d trees", chunk_trees) else ""))
        tv_full <- tree_values_full(r_full, X, chunk_trees)
        rm(r_full); invisible(gc(verbose = FALSE))
        
        ## Reduced file: integrity ----
        red_path <- file.path(red_dir, paste0(stem, ".rds"))
        obj_red  <- readRDS(red_path)
        k_file   <- as.integer(obj_red$model$num.trees)
        red_info <- obj_red$reduction
        if (!is.null(red_info$n_trees_full) && red_info$n_trees_full != N) {
                stop(stem, ": reduced file was cut from a forest with ", red_info$n_trees_full, " trees; full file has ", N)
        }
        tv_red <- tree_values(obj_red$model, X)
        rm(obj_red); invisible(gc(verbose = FALSE))
        ref    <- tv_full[, seq_len(k_file), drop = FALSE]
        integrity <- data.table(model = stem, trees_full = N, trees = k_file, rows = nrow(P),
                                identical = identical(tv_red, ref),
                                max_abs_diff = max(abs(tv_red - ref)))
        rm(ref)
        log_step(sprintf("  reduced file: %d trees | identical to the first %d trees of the full forest: %s",
                         k_file, k_file, integrity$identical))
        
        ## Full-forest reference quantities ----
        S_full   <- sort_rows(tv_full)
        p_lt     <- rowMeans(tv_full < y)
        p_eq     <- rowMeans(tv_full == y)
        pct_full <- p_lt + 0.5 * p_eq                                   # = mid_pct(tv_full, y)
        sv_full  <- pmax(0, p_lt + 0.25 * p_eq - pct_full^2)            # variance of one tree's score
        Q_full   <- q7(S_full, TAUS)
        width90  <- Q_full[, TAUS == 0.95] - Q_full[, TAUS == 0.05]
        cat_full <- cat5(pct_full)
        half     <- N %/% 2
        seed_dev <- 100 * (mid_pct(tv_full[, seq_len(half)], y) - mid_pct(tv_full[, half + seq_len(half)], y)) / sqrt(2)
        
        ## Fit-time check of the rebuilt held-out rows: full-forest median vs. stored median ----
        if (!is.null(fit_check$tp) && "q_median" %in% names(fit_check$tp)) {
                tp  <- fit_check$tp[, .(data.set, scheme_id = as.character(scheme_id), y = value, q_stored = q_median)]
                tp  <- tp[, if (.N == 1) .SD, by = .(data.set, scheme_id, y)]                   # unique keys only
                h   <- P[panel == "heldout", .(row, data.set, scheme_id, y)]
                mm  <- tp[h, on = .(data.set, scheme_id, y), nomatch = 0]
                q50 <- Q_full[mm$row, TAUS == 0.5]
                fit_check$matched <- nrow(mm) / nrow(h)
                fit_check$agree   <- if (nrow(mm)) mean(abs(q50 - mm$q_stored) <= 1e-9 * (1 + abs(mm$q_stored))) else NA_real_
                log_step(sprintf("  held-out rows vs. fit time: %.0f%% matched, %.1f%% of those with identical median",
                                 100 * fit_check$matched, 100 * fit_check$agree))
        }
        if (!is.null(fit_check)) fit_check$tp <- NULL
        
        ## Tree-count ladder ----
        ks <- sort(unique(c(ladder[ladder < N], k_file)))
        ladder_rows <- list(); tau_rows <- list(); perf_rows <- list(); rows_file <- NULL
        ho_rows <- P[panel == "heldout", row]
        
        perf_of <- function(S, k) {
                if (!length(ho_rows)) return(NULL)
                Sh  <- S[ho_rows, , drop = FALSE]; yh <- y[ho_rows]
                cov <- interval_cov(yh, Sh)
                data.table(trees = k, nominal = NOMINAL, coverage = cov,
                           crps = pinball_crps(yh, q7(Sh, EVAL_PROBS), EVAL_PROBS))
        }
        perf_rows[[length(perf_rows) + 1]] <- perf_of(S_full, N)
        
        for (k in ks) {
                tk  <- if (k == k_file) tv_red else tv_full[, seq_len(k), drop = FALSE]
                Sk  <- sort_rows(tk)
                pk  <- mid_pct(tk, y)
                Qk  <- q7(Sk, TAUS)
                d_pct <- 100 * (pk - pct_full)
                sig   <- 100 * sqrt(sv_full * (1 / k - 1 / N))
                z     <- ifelse(sig > 0, d_pct / sig, NA_real_)
                d_tau <- vapply(seq_along(TAUS), function(j) 100 * level_dev(tv_full, Qk[, j], TAUS[j]), numeric(nrow(P)))
                rel   <- abs(Qk - Q_full) / ifelse(width90 > 0, width90, NA_real_)
                ck    <- cat5(pk)
                
                R <- data.table(panel = P$panel, d_pct, z, seed_dev, cat_full, cat_k = ck,
                                tau_absmax = apply(abs(d_tau), 1, max), rel_max = apply(rel, 1, max, na.rm = TRUE))
                ladder_rows[[length(ladder_rows) + 1]] <- R[, .(
                        trees      = k,
                        rows       = .N,
                        pct_rms    = sqrt(mean(d_pct^2)),
                        pct_p95    = quantile(abs(d_pct), 0.95, names = FALSE),
                        pct_max    = max(abs(d_pct)),
                        z_sd       = sd(z, na.rm = TRUE),
                        tau_p95    = quantile(abs(d_tau[.I, ]), 0.95, names = FALSE),
                        tau_max    = max(tau_absmax),
                        rel_median = median(rel_max, na.rm = TRUE),
                        rel_p95    = quantile(rel_max, 0.95, names = FALSE, na.rm = TRUE),
                        cat_changed = mean(cat_k != cat_full),
                        cat_nonadj  = sum(abs(cat_k - cat_full) > 1L),
                        seed_rms   = sqrt(mean(seed_dev^2)),
                        seed_p95   = quantile(abs(seed_dev), 0.95, names = FALSE)
                ), by = panel]
                tau_rows[[length(tau_rows) + 1]] <- rbindlist(lapply(seq_along(TAUS), function(j)
                        data.table(panel = P$panel, tau = TAUS[j], d = d_tau[, j], rel = rel[, j])))[, .(
                                trees = k, tau_p95 = quantile(abs(d), 0.95, names = FALSE), tau_max = max(abs(d)),
                                rel_p95 = quantile(rel, 0.95, names = FALSE, na.rm = TRUE)), by = .(panel, tau)]
                perf_rows[[length(perf_rows) + 1]] <- perf_of(Sk, k)
                
                if (k == k_file) {
                        rows_file <- data.table(P[, .(row, panel, data.set, y)], pct_full, sv = sv_full, pct_red = pk, d_pct, sig, z,
                                                cat_full, cat_red = ck, tau_absmax = R$tau_absmax, rel_max = R$rel_max,
                                                q50_full = Q_full[, TAUS == 0.5], q50_red = Qk[, TAUS == 0.5])
                }
                rm(tk, Sk, Qk, d_tau, rel, R); invisible(gc(verbose = FALSE))
        }
        ladder_tab <- rbindlist(ladder_rows)[, model := stem][]
        tau_tab    <- rbindlist(tau_rows)[, model := stem][]
        perf_tab   <- rbindlist(perf_rows)
        if (nrow(perf_tab)) perf_tab[, model := stem]
        
        ## Footprint of the reduced file ----
        foot <- if (do_foot) footprint(red_path, foot_rows, min(2L, n_threads)) else NULL
        if (!is.null(foot)) {
                foot[, `:=`(model = stem, file_mb = file.size(red_path) / 1e6)]
                log_step(sprintf("  footprint: %.0f MB resident after loading (%.1f s), +%.0f MB transient while predicting %d rows (%.2f s; %s)",
                                 foot$load_mb, foot$load_s, foot$transient_mb, foot$rows, foot$predict_s, foot$method))
        }
        
        at_k <- ladder_tab[trees == k_file]
        for (pn in at_k$panel) {
                s <- at_k[panel == pn]
                log_step(sprintf("  %-10s k = %d: percentile |dev| P95 %.2f pp, max %.2f pp (reseeding P95 %.2f) | categories changed %.1f%% | z SD %.2f",
                                 pn, k_file, s$pct_p95, s$pct_max, s$seed_p95, 100 * s$cat_changed, s$z_sd))
        }
        
        res <- list(model = stem, trees_full = N, trees = k_file, reduction = red_info,
                    settings = list(ladder = ks, taus = TAUS, nominal = NOMINAL, n_heldout = n_heldout,
                                    n_train = n_train, n_offsupport = n_off, seed = seed,
                                    footprint_rows = foot_rows, chunk_trees = chunk_trees,
                                    heldout = if (!is.null(heldout_dir)) heldout_dir else do_heldout,
                                    ranger = as.character(packageVersion("ranger")), platform = R.version$platform),
                    integrity = integrity, fit_check = fit_check, ladder = ladder_tab, taus = tau_tab,
                    performance = perf_tab, rows = rows_file[, model := stem][], footprint = foot)
        saveRDS(res, res_file)
        log_step(sprintf("  saved: %s (%.0f s)", res_file, as.numeric(difftime(Sys.time(), t_model, units = "secs"))))
        rm(tv_full, tv_red, S_full, P, X, res); invisible(gc(verbose = FALSE))
}
log_step("Done")

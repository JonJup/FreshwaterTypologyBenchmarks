### - Reduced QRF benchmark models for the online app - ###
#
# Usage:
#   Rscript 12_reduce_qrf.R --input <dir> --output <dir> [--trees 1000]
#                           [--compress none|gzip|xz] [--pattern <regex>]
#                           [--check-rows 200] [--overwrite]
#
#   --input       folder with the model files of 11_fit_qrf.R (<taxon>__<metric>.rds)
#   --output      folder for the reduced files (same names, so app.R reads them unchanged)
#   --trees k     trees kept per forest (default 1000)
#   --compress    saveRDS compression of the output (default none: fastest to load)
#   --pattern     regular expression on the file stem, e.g. "^invertebrates__"
#   --check-rows  training rows used to verify each reduced forest (default 200)
#   --overwrite   replace existing output files
#
# What is done. Nothing is refitted. The trees of a random forest are grown
# independently, so the first k trees are a random sample of the full forest: the
# same estimator, with more Monte Carlo noise in the per-tree predictive distribution.
#   model     : forest$child.nodeIDs, split.varIDs, split.values (one element per tree)
#               and the columns of random.node.values cut to the first k trees; rows of
#               random.node.values beyond the largest kept tree dropped. Out-of-bag fields
#               of the full forest are removed (predictions, inbag.counts,
#               random.node.values.oob) or set to NA (prediction.error, r.squared),
#               because they describe the full forest, not the reduced one.
#   recipe    : unchanged except recipe$template, cut to the columns app.R reads
#               (predictors, value, data.set).
#   reduction : provenance (source file, kept trees, date, ranger versions).
# Each reduced forest is verified before it is saved: its per-tree predictions for
# --check-rows training rows must be identical to those of the kept trees of the full
# forest. Files are written to a temporary name and renamed when complete.
#
# Memory: the full forest, its reduced copy and ranger's transient copy of the full
# forest during the check (about 1.9 x the forest) coexist. For the largest model
# (2.8 GB in memory) request 12-16 GB.
#
# Output: <output>/<taxon>__<metric>.rds and <output>/reduction_manifest.csv


library(data.table)
library(ranger)

## --- Command-line arguments ----

in_dir     <- "SHINY APP/models/"
out_dir    <- "SHINY APP/small_models"
k_trees    <- 1000
check_rows <- 200
compress = "xz"

if (is.null(in_dir) || is.null(out_dir)) stop("Missing --input or --output argument.")
if (normalizePath(in_dir) == normalizePath(out_dir, mustWork = FALSE)) stop("--output must differ from --input.")
stopifnot(!is.na(k_trees), k_trees >= 1, compress %in% c("none", "gzip", "xz"), check_rows >= 2)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
n_threads <- max(parallel::detectCores(logical = FALSE)-3, 1L, na.rm = TRUE)

log_step <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")
mb       <- function(x) as.numeric(object.size(x)) / 1e6


## --- Model discovery (as in app.R) ----
AUX_SUFFIXES <- c("ale", "datasets", "intervals", "lofo", "metrics",
                  "quantile_importance", "run_info", "split", "tuning")
stems <- sub("\\.rds$", "", list.files(in_dir, pattern = "__.+\\.rds$"))
stems <- stems[!grepl(paste0("_(", paste(AUX_SUFFIXES, collapse = "|"), ")$"), stems)]
if (!length(stems)) stop("No model files (<taxon>__<metric>.rds) in ", in_dir)

cat("\n========================================\n")
cat("Input:      ", in_dir, "\n")
cat("Output:     ", out_dir, "\n")
cat("Trees kept: ", k_trees, "(first k)\n")
cat("Compression:", compress, "\n")
cat("Models:     ", length(stems), "\n")
cat("Threads:    ", n_threads, "\n")
cat("========================================\n\n")


## --- Forest surgery ----
# Elements of ranger's forest with one entry per tree (regression forests). Any other
# per-tree element (e.g. from a different ranger version) stops the script rather
# than being carried over uncut.
PER_TREE   <- c("child.nodeIDs", "split.varIDs", "split.values")
OOB_FIELDS <- c("predictions", "inbag.counts", "random.node.values.oob")

subset_forest <- function(r, keep) {
        n0 <- r$num.trees
        stopifnot(inherits(r, "ranger"), r$treetype == "Regression", !is.null(r$random.node.values),
                  ncol(r$random.node.values) == n0, length(keep) >= 1, max(keep) <= n0, !anyDuplicated(keep))
        per_tree <- names(r$forest)[vapply(r$forest, function(e) is.list(e) && length(e) == n0, NA)]
        if (!setequal(per_tree, PER_TREE)) {
                stop("Per-tree forest elements differ from the expected set. Found: ",
                     paste(per_tree, collapse = ", "), "; expected: ", paste(PER_TREE, collapse = ", "))
        }
        for (nm in per_tree) r$forest[[nm]] <- r$forest[[nm]][keep]
        n_new <- length(keep)
        storage.mode(n_new) <- typeof(n0)                 # ranger stores num.trees as double
        r$forest$num.trees <- n_new
        r$num.trees        <- n_new
        # Node ids of a tree run from 1 to its number of nodes; rows beyond the largest
        # kept tree are padding.
        n_nodes <- max(lengths(r$forest$split.varIDs))
        r$random.node.values <- r$random.node.values[seq_len(n_nodes), keep, drop = FALSE]
        r[intersect(OOB_FIELDS, names(r))] <- NULL
        r$prediction.error <- NA_real_
        r$r.squared        <- NA_real_
        r
}

# Per-tree predictions (n x trees), as tree_values() in app.R
tree_values <- function(model, X) {
        p <- predict(model, data = X, type = "quantiles", what = function(v) v,
                     num.threads = n_threads)$predictions
        p <- as.matrix(p)
        if (nrow(p) != nrow(X)) p <- t(p)
        p
}

# Verification rows: training rows, spread evenly over data sets
check_data <- function(tmpl, preds, n, seed = 1) {
        d <- as.data.table(tmpl)[, c(preds, "data.set"), with = FALSE]
        d <- d[stats::complete.cases(d)]
        set.seed(seed)
        per <- ceiling(n / uniqueN(d$data.set))
        idx <- d[, .(i = .I[sample.int(.N, min(.N, per))]), by = data.set]$i
        as.data.frame(d[sort(head(sample(idx), n)), preds, with = FALSE])
}

ranger_version_fit <- function(stem) {
        f <- file.path(in_dir, paste0(stem, "_run_info.rds"))
        if (!file.exists(f)) return(NA_character_)
        v <- tryCatch(readRDS(f)$versions$ranger, error = function(e) NULL)
        if (is.null(v)) NA_character_ else as.character(v)
}


## --- Reduce each model ----
manifest <- vector("list", length(stems))
for (i in seq_along(stems)) {
        stem     <- stems[i]
        src      <- file.path(in_dir, paste0(stem, ".rds"))
        dst      <- file.path(out_dir, paste0(stem, ".rds"))
        t0       <- Sys.time()
        log_step(sprintf("%2d/%d %s", i, length(stems), stem))
        if (file.exists(dst) && !overwrite) {
                log_step("  exists, skipped (use --overwrite)")
                next
        }
        
        obj <- readRDS(src)
        if (!is.list(obj) || !inherits(obj$model, "ranger")) stop(stem, ": not a model file written by 11_fit_qrf.R.")
        r_full <- obj$model
        preds  <- r_full$forest$independent.variable.names
        n_full <- r_full$num.trees
        k      <- as.integer(min(k_trees, n_full))
        keep   <- seq_len(k)
        mem_full <- mb(r_full)
        
        r_red <- subset_forest(r_full, keep)
        
        ## Verification: identical per-tree predictions for the kept trees
        X  <- check_data(obj$recipe$template, preds, check_rows)
        tv_full <- tree_values(r_full, X)
        rm(r_full); obj$model <- NULL; invisible(gc(verbose = FALSE))
        tv_red  <- tree_values(r_red, X)
        ok <- all(dim(tv_red) == c(nrow(X), k)) && identical(unname(tv_red), unname(tv_full[, keep, drop = FALSE]))
        if (!ok) {
                stop(stem, ": reduced forest does not reproduce the kept trees (max abs difference ",
                     signif(max(abs(tv_red - tv_full[, keep, drop = FALSE])), 3), "). Nothing saved.")
        }
        rm(tv_full, tv_red); invisible(gc(verbose = FALSE))
        
        ## Training data: only the columns app.R uses
        tmpl     <- obj$recipe$template
        mem_tmpl <- mb(tmpl)
        need     <- c(preds, "value", "data.set")
        if (all(need %in% names(tmpl))) {
                obj$recipe$template <- tmpl[, need]
        } else {
                warning(stem, ": template lacks ", paste(setdiff(need, names(tmpl)), collapse = ", "),
                        "; left unchanged.", call. = FALSE)
        }
        rm(tmpl)
        
        obj$model     <- r_red
        obj$reduction <- list(
                source_file   = normalizePath(src),
                source_bytes  = file.size(src),
                n_trees_full  = n_full,
                n_trees       = k,
                kept_trees    = "first k (1:k)",
                created       = Sys.time(),
                ranger_fit    = ranger_version_fit(stem),
                ranger_reduce = as.character(packageVersion("ranger")),
                check_rows    = nrow(X),
                note          = paste("Trees subset from the fitted forest; out-of-bag fields removed.",
                                      "Calibration (_intervals.rds) refers to the full forest.")
        )
        
        tmp <- paste0(dst, ".tmp")
        saveRDS(obj, tmp, compress = if (compress == "none") FALSE else compress)
        if (!file.rename(tmp, dst)) stop("Could not rename ", tmp, " to ", dst)
        Sys.chmod(dst, mode = "0644")                     # readable by the arbitrary UID on Rahti
        
        manifest[[i]] <- data.table(
                model        = stem,
                trees_full   = n_full,
                trees        = k,
                mem_full_mb  = round(mem_full, 1),
                mem_mb       = round(mb(r_red), 1),
                forest_mb    = round(mb(r_red$forest), 1),
                rnv_mb       = round(mb(r_red$random.node.values), 1),
                template_mb  = round(mb(obj$recipe$template), 1),
                template_full_mb = round(mem_tmpl, 1),
                file_full_mb = round(file.size(src) / 1e6, 1),
                file_mb      = round(file.size(dst) / 1e6, 1),
                verified     = ok,
                seconds      = round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
        )
        log_step(sprintf("  %d -> %d trees | memory %.0f -> %.0f MB | file %.0f -> %.0f MB | verified (%d rows)",
                         n_full, k, mem_full, mb(r_red), file.size(src) / 1e6, file.size(dst) / 1e6, nrow(X)))
        # An uncompressed file much larger than its contents means something carries a
        # captured environment (e.g. recipe quosures created inside a function).
        contents <- mb(r_red) + mb(obj$recipe$template)
        if (compress == "none" && file.size(dst) / 1e6 > 1.2 * contents + 5) {
                warning(sprintf("%s: file (%.0f MB) exceeds its contents (%.0f MB); check obj$recipe for captured environments.",
                                stem, file.size(dst) / 1e6, contents), call. = FALSE)
        }
        rm(obj, r_red)
                invisible(gc(verbose = FALSE))
}

manifest <- rbindlist(manifest)
if (nrow(manifest)) {
        mf <- file.path(out_dir, "reduction_manifest.csv")
        if (file.exists(mf)) manifest <- rbind(fread(mf)[!model %in% manifest$model], manifest)
        setorder(manifest, model)
        fwrite(manifest, mf)
        cat("\n")
        print(manifest[, .(model, trees, mem_mb, forest_mb, file_mb, verified)])
        cat(sprintf("\nTotal: memory %.1f -> %.1f GB | files %.1f -> %.1f GB\n",
                    sum(manifest$mem_full_mb) / 1e3, sum(manifest$mem_mb) / 1e3,
                    sum(manifest$file_full_mb) / 1e3, sum(manifest$file_mb) / 1e3))
        cat("  saved:", mf, "\n")
}
log_step("Done")

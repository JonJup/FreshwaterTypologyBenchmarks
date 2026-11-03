################################################################################*
# prepare_qrf_inputs.R
# Shared input preparation for 11_fit_qrf.R and 11a_check_qrf_inputs.R, so the
# check script prepares the data exactly as the fitting script does.
#
# Subsets are identified by FILE NAME, not by the scheme_id stored inside the
# files: the stored variable follows a slightly different convention in some
# sources, which made the old join on the stored variable drop whole subsets.
# Every identifier is mapped to a canonical form "<taxon>_<number, 4 digits>"
# (e.g. "fish_0001"), which absorbs differences in padding, prefixes and
# separators. Where the stored identifier denotes a different subset, it is a
# stale label: the file name is authoritative.
################################################################################*

QRF_TAXA          <- c("diatoms", "fish", "invertebrates", "macrophytes")
QRF_DATE_VARS     <- c("mean_doy", "circ_sd_days")
qrf_predictor_vec <- c(
  "samples", "n_types", "n_vars", "type_evenness",
  "env_r2_after", "env_ch_after", "asw_after", "min_distance",
  "median_distance", "species_rank", "genus_rank", "family_rank",
  "higher_rank", "n_taxa_used",  "vp_env_pure", "vp_shared", "vp_space_pure",
  QRF_DATE_VARS
)
qrf_model_cols <- c("value", qrf_predictor_vec, "scheme_id", "data.set")

## Canonical subset identifier ----
# taxon = first occurrence of a taxon name, number = trailing digit run.
# Anything else in the string (doubled prefixes, data-set names, separators,
# zero padding) is ignored. Returns NA if either part is missing.
canonical_id <- function(x) {
  x   <- tolower(trimws(as.character(x)))
  any_taxon <- paste(QRF_TAXA, collapse = "|")
  tax <- rep(NA_character_, length(x))
  num <- rep(NA_character_, length(x))
  has_tax <- !is.na(x) & grepl(any_taxon, x)
  has_num <- !is.na(x) & grepl("[0-9]+$", x)
  tax[has_tax] <- sub(paste0("^.*?(", any_taxon, ").*$"), "\\1", x[has_tax], perl = TRUE)
  num[has_num] <- sub("^.*?([0-9]+)$", "\\1", x[has_num], perl = TRUE)
  out <- sprintf("%s_%04d", tax, suppressWarnings(as.integer(num)))
  out[is.na(tax) | is.na(num)] <- NA_character_
  out
}

## Read a directory, keyed by file name ----
# scheme_id        : canonical id of the FILE NAME (used for all joins)
# scheme_id_stored : the variable stored in the file (kept for checks only)
read_dir_keyed <- function(path, fill = FALSE) {
  files <- list.files(path, full.names = TRUE)
  if (!length(files)) stop("No files in ", path)
  rbindlist(lapply(files, function(f) {
    d <- as.data.table(readRDS(f))
    # copy(): as.character() on a character column returns the same vector, which the
    # by-reference := below would otherwise overwrite with the file-name id.
    stored <- if ("scheme_id" %in% names(d)) copy(as.character(d$scheme_id)) else NA_character_
    d[, `:=`(scheme_id        = canonical_id(tools::file_path_sans_ext(basename(f))),
             scheme_id_stored = stored,
             source_file      = basename(f))]
  }), fill = fill)
}

## Load and merge ----
# Returns the merged evaluation table (same columns and row order logic as the
# original section 4 of 11_fit_qrf.R) plus two check tables:
#   id_check : one row per file (per id for the date table): canonical id,
#              stored id, whether they differ in format or denote another subset
#   coverage : one row per evaluated subset: found in spatial / taxonomy / dates
load_qrf_inputs <- function(eval_dir = "data/007_evaluations/",
                            ss_dir   = "data/misc/spatial_scale/",
                            tx_dir   = "data/misc/taxonomic_resolution/",
                            dv_file  = "data/misc/date_variables.rds",
                            log      = function(...) cat(..., "\n")) {
  log("Loading evaluation data")
  evl <- read_dir_keyed(eval_dir)
  log("Loading spatial scale data")
  ss  <- read_dir_keyed(ss_dir, fill = TRUE)
  log("Loading taxonomic resolution data")
  tx  <- read_dir_keyed(tx_dir)
  log("Loading sampling-date variables")
  dv  <- as.data.table(readRDS(dv_file))
  setnames(dv, 1L, "scheme_id_stored")
  dv[, `:=`(scheme_id_stored = as.character(scheme_id_stored),
            scheme_id        = canonical_id(scheme_id_stored),
            source_file      = basename(dv_file))]
  
  ## 1. Identifier check ----
  id_cols  <- c("source_file", "scheme_id", "scheme_id_stored")
  id_check <- rbindlist(list(
    evaluation = unique(evl[, ..id_cols]),
    spatial    = unique(ss[,  ..id_cols]),
    taxonomy   = unique(tx[,  ..id_cols]),
    dates      = unique(dv[,  ..id_cols])
  ), idcol = "source")
  id_check[, `:=`(
    unparsed         = is.na(scheme_id),
    format_differs   = !is.na(scheme_id_stored) & scheme_id_stored != scheme_id,
    denotes_other    = !is.na(scheme_id_stored) & !is.na(scheme_id) &
      !is.na(canonical_id(scheme_id_stored)) &
      canonical_id(scheme_id_stored) != scheme_id,
    stored_unparsed  = !is.na(scheme_id_stored) & is.na(canonical_id(scheme_id_stored))
  )]
  # several stored ids in one per-subset file = ambiguous; the date table is one file by design
  id_check[, n_stored_per_file := if (source[1] == "dates") 1L else uniqueN(scheme_id_stored),
           by = .(source, source_file)]
  
  ## 2. One row per subset in the subset-level tables ----
  ss <- unique(ss[, .(scheme_id, min_distance, median_distance)])
  tx <- tx[, setdiff(names(tx), c("taxon", "scheme_id_stored", "source_file")), with = FALSE]
  dv <- dv[, c("scheme_id", QRF_DATE_VARS), with = FALSE]
  for (nm in c("ss", "tx", "dv")) {
    d <- get(nm)
    if (anyDuplicated(d$scheme_id)) {
      stop("Duplicated subset in ", nm, ": ",
           paste(head(unique(d$scheme_id[duplicated(d$scheme_id)]), 10), collapse = ", "))
    }
  }
  
  ## 3. Coverage of every evaluated subset ----
  coverage <- unique(evl[, .(scheme_id, source_file)])
  coverage[, `:=`(taxon       = sub("_[0-9]+$", "", scheme_id),
                  in_spatial  = scheme_id %in% ss$scheme_id,
                  in_taxonomy = scheme_id %in% tx$scheme_id,
                  in_dates    = scheme_id %in% dv$scheme_id)]
  
  ## 4. Merge (as in the original section 4) ----
  log("Merging datasets")
  n_evl <- nrow(evl)
  evl[, c("scheme_id_stored", "source_file") := NULL]
  evl2 <- ss[evl, on = "scheme_id"]
  evl2 <- tx[evl2, on = "scheme_id"]
  evl2 <- dv[evl2, on = "scheme_id"]
  drop_cols <- grep("^i\\.", names(evl2), value = TRUE)
  if (length(drop_cols)) evl2[, (drop_cols) := NULL]
  stopifnot("Row count changed unexpectedly after merge" = nrow(evl2) == n_evl)
  evl2[, taxon_key := sub("_[0-9]+$", "", scheme_id)]
  
  list(evl2 = evl2, id_check = id_check, coverage = coverage)
}

## Structural identifier problems (always fatal) ----
# Stored ids are stale labels; the file name is authoritative. Only problems that
# make the file-name id itself unusable are fatal.
qrf_id_problems <- function(id_check) {
  id_check[unparsed | n_stored_per_file > 1L]
}

## Rows that would be dropped for missing values ----
# One row per subset with at least one incomplete row; `whole_subset` marks
# subsets lost entirely. Uses exactly the columns of the na.omit in 11.
qrf_drop_report <- function(d, model_cols = qrf_model_cols) {
  chk <- setdiff(model_cols, c("scheme_id", "data.set"))
  rep <- d[, {
    inc <- !complete.cases(.SD)
    list(n_rows     = .N,
         n_dropped  = sum(inc),
         na_columns = paste(names(.SD)[vapply(.SD, anyNA, logical(1))], collapse = ";"))
  }, by = .(scheme_id, data.set), .SDcols = chk][n_dropped > 0]
  rep[, whole_subset := n_dropped == n_rows]
  rep[]
}
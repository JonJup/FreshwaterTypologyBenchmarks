### PULSE: benchmark-model predictors from community, environment and typology ###
#
# Computes, for one or more typologies of a set of sites, the predictors of the
# PULSE benchmark models and writes them as a CSV that the batch tab of the app
# reads directly (one row per typology, columns named as the predictors).
# Predictors that cannot be computed from the inputs are left empty; the app
# imputes them.
#
# A second CSV holds the observed type coherence of each typology, computed as
# in the simulation pipeline (10_evaluate_typologies.R). Optionally, a third CSV
# combines both for the app's batch tab: one row per typology x metric with the
# columns typology, metric, observed_value (and taxon, if given) and the predictors.
#
# Inputs (data frames, sf objects or file paths: .csv/.tsv/.txt, .rds, .gpkg/.shp/.geojson)
#   community  wide: one row per site (or sample), a site column and one column per taxon
#              long: one row per site (or sample) x taxon, with site and taxon columns,
#                    optionally abundance, sample id and sampling date
#   env        wide: one row per site, a site column and one column per variable
#              long: site, variable and value columns
#   typology   a) type per site: data frame (site, type), named vector (names = sites),
#                 or unnamed vector in the order of the community sites
#              b) spatial: sf polygons or lines with a type attribute; sites are
#                 matched by location (requires `sites`)
#              A named list of typologies gives one output row per typology.
#   sites      site coordinates: data frame (site, x, y; CRS in `crs`) or sf points.
#              Needed for distances, MEMs (variation partitioning) and spatial typologies.
#   dates      optional sampling dates: data frame (site, date), or a date column in
#              the community data. Needed for mean_doy and circ_sd_days.
#   typology_vars  names of the environmental variables the typology is built on.
#              Needed for n_vars and the variation partitioning.
#   taxon_ranks    optional lookup (taxon, rank) with rank in species / genus /
#              family / higher. Without it, ranks are guessed from the names.
#
# Definitions follow the PULSE simulation pipeline (methods sections 2.4.3, S6):
#   samples           number of sites
#   n_types           number of types
#   n_vars            number of typology variables
#   type_evenness     Pielou evenness of the type sizes
#   env_r2_after      share of the variation of all standardised environmental
#                     variables explained by the types (PERMANOVA R2 on Euclidean distances)
#   env_ch_after      Calinski-Harabasz index, same space
#   asw_after         average silhouette width, same space
#   min_distance,
#   median_distance   smallest and median distance between site locations (km)
#   species_rank ...
#   higher_rank       share of taxa identified to species / genus / family / above family
#   n_taxa_used       taxa occurring at least once
#   vp_env_pure,
#   vp_shared,
#   vp_space_pure     adjusted R2 fractions [a], [b], [c] of the db-RDA variation
#                     partitioning of binary Jaccard dissimilarities (square-rooted)
#                     between the typology variables and up to five MEMs
#   mean_doy,
#   circ_sd_days      circular mean and circular SD (days) of the sampling dates,
#                     sites weighted equally
#
# Coherence metrics of the crisp typology (binary Jaccard dissimilarities of
# presence-absence data); fuzzy typologies are not supported:
#   Classification Strength   vegan::meandist(), summary()$CS
#   ANOSIM R                  vegan::anosim()$statistic
#   PERMANOVA R2              vegan::adonis2(d ~ type)$R2
#   AucZeta                   mean over types (>= 2 sites) of the area under the
#                             within-type zeta decline minus that of all sites
#
# Usage: edit the RUN block at the end and run `Rscript pulse_predictors.R`,
# or source() the file and call pulse_predictors() yourself.

suppressPackageStartupMessages(library(data.table))

# ==============================================================================*
# 1. SETTINGS ----
# ==============================================================================*

PROJ_CRS        <- 3035    # projected CRS for distances, MEMs and spatial matching (ETRS89-LAEA Europe)
MIN_OCCURRENCES <- 1L      # taxa with fewer occurrences are dropped (training communities contained
# only the taxa kept for the HMSC models; set to that threshold to match)
LINE_SNAP_MAX_M <- 500     # line typologies: max distance (m) from a site to its nearest line
POLY_SNAP_MAX_M <- 0       # polygon typologies: sites outside all polygons take the nearest
# polygon within this distance (m); 0 = such sites are dropped

# Variation partitioning (as in the simulation pipeline)
MEM_PREVALENCE  <- 0.05    # taxa used for MEM selection: prevalence >= this
MEM_ALPHA       <- 0.005   # Holm-adjusted p-value threshold
MEM_MAX         <- 5L      # at most this many MEMs

# Coherence (as in 10_evaluate_typologies.R)
ZETA_MAX <- 10    # highest zeta order

METRICS <- c("Classification Strength", "ANOSIM R", "AucZeta", "PERMANOVA R2")

PREDICTORS <- c("samples", "n_types", "n_vars", "type_evenness",
                "env_r2_after", "env_ch_after", "asw_after",
                "min_distance", "median_distance",
                "species_rank", "genus_rank", "family_rank", "higher_rank",
                "n_taxa_used", "vp_env_pure", "vp_shared", "vp_space_pure",
                "mean_doy", "circ_sd_days")

has_pkg <- function(p) suppressMessages(requireNamespace(p, quietly = TRUE))
say     <- function(...) cat(sprintf(...), "\n", sep = "")


# ==============================================================================*
# 2. INPUT HANDLING ----
# ==============================================================================*

read_any <- function(x) {
        if (is.null(x) || !is.character(x) || length(x) != 1L) return(x)
        if (!file.exists(x)) stop("File not found: ", x)
        ext <- tolower(tools::file_ext(x))
        switch(ext,
               csv = , tsv = , txt = fread(x),
               rds = readRDS(x),
               gpkg = , shp = , geojson = , json = , kml = {
                       if (!has_pkg("sf")) stop("Package 'sf' is needed to read ", x)
                       sf::st_read(x, quiet = TRUE)
               },
               stop("Unsupported file type: ", x))
}

drop_geometry <- function(x) if (inherits(x, "sf")) sf::st_drop_geometry(x) else x

check_cols <- function(x, cols, what) {
        miss <- setdiff(cols, names(x))
        if (length(miss)) stop(what, ": column(s) not found: ", paste(miss, collapse = ", "))
}

# Community -> list(Y = 0/1 site x taxon matrix, dates = data.table(site, date) or NULL)
read_community <- function(x, site_col = "site", taxon_col = "taxon", value_col = NULL,
                           date_col = "date", format = c("auto", "wide", "long"), exclude = NULL) {
        format <- match.arg(format)
        x <- as.data.table(drop_geometry(read_any(x)))
        check_cols(x, site_col, "community")
        if (format == "auto") format <- if (taxon_col %in% names(x)) "long" else "wide"
        has_date <- !is.null(date_col) && date_col %in% names(x)
        
        if (format == "long") {
                check_cols(x, taxon_col, "community")
                v <- if (!is.null(value_col)) { check_cols(x, value_col, "community"); as.numeric(x[[value_col]]) } else 1
                d <- data.table(site = as.character(x[[site_col]]), taxon = as.character(x[[taxon_col]]), v = v)
                d <- d[!is.na(site) & !is.na(taxon) & nzchar(taxon)]
                d[is.na(v), v := 0]
                Y <- dcast(d, site ~ taxon, value.var = "v", fill = 0,
                           fun.aggregate = function(z) as.numeric(any(z > 0)))
        } else {
                taxa <- setdiff(names(x), c(site_col, date_col, exclude))
                num  <- vapply(x[, ..taxa], is.numeric, TRUE)
                if (any(!num)) say("  community: ignoring non-numeric columns: %s", paste(taxa[!num], collapse = ", "))
                taxa <- taxa[num]
                d <- x[, c(site_col, taxa), with = FALSE]
                setnames(d, site_col, "site")
                d[, site := as.character(site)]
                for (tx in taxa) set(d, j = tx, value = as.numeric(d[[tx]] > 0 & !is.na(d[[tx]])))
                Y <- d[, lapply(.SD, max), by = site]          # several samples per site -> presence in any
        }
        sites <- Y$site
        Y <- as.matrix(Y[, -"site"]); rownames(Y) <- sites
        dates <- if (has_date) unique(data.table(site = as.character(x[[site_col]]), date = x[[date_col]])) else NULL
        list(Y = Y, dates = dates)
}

# Environment -> numeric matrix, rows = sites
read_env <- function(x, site_col = "site", var_col = "variable", value_col = "value",
                     format = c("auto", "wide", "long"), exclude = NULL) {
        format <- match.arg(format)
        x <- as.data.table(drop_geometry(read_any(x)))
        check_cols(x, site_col, "env")
        if (format == "auto") format <- if (all(c(var_col, value_col) %in% names(x))) "long" else "wide"
        if (format == "long") {
                d <- data.table(site = as.character(x[[site_col]]), variable = as.character(x[[var_col]]),
                                value = as.numeric(x[[value_col]]))
                d <- dcast(d, site ~ variable, value.var = "value", fun.aggregate = mean, fill = NA)
        } else {
                vars <- setdiff(names(x), c(site_col, exclude))
                num  <- vapply(x[, ..vars], is.numeric, TRUE)
                if (any(!num)) say("  env: ignoring non-numeric columns: %s", paste(vars[!num], collapse = ", "))
                d <- x[, c(site_col, vars[num]), with = FALSE]
                setnames(d, site_col, "site")
                d[, site := as.character(site)]
                if (anyDuplicated(d$site)) {
                        say("  env: several rows per site, averaging them")
                        d <- d[, lapply(.SD, mean, na.rm = TRUE), by = site]
                }
        }
        E <- as.matrix(d[, -"site"]); rownames(E) <- d$site
        E[is.nan(E)] <- NA
        E
}

# Sites -> sf points in PROJ_CRS with a `site` column
read_sites <- function(x, site_col = "site", x_col = "x", y_col = "y", crs = 4326) {
        if (is.null(x)) return(NULL)
        if (!has_pkg("sf")) stop("Package 'sf' is needed for site coordinates.")
        x <- read_any(x)
        if (inherits(x, "sf")) {
                check_cols(x, site_col, "sites")
                pts <- sf::st_sf(site = as.character(x[[site_col]]), geometry = sf::st_geometry(x))
                if (!all(sf::st_geometry_type(pts) == "POINT")) pts <- suppressWarnings(sf::st_centroid(pts))
        } else {
                x <- as.data.table(x)
                check_cols(x, c(site_col, x_col, y_col), "sites")
                x <- unique(x[, c(site_col, x_col, y_col), with = FALSE])
                x <- x[!is.na(get(x_col)) & !is.na(get(y_col))]
                pts <- sf::st_as_sf(data.frame(site = as.character(x[[site_col]]),
                                               X = as.numeric(x[[x_col]]), Y = as.numeric(x[[y_col]])),
                                    coords = c("X", "Y"), crs = crs)
        }
        if (anyDuplicated(pts$site)) {
                say("  sites: several coordinates per site, keeping the first")
                pts <- pts[!duplicated(pts$site), ]
        }
        sf::st_transform(pts, PROJ_CRS)
}

# Typology -> named character vector (names = sites), NA where no type
assign_types <- function(typology, site_ids, pts = NULL, site_col = "site", type_col = "type") {
        typology <- read_any(typology)
        
        if (inherits(typology, "sf")) {
                if (is.null(pts)) stop("A spatial typology needs site coordinates (`sites`).")
                check_cols(typology, type_col, "typology")
                typ  <- sf::st_transform(typology[, type_col], PROJ_CRS)
                pts  <- pts[pts$site %in% site_ids, ]
                gt   <- unique(as.character(sf::st_geometry_type(typ)))
                out  <- setNames(rep(NA_character_, length(site_ids)), site_ids)
                
                if (all(gt %in% c("POLYGON", "MULTIPOLYGON"))) {
                        typ <- sf::st_make_valid(typ)
                        hit <- sf::st_intersects(pts, typ)
                        n_hit <- lengths(hit)
                        if (any(n_hit > 1)) say("  typology: %d sites fall in several polygons, using the first", sum(n_hit > 1))
                        idx <- vapply(hit, function(h) if (length(h)) h[1] else NA_integer_, 1L)
                        if (POLY_SNAP_MAX_M > 0 && anyNA(idx)) {
                                miss <- which(is.na(idx))
                                near <- sf::st_nearest_feature(pts[miss, ], typ)
                                dd   <- as.numeric(sf::st_distance(pts[miss, ], typ[near, ], by_element = TRUE))
                                idx[miss[dd <= POLY_SNAP_MAX_M]] <- near[dd <= POLY_SNAP_MAX_M]
                        }
                        out[pts$site] <- as.character(typ[[type_col]][idx])
                } else if (all(gt %in% c("LINESTRING", "MULTILINESTRING"))) {
                        near <- sf::st_nearest_feature(pts, typ)
                        dd   <- as.numeric(sf::st_distance(pts, typ[near, ], by_element = TRUE))
                        if (any(dd > LINE_SNAP_MAX_M)) {
                                say("  typology: %d sites are more than %g m from any line and get no type",
                                    sum(dd > LINE_SNAP_MAX_M), LINE_SNAP_MAX_M)
                        }
                        val <- as.character(typ[[type_col]][near]); val[dd > LINE_SNAP_MAX_M] <- NA
                        out[pts$site] <- val
                } else {
                        stop("Spatial typology must contain only polygons or only lines (found: ", paste(gt, collapse = ", "), ").")
                }
                return(out)
        }
        
        if (is.data.frame(typology)) {
                typology <- as.data.table(typology)
                check_cols(typology, c(site_col, type_col), "typology")
                m <- setNames(as.character(typology[[type_col]]), as.character(typology[[site_col]]))
                return(setNames(unname(m[site_ids]), site_ids))
        }
        if (is.atomic(typology)) {
                if (!is.null(names(typology))) return(setNames(as.character(typology[site_ids]), site_ids))
                if (length(typology) != length(site_ids)) {
                        stop("An unnamed typology vector must have one entry per community site (", length(site_ids), ").")
                }
                say("  typology: unnamed vector, matched to the community sites by position")
                return(setNames(as.character(typology), site_ids))
        }
        stop("Unsupported typology input.")
}


# ==============================================================================*
# 3. PREDICTORS ----
# ==============================================================================*

## 3.1 Typology descriptors in the space of all standardised environmental variables
typology_descriptors <- function(E, types) {
        out <- list(n_types = NA_real_, type_evenness = NA_real_,
                    env_r2_after = NA_real_, env_ch_after = NA_real_, asw_after = NA_real_)
        f <- factor(types); k <- nlevels(f); n <- length(f)
        out$n_types <- k
        p <- as.numeric(table(f)) / n
        if (k > 1) out$type_evenness <- -sum(p * log(p)) / log(k)
        
        keep <- apply(E, 2, function(z) isTRUE(sd(z) > 0))
        if (!any(keep) || k < 2 || n <= k) return(out)
        Z <- scale(E[, keep, drop = FALSE])
        
        ss_tot <- sum(Z^2)                                              # columns are centred
        ss_w   <- sum(vapply(split(seq_len(n), f), function(i)
                sum(sweep(Z[i, , drop = FALSE], 2, colMeans(Z[i, , drop = FALSE]))^2), 0))
        out$env_r2_after <- 1 - ss_w / ss_tot
        if (ss_w > 0) out$env_ch_after <- ((ss_tot - ss_w) / (k - 1)) / (ss_w / (n - k))
        if (has_pkg("cluster")) {
                sil <- cluster::silhouette(as.integer(f), dist(Z))
                out$asw_after <- mean(sil[, "sil_width"])
        }
        out
}

## 3.2 Distances between site locations (km)
site_distances <- function(pts) {
        xy <- unique(round(sf::st_coordinates(pts), 3))
        if (nrow(xy) < 2) return(list(min_distance = NA_real_, median_distance = NA_real_))
        d <- dist(xy) / 1000
        list(min_distance = min(d), median_distance = median(d))
}

## 3.3 Taxonomic resolution
# Guess the rank from a taxon name: two words with a lowercase epithet -> species;
# family-group suffixes -> family; names of higher taxa -> higher; else genus.
QUALIFIERS <- c("sp", "spp", "ssp", "cf", "aff", "agg", "gr", "group", "gen", "indet", "sect",
                "subg", "subgen", "juv", "lv", "ad", "larvae", "s", "l", "str", "sl", "complex", "type")
HIGHER_SUFFIX <- "(oidea|ptera|poda|ales|phyceae|phyta|phytina|opsida|mycetes|aria|acea)$"
HIGHER_NAMES  <- c("oligochaeta", "hirudinea", "turbellaria", "nematoda", "nematomorpha", "hydrachnidia",
                   "hydracarina", "acari", "acarina", "bivalvia", "gastropoda", "ostracoda", "copepoda",
                   "cladocera", "bryozoa", "porifera", "hydrozoa", "cnidaria", "tricladida", "odonata",
                   "annelida", "mollusca", "crustacea", "insecta", "arachnida", "collembola", "isopoda",
                   "amphipoda", "decapoda", "polychaeta", "nemertea", "tardigrada", "rotifera",
                   "bryophyta", "charophyceae", "chlorophyta", "cyanobacteria", "bacillariophyta")

guess_rank <- function(name) {
        w <- strsplit(trimws(gsub("[._]+", " ", name)), "\\s+")
        vapply(w, function(x) {
                x <- x[!tolower(x) %in% QUALIFIERS & !grepl("^\\(|\\)$|^[0-9]", x)]
                if (length(x) >= 2 && grepl("^[a-z-]+$", x[2])) return("species")   # lowercase epithet
                if (!length(x)) return(NA_character_)
                g <- tolower(x[1])
                if (grepl("(idae|aceae|inae|ini)$", g)) return("family")
                if (g %in% HIGHER_NAMES || grepl(HIGHER_SUFFIX, g)) return("higher")
                "genus"
        }, "")
}

taxonomic_shares <- function(taxa, taxon_ranks = NULL) {
        rk <- setNames(rep(NA_character_, length(taxa)), taxa)
        if (!is.null(taxon_ranks)) {
                tr <- as.data.table(read_any(taxon_ranks))
                check_cols(tr, c("taxon", "rank"), "taxon_ranks")
                m  <- setNames(tolower(as.character(tr$rank)), as.character(tr$taxon))
                rk[] <- m[taxa]
        }
        guess <- is.na(rk)
        rk[guess] <- guess_rank(taxa[guess])
        rk[rk %in% c("subspecies", "variety", "form", "species group", "aggregate")] <- "species"
        rk[rk %in% c("subfamily", "tribe", "subtribe")] <- "family"
        rk[!rk %in% c("species", "genus", "family") & !is.na(rk)] <- "higher"
        list(shares = vapply(c(species = "species", genus = "genus", family = "family", higher = "higher"),
                             function(r) mean(rk == r, na.rm = TRUE), 0),
             ranks = rk, guessed = sum(guess))
}

## 3.4 Variation partitioning with MEMs (as community_predictors() in the simulation pipeline)

# Candidate MEMs: distance-based MEMs with non-null autocorrelation. Uses
# adespatial::dbmem() when installed; otherwise the same construction directly
# (truncation at the longest edge of the minimum spanning tree, weights
# 1 - (d / 4t)^2, eigenvectors of the doubly centred weight matrix).
candidate_mems <- function(xy) {
        if (has_pkg("adespatial")) {
                m <- adespatial::dbmem(xy, MEM.autocor = "non-null", store.listw = FALSE, silent = TRUE)
                return(as.matrix(m))
        }
        if (!has_pkg("vegan")) stop("Package 'vegan' (or 'adespatial') is needed for MEMs.")
        D <- as.matrix(dist(xy)); n <- nrow(D)
        t <- max(vegan::spantree(as.dist(D))$dist)
        W <- (1 - (D / (4 * t))^2) * (D > 0 & D <= t * (1 + 1e-6))
        Wc <- W - outer(rowMeans(W), colMeans(W), "+") + mean(W)
        e <- eigen(Wc, symmetric = TRUE)
        keep <- abs(e$values) > 1e-8 * max(abs(e$values))
        M <- e$vectors[, keep, drop = FALSE] * sqrt(n)
        colnames(M) <- paste0("MEM", seq_len(ncol(M)))
        M
}

# Up to MEM_MAX MEMs that correlate most often with the residuals of taxon-wise
# binomial GLMs on the typology variables (Holm-adjusted p < MEM_ALPHA).
select_mems <- function(Y, E, M) {
        counts <- integer(ncol(M)); best_p <- rep(Inf, ncol(M))
        Ed <- as.data.frame(E)
        for (j in seq_len(ncol(Y))) {
                y <- Y[, j]
                if (var(y) == 0 || mean(y) < MEM_PREVALENCE) next
                mod <- tryCatch(suppressWarnings(glm(y ~ ., data = cbind(data.frame(y = y), Ed), family = binomial)),
                                error = function(e) NULL)
                if (is.null(mod)) next
                r  <- residuals(mod)                                    # deviance residuals
                pa <- p.adjust(apply(M, 2, function(m) cor.test(m, r)$p.value), method = "holm")
                counts <- counts + (pa < MEM_ALPHA)
                best_p <- pmin(best_p, pa)
        }
        n_max <- min(MEM_MAX, nrow(Y) - 2L - qr(as.matrix(E))$rank)
        if (n_max < 1L) return(NULL)
        ord <- order(-counts, best_p)
        sel <- if (any(counts > 0)) head(ord[counts[ord] > 0], n_max) else which.min(best_p)
        M[, sel, drop = FALSE]
}

variation_partitioning <- function(Y, E, xy) {
        out <- list(n_taxa_used = NA_real_, vp_env_pure = NA_real_, vp_shared = NA_real_,
                    vp_space_pure = NA_real_, n_mem = NA_real_)
        sp <- colSums(Y) > 0
        si <- rowSums(Y[, sp, drop = FALSE]) > 0
        out$n_taxa_used <- sum(sp)
        if (is.null(E) || is.null(xy) || !has_pkg("vegan")) return(out)
        M  <- candidate_mems(xy)
        Y  <- Y[si, sp, drop = FALSE]
        E  <- E[si, , drop = FALSE]
        M  <- M[si, , drop = FALSE]
        if (nrow(Y) <= ncol(E) + 2L || ncol(Y) < 2L) return(out)
        Ms <- select_mems(Y, E, M)
        if (is.null(Ms)) return(out)
        out$n_mem <- ncol(Ms)
        d  <- vegan::vegdist(Y, method = "jaccard", binary = TRUE)
        vp <- tryCatch(suppressWarnings(vegan::varpart(d, as.data.frame(E), as.data.frame(Ms), sqrt.dist = TRUE)),
                       error = function(e) NULL)
        if (is.null(vp)) return(out)
        f <- vp$part$indfract$Adj.R.squared                              # [a], [b], [c], [d]
        out$vp_env_pure <- f[1]; out$vp_shared <- f[2]; out$vp_space_pure <- f[3]
        out
}

## 3.5 Sampling dates (as sampling_time() in the pipeline)
sampling_season <- function(dates) {
        out <- list(mean_doy = NA_real_, circ_sd_days = NA_real_)
        if (is.null(dates) || !nrow(dates)) return(out)
        d <- copy(dates)
        d[, date := as.Date(date)]
        d <- d[!is.na(date)]
        if (!nrow(d)) return(out)
        yr <- as.integer(format(d$date, "%Y"))
        n_days <- ifelse((yr %% 4 == 0 & yr %% 100 != 0) | yr %% 400 == 0, 366, 365)
        d[, theta := 2 * pi * (as.integer(format(date, "%j")) - 0.5) / n_days]
        s <- d[, .(s_sin = mean(sin(theta)), s_cos = mean(cos(theta))), by = site]    # sites weigh equally
        ts <- mean(s$s_sin); tc <- mean(s$s_cos)
        R  <- sqrt(ts^2 + tc^2)
        if (R > 1e-8) out$mean_doy <- (atan2(ts, tc) %% (2 * pi)) / (2 * pi) * 365.25 + 0.5
        out$circ_sd_days <- sqrt(-2 * log(min(R, 1))) * 365.25 / (2 * pi)
        out
}


## 3.6 Coherence metrics

# Exact (expected) zeta diversity of orders 1..k from species occupancies
# (zetadiv::Zeta.order.ex), rescaled by zeta_1. From exact_zeta.R.
exact_zeta <- function(pa, orders, rescale = TRUE) {
        pa  <- as.matrix(pa) > 0
        N   <- nrow(pa)
        if (max(orders) > N) stop("Zeta order ", max(orders), " exceeds number of sites (", N, ").")
        occ <- colSums(pa)
        occ <- occ[occ > 0]
        if (length(occ) == 0) return(rep(NA_real_, length(orders)))
        z <- vapply(orders, function(i) sum(exp(lchoose(occ, i) - lchoose(N, i))), numeric(1))
        if (rescale) z <- z / z[1]
        z
}

# Area under a curve by the trapezoidal rule. From calculate_auc.R.
calculate_auc <- function(x, y) {
        # Sort x and y together (in case they're not already ordered)
        ord <- order(x)
        x <- x[ord]
        y <- y[ord]
        
        # Calculate width of each trapezoid
        dx <- diff(x)
        
        # Calculate mean height of each trapezoid
        mean_height <- (y[-1] + y[-length(y)]) / 2
        
        # Calculate area as sum of trapezoid areas
        area <- sum(dx * mean_height)
        
        return(area)
}

# AucZeta: for every type with >= 2 sites, AUC of its zeta decline (orders
# 1..min(n_t, ZETA_MAX)) minus the AUC of all sites over the same orders; mean over types.
auc_zeta <- function(Y, types) {
        sbt <- split(seq_along(types), types)
        sbt <- sbt[lengths(sbt) >= 2]
        if (!length(sbt)) return(NA_real_)
        z_base <- exact_zeta(Y, orders = seq_len(min(ZETA_MAX, max(lengths(sbt)))))
        rel <- vapply(sbt, function(idx) {
                ord <- seq_len(min(length(idx), ZETA_MAX))
                z_w <- exact_zeta(Y[idx, , drop = FALSE], orders = ord)
                calculate_auc(ord, z_w) - calculate_auc(ord, z_base[ord])
        }, 0)
        mean(rel, na.rm = TRUE)
}

coherence_metrics <- function(Y, types) {
        out <- setNames(as.list(rep(NA_real_, length(METRICS))), METRICS)
        n_empty <- sum(rowSums(Y) == 0)
        if (n_empty) say("  coherence: %d sites without taxa (kept, as in the simulations)", n_empty)
        d <- dist(Y > 0, method = "binary")                    # binary Jaccard, as parallelDist(method = "binary")
        g <- factor(types)
        if (nlevels(g) < 2) return(out)
        out[["Classification Strength"]] <- unname(unlist(summary(vegan::meandist(d, g)))["CS"])
        out[["ANOSIM R"]]     <- unname(vegan::anosim(d, g, permutations = 0)$statistic)
        out[["AucZeta"]]      <- auc_zeta(Y, as.character(types))
        out[["PERMANOVA R2"]] <- vegan::adonis2(d ~ g, permutations = 0)$R2[1]
        out
}


# ==============================================================================*
# 4. MAIN FUNCTION ----
# ==============================================================================*

pulse_predictors <- function(community, env, typology, sites = NULL, dates = NULL,
                             typology_vars = NULL, taxon_ranks = NULL,
                             site_col = "site", taxon_col = "taxon", value_col = NULL, date_col = "date",
                             community_format = "auto", env_format = "auto",
                             env_var_col = "variable", env_value_col = "value",
                             x_col = "x", y_col = "y", crs = 4326, type_col = "type",
                             exclude = NULL, output = NULL, coherence_output = NULL,
                             batch_output = NULL, taxon = NULL) {
        
        say("Reading inputs")
        comm <- read_community(community, site_col, taxon_col, value_col, date_col, community_format,
                               exclude = c(exclude, x_col, y_col))
        Eall <- read_env(env, site_col, env_var_col, env_value_col, env_format,
                         exclude = c(exclude, x_col, y_col, date_col))
        pts  <- read_sites(sites, site_col, x_col, y_col, crs)
        if (!is.null(dates)) {
                dates <- as.data.table(drop_geometry(read_any(dates)))
                check_cols(dates, c(site_col, date_col), "dates")
                dates <- data.table(site = as.character(dates[[site_col]]), date = dates[[date_col]])
        } else dates <- comm$dates
        say("  community: %d sites x %d taxa | env: %d sites x %d variables | coordinates: %s | dates: %s",
            nrow(comm$Y), ncol(comm$Y), nrow(Eall), ncol(Eall),
            if (is.null(pts)) "none" else paste(nrow(pts), "sites"),
            if (is.null(dates)) "none" else paste(uniqueN(dates$site), "sites"))
        
        # Missing environmental values: median of the variable (training data were imputed)
        na_env <- colSums(is.na(Eall))
        if (any(na_env > 0)) {
                say("  env: %d missing values in %d variables, replaced by the variable median",
                    sum(na_env), sum(na_env > 0))
                for (j in which(na_env > 0)) Eall[is.na(Eall[, j]), j] <- median(Eall[, j], na.rm = TRUE)
        }
        
        # Several typologies -> one row each
        typologies <- if (is.list(typology) && !is.data.frame(typology)) typology else list(typology = typology)
        if (is.null(names(typologies)) || any(!nzchar(names(typologies)))) {
                names(typologies) <- paste0("typology_", seq_along(typologies))
        }
        tvars_list <- if (is.list(typology_vars)) typology_vars else rep(list(typology_vars), length(typologies))
        if (is.null(names(tvars_list))) names(tvars_list) <- names(typologies)
        
        rows <- lapply(names(typologies), function(nm) {
                say("\nTypology: %s", nm)
                types <- assign_types(typologies[[nm]], rownames(comm$Y), pts, site_col, type_col)
                use <- intersect(names(types)[!is.na(types)], rownames(Eall))
                if (!is.null(pts)) {
                        no_xy <- setdiff(use, pts$site)
                        if (length(no_xy)) say("  %d sites have no coordinates: used without distances and MEMs", length(no_xy))
                }
                dropped <- setdiff(rownames(comm$Y), use)
                if (length(dropped)) {
                        say("  %d of %d community sites dropped (no type or no environmental data)",
                            length(dropped), nrow(comm$Y))
                }
                if (length(use) < 3) stop("Typology ", nm, ": fewer than three sites with type and environment.")
                
                ty <- types[use]
                Y  <- comm$Y[use, , drop = FALSE]
                Y  <- Y[, colSums(Y) >= MIN_OCCURRENCES, drop = FALSE]
                E  <- Eall[use, , drop = FALSE]
                sizes <- table(ty)
                say("  sites: %d | types: %d (sizes %d-%d)", length(use), length(sizes), min(sizes), max(sizes))
                if (any(sizes < 3)) say("  NOTE: %d types have fewer than 3 sites; training typologies had at least 3", sum(sizes < 3))
                
                res <- list(samples = length(use))
                res <- c(res, typology_descriptors(E, ty))
                
                # Typology variables
                tv <- tvars_list[[nm]]
                Et <- NULL
                if (!is.null(tv)) {
                        miss <- setdiff(tv, colnames(Eall))
                        if (length(miss)) stop("Typology ", nm, ": typology_vars not in env: ", paste(miss, collapse = ", "))
                        Et <- E[, tv, drop = FALSE]
                        res$n_vars <- length(tv)
                } else {
                        res$n_vars <- NA_real_
                        say("  no typology_vars given: n_vars and the variation partitioning are left empty")
                }
                
                # Space
                xy <- NULL
                if (!is.null(pts) && all(use %in% pts$site)) {
                        p  <- pts[match(use, pts$site), ]
                        res <- c(res, site_distances(p))
                        xy <- sf::st_coordinates(p) / 1000
                } else {
                        res$min_distance <- res$median_distance <- NA_real_
                }
                
                # Taxa
                tx <- taxonomic_shares(colnames(Y)[colSums(Y) > 0], taxon_ranks)
                res$species_rank <- tx$shares[["species"]]; res$genus_rank <- tx$shares[["genus"]]
                res$family_rank  <- tx$shares[["family"]];  res$higher_rank <- tx$shares[["higher"]]
                if (tx$guessed) say("  taxonomic ranks guessed from names for %d taxa (check with taxon_ranks)", tx$guessed)
                
                vp <- variation_partitioning(Y, Et, xy)
                res$n_taxa_used <- vp$n_taxa_used
                res$vp_env_pure <- vp$vp_env_pure; res$vp_shared <- vp$vp_shared; res$vp_space_pure <- vp$vp_space_pure
                if (!is.na(vp$n_mem)) say("  variation partitioning: %d MEMs selected", vp$n_mem)
                
                res <- c(res, sampling_season(if (is.null(dates)) NULL else dates[site %in% use]))
                
                row <- as.data.table(c(list(typology = nm), res[PREDICTORS]))
                n_na <- sum(is.na(unlist(res[PREDICTORS])))
                say("  %d of %d predictors computed%s", length(PREDICTORS) - n_na, length(PREDICTORS),
                    if (n_na) paste0(" (empty, imputed by the app: ", paste(PREDICTORS[is.na(unlist(res[PREDICTORS]))], collapse = ", "), ")") else "")
                
                # Observed type coherence (same sites and taxa as the predictors)
                coh <- as.data.table(c(list(typology = nm), coherence_metrics(Y, ty)))
                list(pred = row, coh = coh)
        })
        
        out <- rbindlist(lapply(rows, `[[`, "pred"))
        coh <- rbindlist(lapply(rows, `[[`, "coh"))
        if (!is.null(output)) {
                fwrite(out, output, na = "")
                say("\nWritten: %s", output)
        }
        if (!is.null(coherence_output)) {
                fwrite(coh, coherence_output, na = "")
                say("Written: %s", coherence_output)
        }
        # One batch file for the app: a row per typology x metric, the metric in a column
        batch <- NULL
        if (!is.null(batch_output)) {
                batch <- rbindlist(lapply(METRICS, function(m) {
                        b <- data.table(typology = out$typology, metric = m, observed_value = coh[[m]])
                        if (!is.null(taxon)) b[, taxon := taxon]
                        cbind(b, out[, -"typology"])
                }))
                batch <- batch[order(match(typology, out$typology))]    # input order of typologies, metrics within
                fwrite(batch, batch_output, na = "")
                say("Written: %s (%d typologies x %d metrics)", batch_output, nrow(out), length(METRICS))
        }
        invisible(list(predictors = out, coherence = coh, batch = batch))
}


# ==============================================================================*
# 5. RUN ----  (only when the file is run with Rscript, not when sourced)
# ==============================================================================*

# if (sys.nframe() == 0L) {
#         pulse_predictors(
#                 community     = "data/community.csv",      # wide or long
#                 env           = "data/environment.csv",    # wide or long
#                 typology      = "data/typology.gpkg",      # polygons/lines with column `type`, or a site-type table
#                 sites         = "data/sites.csv",          # site, x, y (longitude, latitude)
#                 typology_vars = c("altitude", "catchment_area", "geology_calcareous",
#                                   "geology_siliceous", "geology_sedimentary"),
#                 taxon_ranks   = NULL,                      # optional table: taxon, rank
#                 site_col      = "site",
#                 taxon_col     = "taxon",                   # long community data only
#                 value_col     = NULL,                      # abundance column (long data); NULL = presence rows
#                 date_col      = "date",                    # sampling dates in the community data, if any
#                 x_col         = "x", y_col = "y", crs = 4326,
#                 type_col      = "type",
#                 output           = "pulse_predictors.csv",
#                 coherence_output = "pulse_coherence.csv",
#                 batch_output     = "pulse_batch.csv",      # optional: all metrics in one file for the app's batch tab
#                 taxon            = NULL                    # optional: taxon of the batch rows, e.g. "macroinvertebrates";
#                 # NULL = the taxon selected in the app
#         )
# }
# ==============================================================================*
# 5. RUN ----  (only when the file is run with Rscript, not when sourced)
# ==============================================================================*


        pulse_predictors(
                community = "data/application_example/brt_fish_community.csv", # wide or long
                env = "data/application_example/brt_fish_env.csv", # wide or long
                typology = "E://Arbeit/data/typology_systems/broad_river_types/rawData/m_river_fec_broad_type.shp", # polygons/lines with column `type`, or a site-type table
                sites = "data/application_example/brt_fish_sites.csv", # site, x, y (longitude, latitude)
                typology_vars = c(
                        "elevation",
                        "upstream_catchment_area",
                        "area_calcareous",
                        "area_siliceous",
                        "area_sediment"
                ),
                taxon_ranks = NULL, # optional table: taxon, rank
                site_col = "site",
                taxon_col = "taxon", # long community data only
                value_col = "value", # abundance column (long data); NULL = presence rows
                date_col = "date", # sampling dates in the community data, if any
                x_col = "x", y_col = "y", crs = 3035,
                type_col = "m_btype12",
                output = "data/application_example/pulse_predictors.csv",
                coherence_output = "data/application_example/pulse_coherence.csv",
                batch_output = "data/application_example/pulse_batch.csv" # optional: one batch file per metric; NULL = none
        )


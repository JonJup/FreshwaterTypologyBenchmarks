### - Sampling time of a subset on the annual cycle - ###


library(data.table)
library(lubridate)


sampling_time <- function(samples, id = "scheme_id", date = "date", site = NULL) {
        d <- data.table(id   = samples[[id]],
                        date = as_date(samples[[date]]),
                        site = if (is.null(site)) seq_len(nrow(samples)) else samples[[site]])
        missing <- d[, .(n_missing = sum(is.na(date))), by = id]
        d <- d[!is.na(date)]
        d[, theta := 2 * pi * (yday(date) - 0.5) / fifelse(leap_year(date), 366, 365)]
        
        by_site <- d[, .(s_sin = mean(sin(theta)), s_cos = mean(cos(theta)), n_dates = .N),
                     by = .(id, site)]
        out <- by_site[, .(time_sin = mean(s_sin), time_cos = mean(s_cos),
                           n_sites = .N, n_dates = sum(n_dates)), by = id]
        
        out[, concentration := sqrt(time_sin^2 + time_cos^2)]
        out[, mean_doy := fifelse(concentration > 1e-8,
                                  (atan2(time_sin, time_cos) %% (2 * pi)) / (2 * pi) * 365.25 + 0.5,
                                  NA_real_)]
        out[, circ_sd_days := sqrt(-2 * log(pmin(concentration, 1))) * 365.25 / (2 * pi)]
        out[, mean_month := month(as_date("2001-01-01") + round(mean_doy) - 1)]
        
        out <- merge(out, missing, by = "id", all = TRUE)        # subsets without any date stay, as NA
        setnames(out, "id", id)
        out[]
}

# rebuild_env_for_coverage_check
# 1. Setup ----------------------------------------------------------------------


suppressPackageStartupMessages({
        library(data.table)
        library(lubridate)
        library(sf)
        library(Hmsc)
        library(adespatial)
        library(jsonify)
        library(mice)
        library(magrittr)
})


# custom function ---------------------------------------------------------
prep_data <- function(x) {return(x[sample_type == "full"])}


# Each sampling date becomes an angle on the annual cycle (day of year, leap
# years respected). Per subset, the predictors are the mean sine and cosine:
#   time_sin, time_cos   -> QRF predictors
# Their mean vector also gives, for reporting only:
#   mean_doy        direction = typical sampling day (circular mean)
#   concentration   length    = 1 if all samples share one day, -> 0 if spread over the year
#   circ_sd_days    circular SD in days, sqrt(-2 log(concentration))
#   mean_month      month of mean_doy (labels, season comparison)
# With site given, angles are averaged within site first, so repeatedly visited
# sites count once (like in the typology evaluation).

sampling_time <- function(samples, ID =  o.scheme.number, date = "eventDate", site = "siteID") {
        d <- data.table(id   = ID,
                        date = as_date(samples[[date]]),
                        site = samples[[site]])
        missing <- d[, .(n_missing = sum(is.na(date))), by = id]
        d <- d[!is.na(date)]
        d <- unique(d, by = c("date", "site"))
        d[, theta := 2 * pi * (yday(date) - 0.5) / fifelse(leap_year(date), 366, 365)]
        
        by_site <- d[, .(
                s_sin = mean(sin(theta)), 
                s_cos = mean(cos(theta)),
                n_dates = .N),
                by = .(id, site)]
        out <- by_site[, .(time_sin = mean(s_sin), time_cos = mean(s_cos),
                           n_sites = .N, n_dates = sum(n_dates)), by = id]
        
        out[, concentration := sqrt(time_sin^2 + time_cos^2)]
        out[, mean_doy := fifelse(concentration > 1e-8,
                                  (atan2(time_sin, time_cos) %% (2 * pi)) / (2 * pi) * 365.25 + 0.5,
                                  NA_real_)]
        out[, circ_sd_days := sqrt(-2 * log(pmin(concentration, 1))) * 365.25 / (2 * pi)]
        out[, mean_month := month(as_date("2001-01-01") + round(mean_doy) - 1)]
        
        out <- merge(out, missing, by = "id", all = TRUE)        # subsets without any date stay, as NA
        setnames(out, "id", id)
        out[]
}



# load data ---------------------------------------------------------------
schemesD <- readRDS("data/000_biota/03_diatoms_scheme.rds")
schemesF <- readRDS("data/000_biota/03_fish_scheme.rds")
schemesI <- readRDS("data/000_biota/03_invertebrates_scheme.rds")
schemesM <- readRDS("data/000_biota/03_macrophytes_scheme.rds")

vpD %<>% sprintf("%04d", .) %>% paste0("diatoms_", .)
vpF %<>% sprintf("%04d", .) %>% paste0("fish_", .)
vpI %<>% sprintf("%04d", .) %>% paste0("invertebrates_", .)
vpM %<>% sprintf("%04d", .) %>% paste0("macrophytes_", .)

# vpD <- list.files(paste0(path_dir, "pulseD/data/005_variation_partitioning/")) |> gsub(".rds", "", x = _)
# vpF <- list.files(paste0(path_dir, "pulseF/data/005_variation_partitioning/")) |> gsub(".rds", "", x = _)
# vpI <- list.files(paste0(path_dir, "pulseI/data/005_variation_partitioning/")) |> gsub(".rds", "", x = _)
# vpM <- list.files(paste0(path_dir, "pulseM/data/005_variation_partitioning/")) |> gsub(".rds", "", x = _)


Di <- readRDS("data/000_biota/02_diatoms_w_environment.rds")
Fi <- readRDS("data/000_biota/02_fish_w_environment.rds")
In <- readRDS("data/000_biota/02_invertebrates_w_environment.rds")
Ma <- readRDS("data/000_biota/02_macrophytes_w_environment.rds")

biota <- list(
        Di,
        Fi,
        In,
        Ma
)

scheme <- list(
        schemesD,
        schemesF,
        schemesI,
        schemesM
)

vp <- 
        list(
                vpD,
                vpF,
                vpI,
                vpM 
        )

out <- 
        list(
                diatoms = vector(mode = "list"),
                fish    = vector(mode = "list"),
                invertebrates = vector(mode = "list"),
                macrophytes   = vector(mode = "list")
        )


for (i in seq_along(scheme)){
        cat("taxon:" ,i, "\n")
        ss <- scheme[[i]]
        b.bio <- biota[[i]]
        i.vp <- vp[[i]]
        
        for (o in seq_len(nrow(ss))){
                cat("scheme:" ,o, "\r")
                o.scheme <- ss[o, ]
                
                #if (o.scheme$sample_type != "full") next
                if (!o.scheme$scheme_id %in% i.vp) next
                
                o.scheme.number <- o.scheme$scheme_id
                
                set.seed(o)
                #cat(paste("Processing scheme", o, "of", nrow(b.scheme), "-", bio.names, "\n"))
                # Filter by dataset, year, and focal months
                o.data <- b.bio[data.set == o.scheme$data.set]
                
                # If focal_months is not NA/Empty, filter by it
                # If month and year are NA just use all the data
                focal_months <- as.numeric(o.scheme$focal_months[[1]])
                if (!all(is.na(o.data$eventDate)) && length(focal_months) > 0 && !all(is.na(focal_months))) {
                        o.data <- o.data[eventYear == o.scheme$eventYear & month(eventDate) %in% focal_months]
                } else if (!all(is.na(o.data$eventYear))){
                        o.data <- o.data[eventYear == o.scheme$eventYear]
                } 
                
                out[[i]][[length(out[[i]]) + 1]] <- sampling_time(o.data)
        }
        
        
}
out2 <- list()
out2[[1]] <- rbindlist(out[[1]])
out2[[2]] <- rbindlist(out[[2]])
out2[[3]] <- rbindlist(out[[3]])
out2[[4]] <- rbindlist(out[[4]])

names(out2) <- names(out)
out2$diatoms

lapply(1:4, function(x){
        saveRDS(out2[[x]], file = paste0("data/", names(out)[x],"_date_variables.rds"))
})

# saveRDS(object = out2, file = "data/date_variables.rds")      

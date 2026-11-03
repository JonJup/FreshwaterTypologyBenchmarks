library(data.table)
library(flextable)
library(officer)

setwd("//turso-fs.cs.helsinki.fi/turso/home/jonjupke/")


files <- vector(mode="list", length = 4)
tax <- c("D", "F", "I", "M")
for (i in tax){
        
        root <- paste0("pulse",i,"/data/008_qrf/")
        files[[which(tax == i)]] <- list.files(root, 
                                               pattern = "_tuning\\.rds$",
                                               recursive = TRUE, 
                                               full.names = TRUE)
}


files <- unlist(files)

tab_s8 <- rbindlist(lapply(files, function(f) {
        x <- readRDS(f); s <- as.data.table(x$summary)          # already sorted by CRPS
        sel   <- s[1]
        other <- s[bagging != sel$bagging][1]                   # best candidate of the other unit
        stopifnot(sel$mtry == x$best$mtry, sel$min_n == x$best$min_n)
        data.table(run              = sub("_tuning\\.rds$", "", basename(f)),
                   mtry             = sel$mtry, 
                   min_n            = sel$min_n, bagging = sel$bagging,
                   crps_cv          = sel$crps, 
                   crps_se          = sel$crps_se,
                   delta_other_unit = other$crps / sel$crps - 1,
                   n_candidates     = nrow(s)
                   )
}))

tab_s8[, c("taxon", "metric") := tstrsplit(run, "__", fixed = TRUE)]
tab_s8[, run := NULL]
tab_s8 <- tab_s8[, c(8,9,1,2,3)]
tab_s8[, metric := gsub("_", "\\ ", metric)]

ft <- flextable(tab_s8) |>
        autofit() |>
        theme_booktabs()


# save to file  -----------------------------------------------------------

setwd(rstudioapi::getActiveProject())
save_as_docx(ft, path = "data/supplement/S7_QRF_TUNING_table.docx")


gelman_check <- function(posterior, skip = F){

        if (!skip) {
                o.ge    <- gelman.diag(x = posterior$Beta, multivariate = FALSE)
        } else {
                o.ge <- posterior
        }
        dt <- data.table(rowname = rownames(o.ge$psrf))
        dt[, taxon :=  regmatches(rowname, regexpr(",\\s.*\\(S", rowname))]
        dt[, taxon :=   gsub("^,\\s*", "", taxon)]
        dt[, taxon :=   gsub("\\ \\(S", "", taxon)]
        dt[, variable := gsub("^B\\[(.+?) \\(C[0-9]+\\).*", "\\1",  rowname)]
        dt[, rowname := NULL]
        dt[, psrf := o.ge[[1]][,1]]
        dt[, taxon_mean := mean(psrf), by = "taxon"]
        dt[, var_mean := mean(psrf), by = "variable"]
        dt <- dt[! grep("MEM", variable)]
        dt.tax <- unique(dt, by = "taxon")
        dt.var <- unique(dt, by = "variable")
        dt <- list(taxon = dt.tax, variable = dt.var)
}
        
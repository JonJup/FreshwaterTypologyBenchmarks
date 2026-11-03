
library(data.table)
library(flextable)

# co <- fread(file.path("data/QRF_RESULTS/tables", "dataset_concentration.csv"))

co <- fread(file.path("data/QRF_RESULTS/tables", "dataset_concentration.csv"))[data_sets %in% c("all, before cap", "all")]

# One value per taxon: check before collapsing over metrics
stopifnot(co[, uniqueN(n_eff_rows), by = .(taxon, data_sets)][, all(V1 == 1)])

b <- unique(co[data_sets == "all, before cap",
               .(taxon, n_datasets, neff_b = n_eff_rows, top3_b = top3_share, max_b = largest_share)])
a <- unique(co[data_sets == "all",
               .(taxon, neff_a = n_eff_rows, top3_a = top3_share, max_a = largest_share)])
tab <- b[a, on = "taxon"][, .(
        Taxon                = taxon,
        `Data sets`          = n_datasets,
        `Effective number`   = sprintf("%.1f → %.1f", neff_b, neff_a),
        `Top 3 (%)`          = sprintf("%.1f → %.1f", 100 * top3_b, 100 * top3_a),
        `Largest (%)`        = sprintf("%.1f → %.1f", 100 * max_b, 100 * max_a))]

ft <- flextable(head(tab)) |>
        autofit() |>
        theme_booktabs()

# Save straight to a .docx
save_as_docx(ft, path = "data/supplement/S7_DataSet_effectiveSize_table.docx")


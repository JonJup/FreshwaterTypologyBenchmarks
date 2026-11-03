### - Table: predictors of the quantile random forests - ###
#
# Builds a description table of the QRF predictors (names, classes and figure
# labels as in qrf_summary.R) and writes it to a Word document.
# Edit the `Description` and `Unit` entries below as needed; the order of the
# rows within a class is the order in which they are listed here.

suppressPackageStartupMessages({
        library(data.table)
        library(flextable)
        library(officer)
})

out_file <- "data/supplement/S7_qrf_predictor_table.docx"

#===============================================================================*
# 1. Table content ----
#===============================================================================*
pred_tab <- rbindlist(list(
        ## Sampling effort
        list("Sampling effort", "samples", "# Samples",
             "Number of samples in the data set.",
             "count"),
        ## Typology structure
        list("Typology structure", "n_types", "# Types",
             "Number of types of the typology scheme.",
             "count"),
        list("Typology structure", "n_vars", "# Variables",
             "Number of environmental variables used to delineate the types.",
             "count"),
        list("Typology structure", "type_evenness", "Type Evenness",
             "Evenness of the distribution of samples among types (1 = all types equally frequent).",
             "0\u20131"),
        ## Environmental separation
        list("Environmental separation", "env_r2_after", "Environmental R2",
             "PERMANOVA R\u00b2 of the environmental variables explained by type membership.",
             "0\u20131"),
        list("Environmental separation", "env_ch_after", "Calinski\u2013Harabasz",
             "Calinski\u2013Harabasz index of the types in environmental space (ratio of between- to within-type dispersion).",
             "\u2265 0"),
        list("Environmental separation", "asw_after", "Silhouette Width",
             "Average silhouette width of the types in environmental space.",
             "\u22121 to 1"),
        ## Spatial extent
        list("Spatial extent", "min_distance", "Spatial Distance (Minimum)",
             "Minimum geographic distance between samples of the data set.",
             "km"),
        list("Spatial extent", "median_distance", "Spatial Distance (Median)",
             "Median pairwise geographic distance between samples of the data set.",
             "km"),
        ## Taxonomic resolution
        list("Taxonomic resolution", "species_rank", "Species Rank (%)",
             "Percentage of taxa identified to species level.",
             "%"),
        list("Taxonomic resolution", "genus_rank", "Genus Rank (%)",
             "Percentage of taxa identified to genus level.",
             "%"),
        list("Taxonomic resolution", "family_rank", "Family Rank (%)",
             "Percentage of taxa identified to family level.",
             "%"),
        list("Taxonomic resolution", "higher_rank", "Higher Rank (%)",
             "Percentage of taxa identified above family level.",
             "%"),
        ## Taxon richness
        list("Taxon richness", "n_taxa_used", "# Taxa",
             "Number of taxa in the community data used for the evaluation.",
             "count"),
        ## Variation partitioning
        list("Variation partitioning", "vp_env_pure", "Environment fraction",
             "Fraction of community variation explained purely by the environment (variation partitioning).",
             "0\u20131"),
        list("Variation partitioning", "vp_space_pure", "Space fraction",
             "Fraction of community variation explained purely by space.",
             "0\u20131"),
        list("Variation partitioning", "vp_shared", "Shared",
             "Fraction of community variation jointly explained by environment and space.",
             "0\u20131"),
        ## Sampling season
        list("Sampling season", "mean_doy", "Mean sampling date",
             "Mean day of the year on which the samples were taken.",
             "day of year"),
        list("Sampling season", "circ_sd_days", "Seasonal spread",
             "Circular standard deviation of the sampling dates.",
             "days")
))
setnames(pred_tab, c("Class", "Predictor", "Label", "Description", "Unit"))

# Keep classes in the order listed above
pred_tab[, Class := factor(Class, levels = unique(Class))]
setorder(pred_tab, Class)
pred_tab[, Class := as.character(Class)]

#===============================================================================*
# 2. Flextable ----
#===============================================================================*
font_name <- "Arial"
font_size <- 9

ft <- flextable(pred_tab) |>
        set_header_labels(Class       = "Predictor class",
                          Predictor   = "Variable",
                          Label       = "Label in figures",
                          Description = "Description",
                          Unit        = "Unit / range") |>
        merge_v(j = "Class") |>
        valign(j = "Class", valign = "top", part = "body") |>
        theme_booktabs() |>
        font(fontname = font_name, part = "all") |>
        fontsize(size = font_size, part = "all") |>
        bold(part = "header") |>
        italic(j = "Predictor", part = "body") |>      # variable names as in the code
        align(j = "Unit", align = "center", part = "all") |>
        padding(padding.top = 1.5, padding.bottom = 1.5, part = "all") |>
        width(j = c("Class", "Predictor", "Label", "Description", "Unit"),
              width = c(1.5, 1.2, 1.7, 3.9, 0.9)) |>
        # thin line between predictor classes
        hline(i = which(!duplicated(pred_tab$Class, fromLast = TRUE))[-uniqueN(pred_tab$Class)],
              border = fp_border(color = "grey60", width = 0.5), part = "body") |>
        set_caption(caption = as_paragraph(
                as_chunk("Table S1. ", props = fp_text(bold = TRUE, font.family = font_name,
                                                       font.size = font_size + 1)),
                as_chunk("Predictors of the quantile random forests.",
                         props = fp_text(font.family = font_name, font.size = font_size + 1))
        ), word_stylename = "Table Caption") |>
        add_footer_lines(paste0(
                "Each row of the training data is one typology scheme evaluated on one data set; ",
                "predictors describe the typology scheme, the environmental separation of its types, ",
                "and the biological data set it is evaluated with. ",
                "Predictor classes are those used to aggregate importance in the results.")) |>
        font(fontname = font_name, part = "footer") |>
        fontsize(size = font_size - 1, part = "footer")

#===============================================================================*
# 3. Word document ----
#===============================================================================*
landscape <- prop_section(page_size = page_size(orient = "landscape"),
                          page_margins = page_mar(top = 0.6, bottom = 0.6, left = 0.8, right = 0.8))

read_docx() |>
        body_add_flextable(ft, align = "left") |>
        body_set_default_section(landscape) |>
        print(target = out_file)

cat("saved:", normalizePath(out_file), "\n")
pac
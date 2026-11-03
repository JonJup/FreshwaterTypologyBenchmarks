
# setup -------------------------------------------------------------------


library(data.table)
library(dplyr)
library(officer)
library(flextable)


# load data ---------------------------------------------------------------
files      <- 
        list.files("E://Arbeit/data/biota/MIDFIRE/data/",
                   pattern = "\\.rds$",
                   full.names = T)
stopifnot(length(files) == 4)
bio.list <- lapply(files, readRDS)
bio      <- rbindlist(bio.list, use.name = TRUE)


rank_levels <- unique(bio$taxonRank)

# collapse everything coarser than order into "higher"
bio[, rank_cat := fifelse(taxonRank %in% rank_levels[1:4], taxonRank, "higher")]
bio[, rank_cat := factor(rank_cat, levels = rank_levels)]

# long format: one row per dataset x rank, fraction of that dataset's observations
frac_long <- bio[, .N, by = .(data.set, rank_cat)][CJ(data.set = unique(bio$data.set), 
                                                      rank_cat = factor(rank_levels, levels = rank_levels)),
                                                   on = .(data.set, rank_cat)][is.na(N), N := 0L][, frac := N / sum(N), by = data.set]

# wide format: one row per dataset, one column per rank
frac_wide <- dcast(frac_long, data.set ~ rank_cat, value.var = "frac")
names(frac_wide)[2] <- "higher_rank"
frac_wide <- frac_wide[, c(1,3,4,5,2)] 
frac_wide[, species := round(species*100,0)] 
frac_wide[, genus := round(genus*100,0)] 
frac_wide[, family := round(family*100,0)] 
frac_wide[, higher_rank := round(higher_rank*100,0)] 
frac_wide[, 'taxonomic_resolution' := paste(species, genus, family, higher_rank, sep = "/")]


bio2 <- bio |>
        mutate(eventMonth = lubridate::month(eventDate)) |> 
        group_by(data.set) |> 
        summarize(
                sampling_start  = min(eventYear, na.rm = T),
                sampling_end    = max(eventYear, na.rm = T),
                Months_min          = min(sort(unique(eventMonth)), na.rm = T),
                Months_max          = max(sort(unique(eventMonth)), na.rm = T),
                'number of sites' = uniqueN(siteID),
                'number of samples' = uniqueN(eventID),
                protocol = unique(methodName),
                protocol2 = unique(methodReference),
                citation = unique(resourceReference)
                  )

bio3 <- bio2 %>%
        mutate(
                `sampling years` = if_else(
                        sampling_start == sampling_end,
                        as.character(sampling_start),
                        paste(sampling_start, sampling_end, sep = " - ")
                ),
                `sampling month` = if_else(
                        Months_min == Months_max,
                        as.character(Months_min),
                        paste(Months_min, Months_max, sep = " - ")
                )
        ) |> 
        select(data.set, 'sampling years', 'sampling month', 'number of sites', 'number of samples', protocol, protocol2, citation) |> 
        left_join(
                select(frac_wide, data.set, taxonomic_resolution),
                by = "data.set"
                ) |> 
        mutate(
                taxon := gsub("_.*", "", data.set)
        ) |> 
        mutate(
                data.set = gsub("diatoms_|fish_|invertebrates_|macrophytes_", "", data.set)
        ) |> 
        mutate(
                country := gsub("_.*", "", data.set)
        ) |>
        mutate(
                data.set = gsub("^.*_", "", data.set)
        )

setDT(bio3)

bio3[data.set == "aqem", data.set := "Project AQEM"]
bio3[data.set == "basin", data.set := "WFD monitoring (Po River Basin)"]
bio3[data.set == "BDM", data.set := "Monitoring BDM"]
bio3[data.set == "ecosurv",    data.set := "Project Ecosurv"]
bio3[data.set == "STAR",       data.set := "Project STAR"]
bio3[data.set == "STARS",       data.set := "Project STAR"]
bio3[data.set == "WISER",       data.set := "Project WISER"]
bio3[data.set == "drypond",    data.set := "Project Drypond"]
bio3[data.set == "biodrought", data.set := "Project Biodrought"]
bio3[data.set == "dryfish", data.set := "Project Dryfish"]
bio3[data.set == "krknose", data.set := "Project Krknose"]
bio3[data.set == "salamander", data.set := "Project Salamander"]
bio3[data.set == "rcs", data.set := "WFD monitoring (RCS)"]
bio3[data.set == "NAWA", data.set := "Monitoring NAWA"]
bio3[data.set == "nawa", data.set := "Monitoring NAWA"]
bio3[data.set == "naiades", data.set := "WFD monitoring (Naïades)"]
bio3[data.set == "monitoring2", data.set := "WFD monitoring (2)"]
bio3[data.set == "monitoring", data.set := "WFD monitoring"]
bio3[data.set == "intercalibration", data.set := "WFD intercalibration"]
bio3[data.set == "INRAE", data.set := "WFD monitoring (INRAE)"]
bio3[data.set == "ebro", data.set := "WFD monitoring (Ebro Basin)"]
bio3[data.set == "duero", data.set := "WFD monitoring (Duero Basin)"]
bio3[data.set == "BMST", data.set := "Biodiversity Monitoring South Tyrol"]
bio3[data.set == "chmi", data.set := "WFD monitoring"]
bio3[data.set == "cimpean", data.set := "study01"]




bio3[,country := stringr::str_to_title(country)]
bio3[,taxon := stringr::str_to_title(taxon)]

names(bio3)
names(bio3)[1] <- "Data Set"
names(bio3)[2] <- "Years"
names(bio3)[3] <- "Months"
names(bio3)[4] <- "No. Sites"
names(bio3)[5] <- "No. Samples"
names(bio3)[6] <- "Protocol"
names(bio3)[7] <- "Method Reference"
names(bio3)[8] <- "Reference"
names(bio3)[9] <- "Resolution"
names(bio3)[10] <- "Taxon"
names(bio3)[11] <- "Country"

bio3 <- bio3[, c(10,11,1,2,3,4,5,8,6,7)] 


bio3[Months == "Inf - -Inf", Months := NA]
bio3[Years == "Inf - -Inf", Years := NA]
bio3[Country == "Czech", Country := "Czech Republic"]

fwrite(bio3, "data/supplement/SM2_Bio_Table.csv")


# ft <- flextable(bio3) |>
#         theme_booktabs() |>
#         fontsize(size = 7, part = "all") |>
#         padding(padding = 1, part = "all") |>
#         set_table_properties(layout = "autofit", width = 1) |>  # width = 1 -> 100% of page width
#         autofit()
# 
# # Landscape section with narrow margins
# landscape <- prop_section(
#         page_size = page_size(orient = "landscape"),
#         page_margins = page_mar(left = 0.5, right = 0.5, top = 0.5, bottom = 0.5)
# )
# 
# save_as_docx(
#         ft,
#         path = "data/supplement/S1_DataSet_Bio_table.docx",
#         pr_section = landscape
# )

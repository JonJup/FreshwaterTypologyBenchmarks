## Data for SM: CHANGES INCURRED THROUGH FILTERING

# setup -------------------------------------------------------------------
library(data.table)
library(tidyverse)
library(magrittr)
library(sf)
library(maptiles)
library(tmap)
library(scales)

# working directory
setwd(rstudioapi::getActiveProject())

# custom function 
source("code/functions/find_max_consequtive_sum.R")

# MIDFIRE data ------------------------------------------------------------
d <- fread("E://Arbeit/data/biota/MIDFIRE/data/csv/midfire_diatoms.csv")
f <- fread("E://Arbeit/data/biota/MIDFIRE/data/csv/midfire_fish.csv")
i <- readRDS("E://Arbeit/data/biota/invertebrate data/00_combine_data/midfire_invertebrates.rds")
m <- fread("E://Arbeit/data/biota/MIDFIRE/data/csv/midfire_macrophytes.csv")

# Combine midfire data
d[, taxon := "diatoms"]
f[, taxon := "fish"]
i[, taxon := "invertebrates"]
m[, taxon := "macrophytes"]


# combine all four taxa into one data set 
data <- rbindlist(list(d,f,i,m), use.names=T)
# one row per sample 
data2 <- unique(data, by = c("taxon", "eventID"))

unique(data2$data.set)
data2[!is.na(eventYear), uniqueN(eventYear)]

# original number of taxa
m[!is.na(species), uniqueN(species)]
m[!is.na(genus)  , uniqueN(genus)]
m[!is.na(family) , uniqueN(family)]
m[!is.na(order)  , uniqueN(order)]


# Now we go through the filtering process 
n.data.sets      <- uniqueN(data2$data.set)
ud               <- unique(data2$data.set)

min_samples <- 100

# list to store results of loop 
result_list <- vector(mode = "list")

## 2.2 run Loop ----
for (d in 1:length(ud)) {
        print(d)
        ds.result.list <- list()
        #- select data set
        ds.data.set <- ud[d]
        ds.bio      <- data2[data.set == ds.data.set]
        #- how many years?
        if (class(ds.bio$eventYear) == "list") {
                ds.bio[, eventYear := unlist(eventYear)]
        }
        ds.n.year <- uniqueN(ds.bio$eventYear)
        #- how many samples per year?
        ds.n.samples <- ds.bio[, uniqueN(eventID), by = "eventYear"]
        
        if (all(ds.n.samples$V1 < min_samples)) {
                rm(list = ls()[grepl(pattern = "^ds\\.", x = ls())])
                next()
        }
        #- select years with more than 100 samples
        ds.n.samples <- ds.n.samples[V1 >= min_samples]
        #- remove NA row, if one of the entries in the year table is NA
        if (!all(is.na(ds.n.samples$eventYear)) &
            any(is.na(ds.n.samples$eventYear))) {
                ds.n.samples <- ds.n.samples[-which(is.na(ds.n.samples$eventYear))]
        }
        
        # START LOOP i OVER years in data set
        for (i in 1:nrow(ds.n.samples)) {
                # ONLY NA as YEAR
                if (all(is.na(ds.bio$eventYear))) {
                        i.data    <- copy(ds.bio)
                        # ONLY NA as DATE
                } else if (all(is.na(ds.bio$eventDate))) {
                        i.data <- ds.bio[eventYear == ds.n.samples$eventYear[i]]
                } else {
                        #- Create subset of focal year.
                        i.data <- ds.bio[eventYear == ds.n.samples$eventYear[i]]
                        #--- Check seasons. To prevent strong seasonal changes from influencing
                        #--- the community composition, we identify the three consecutive month
                        #--- with the most samples.
                        #- create month variable
                        i.data[, month := month(eventDate)]
                        #- count samples per month
                        i.month_table <- unique(i.data, by = "eventID")
                        i.month_table <- i.month_table$month
                        i.month_table <- table(i.month_table)
                        #- Which three consecutive months have the most samples?
                        #- This function is loaded in the beginning of the script
                        if (length(i.month_table) > 2) {
                                i.max_month    <- find_max_consecutive_sum(i.month_table)
                                i.focal.months <- names(i.max_month$values)
                                
                                #- create subset of i.data only containing the focal months
                                i.data <- i.data[month %in% i.focal.months]
                                
                        } else {
                                #- Are the two months consecutive?
                                i.monthdiff <- diff(as.numeric(names(
                                        i.month_table
                                )))
                                #- No or only one month
                                if (length(i.monthdiff) == 0) {
                                        # Does one of the month have more than 100 samples?
                                        if (any(i.month_table > min_samples)) {
                                                i.focal.months <- names(which.max(
                                                        i.month_table
                                                ))
                                                i.data <- i.data[month == i.focal.months]
                                        } else {
                                                counter[[b]] <- counter[[b]] + 1
                                                rm(list = ls()[grepl(
                                                        pattern = "^i\\.",
                                                        x = ls()
                                                )])
                                                next()
                                        }
                                } else if (i.monthdiff > 2) {
                                        # Does one of the month have more than 100 samples?
                                        if (any(i.month_table > min_samples)) {
                                                i.focal.months <- names(which.max(
                                                        i.month_table
                                                ))
                                                i.data <- i.data[month == i.focal.months]
                                        } else {
                                                counter[[b]] <- counter[[b]] + 1
                                                rm(list = ls()[grepl(
                                                        pattern = "^i\\.",
                                                        x = ls()
                                                )])
                                                next()
                                        }
                                } else {
                                        #- yes
                                        i.focal.months <- names(i.month_table)
                                        i.data <- i.data[month %in% i.focal.months]
                                }
                        }
                }
                
                #- how many samples are left?
                i.samples <- uniqueN(i.data$eventID)
                if (i.samples < min_samples) {
                        rm(list = ls()[grepl("^i\\.", x = ls())])
                        next()
                }
                i.out <- data.table(
                        taxon    = unique(i.data$taxon),
                        data.set = ud[d],
                        eventYear = unique(i.data$eventYear)
                )
                if (all(is.na(ds.bio$eventYear)) |
                    all(is.na(ds.bio$eventDate))) {
                        i.out[, focal_months := list(list(NA))]
                } else {
                        i.out[, focal_months := list(list(i.focal.months))]
                }
                ds.result.list[[length(ds.result.list) + 1]] <- i.out
                rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
        } # END of loop i over years in data.set d
        
        # combine results of all years.
        ds.results.list <- rbindlist(ds.result.list)
        result_list[[length(result_list) + 1]] <- ds.results.list
        rm(list = ls()[grepl(pattern = "^ds\\.", x = ls())])
        
} # END d of loop over data sets in taxonomic group b
result.data <- rbindlist(result_list, fill = T)


# Filter the full data set to the reduced one 
data2[,included := FALSE]
data2[, month := month(eventDate)]
for (i in 1:nrow(result.data)){
        print(i)
        data2[data.set == result.data[i, data.set] & 
                      eventYear == result.data[i, eventYear] & 
                      month %in% unlist(result.data[i, focal_months]),
              included := TRUE
              ]
        
        
}
# only samples which are included  
data3 <- data2[included == TRUE]
# all included data
data4 <- data[eventID %in% data3$eventID]
        
uniqueN(data3$data.set)
data2[!is.na(eventYear), uniqueN(eventYear)]

sort(data3[!is.na(eventYear), unique(eventYear)])

cat("The complete (unfiltered) data set comprises", n.data.sets, "unique data sets with a total of", nrow(data2),"samples, spanning", data2[!is.na(eventYear), uniqueN(eventYear)], "years")

cat("After filtering", uniqueN(data3$data.set), "data sets remain, totaling ", nrow(data3)," samples across ",data3[!is.na(eventYear), uniqueN(eventYear)]," years") 
cat("This corresponds to a reduction of",round(100- uniqueN(data3$data.set)/n.data.sets *100,2),"% of data sest")
cat(" and",round(100- nrow(data3)/nrow(data2) *100,2),"% of individal samples")

# spatial  ----------------------------------------------------------------
d2j <- select(data2, taxon, eventID, included)
data <- left_join(data, d2j, by = c("taxon", "eventID"))
data[, all_missing := sum(included), by = c("eventID", "taxon")]
data[, all_missing_bool := case_when(all_missing == 0 ~ TRUE, TRUE ~FALSE)]
data6 <- data[(all_missing_bool == FALSE & all_missing != 0) | all_missing_bool == TRUE, ]
data.sf <- unique(data6, by = c("taxon", "siteID"))
sites <- st_as_sf(data.sf, coords = c("x.coord", "y.coord"), crs = 3035)

custom.color.palette <- c("#EC6B4F", "#65F78D")
basemap.tile <-
        get_tiles(
                x = st_union(sites),
                provider = "Esri.WorldTerrain",
                zoom = 4, crop = TRUE)

# 1. Refactor the labels
# 2. Re-order so 'No' (Missing) is drawn on top of 'Yes'
sites <- sites %>%
        mutate(all_missing_bool = factor(all_missing_bool, 
                                         levels = c(FALSE, TRUE), 
                                         labels = c("Yes", "No"))) %>%
        arrange(all_missing_bool)
sites_ordered <- sites[order(sites$all_missing_bool, decreasing = T), ]

map.static <-
        tm_shape(basemap.tile) +
        tm_rgb() +
        tm_shape(
                sites_ordered
        ) + 
        tm_dots(
                fill = "all_missing_bool",
                size = 0.3, 
                fill_alpha = 0.2,
                # Define the title and colors here
                fill.scale = tm_scale_categorical(
                        values = c("Yes" = "#377eb8", "No" = "#e41a1c")
                ),
                fill.legend = tm_legend(
                        title = "Are data included in the analysis?", 
                        position = tm_pos_out(cell.h = "center", cell.v = "bottom")
                )
        ) + 
        tm_facets(by = "taxon", free.coords = FALSE, nrow = 2) +
        # This makes the legend symbols larger and moves it outside
        tm_layout(
                legend.title.size = 1.1,
                legend.text.size = 0.9
        )

tmap_save(tm       = map.static, 
          filename = "output/figures/supplement/map.png", 
          dpi      = 300,
          height   = 10,
          width    = 10)


# taxonomic ---------------------------------------------------------------

cat("In the filtered data set, diatoms are represented by ", data4[taxon == "diatoms" & !is.na(species) , uniqueN(species)] ,"species (or species complexes) belonging to ", data4[taxon == "diatoms" & !is.na(genus) , uniqueN(genus)],"genera, ", data4[taxon == "diatoms" & !is.na(family) , uniqueN(family)], "families, and ", data4[taxon == "diatoms" & !is.na(order)  , uniqueN(order)]," orders. ")
cat("Fish are represented by" ,data4[taxon == "fish"          & !is.na(species) , uniqueN(species)], "species in" ,data4[taxon == "fish"          & !is.na(genus) , uniqueN(genus)], "genera," ,data4[taxon == "fish"          & !is.na(family) , uniqueN(family)], "families, and ",data4[taxon == "fish"          & !is.na(order) , uniqueN(order)],"orders. ")
cat("Invertebrates comprise " ,data4[taxon == "invertebrates" & !is.na(species) , uniqueN(species)], "species in" ,data4[taxon == "invertebrates" & !is.na(genus) , uniqueN(genus)], "genera," ,data4[taxon == "invertebrates" & !is.na(family) , uniqueN(family)], "families, and ",data4[taxon == "invertebrates" & !is.na(order) , uniqueN(order)],"orders. ")
cat(   "Macrophytes include " ,data4[taxon == "macrophytes"   & !is.na(species) , uniqueN(species)], "species in" ,data4[taxon == "macrophytes"   & !is.na(genus) , uniqueN(genus)], "genera," ,data4[taxon == "macrophytes"   & !is.na(family) , uniqueN(family)], "families, and ",data4[taxon == "macrophytes"   & !is.na(order) , uniqueN(order)],"orders. ")


# LOST taxa
txn <- "invertebrates"
lost_orders  <- setdiff(data[taxon == txn & !is.na(order), unique(order)], data4[taxon == txn & !is.na(order) , unique(order)])
check_data   <- data[taxon == txn & !order %in% lost_orders]
lost_families <- setdiff(check_data[!is.na(family),  unique(family) ], data4[taxon == txn & !is.na(family) , unique(family)])
check_data   <- check_data[!family %in% lost_families]
lost_genera <- setdiff(check_data[!is.na(genus),   unique(genus)  ], data4[taxon == txn & !is.na(genus) , unique(genus)])
check_data   <- check_data[!genus %in% lost_genera]
lost_species <- setdiff(check_data[!is.na(species), unique(species)], data4[taxon == txn & !is.na(species) , unique(species)])

cat("In", txn, ", we lost", length(lost_species)," species", length(lost_genera),"genera", length(lost_families),"families and", length(lost_orders),"orders")
cat("Together, these Families accounted for",round((nrow(data[taxon == txn & family  %in% lost_families])/nrow(data))*100, 4), "of the data")
cat("Together, these Genera accounted for"  ,round((nrow(data[taxon == txn & genus   %in% lost_genera  ])/nrow(data))*100, 4), "of the data")
cat("Together, these Species accounted for", round((nrow(data[taxon == txn & species %in% lost_species ])/nrow(data))*100, 4), "of the data")
paste(lost_species, collapse = "*, *")

# temporal  ----------------------------------------------------------------

# 1. Create a clear label for the legend
# Setting "Kept" as the first level puts it at the bottom of the stacked bar
data6[, subset_status := factor(ifelse(all_missing_bool, "Discarded", "Kept"), 
                             levels = c("Kept", "Discarded"))]
data6[, month := month(eventDate)]

data6 %>% 
        filter(eventYear > 1990 & eventYear < 2024) %>%
        unique(by = c("eventID", "taxon")) %>% 
        ggplot(aes(x = eventYear, fill = subset_status)) +
        geom_bar(position = "stack") +
        scale_fill_manual(values = c("Kept" = "#2c7bb6", "Discarded" = "#d7191c")) +
        scale_x_continuous(breaks = unique(data6$eventYear)) + # Ensures integer years on x-axis
        labs(title = "Effect of Subsetting on Yearly Distribution",
             x = "Year",
             y = "Number of Data Samples",
             fill = "Status") +
        theme_minimal() +
        theme(legend.position = "top",
              axis.text.x = element_text(angle = -45))

ggsave("output/figures/supplement/histogram_years.tiff", dpi = 600)
## removed fraction 
data6[, kept_fraction := sum(subset_status == "Kept")/.N, by = "eventYear"]
data6[, n_year := uniqueN(eventID), by = "eventYear"]

data7 <- copy(data6)
data7 <- unique(data7, by = "eventYear")
data7[eventYear > 1990 & eventYear < 2024] %>% 
        ggplot(aes(x = eventYear, y = kept_fraction)) + 
        
        # 1. Line aesthetics: Make the line subdued so the points are the star of the show
        geom_line(color = "grey70", linewidth = 0.8) + 
        
        # 2. Point aesthetics: Add a nice pop of color
        geom_point(aes(size = n_year, alpha = n_year), color = "#2c7bb6") +
        
        # 3. Size scale: scale_size_area ensures a value of 0 maps to 0 size, 
        # making the visual weight strictly proportional to your counts.
        scale_size_area(max_size = 10, name = "Total # Samples") + 
        
        # 4. Alpha scale: Restrict the range so small points don't vanish entirely,
        # and remove its legend so it doesn't duplicate the size legend.
        scale_alpha_continuous(range = c(0.4, 0.9), guide = "none") +
        
        # 5. Y-axis: Format as percentages (e.g., 80% instead of 0.8)
        # Consider adding limits = c(0, 1) inside this function if you want the axis anchored to 0-100%
        scale_y_continuous(labels = scales::percent_format()) + 
        
        # 6. Labels: Add clear titles and axis names
        labs(
                title = "Proportion of Data Kept Post-Subsetting",
                subtitle = "Larger, darker points indicate years with higher initial data volume",
                x = "Event Year",
                y = "Fraction Kept"
        ) +
        
        # 7. Theme: A clean, modern theme with a well-placed legend
        theme_minimal() +
        theme(
                legend.position = "bottom",
                plot.title = element_text(face = "bold", size = 14),
                panel.grid.minor = element_blank() # Removes in-between grid lines for less clutter
        )
ggsave("output/figures/supplement/fraction_years_kept.tiff", dpi = 600)

data6 %>%
        filter(eventYear > 1990 & eventYear < 2024) %>%
        filter(!is.na(month)) %>% 
        unique(by = c("eventID", "taxon")) %>%
        ggplot(aes(x = month, fill = subset_status)) +
        geom_bar()+
        scale_fill_manual(values = c("Kept" = "#2c7bb6", "Discarded" = "#d7191c")) +
        scale_x_continuous(breaks = 1:12) + 
        
        labs(title = "Effect of Subsetting on Monthly Distribution",
             x = "Month",
             y = "Number of Samples",
             fill = "Status") +
        theme_minimal() + 
        theme(
                legend.position = "top"
        )
ggsave("output/figures/supplement/histogram_month.tiff", dpi = 600)

data6[, kept_fraction := sum(subset_status == "Kept")/.N, by = "month"]
data6[, n_month := uniqueN(eventID), by = "month"]

data7 <- copy(data6)
data7 <- unique(data7, by = "month")
data7[eventYear > 1990 & eventYear < 2024] %>% 
        ggplot(aes(x = month, y = kept_fraction)) + 
        
        # 1. Line aesthetics: Make the line subdued so the points are the star of the show
        geom_line(color = "grey70", linewidth = 0.8) + 
        
        # 2. Point aesthetics: Add a nice pop of color
        geom_point(aes(size = n_year, alpha = n_year), color = "#2c7bb6") +
        
        # 3. Size scale: scale_size_area ensures a value of 0 maps to 0 size, 
        # making the visual weight strictly proportional to your counts.
        scale_size_area(max_size = 10, name = "Total # Samples") + 
        
        # 4. Alpha scale: Restrict the range so small points don't vanish entirely,
        # and remove its legend so it doesn't duplicate the size legend.
        scale_alpha_continuous(range = c(0.4, 0.9), guide = "none") +
        
        # 5. Y-axis: Format as percentages
        scale_y_continuous(labels = scales::percent_format()) + 
        
        # --- NEW: X-axis: Force breaks at every integer from 1 to 12 ---
        scale_x_continuous(breaks = 1:12) + 
        
        # 6. Labels: Add clear titles and axis names
        labs(
                title = "Proportion of Data Kept Post-Subsetting",
                subtitle = "Larger, darker points indicate years with higher initial data volume",
                x = "Event Month",
                y = "Fraction Kept"
        ) +
        
        # 7. Theme: A clean, modern theme with a well-placed legend
        theme_minimal() +
        theme(
                legend.position = "bottom",
                plot.title = element_text(face = "bold", size = 14),
                panel.grid.minor = element_blank() # Removes in-between grid lines for less clutter
        )
ggsave("output/figures/supplement/fraction_month_kept.tiff", dpi = 600)

## Pressure Gradient Representation ## 
## What parts of the pressure gradient do the midfire data represent? 
## For this we take the sites used in the models that pass all stages, and look
## at their stressor values. 
## 

# library -----------------------------------------------------------------
library(data.table)
library(ggplot2)
library(ggdist)
library(dplyr)
library(tmap)
library(sf)

# custom functions --------------------------------------------------------
get_sites <- function(x){
        x2 <- x[1:1000, ]
        y <- unique(x2, by = "siteID")
}

# Paths ---------------------------------------------------------------
pressures    <- "E://Arbeit/data/my_stressor_db/processedData/PCA_stressors.rds"
bio          <- list.files(path = "data/000_biota", pattern = "^02", full.names = TRUE)

# load data ---------------------------------------------------------------
pressures <- readRDS(pressures)
bio       <- lapply(bio, readRDS)

# prepare data ------------------------------------------------------------
sites <- lapply(bio, function(x) unique(x, by = "siteID"))
pressures_sub <- pressures[, c("ID", "dist_ref")]
sitesXp <- lapply(sites, function(x) pressures_sub[x, on = "ID"])

hist(sitesXp[[1]]$dist_ref)
hist(sitesXp[[2]]$dist_ref)
hist(sitesXp[[3]]$dist_ref)
hist(sitesXp[[4]]$dist_ref)

# build_combined_data_set -------------------------------------------------

dt1 <- data.table(distance = sitesXp[[1]]$dist_ref, group = "diatoms")
dt2 <- data.table(distance = sitesXp[[2]]$dist_ref, group = "fish")
dt3 <- data.table(distance = sitesXp[[3]]$dist_ref, group = "invertebrates")
dt4 <- data.table(distance = sitesXp[[4]]$dist_ref, group = "macrophytes")
dt5 <- data.table(distance = pressures$dist_ref, group = "Europe")

dt <- rbindlist(list(dt1,dt2, dt3,dt4,dt5))

dt$group <- factor(dt$group, levels  = c("macrophytes", "invertebrates", "fish", "diatoms", "Europe"))

LD1 <- quantile(dt5$distance, .01)
LD10 <- quantile(dt5$distance, .1)
LD25 <- quantile(dt5$distance, .25)

# install.packages(c("ggplot2", "ggdist", "dplyr"))

p <- ggplot(dt, aes(x = group, y = distance, fill = group, colour = group)) +
        # Cloud: half-violin (density)
        stat_halfeye(
                adjust = 0.6,
                width = 0.6,
                justification = -0.2,
                .width = 0,
                point_colour = NA,
                alpha = 0.7
        ) +
        # Box: slim boxplot between cloud and rain
        geom_boxplot(
                width = 0.12,
                outlier.shape = NA,
                alpha = 0.4,
                colour = "grey20"
        ) +
        # Rain: jittered raw points
        geom_point(
                position = position_jitter(width = 0.06, height = 0, seed = 1),
                size = 1.6,
                alpha = 0.5,
                shape = 16,
                # shift points to the left of the box
                aes(x = as.numeric(factor(group)) - 0.22)
        ) +
        scale_fill_brewer(palette = "Set2") +
        scale_colour_brewer(palette = "Set2") +
        coord_flip() +   # remove for vertical orientation
        labs(x = NULL, y = "Latent distance to reference condition") +
        theme_minimal(base_size = 13) +
        theme(legend.position = "none",
              panel.grid.minor = element_blank()) + 
        # geom_hline(yintercept = line1) + 
        # geom_hline(yintercept = line10) + 
        geom_hline(yintercept = LD25) 

ggsave("figures/supplement/S2/pressure_raincloud.png")

# Maps  -------------------------------------------------------------------
# How would dropping sites affect distribution patterns 

stadia_key <- Sys.getenv("STADIA_KEY")
tile_url <- paste0("https://tiles.stadiamaps.com/tiles/stamen_terrain_background/",
                   "{z}/{x}/{y}.png?api_key=", stadia_key)

pad        <- 0

pad_bbox <- function(x, frac) {
        bb <- st_bbox(x)
        dx <- (bb["xmax"] - bb["xmin"]) * frac
        dy <- (bb["ymax"] - bb["ymin"]) * frac
        bb + c(-dx, -dy, dx, dy)
}

make_map <- function(x,title){
        x <- st_as_sf(x, coords = c("x.coord", "y.coord"), crs = 3035)
        bb <- pad_bbox(x, pad)
        p <- tm_basemap(server = tile_url) +
                tm_shape(x, bbox = bb) +
                tm_dots(col = "black", lwd = 0.5, size = 0.6, alpha = .1) +
                #tm_title(paste0(focal_data_set, " (", focal_year, ")")) +
                tm_layout(legend.outside = TRUE) + 
                tm_title(title)
        p
}
# tes <- make_map(focal_sites)
taxa = c("diatoms", "fish", "invertebrates", "macrophytes")
map_path = "figures/supplement/S2/"
# Loop over taxa 
for (i in 1:4){
        
        # first the data we would keep if we subset to just 10%LD panEuropean
        focal_sites <- sitesXp[[i]]
        ld10_sites  <- focal_sites[dist_ref <= LD10]
        
        # Now take 10 of data based on what is in the data 
        LD10_focal <- quantile(focal_sites$dist_ref, .1, na.rm = T)
        LD10F_sites <- focal_sites[dist_ref <= LD10_focal]
        
        # Now do this by data set 
        LD10FDS_sites <- focal_sites[, LD10DS := quantile(dist_ref, .1, na.rm = TRUE), by = "data.set"]
        LD10FDS_sites <- LD10FDS_sites[dist_ref <= LD10DS]
        
        m1 <- make_map(focal_sites,   title = paste(taxa[i], "- all sites (", nrow(focal_sites), ")"))     
        m2 <- make_map(ld10_sites,    title = paste(taxa[i], "- 10% Europe(", nrow(ld10_sites), ")"))      
        m3 <- make_map(LD10F_sites,   title = paste(taxa[i], "- 10% Sites(", nrow(LD10F_sites), ")"))         
        m4 <- make_map(LD10FDS_sites, title = paste(taxa[i], "- 10% per data.set(", nrow(LD10FDS_sites), ")"))            
        m5 <- tmap_arrange(m1,m2,m3,m4, ncol = 2)

        file_name <- paste0(map_path, taxa[i], "_stressor_map.png")
        tmap_save(m5, filename = file_name)
        
}



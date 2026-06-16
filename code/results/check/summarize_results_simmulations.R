# ==============================================================================
# setup
# ==============================================================================

library(data.table)
library(ggplot2)
library(tidyverse)

setwd("E://Arbeit/Projekte/02_ongoing/PULSE/wp1/data/")

# ==============================================================================
# How many simulations? 
# ==============================================================================

wd <- c("E://Arbeit/Projekte/02_ongoing/PULSE/wp1/data/pulseD",
        "E://Arbeit/Projekte/02_ongoing/PULSE/wp1/data/pulseF",
        "E://Arbeit/Projekte/02_ongoing/PULSE/wp1/data/pulseI",
        "E://Arbeit/Projekte/02_ongoing/PULSE/wp1/data/pulseM")
vp <- c()
sm <- c()
ms <- list()

for (i in 1:4){
        setwd(wd[i])
        vp[i] <- length(list.files("data/005_variation_partitioning/"))
        sm[i] <- length(list.files("data/006_simulated_data/"))
        ms[[i]] <- setdiff(list.files("data/005_variation_partitioning/"), list.files("data/006_simulated_data/"))
        ms[[i]] <- gsub("\\.rds", "", ms[[i]])
        ms[[i]] <- gsub(".*_", "", ms[[i]])
        ms[[i]] <- as.numeric(ms[[i]])
}

cat("We have a total of", sum(vp), "models")
cat("We have a total of", sum(sm), "simulations")
cat("We lost ", length(ms[[1]]), "diatom schemes")
cat("We lost ", length(ms[[2]]), "fish schemes")
cat("We lost ", length(ms[[3]]), "invertebrates schemes")
cat("We lost ", length(ms[[4]]), "macrophyte schemes")

# Make sure that all schemes are accounted for
sum(sapply(ms, length)) + sum(sm) == sum(vp)
# Simulations failed for: 
round(sum(sapply(ms, length)) / sum(vp) * 100,2)

# ==============================================================================
# Detailed output 
# ==============================================================================
setwd("E://Arbeit/Projekte/02_ongoing/PULSE/wp1/data/")
filesD <- list.files("pulseD/data/misc/simulation_diagnostics/", full.names = T)
filesF <- list.files("pulseF/data/misc/simulation_diagnostics/", full.names = T)
filesI <- list.files("pulseI/data/misc/simulation_diagnostics/", full.names = T)
filesM <- list.files("pulseM/data/misc/simulation_diagnostics/", full.names = T)

dataD <- lapply(filesD, readRDS)
dataF <- lapply(filesF, readRDS)
dataI <- lapply(filesI, readRDS)
dataM <- lapply(filesM, readRDS)

for (i in seq_along(dataD)) dataD[[i]]$model = gsub("pulseD/data/misc/simulation_diagnostics/","",filesD[i])
for (i in seq_along(dataF)) dataF[[i]]$model = gsub("pulseD/data/misc/simulation_diagnostics/","",filesF[i])
for (i in seq_along(dataI)) dataI[[i]]$model = gsub("pulseD/data/misc/simulation_diagnostics/","",filesI[i])
for (i in seq_along(dataM)) dataM[[i]]$model = gsub("pulseD/data/misc/simulation_diagnostics/","",filesM[i])


diag_dtD <- rbindlist(dataD, use.names =T)
diag_dtF <- rbindlist(dataF, use.names =T)
diag_dtI <- rbindlist(dataI, use.names =T)
diag_dtM <- rbindlist(dataM, use.names =T)
diag_dtD$taxon <- "diatoms"
diag_dtF$taxon <- "fish"
diag_dtI$taxon <- "invertebrates"
diag_dtM$taxon <- "macrophytes"

diag_dt <- rbindlist(list(
        diag_dtD,
        diag_dtF,
        diag_dtI,
        diag_dtM  
), use.names =T)

cat("\n========================================\n")
cat("DIAGNOSTIC SUMMARY\n")
cat("========================================\n\n")

cat(sprintf("Total iterations attempted:       %d\n", max(diag_dt$iteration)))
cat(sprintf("Total parameter combos evaluated: %d\n", nrow(diag_dt)))
cat(sprintf("Passed Isolation Forest:          %d (%.1f%%)\n",
            sum(diag_dt$passed_isoforest, na.rm = TRUE),
            100 * mean(diag_dt$passed_isoforest, na.rm = TRUE)))
cat(sprintf("Passed Mahalanobis:               %d (%.1f%%)\n",
            sum(diag_dt$passed_mahalanobis == TRUE, na.rm = TRUE),
            100 * mean(diag_dt$passed_mahalanobis == TRUE, na.rm = TRUE)))
cat(sprintf("Final PASS:                       %d (%.1f%%)\n",
            sum(diag_dt$final_status == "PASS"),
            100 * mean(diag_dt$final_status == "PASS")))

cat("\n--- Status Breakdown ---\n")
print(diag_dt[, .N, by = final_status][order(-N)])

cat("\n--- Separation: Original vs Effective (PASS only) ---\n")
pass_dt <- diag_dt[final_status == "PASS"]
if (nrow(pass_dt) > 0) {
        cat(sprintf("  Original  — mean: %.3f, sd: %.3f, range: [%.3f, %.3f]\n",
                    mean(pass_dt$quality_original), sd(pass_dt$quality_original),
                    min(pass_dt$quality_original), max(pass_dt$quality_original)))
        cat(sprintf("  Effective — mean: %.3f, sd: %.3f, range: [%.3f, %.3f]\n",
                    mean(pass_dt$quality_effective), sd(pass_dt$quality_effective),
                    min(pass_dt$quality_effective), max(pass_dt$quality_effective)))
        cat(sprintf("  Mean shift (effective - original): %.4f\n",
                    mean(pass_dt$quality_effective - pass_dt$quality_original)))
}


cat("\n--- ASW: Before vs After Manipulation (PASS only) ---\n")
if (nrow(pass_dt) > 0) {
        cat("  [Subspace = selected variables only | Full = all variables]\n\n")
        cat(sprintf("  Before (subspace) — mean: %.4f, range: [%.4f, %.4f]\n",
                    mean(pass_dt$asw_orig_sub), min(pass_dt$asw_orig_sub), max(pass_dt$asw_orig_sub)))
        cat(sprintf("  After  (subspace) — mean: %.4f, range: [%.4f, %.4f]\n",
                    mean(pass_dt$asw_after_sub, na.rm = TRUE),
                    min(pass_dt$asw_after_sub, na.rm = TRUE),
                    max(pass_dt$asw_after_sub, na.rm = TRUE)))
        cat(sprintf("  Subspace change:    %.4f\n\n",
                    mean(pass_dt$asw_after_sub - pass_dt$asw_orig_sub, na.rm = TRUE)))
        cat(sprintf("  Before (full)     — mean: %.4f, range: [%.4f, %.4f]\n",
                    mean(pass_dt$asw_orig_full), min(pass_dt$asw_orig_full), max(pass_dt$asw_orig_full)))
        cat(sprintf("  After  (full)     — mean: %.4f, range: [%.4f, %.4f]\n",
                    mean(pass_dt$asw_after_full, na.rm = TRUE),
                    min(pass_dt$asw_after_full, na.rm = TRUE),
                    max(pass_dt$asw_after_full, na.rm = TRUE)))
        cat(sprintf("  Full-space change:  %.4f\n",
                    mean(pass_dt$asw_after_full - pass_dt$asw_orig_full, na.rm = TRUE)))
}

# ==============================================================================
# Which ones failed the simulation?
# ==============================================================================

setwd("E://Arbeit/Projekte/02_ongoing/PULSE/wp1/data/")
sd <- readRDS("pulseD/data/000_biota/03_diatoms_scheme.rds")
sf <- readRDS("pulseF/data/000_biota/03_fish_scheme.rds")
si <- readRDS("pulseI/data/000_biota/03_invertebrates_scheme.rds")
sm <- readRDS("pulseM/data/000_biota/03_macrophytes_scheme.rds")

sd$included <- TRUE
sf$included <- TRUE
si$included <- TRUE
sm$included <- TRUE
s <- rbindlist(list(sd,sf,si,sm))

for (j in 1:4){
        taxon <- c("diatoms", "fish", "invertebrates", "macrophytes")[j]
        for (i in ms[[j]]) {
                id <- case_when(
                        i<10 ~ paste0(taxon,"_000",i),
                        i>=10 & i<100 ~ paste0(taxon,"_00",i),
                        i>=100 & i<1000 ~ paste0(taxon,"_0",i)
                )
                
                s[scheme_id == id, included := FALSE]
        }    
}

# "diatoms_hungary_ecosurv"  "diatoms_switzerland_NAWA"
s %>%  filter(included == TRUE) %>% pull(data.set) %>%  unique() %>% setdiff(unique(s[taxon == "diatoms", data.set]), .)
s %>%  filter(included == TRUE) %>% pull(data.set) %>%  unique() %>% setdiff(unique(s[taxon == "fish", data.set]), .)
s %>%  filter(included == TRUE) %>% pull(data.set) %>%  unique() %>% setdiff(unique(s[taxon == "invertebrates", data.set]), .)
s %>%  filter(included == TRUE) %>% pull(data.set) %>%  unique() %>% setdiff(unique(s[taxon == "macrophytes", data.set]), .)




for (i in 1:nrow(sf)) if(i %in% ff) sf$included[i]<-FALSE
for (i in 1:nrow(si)) if(i %in% fi) si$included[i]<-FALSE
for (i in 1:nrow(sm)) if(i %in% fm) sm$included[i]<-FALSE


sd %>% ggplot(aes(y = samples, x = included)) + geom_boxplot()
sf %>% ggplot(aes(y = samples, x = included)) + geom_boxplot()
si %>% ggplot(aes(y = samples, x = included)) + geom_boxplot()
sm %>% ggplot(aes(y = samples, x = included)) + geom_boxplot()

sd %>% ggplot(aes(y = eventYear  , x = included)) + geom_boxplot()
sf %>% ggplot(aes(y = eventYear  , x = included)) + geom_boxplot()
si %>% ggplot(aes(y = eventYear  , x = included)) + geom_boxplot()
sm %>% ggplot(aes(y = eventYear  , x = included)) + geom_boxplot()

sd %>% ggplot(aes(x = included)) + geom_bar(aes(fill = sample_type))
sf %>% ggplot(aes(x = included)) + geom_bar(aes(fill = sample_type))
si %>% ggplot(aes(x = included)) + geom_bar(aes(fill = sample_type))
sm %>% ggplot(aes(x = included)) + geom_bar(aes(fill = sample_type))


# 1. Unnest the list column so each month is a separate row
pd <- sd %>% unnest(focal_months) %>% mutate(focal_months = as.factor(focal_months))
pf <- sf %>% unnest(focal_months) %>% mutate(focal_months = as.factor(focal_months))
pi <- si %>% unnest(focal_months) %>% mutate(focal_months = as.factor(focal_months))
pm <- sm %>% unnest(focal_months) %>% mutate(focal_months = as.factor(focal_months))

monthplot <- function(x){
        ggplot(x, aes(x = focal_months, fill = included)) +
                geom_bar(position = "fill") + # "fill" shows proportion (0 to 1)
                scale_y_continuous(labels = scales::percent) +
                labs(
                        title = "Prevalence of 'Included' Status by Month",
                        x = "Month",
                        y = "Percentage",
                        fill = "Included"
                ) +
                theme_minimal()
}

monthplot(pd)
monthplot(pf)
monthplot(pi)
monthplot(pm)

## Excluded data sets 
# "diatoms_hungary_ecosurv"  "diatoms_switzerland_NAWA"
sd %>%  filter(included == TRUE) %>% pull(data.set) %>%  unique() %>% setdiff(unique(sd$data.set), .)
# "invertebrates_italy_po"         "invertebrates_switzerland_BDM"  "invertebrates_switzerland_NAWA"
si %>%  filter(included == TRUE) %>% pull(data.set) %>%  unique() %>% setdiff(unique(si$data.set), .)
# "fish_hungary_ecosurv"
sf %>%  filter(included == TRUE) %>% pull(data.set) %>%  unique() %>% setdiff(unique(sf$data.set), .)
# none
sm %>% filter(included == TRUE) %>% pull(data.set) %>%  unique() %>% setdiff(unique(sm$data.set), .)


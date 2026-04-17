## Investigate simulation failures 



# setup -------------------------------------------------------------------
library(data.table)
library(tidyverse)

setwd("E://Arbeit/Projekte/02_ongoing/PULSE/wp1/data/")
sd <- readRDS("pulseD/data/000_biota/03_diatoms_scheme.rds")
sf <- readRDS("pulseF/data/000_biota/03_fish_scheme.rds")
si <- readRDS("pulseI/data/000_biota/03_invertebrates_scheme.rds")
sm <- readRDS("pulseM/data/000_biota/03_macrophytes_scheme.rds")


# plots -------------------------------------------------------------------
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


# monthplots --------------------------------------------------------------
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
 
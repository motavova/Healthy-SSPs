# =============================================================================
# Regional Summary Measures
# Updated: 14/08/2026
# Author: Martina Otavova
# =============================================================================



setwd("/home/otavova/Healthy-SSPs/SSPs_CVD_mortality")

library(dplyr)

out_dir <- "Results/RegionalSummaryMeasures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source("Transcript/Sellers_regions.R")

region_lookup <- bind_rows(
  sellers_regions,
  tibble(iso3 = c("CMR", "KIR"), region = c("Sub-Saharan Africa", "East Asia & Pacific"))
)

# =============================================================================
# 1. Load observed + projected series, attach region
# =============================================================================

obs_proj <- readRDS("Results/Projections/obs_and_proj.rds") %>%
  mutate(iso3 = as.character(iso3), sex = as.character(sex)) %>%
  select(scenario, iso3, age, sex, year, pop,
        pred_cardio_deaths = cardio_deaths,
        pred_allcause_deaths = allcause_deaths) %>%
  left_join(region_lookup, by = "iso3")



# =============================================================================
# 2. Regional summaries (pooled-sex, and by sex)
# =============================================================================

regional_summary <- obs_proj %>%
  filter(!is.na(region)) %>%
  group_by(scenario, region, year) %>%
  summarise(
    cardio_deaths   = sum(pred_cardio_deaths,   na.rm = TRUE),
    allcause_deaths = sum(pred_allcause_deaths, na.rm = TRUE),
    n_countries     = n_distinct(iso3),
    .groups = "drop"
  ) %>%
  mutate(
    cardio_share           = cardio_deaths / allcause_deaths,
    cardio_deaths_millions = cardio_deaths / 1e6
  )

regional_summary_sex <- obs_proj %>%
  filter(!is.na(region)) %>%
  group_by(scenario, region, sex, year) %>%
  summarise(
    cardio_deaths   = sum(pred_cardio_deaths,   na.rm = TRUE),
    allcause_deaths = sum(pred_allcause_deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cardio_share = cardio_deaths / allcause_deaths)

saveRDS(regional_summary,     file.path(out_dir, "regional_summary_sellers.rds"))
saveRDS(regional_summary_sex, file.path(out_dir, "regional_summary_sellers_sex.rds"))



# =============================================================================
# 3. Repeat using World Bank regions instead of Sellers regions
# =============================================================================

source("Transcript/WorldBank_regions.R")

obs_proj_wb <- readRDS("Results/Projections/obs_and_proj.rds") %>%
  mutate(iso3 = as.character(iso3), sex = as.character(sex)) %>%
  select(scenario, iso3, age, sex, year, pop,
        pred_cardio_deaths = cardio_deaths,
        pred_allcause_deaths = allcause_deaths) %>%
  left_join(worldbank_regions, by = "iso3")

regional_summary_wb <- obs_proj_wb %>%
  filter(!is.na(region)) %>%
  group_by(scenario, region, year) %>%
  summarise(
    cardio_deaths   = sum(pred_cardio_deaths,   na.rm = TRUE),
    allcause_deaths = sum(pred_allcause_deaths, na.rm = TRUE),
    n_countries     = n_distinct(iso3),
    .groups = "drop"
  ) %>%
  mutate(
    cardio_share           = cardio_deaths / allcause_deaths,
    cardio_deaths_millions = cardio_deaths / 1e6
  )

regional_summary_wb_sex <- obs_proj_wb %>%
  filter(!is.na(region)) %>%
  group_by(scenario, region, sex, year) %>%
  summarise(
    cardio_deaths   = sum(pred_cardio_deaths,   na.rm = TRUE),
    allcause_deaths = sum(pred_allcause_deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cardio_share = cardio_deaths / allcause_deaths)

saveRDS(regional_summary_wb,     file.path(out_dir, "regional_summary_worldbank.rds"))
saveRDS(regional_summary_wb_sex, file.path(out_dir, "regional_summary_worldbank_sex.rds"))


# =============================================================================
# Regional Summary Measures
# Updated: 14/08/2026
# Author: Martina Otavova
# =============================================================================



setwd("/home/otavova/Healthy-SSPs/SSPs_CVD_mortality")

library(dplyr)

out_dir <- "Results/RegionalSummaryMeasures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source("Transcript/WorldBankClassification.R")

region_lookup <- world_bank_classification  # CMR/KIR are now covered directly in the table above; no manual patch needed

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

saveRDS(regional_summary,     file.path(out_dir, "regional_summary_worldbankupdate.rds"))
saveRDS(regional_summary_sex, file.path(out_dir, "regional_summary_worldbankupdate_sex.rds"))

# Note: this used to also produce a second scheme here ("World Bank regions",
# via WorldBank_regions.R / countrycode's region field, no High-Income
# carve-out). That scheme's country membership was stale (didn't reflect the
# World Bank's current classification -- see WorldBankClassification.R's
# header) and has been retired; its last output is archived at
# Results/History/RegionalSummaries_old_worldbank/. WorldBankClassification.R
# above is now the single, current-vintage regional scheme used everywhere.


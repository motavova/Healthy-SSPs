# =============================================================================
# Age- and Sex-Specific Regional CVD Mortality Rates
#
# Builds a region x scenario x year x age x sex table of cardiovascular
# mortality rates (deaths per capita and per 100,000), aggregating the
# country-level obs_and_proj.rds series up to the World Bank's current
# regional classification (see WorldBankClassification.R).
# Used to support age/sex decomposition analyses (e.g. the sex crossover
# and age-standardization checks in the report's "Differences Across Sex"
# section) without recomputing country-level joins each time.
# =============================================================================

setwd("/home/otavova/Healthy-SSPs/SSPs_CVD_mortality")

.libPaths(c("/home/otavova/R/library"))

library(dplyr)

out_dir <- "Results/RegionalSummaryMeasures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source("Transcript/WorldBankClassification.R")

region_lookup <- world_bank_classification  # CMR/KIR are now covered directly in the table above; no manual patch needed

obs_proj <- readRDS("Results/Projections/obs_and_proj.rds") %>%
  mutate(iso3 = as.character(iso3), sex = as.character(sex)) %>%
  left_join(region_lookup, by = "iso3") %>%
  filter(!is.na(region))                    # drop any territories with no region assignment

age_levels <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44",
                "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84",
                "85-89", "90-94", "95+")

age_sex_region_rates <- obs_proj %>%
  mutate(age = factor(age, levels = age_levels)) %>%
  group_by(scenario, region, year, age, sex) %>%
  summarise(
    cardio_deaths   = sum(cardio_deaths,   na.rm = TRUE),   # regional total cardio deaths
    allcause_deaths = sum(allcause_deaths, na.rm = TRUE),   # regional total all-cause deaths
    pop             = sum(pop,             na.rm = TRUE),   # regional total population
    n_countries     = n_distinct(iso3),
    .groups = "drop"
  ) %>%
  mutate(
    cardio_mortality_rate       = cardio_deaths / pop,             # deaths per capita
    cardio_mortality_rate_100k  = cardio_mortality_rate * 1e5,       # deaths per 100,000
    allcause_mortality_rate     = allcause_deaths / pop,
    allcause_mortality_rate_100k = allcause_mortality_rate * 1e5,
    cardio_share                = cardio_deaths / allcause_deaths
  ) %>%
  arrange(region, scenario, year, age, sex)

saveRDS(age_sex_region_rates, file.path(out_dir, "age_sex_region_mortality_rates.rds"))
write.csv(age_sex_region_rates, file.path(out_dir, "age_sex_region_mortality_rates.csv"), row.names = FALSE)

cat("Saved:\n")
cat(" -", file.path(out_dir, "age_sex_region_mortality_rates.rds"), "\n")
cat(" -", file.path(out_dir, "age_sex_region_mortality_rates.csv"), "\n")
cat("Rows:", nrow(age_sex_region_rates), "\n")
cat("Regions:", paste(sort(unique(age_sex_region_rates$region)), collapse = ", "), "\n")
cat("Scenarios:", paste(sort(unique(age_sex_region_rates$scenario)), collapse = ", "), "\n")
cat("Years:", paste(range(age_sex_region_rates$year), collapse = "-"), "\n")
cat("Age groups:", length(levels(age_sex_region_rates$age)), "\n")

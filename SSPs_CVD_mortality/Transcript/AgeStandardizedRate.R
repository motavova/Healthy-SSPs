# =============================================================================
# AGE-STANDARDIZED MORTALITY RATE (CVD + ALL-CAUSE)
# Updated: 08/07/2026
# Author: Martina Otavova
# =============================================================================
#
# Direct age standardization of the observed + projected death series
# (CVD_projection.R output) onto the WHO World Standard Population
# (Ahmad et al. 2001), producing age-standardized CVD and all-cause
# mortality rates per 100,000 population for every scenario x country x
# year, both by sex and with sexes pooled.
#
# Output (Results/AgeStandardizedRate/):
#   asmr_by_sex.rds     -- scenario x iso3 x sex x year x cause (CVD / All-cause)
#   asmr_both_sexes.rds -- scenario x iso3 x year x cause (CVD / All-cause), sexes pooled
# =============================================================================

.libPaths(c("/home/otavova/R/library"))  # point R at the project-specific package library

library(dplyr)  # data wrangling (mutate, filter, joins, group_by, summarise)
library(tibble) # tribble() for the standard-population weight table

out_dir <- "Results/AgeStandardizedRate"  # where the standardized-rate tables get saved
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)  # create out_dir if it doesn't exist yet (silently)

# =============================================================================
# 1. WHO World Standard Population weights (Ahmad et al. 2001)
# =============================================================================
# 5-year age-band weights per 100 population, used to reweight each age
# group's observed rate onto a common reference population before summing.
# The 95-99 and 100+ bands are combined into "95+" to match the age bands
# present in the mortality data.

who_std_pop <- tribble(
  ~age,    ~std_weight,
  "0-4",     8.86,
  "5-9",     8.69,
  "10-14",   8.60,
  "15-19",   8.47,
  "20-24",   8.22,
  "25-29",   7.93,
  "30-34",   7.61,
  "35-39",   7.15,
  "40-44",   6.59,
  "45-49",   6.04,
  "50-54",   5.37,
  "55-59",   4.55,
  "60-64",   3.72,
  "65-69",   2.96,
  "70-74",   2.21,
  "75-79",   1.52,
  "80-84",   0.91,
  "85-89",   0.44,
  "90-94",   0.15,
  "95+",     0.045
)  # weights sum to ~100; any rounding is immaterial since standardization renormalizes by the summed weight actually used

# =============================================================================
# 2. Load observed history + projections (CVD_projection.R output)
# =============================================================================

obs_and_proj <- readRDS("Results/Projections/obs_and_proj.rds") %>%  # load the combined observed + projected dataset
  mutate(iso3 = as.character(iso3), sex = as.character(sex))  # ensure iso3/sex are plain character vectors, not factors

message("Rows loaded: ", nrow(obs_and_proj))
message("Scenarios: ", paste(sort(unique(obs_and_proj$scenario)), collapse = ", "))
message("Age groups: ", n_distinct(obs_and_proj$age), " (expected 20)")

# =============================================================================
# 3. Death + population counts by age: by sex, and with sexes pooled
# =============================================================================

counts_by_sex <- obs_and_proj %>%
  select(scenario, iso3, age, sex, year, cardio_deaths, allcause_deaths, pop)

counts_both_sexes <- obs_and_proj %>%
  group_by(scenario, iso3, age, year) %>%  # pool male + female counts before computing any rate (not an average of two rates)
  summarise(cardio_deaths = sum(cardio_deaths, na.rm = TRUE),
           allcause_deaths = sum(allcause_deaths, na.rm = TRUE),
           pop = sum(pop, na.rm = TRUE), .groups = "drop")

# =============================================================================
# 4. Direct standardization onto the WHO World Standard Population
# =============================================================================
# ASR = sum(age-specific rate * standard weight) / sum(standard weight used)
# Dividing by the summed weight actually used (rather than a fixed 100)
# re-normalizes automatically when an age group is missing for a given
# country/year/scenario, instead of silently understating the rate.

standardize <- function(df, deaths_col, cause_label, group_cols) {
  df %>%
    mutate(rate_per_100k = .data[[deaths_col]] / pop * 1e5) %>%  # age-specific mortality rate per 100,000 population for this cause
    inner_join(who_std_pop, by = "age") %>%  # attach each age band's standard-population weight
    group_by(across(all_of(group_cols))) %>%
    summarise(
      asmr_per_100k   = sum(rate_per_100k * std_weight) / sum(std_weight),  # weighted average = the age-standardized rate
      deaths          = sum(.data[[deaths_col]]),                            # total deaths (this cause) across the age groups used
      pop             = sum(pop),                                               # total population across the age groups used
      n_age_groups    = n(),                                                      # how many of the 20 age bands contributed
      weight_coverage = sum(std_weight) / sum(who_std_pop$std_weight),              # fraction of the standard population covered
      .groups = "drop"
    ) %>%
    mutate(cause = cause_label)
}

asmr_by_sex <- bind_rows(
  standardize(counts_by_sex, "cardio_deaths",   "CVD",       c("scenario", "iso3", "sex", "year")),
  standardize(counts_by_sex, "allcause_deaths", "All-cause", c("scenario", "iso3", "sex", "year"))
)  # age-standardized rate, by sex, both causes

asmr_both_sexes <- bind_rows(
  standardize(counts_both_sexes, "cardio_deaths",   "CVD",       c("scenario", "iso3", "year")),
  standardize(counts_both_sexes, "allcause_deaths", "All-cause", c("scenario", "iso3", "year"))
)  # age-standardized rate, sexes pooled, both causes

n_incomplete_sex   <- sum(asmr_by_sex$n_age_groups < 20)  # groups standardized on fewer than all 20 age bands
n_incomplete_both  <- sum(asmr_both_sexes$n_age_groups < 20)
message("asmr_by_sex: ", nrow(asmr_by_sex), " rows, ", n_incomplete_sex, " with <20 age groups")
message("asmr_both_sexes: ", nrow(asmr_both_sexes), " rows, ", n_incomplete_both, " with <20 age groups")

# =============================================================================
# 5. Save results
# =============================================================================

saveRDS(asmr_by_sex, file.path(out_dir, "asmr_by_sex.rds"))  # persist the by-sex age-standardized rates
saveRDS(asmr_both_sexes, file.path(out_dir, "asmr_both_sexes.rds"))  # persist the pooled-sex age-standardized rates

cat("Saved age-standardized CVD and all-cause mortality rates to", out_dir, "\n")

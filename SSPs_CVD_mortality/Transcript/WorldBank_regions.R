# =============================================================================
# WorldBank_regions.R
# =============================================================================
# Purpose : Region lookup table using the World Bank's traditional scheme --
#           every country keeps its own geographic region (7 regions,
#           including North America) regardless of income level. This is
#           unlike Sellers_regions.R, which pulls high-income countries out
#           of their geography into a single "High-Income" bucket.
#           Regions sourced from the `countrycode` package's `region` field,
#           which mirrors the World Bank's current regional classification.
# Usage   : source("WorldBank_regions.R")  ->  adds `worldbank_regions` data frame
#           then: left_join(your_data, worldbank_regions, by = "iso3")
# =============================================================================

setwd("/home/otavova/Healthy-SSPs/SSPs_CVD_mortality")

library(dplyr)
library(countrycode)

obs_proj <- readRDS("Results/Projections/obs_and_proj.rds")

worldbank_regions <- tibble(iso3 = sort(unique(as.character(obs_proj$iso3)))) %>%
  mutate(region = countrycode(iso3, "iso3c", "region",
                              custom_match = c(PRI = "Latin America & Caribbean")))


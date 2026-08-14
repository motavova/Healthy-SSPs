# =============================================================================
# Global Plots
# Updated: 14/08/2026
# Author: Martina Otavova
#
# Global (all-country) CVD summary figures by SSP, 2030-2100:
#   fig_global_ssp_deaths.png    -- CVD share of all-cause deaths
#   fig_global_cvd_deaths.png    -- number of CVD deaths (millions)
#   fig_global_asmr.png          -- age-standardized CVD mortality rate
#   fig_global_asmr_premature.png -- age-standardized CVD mortality rate, ages 0-69
# =============================================================================

setwd("/home/otavova/Healthy-SSPs/SSPs_CVD_mortality")

.libPaths(c("/home/otavova/R/library"))

library(dplyr)
library(ggplot2)
library(scales)
library(tibble)

out_dir <- "Results/Plots/Global"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ssp_colours <- c("SSP1" = "#1B9E77", "SSP2" = "#D95F02", "SSP3" = "#7570B3",
                "SSP4" = "#E7298A", "SSP5" = "#66A61E")

plot_from_year <- 2030

obs_and_proj <- readRDS("Results/Projections/obs_and_proj.rds") %>%
  mutate(iso3 = as.character(iso3), sex = as.character(sex)) %>%
  filter(scenario != "historical", year >= plot_from_year)

# =============================================================================
# 1. fig_global_ssp_deaths: CVD share of all-cause deaths, by SSP
# =============================================================================

global_summary <- obs_and_proj %>%
  group_by(scenario, year) %>%
  summarise(
    cardio_deaths   = sum(cardio_deaths,   na.rm = TRUE),
    allcause_deaths = sum(allcause_deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cardio_share           = cardio_deaths / allcause_deaths,
    cardio_deaths_millions = cardio_deaths / 1e6
  )

p_share <- ggplot(global_summary, aes(x = year, y = cardio_share, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = ssp_colours, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = c(2030, 2050, 2075, 2100)) +
  labs(title = "Cardiovascular share of all-cause deaths",
      x = NULL, y = "Proportion of deaths") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "fig_global_ssp_deaths.png"), p_share, width = 10, height = 6, dpi = 300)

# =============================================================================
# 2. fig_global_cvd_deaths: number of CVD deaths, by SSP
# =============================================================================

p_deaths <- ggplot(global_summary, aes(x = year, y = cardio_deaths_millions, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = ssp_colours, name = NULL) +
  scale_x_continuous(breaks = c(2030, 2050, 2075, 2100)) +
  labs(title = "Global cardiovascular deaths",
      x = NULL, y = "CVD deaths (millions)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "fig_global_cvd_deaths.png"), p_deaths, width = 10, height = 6, dpi = 300)

# =============================================================================
# 3. fig_global_asmr: age-standardized CVD mortality rate, by SSP
# =============================================================================
# Direct standardization onto the WHO World Standard Population (Ahmad et
# al. 2001), same reference weights as AgeStandardizedRate.R, applied here
# to global (all-country, both-sex) age-specific death and population
# totals instead of country-level ones.

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
)

premature_ages <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34",
                    "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69")  # ages 0-69

global_age <- obs_and_proj %>%
  group_by(scenario, age, year) %>%
  summarise(
    cardio_deaths = sum(cardio_deaths, na.rm = TRUE),
    pop           = sum(pop,           na.rm = TRUE),
    .groups = "drop"
  )

global_asmr <- global_age %>%
  mutate(rate_per_100k = cardio_deaths / pop * 1e5) %>%
  inner_join(who_std_pop, by = "age") %>%
  group_by(scenario, year) %>%
  summarise(asmr_per_100k = sum(rate_per_100k * std_weight) / sum(std_weight), .groups = "drop")

p_asmr <- ggplot(global_asmr, aes(x = year, y = asmr_per_100k, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.3) +
  scale_colour_manual(values = ssp_colours, name = NULL) +
  scale_x_continuous(breaks = c(2030, 2050, 2075, 2100)) +
  labs(title = "Global age-standardized cardiovascular mortality rate",
      x = NULL, y = "CVD ASMR (per 100,000)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "fig_global_asmr.png"), p_asmr, width = 10, height = 6, dpi = 300)

# =============================================================================
# 4. fig_global_asmr_premature: age-standardized CVD mortality rate, ages 0-69
# =============================================================================
# Same standardization as above, but restricted to the under-70 age bands, so
# the rate reflects premature (pre-70) mortality risk rather than the full
# age range.

global_asmr_premature <- global_age %>%
  filter(age %in% premature_ages) %>%
  mutate(rate_per_100k = cardio_deaths / pop * 1e5) %>%
  inner_join(who_std_pop, by = "age") %>%
  group_by(scenario, year) %>%
  summarise(asmr_per_100k = sum(rate_per_100k * std_weight) / sum(std_weight), .groups = "drop")

p_asmr_premature <- ggplot(global_asmr_premature, aes(x = year, y = asmr_per_100k, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.3) +
  scale_colour_manual(values = ssp_colours, name = NULL) +
  scale_x_continuous(breaks = c(2030, 2050, 2075, 2100)) +
  labs(title = "Global age-standardized premature cardiovascular mortality rate",
      x = NULL, y = "CVD ASMR, ages 0-69 (per 100,000)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "fig_global_asmr_premature.png"), p_asmr_premature, width = 10, height = 6, dpi = 300)

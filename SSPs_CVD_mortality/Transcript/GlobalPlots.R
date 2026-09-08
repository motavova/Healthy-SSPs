# =============================================================================
# Global Plots
# Updated: 24/08/2026
# Author: Martina Otavova
#
# Global (all-country) CVD summary figures by SSP, 1990-2100, historical +
# SSP1-5, with each SSP line connected seamlessly to the historical line
# (no visual gap at the 2023-2025 join, even though the underlying data has
# no 2024 point -- see CVD_projection.R for why):
#   fig_global_ssp_deaths.png/.svg    -- CVD share of all-cause deaths
#   fig_global_cvd_deaths.png/.svg    -- number of CVD deaths (millions)
#   fig_global_asmr.png/.svg          -- age-standardized CVD mortality rate
#   fig_global_asmr_premature.png/.svg -- age-standardized CVD mortality rate, ages 0-69
#
# SVG output uses grDevices::svg() (Cairo-backed) rather than the svglite
# package, which cannot be installed in this environment (missing system
# fontconfig/freetype headers). Both produce valid vector SVG; svglite
# tends to yield tidier markup, but the rendered output is equivalent for
# presentation purposes.
# =============================================================================

setwd("/home/otavova/Healthy-SSPs/SSPs_CVD_mortality")

.libPaths(c("/home/otavova/R/library"))

library(dplyr)
library(ggplot2)
library(scales)
library(tibble)

out_dir <- "Results/Plots/Global"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ssp_colours <- c("historical" = "grey35", "SSP1" = "#1B9E77", "SSP2" = "#D95F02",
                 "SSP3" = "#7570B3", "SSP4" = "#E7298A", "SSP5" = "#66A61E")
ssp_breaks  <- c("historical", "SSP1", "SSP2", "SSP3", "SSP4", "SSP5")
ssp_labels  <- c("Historical", "SSP1", "SSP2", "SSP3", "SSP4", "SSP5")

last_hist_year <- 2023   # last year of observed data
proj_start_year <- 2025  # earliest SSP projection year -- there is no 2024 row (see CVD_projection.R)

obs_and_proj <- readRDS("Results/Projections/obs_and_proj.rds") %>%
  mutate(iso3 = as.character(iso3), sex = as.character(sex)) %>%
  filter(scenario == "historical" | year >= proj_start_year)

# Duplicates each SSP's row at last_hist_year with the historical value, so
# every colored SSP line starts from the exact point where the grey
# historical line ends -- geom_line then draws an unbroken line through the
# 2023-2025 gap instead of a visual jump.
connect_to_historical <- function(df, value_cols) {
  hist_anchor <- df %>% filter(scenario == "historical", year == last_hist_year)
  ssps <- setdiff(unique(df$scenario), "historical")
  anchors <- lapply(ssps, function(s) hist_anchor %>% mutate(scenario = s)) %>% bind_rows()
  bind_rows(df, anchors) %>% arrange(scenario, year)
}

save_both <- function(plot, filename_stem, width, height, dpi = 300) {
  ggsave(file.path(out_dir, paste0(filename_stem, ".png")), plot, width = width, height = height, dpi = dpi)
  ggsave(file.path(out_dir, paste0(filename_stem, ".svg")), plot, width = width, height = height, device = grDevices::svg)
}

base_theme <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

x_breaks <- c(1990, 2010, proj_start_year, 2050, 2075, 2100)

# =============================================================================
# 1. fig_global_ssp_deaths: CVD share of all-cause deaths, by SSP
# 2. fig_global_cvd_deaths: number of CVD deaths, by SSP
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
  ) %>%
  connect_to_historical()

p_share <- ggplot(global_summary, aes(x = year, y = cardio_share, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = ssp_colours, breaks = ssp_breaks, labels = ssp_labels, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = x_breaks) +
  labs(title = "Cardiovascular share of all-cause deaths",
      subtitle = "Historical (1990-2023) and SSP projections (2024-2100)",  # label convention: year after last historical point; underlying data still starts 2025 (no 2024 row exists)
      x = NULL, y = "Proportion of deaths") +
  base_theme +
  theme(plot.title    = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(size = 13),
        axis.text     = element_text(size = 13), axis.title = element_text(size = 15),
        legend.text   = element_text(size = 13))

save_both(p_share, "fig_global_ssp_deaths", width = 10, height = 6)

p_deaths <- ggplot(global_summary, aes(x = year, y = cardio_deaths_millions, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = ssp_colours, breaks = ssp_breaks, labels = ssp_labels, name = NULL) +
  scale_x_continuous(breaks = x_breaks) +
  labs(title = "Global cardiovascular deaths",
      subtitle = "Historical (1990-2023) and SSP projections (2025-2100)",
      x = NULL, y = "CVD deaths (millions)") +
  base_theme

save_both(p_deaths, "fig_global_cvd_deaths", width = 10, height = 6)

# =============================================================================
# 3. fig_global_asmr: age-standardized CVD mortality rate, by SSP
# 4. fig_global_asmr_premature: age-standardized CVD mortality rate, ages 0-69
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

standardize_global <- function(df) {
  df %>%
    mutate(rate_per_100k = cardio_deaths / pop * 1e5) %>%
    inner_join(who_std_pop, by = "age") %>%
    group_by(scenario, year) %>%
    summarise(asmr_per_100k = sum(rate_per_100k * std_weight) / sum(std_weight), .groups = "drop")
}

global_asmr <- standardize_global(global_age) %>% connect_to_historical()

p_asmr <- ggplot(global_asmr, aes(x = year, y = asmr_per_100k, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = ssp_colours, breaks = ssp_breaks, labels = ssp_labels, name = NULL) +
  scale_x_continuous(breaks = x_breaks) +
  labs(title = "Global age-standardized cardiovascular mortality rate",
      subtitle = "Historical (1990-2023) and SSP projections (2025-2100)",
      x = NULL, y = "CVD ASMR (per 100,000)") +
  base_theme

save_both(p_asmr, "fig_global_asmr", width = 10, height = 6)

global_asmr_premature <- global_age %>%
  filter(age %in% premature_ages) %>%
  standardize_global() %>%
  connect_to_historical()

p_asmr_premature <- ggplot(global_asmr_premature, aes(x = year, y = asmr_per_100k, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = ssp_colours, breaks = ssp_breaks, labels = ssp_labels, name = NULL) +
  scale_x_continuous(breaks = x_breaks) +
  labs(title = "Global age-standardized premature cardiovascular mortality rate",
      subtitle = "Historical (1990-2023) and SSP projections (2025-2100)",
      x = NULL, y = "CVD ASMR, ages 0-69 (per 100,000)") +
  base_theme

save_both(p_asmr_premature, "fig_global_asmr_premature", width = 10, height = 6)

cat("Saved (PNG + SVG) to", out_dir, ":\n")
cat(" - fig_global_ssp_deaths\n - fig_global_cvd_deaths\n - fig_global_asmr\n - fig_global_asmr_premature\n")

# =============================================================================
# 5. Projections-only variants (dots): drop the observed part, start at 2030,
#    mark every 5-year data point, SVG only. New files -- the four originals
#    above are left untouched.
# =============================================================================

dots_start_year <- 2030

save_svg <- function(plot, filename_stem, width, height) {
  ggsave(file.path(out_dir, paste0(filename_stem, ".svg")), plot, width = width, height = height, device = grDevices::svg)
}

dots_theme <- theme_bw(base_size = 15) +
  theme(
    legend.position    = "bottom",
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 18),
    plot.subtitle      = element_text(size = 13),
    axis.title         = element_text(size = 15),
    axis.text          = element_text(size = 13),
    legend.text        = element_text(size = 13)
  )

proj_only <- obs_and_proj %>% filter(scenario != "historical", year >= dots_start_year)

global_summary_dots <- proj_only %>%
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

dots_x_breaks <- seq(dots_start_year, 2100, by = 10)

p_share_dots <- ggplot(global_summary_dots, aes(x = year, y = cardio_share, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = ssp_colours, breaks = ssp_breaks[-1], labels = ssp_labels[-1], name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = dots_x_breaks) +
  labs(title = "Cardiovascular share of all-cause deaths",
      subtitle = paste0("SSP projections, ", dots_start_year, "-2100"),
      x = NULL, y = "Proportion of deaths") +
  dots_theme

save_svg(p_share_dots, "fig_global_ssp_deaths_dots", width = 10, height = 6)

p_deaths_dots <- ggplot(global_summary_dots, aes(x = year, y = cardio_deaths_millions, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = ssp_colours, breaks = ssp_breaks[-1], labels = ssp_labels[-1], name = NULL) +
  scale_x_continuous(breaks = dots_x_breaks) +
  labs(title = "Global cardiovascular deaths",
      subtitle = paste0("SSP projections, ", dots_start_year, "-2100"),
      x = NULL, y = "CVD deaths (millions)") +
  dots_theme

save_svg(p_deaths_dots, "fig_global_cvd_deaths_dots", width = 10, height = 6)

global_age_dots <- proj_only %>%
  group_by(scenario, age, year) %>%
  summarise(
    cardio_deaths = sum(cardio_deaths, na.rm = TRUE),
    pop           = sum(pop,           na.rm = TRUE),
    .groups = "drop"
  )

global_asmr_dots <- standardize_global(global_age_dots)

p_asmr_dots <- ggplot(global_asmr_dots, aes(x = year, y = asmr_per_100k, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = ssp_colours, breaks = ssp_breaks[-1], labels = ssp_labels[-1], name = NULL) +
  scale_x_continuous(breaks = dots_x_breaks) +
  labs(title = "Global age-standardized cardiovascular mortality rate",
      subtitle = paste0("SSP projections, ", dots_start_year, "-2100"),
      x = NULL, y = "CVD ASMR (per 100,000)") +
  dots_theme

save_svg(p_asmr_dots, "fig_global_asmr_dots", width = 10, height = 6)

global_asmr_premature_dots <- global_age_dots %>%
  filter(age %in% premature_ages) %>%
  standardize_global()

p_asmr_premature_dots <- ggplot(global_asmr_premature_dots, aes(x = year, y = asmr_per_100k, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = ssp_colours, breaks = ssp_breaks[-1], labels = ssp_labels[-1], name = NULL) +
  scale_x_continuous(breaks = dots_x_breaks) +
  labs(title = "Global age-standardized premature cardiovascular mortality rate",
      subtitle = paste0("SSP projections, ", dots_start_year, "-2100"),
      x = NULL, y = "CVD ASMR, ages 0-69 (per 100,000)") +
  dots_theme

save_svg(p_asmr_premature_dots, "fig_global_asmr_premature_dots", width = 10, height = 6)

cat("\nSaved (SVG only) to", out_dir, ":\n")
cat(" - fig_global_ssp_deaths_dots\n - fig_global_cvd_deaths_dots\n - fig_global_asmr_dots\n - fig_global_asmr_premature_dots\n")

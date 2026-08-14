# =============================================================================
# LifeExpectancy.R
#
# Life expectancy with and without cardiovascular disease (CVD), computed
# from this project's projected series (Results/Projections/obs_and_proj.rds,
# SSP1-5 scenarios), Chiang's method (standard multiple-decrement /
# "cause-deleted" life table, 5-year age groups). Restricted to a 5-year step
# grid (year %% 5 == 0), which is already the native step for the SSP data.
#
# For each scenario x country x sex x year, two life tables are built from
# the same all-cause death rates:
#   (a) "all-cause"  -- includes every cause of death, including CVD
#   (b) "non-cardio" -- CVD's share of the death rate is removed from every
#                       age group before building the table
# The life-expectancy-at-birth from each is e0_allcause and e0_noncardio.
# Because removing a cause of death can only lengthen life expectancy,
# e0_noncardio >= e0_allcause always. The gap between them,
#   e0_gap_cardio = e0_noncardio - e0_allcause
# is the number of years of life expectancy CVD "costs" that country/sex/
# scenario/year -- i.e. how much longer people would be expected to live if
# nobody died of cardiovascular disease.
#
# Outputs (Results/LifeExpectancy/):
#   e0_projected.rds                    - one row per scenario x country x sex x year
#   e0_projected_both_sexes.rds         - as above, sexes pooled
#   e0_gap_by_ssp_year.rds              - e0_gap_cardio averaged across
#                                          countries, by SSP x year ("across
#                                          SSPs" comparison)
#   e0_projected_summary_by_ssp_year.rds - summary stats (mean/sd/median/
#                                          p10/p90/min/max) of e0_allcause,
#                                          e0_noncardio, e0_gap_cardio, by
#                                          SSP x year x sex
#   e0_gap_ssp_spread_by_country.rds    - for each country x year, how much
#                                          e0_gap_cardio varies across the 5
#                                          SSPs (max - min) -- which countries'
#                                          CVD life-expectancy toll is most
#                                          sensitive to SSP scenario
#   fig2_gap_by_ssp_over_time.png
#   fig3_gap_spread_across_countries.png
#   country_plots/<Country>.png         - e0_gap_cardio by SSP over time,
#                                          one PNG per country (sexes pooled)
# =============================================================================

# =============================================================================
# 0. Setup
# =============================================================================

setwd("/home/otavova/Healthy-SSPs/SSPs_CVD_mortality")

.libPaths(c("/home/otavova/R/library"))

library(dplyr)
library(ggplot2)
library(countrycode)

out_dir <- "Results/LifeExpectancy"
country_plots_dir <- file.path(out_dir, "country_plots")
dir.create(country_plots_dir, showWarnings = FALSE, recursive = TRUE)

ssp_colours <- c("SSP1" = "#1b9e77", "SSP2" = "#d95f02", "SSP3" = "#7570b3",
                 "SSP4" = "#e7298a", "SSP5" = "#66a61e")

eps <- 1e-9

age_order <- c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
              "40-44","45-49","50-54","55-59","60-64","65-69","70-74",
              "75-79","80-84","85-89","90-94","95+")

sanitize_filename <- function(iso) {
  name <- countrycode(iso, "iso3c", "country.name", warn = FALSE)  # look up the full country name from its ISO3 code
  if (is.na(name) || iso == "TWN") name <- if (iso == "TWN") "Taiwan" else iso  # fall back to "Taiwan" for TWN or to the raw code if lookup failed
  gsub("[^A-Za-z0-9]+", "_", name)  # replace any non-alphanumeric characters with underscores to make a safe filename
}

# =============================================================================
# 1. Life table function
#
# Builds two life tables at once (all-cause, and cause-deleted with CVD
# removed) for every group in `df` (a group is one scenario x country x sex
# x year, or the sex-pooled equivalent), and returns life expectancy at every
# age for both. Vectorized across all groups simultaneously (instead of
# looping group-by-group), because the projected panel has a very large
# number of groups -- a per-group loop would be far too slow.
#
# Required columns in `df`:
#   the columns named in `group_vars` (e.g. c("scenario","iso3","sex","year"))
#   age       : one of the 20 age_order labels
#   mx_all    : annual all-cause central death rate  (all-cause deaths / population)
#   cvd_share : fraction of all-cause deaths that are cardiovascular, at that age
#
# Adds two columns:
#   ex          : remaining life expectancy at the START of each age group,
#                 all-cause. ex for age "0-4" is life expectancy at birth (e0).
#   ex_noncardio: same, but with cardiovascular deaths removed from every age
# =============================================================================

add_life_expectancy <- function(df, group_vars) {

  # Standard 5-year life table conventions:
  #   n  = width of each age interval, in years (open-ended for the last group)
  #   ax = average number of years lived in the interval by those who die in
  #        it. 1.5 for "0-4" (infant/child deaths cluster early in the
  #        interval); 2.5 (the interval midpoint) for every other age group.
  n_by_age  <- setNames(c(rep(5, 19), Inf), age_order)
  ax_by_age <- setNames(c(1.5, rep(2.5, 19)), age_order)

  df <- df %>%
    mutate(
      age   = factor(as.character(age), levels = age_order),
      n_yrs = n_by_age[as.character(age)],
      ax    = ax_by_age[as.character(age)],

      # --- Cause-deleted death rate: remove CVD's share of all-cause
      #     mortality at this age, leaving only non-cardiovascular causes ---
      mx_noncardio = mx_all * (1 - pmin(pmax(cvd_share, 0), 1)),

      # --- Chiang's formula: probability of dying within the interval,
      #     given the central death rate. Same formula applied twice: once
      #     to the all-cause rate, once to the cause-deleted rate. ---
      q_all        = n_yrs * mx_all        / (1 + (n_yrs - ax) * mx_all),
      q_noncardio  = n_yrs * mx_noncardio  / (1 + (n_yrs - ax) * mx_noncardio),

      # Clamp to a valid probability, and force q = 1 in the open-ended last
      # interval (everyone who reaches it eventually dies within it).
      q_all       = if_else(is.infinite(n_yrs), 1, pmin(pmax(q_all, 0), 1)),
      q_noncardio = if_else(is.infinite(n_yrs), 1, pmin(pmax(q_noncardio, 0), 1))
    )

  df %>%
    group_by(across(all_of(group_vars))) %>%
    arrange(age, .by_group = TRUE) %>%
    mutate(
      # lx = survivors to the START of each interval, per 100,000 births.
      # cumprod(1 - q) is the running probability of surviving THROUGH each
      # interval; lag(...) shifts it down one row so lx for "0-4" is exactly
      # 100,000 (nobody has died yet at birth).
      lx           = 100000 * lag(cumprod(1 - q_all),       default = 1),
      lx_noncardio = 100000 * lag(cumprod(1 - q_noncardio), default = 1),

      dx           = lx           * q_all,         # deaths within the interval, all-cause
      dx_noncardio = lx_noncardio * q_noncardio,    # deaths within the interval, cause-deleted

      # Lx = person-years lived within the interval by the cohort.
      # Open-ended last interval: assume a constant death rate from the start
      # of the interval onward, so Lx = lx / mx.
      Lx           = if_else(is.finite(n_yrs),
                              n_yrs * lx           - (n_yrs - ax) * dx,
                              if_else(mx_all       > 0, lx           / mx_all,       0)),
      Lx_noncardio = if_else(is.finite(n_yrs),
                              n_yrs * lx_noncardio - (n_yrs - ax) * dx_noncardio,
                              if_else(mx_noncardio > 0, lx_noncardio / mx_noncardio, 0)),

      # Tx = total remaining person-years from this age onward = this
      # interval's Lx plus every older interval's Lx. rev(cumsum(rev(.)))
      # computes that "sum from here to the end" for every row at once.
      Tx           = rev(cumsum(rev(Lx))),
      Tx_noncardio = rev(cumsum(rev(Lx_noncardio))),

      ex           = Tx           / lx,
      ex_noncardio = Tx_noncardio / lx_noncardio
    ) %>%
    ungroup()
}

# =============================================================================
# 2. Load projected series (SSP1-5, 5-year steps)
# =============================================================================

obs_and_proj <- readRDS("Results/Projections/obs_and_proj.rds") %>%
  mutate(iso3 = as.character(iso3), sex = as.character(sex)) %>%
  filter(scenario != "historical", year %% 5 == 0)

# =============================================================================
# 3. Life tables: by sex, and with sexes pooled
# =============================================================================

proj_input <- obs_and_proj %>%
  mutate(
    age       = factor(as.character(age), levels = age_order),
    mx_all    = allcause_deaths / pmax(pop, 1),
    cvd_share = cardio_deaths   / pmax(allcause_deaths, eps)
  ) %>%
  filter(!is.na(age))

proj_lt <- proj_input %>%
  add_life_expectancy(group_vars = c("scenario", "iso3", "sex", "year")) %>%
  filter(age == "0-4") %>%
  transmute(
    scenario, iso3, sex, year,
    e0_allcause  = ex,
    e0_noncardio = ex_noncardio,
    e0_gap_cardio = ex_noncardio - ex
  )

saveRDS(proj_lt, file.path(out_dir, "e0_projected.rds"))

proj_input_both_sexes <- obs_and_proj %>%
  group_by(scenario, iso3, age, year) %>%  # pool male + female counts before computing any rate (not an average of two rates)
  summarise(
    cardio_deaths   = sum(cardio_deaths,   na.rm = TRUE),
    allcause_deaths = sum(allcause_deaths, na.rm = TRUE),
    pop             = sum(pop,             na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    age       = factor(as.character(age), levels = age_order),
    mx_all    = allcause_deaths / pmax(pop, 1),
    cvd_share = cardio_deaths   / pmax(allcause_deaths, eps)
  ) %>%
  filter(!is.na(age))

proj_lt_both_sexes <- proj_input_both_sexes %>%
  add_life_expectancy(group_vars = c("scenario", "iso3", "year")) %>%
  filter(age == "0-4") %>%
  transmute(
    scenario, iso3, year,
    e0_allcause  = ex,
    e0_noncardio = ex_noncardio,
    e0_gap_cardio = ex_noncardio - ex
  )

saveRDS(proj_lt_both_sexes, file.path(out_dir, "e0_projected_both_sexes.rds"))

rm(proj_input, proj_input_both_sexes, obs_and_proj); invisible(gc())

# =============================================================================
# 4. "Across SSPs" comparison #1: average e0_gap_cardio by SSP and year,
#    across all projected countries -- shows how the SSP scenario itself
#    changes the overall CVD life-expectancy toll over time.
# =============================================================================

gap_by_ssp_year <- proj_lt %>%
  group_by(scenario, sex, year) %>%
  summarise(
    n_countries      = n(),
    mean_gap_cardio  = mean(e0_gap_cardio),
    mean_e0_allcause = mean(e0_allcause),
    .groups = "drop"
  )

saveRDS(gap_by_ssp_year, file.path(out_dir, "e0_gap_by_ssp_year.rds"))

# --- Summary statistics of life expectancy at birth (all-cause and
#     non-cardio) across countries, by SSP x year x sex ---
proj_e0_summary_by_ssp_year <- proj_lt %>%
  group_by(scenario, sex, year) %>%
  summarise(
    n_countries        = n(),
    mean_e0_allcause   = mean(e0_allcause),   sd_e0_allcause   = sd(e0_allcause),
    median_e0_allcause = median(e0_allcause),
    p10_e0_allcause    = quantile(e0_allcause, 0.10), p90_e0_allcause = quantile(e0_allcause, 0.90),
    mean_e0_noncardio   = mean(e0_noncardio),   sd_e0_noncardio   = sd(e0_noncardio),
    median_e0_noncardio = median(e0_noncardio),
    p10_e0_noncardio    = quantile(e0_noncardio, 0.10), p90_e0_noncardio = quantile(e0_noncardio, 0.90),
    mean_gap_cardio     = mean(e0_gap_cardio),  sd_gap_cardio     = sd(e0_gap_cardio),
    .groups = "drop"
  )

saveRDS(proj_e0_summary_by_ssp_year, file.path(out_dir, "e0_projected_summary_by_ssp_year.rds"))

p_gap_ssp <- ggplot(gap_by_ssp_year, aes(x = year, y = mean_gap_cardio, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ sex) +
  scale_colour_manual(values = ssp_colours, name = "SSP") +
  labs(title = "Years of life expectancy lost to CVD, by SSP (mean across countries)",
      x = "Year", y = "e0_gap_cardio (years)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "fig2_gap_by_ssp_over_time.png"), p_gap_ssp,
      width = 10, height = 5.5, dpi = 300)

# =============================================================================
# 5. "Across SSPs" comparison #2: for each country x sex x year, how much
#    does e0_gap_cardio vary across the 5 SSPs? (max - min). This is the
#    "not just within the country, but across SSPs" comparison -- it
#    quantifies, for every country, how sensitive its CVD life-expectancy
#    toll is to which SSP scenario is assumed.
# =============================================================================

gap_ssp_spread <- proj_lt %>%
  group_by(iso3, sex, year) %>%
  summarise(
    min_gap_cardio   = min(e0_gap_cardio),
    max_gap_cardio   = max(e0_gap_cardio),
    ssp_spread_gap   = max(e0_gap_cardio) - min(e0_gap_cardio),   # "across SSPs" difference
    .groups = "drop"
  )

saveRDS(gap_ssp_spread, file.path(out_dir, "e0_gap_ssp_spread_by_country.rds"))

p_gap_spread <- gap_ssp_spread %>%
  filter(year %in% c(2030, 2050, 2080, 2100)) %>%
  ggplot(aes(x = factor(year), y = ssp_spread_gap, fill = sex)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Spread of e0_gap_cardio across the 5 SSPs, by country",
      subtitle = "One point per country: max(e0_gap_cardio across SSPs) - min(...)",
      x = "Year", y = "Cross-SSP spread (years)", fill = "Sex") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "fig3_gap_spread_across_countries.png"), p_gap_spread,
      width = 9, height = 5.5, dpi = 300)

# =============================================================================
# 6. Per-country plots: e0_gap_cardio by SSP over time (sexes pooled)
# =============================================================================

countries_plot <- sort(unique(proj_lt_both_sexes$iso3))

make_gap_plot <- function(iso, out_path) {
  df_i <- proj_lt_both_sexes %>% filter(iso3 == iso) %>% arrange(scenario, year)
  if (nrow(df_i) == 0) return(invisible(NULL))

  country_name <- countrycode(iso, "iso3c", "country.name", warn = FALSE)
  if (is.na(country_name)) country_name <- iso

  p <- ggplot(df_i, aes(x = year, y = e0_gap_cardio, colour = scenario)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.3) +
    scale_colour_manual(values = ssp_colours, name = "SSP") +
    labs(title = country_name, x = "Year", y = "e0_gap_cardio (years)") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom", plot.title = element_text(face = "plain"))

  ggsave(out_path, p, width = 7, height = 5, dpi = 150)
}

for (iso in countries_plot) {
  make_gap_plot(iso, file.path(country_plots_dir, paste0(sanitize_filename(iso), ".png")))
}

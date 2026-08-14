# =============================================================================
# COUNTRY PLOTS
# Updated: 14/08/2026
# Author: Martina Otavova
# =============================================================================
#
#   country_plots/                   -- CVD deaths (count), all ages, both
#                                        sexes combined, one PNG per country
#   country_plots_share/             -- CVD share (%), all ages, both sexes
#                                        combined, one PNG per country
#   country_age_sex_plots/           -- CVD deaths (count), one PNG per
#                                        country x age band x sex
#   country_age_sex_plots_share/     -- CVD share (%), one PNG per
#                                        country x age band x sex
#
# =============================================================================

.libPaths(c("/home/otavova/R/library"))  # point R at the custom library location for package lookups

library(dplyr)       # data wrangling (filter/mutate/group_by/summarise/joins)
library(tidyr)        # data reshaping helpers (not directly used below but loaded for pipeline consistency)
library(ggplot2)      # plotting
library(countrycode)  # convert ISO3 country codes to human-readable country names

out_dir <- "Results/CountryPlots"  # root folder where all plot subfolders will be created

country_plots_dir       <- file.path(out_dir, "country_plots")             # path for country-level death-count plots
country_plots_share_dir <- file.path(out_dir, "country_plots_share")       # path for country-level share plots
age_sex_dir             <- file.path(out_dir, "country_age_sex_plots")     # path for country x age x sex death plots
age_sex_share_dir       <- file.path(out_dir, "country_age_sex_plots_share") # path for country x age x sex share plots
for (d in c(country_plots_dir, country_plots_share_dir, age_sex_dir, age_sex_share_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)  # create each output directory (and parents) if missing, silently
}

age_order <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34",
              "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69",
              "70-74", "75-79", "80-84", "85-89", "90-94", "95+")  # canonical ordering of age bands used to sort facets/filenames

age_to_token <- function(age) {
  if (grepl("\\+$", age)) sub("\\+$", "_", age) else gsub("-", "_", age)  # turn an age band like "95+" -> "95_" or "0-4" -> "0_4" for safe filenames
}

ssp_colors_obs <- c(
  "Observed" = "black",
  "SSP1" = "#1B9E77", "SSP2" = "#D95F02", "SSP3" = "#7570B3",
  "SSP4" = "#E7298A", "SSP5" = "#66A61E"
)  # fixed color mapping so each scenario/series always renders in the same color across plots
legend_order <- c("Observed", "SSP1", "SSP2", "SSP3", "SSP4", "SSP5")  # controls the order series appear in plot legends

sanitize_filename <- function(iso) {
  name <- countrycode(iso, "iso3c", "country.name", warn = FALSE)  # look up the full country name from its ISO3 code
  if (is.na(name) || iso == "TWN") name <- if (iso == "TWN") "Taiwan" else iso  # fall back to "Taiwan" for TWN (not in countrycode's standard mapping) or to the raw code if lookup failed
  gsub("[^A-Za-z0-9]+", "_", name)  # replace any non-alphanumeric characters with underscores to make a safe filename
}

# =============================================================================
# 1. Load observed history + our projections (CVD_projection.R output)
# =============================================================================

obs_and_proj <- readRDS("Results/Projections/obs_and_proj.rds") %>%  # load the combined observed + projected dataset produced upstream
  mutate(iso3 = as.character(iso3), sex = as.character(sex))  # ensure iso3/sex are plain character vectors, not factors

obs_raw <- obs_and_proj %>%
  filter(scenario == "historical") %>%  # keep only the observed/historical rows
  select(iso3, age, sex, year, death = cardio_deaths, allcause_death = allcause_deaths)  # keep and rename the relevant columns

proj_raw <- obs_and_proj %>%
  filter(scenario != "historical") %>%  # keep only the projected (SSP scenario) rows
  select(scenario, iso3, age, sex, year, pop,
        pred_cardio_deaths = cardio_deaths,
        pred_allcause_deaths = allcause_deaths)  # keep and rename the relevant projection columns

rm(obs_and_proj)  # free memory now that obs_raw/proj_raw have been split out

countries_plot <- sort(unique(proj_raw$iso3))  # alphabetically sorted list of all countries that have projections

last_obs_year   <- max(obs_raw$year)   # most recent year with observed data (used to anchor projections visually)
first_proj_year <- min(proj_raw$year)  # first year with projected data (used to anchor projections visually)

# =============================================================================
# 2. Country-level (both sexes combined): deaths + share
# =============================================================================

obs_country <- obs_raw %>%
  group_by(iso3, year) %>%  # collapse across age and sex, keeping only country x year
  summarise(cardio_deaths = sum(death, na.rm = TRUE),
           allcause_deaths = sum(allcause_death, na.rm = TRUE), .groups = "drop") %>%  # sum deaths across all ages/sexes
  mutate(cardio_share_pct = 100 * cardio_deaths / allcause_deaths)  # compute CVD deaths as a % of all-cause deaths

proj_country <- proj_raw %>%
  group_by(iso3, scenario, year) %>%  # collapse across age and sex, keeping country x scenario x year
  summarise(cardio_deaths = sum(pred_cardio_deaths, na.rm = TRUE),
           allcause_deaths = sum(pred_allcause_deaths, na.rm = TRUE), .groups = "drop") %>%  # sum predicted deaths across all ages/sexes
  mutate(cardio_share_pct = 100 * cardio_deaths / allcause_deaths)  # compute predicted CVD share of all-cause deaths

make_country_plot <- function(iso, value_col, y_label, out_path) {
  obs_i  <- obs_country  %>% filter(iso3 == iso) %>% arrange(year)  # observed series for this country, sorted by year
  proj_i <- proj_country %>% filter(iso3 == iso) %>% arrange(scenario, year)  # projected series for this country, sorted by scenario/year
  if (nrow(obs_i) == 0 || nrow(proj_i) == 0) return(invisible(NULL))  # skip plotting if either series is missing

  anchor_obs  <- obs_i[[value_col]][obs_i$year == last_obs_year]  # observed value in the last historical year
  anchor_proj <- mean(proj_i[[value_col]][proj_i$year == first_proj_year], na.rm = TRUE)  # average projected value across scenarios in the first projection year
  offset <- anchor_proj - anchor_obs  # gap between observed and projected at the join point

  obs_plot  <- obs_i  %>% transmute(year, series = "Observed", value = .data[[value_col]] + offset)  # shift observed series by the offset so it visually connects to the projections
  proj_plot <- proj_i %>% transmute(year, series = scenario,   value = .data[[value_col]])  # projected series, labeled by scenario name

  plot_df <- bind_rows(obs_plot, proj_plot) %>% mutate(series = factor(series, levels = legend_order))  # combine observed + projected and fix factor levels for consistent legend order

  country_name <- countrycode(iso, "iso3c", "country.name", warn = FALSE)  # get full country name for the plot title
  if (is.na(country_name)) country_name <- iso  # fall back to the ISO code if the name lookup failed

  proj_plot <- proj_plot %>% mutate(series = factor(series, levels = legend_order))  # apply the same factor levels to the projection-only subset (used for point markers)

  p <- ggplot(plot_df, aes(x = year, y = value, color = series)) +
    geom_line(linewidth = 0.9) +  # draw a line for every series (observed + each scenario)
    geom_point(data = proj_plot, size = 1.3) +  # add point markers only on the projected portion
    scale_color_manual(values = ssp_colors_obs, breaks = legend_order) +  # apply the fixed color scheme and legend ordering
    labs(title = country_name, x = "Year", y = y_label, color = NULL) +  # set title/axis labels and hide the legend title
    theme_bw(base_size = 14) +  # use a clean white-background theme with larger base font
    theme(legend.position = "bottom", plot.title = element_text(face = "plain"))  # put legend below the plot, non-bold title

  ggsave(out_path, p, width = 7, height = 5, dpi = 150)  # write the plot to disk as a PNG
}

for (iso in countries_plot) {
  make_country_plot(iso, "cardio_deaths", "CVD deaths",
                    file.path(country_plots_dir, paste0(sanitize_filename(iso), ".png")))  # save the death-count plot for this country
  make_country_plot(iso, "cardio_share_pct", "CVD share of all-cause deaths (%)",
                    file.path(country_plots_share_dir, paste0(sanitize_filename(iso), ".png")))  # save the CVD-share plot for this country
}

# =============================================================================
# 3. Country x age x sex: deaths
# =============================================================================

obs_as <- obs_raw %>%
  select(iso3, age, sex, year, death) %>%  # keep only the columns needed for age/sex-level plots
  rename(cardio_deaths = death)  # rename for consistency with the projection data below

proj_as <- proj_raw %>%
  select(scenario, iso3, age, sex, year, pred_cardio_deaths) %>%  # keep only the columns needed for age/sex-level plots
  rename(cardio_deaths = pred_cardio_deaths)  # rename for consistency with the observed data above

group_keys_as <- c("iso3", "age", "sex")  # grouping key used repeatedly to join/aggregate at the country x age x sex level

anchor_obs_as <- obs_as %>%
  filter(year == last_obs_year) %>%  # observed value in the last historical year, per group
  select(all_of(group_keys_as), anchor_obs = cardio_deaths)  # keep the group keys plus the anchor value

anchor_proj_as <- proj_as %>%
  filter(year == first_proj_year) %>%  # projected values in the first projection year, per group
  group_by(across(all_of(group_keys_as))) %>%
  summarise(anchor_proj = mean(cardio_deaths, na.rm = TRUE), .groups = "drop")  # average across scenarios to get one anchor value per group

offsets_as <- anchor_obs_as %>%
  inner_join(anchor_proj_as, by = group_keys_as) %>%  # match observed and projected anchors by group
  mutate(offset = anchor_proj - anchor_obs) %>%  # compute the shift needed to connect observed to projected
  select(all_of(group_keys_as), offset)  # keep just the group keys and the offset

obs_shifted_as <- obs_as %>%
  inner_join(offsets_as, by = group_keys_as) %>%  # attach the group-specific offset to every observed row
  mutate(series = "Observed", value = cardio_deaths + offset) %>%  # shift observed values and label the series
  select(all_of(group_keys_as), year, series, value)  # keep only the columns needed for plotting

proj_series_as <- proj_as %>%
  mutate(series = scenario, value = cardio_deaths) %>%  # label projected rows by their scenario name
  select(all_of(group_keys_as), year, series, value)  # keep only the columns needed for plotting

plot_df_as <- bind_rows(obs_shifted_as, proj_series_as) %>%
  mutate(series = factor(series, levels = legend_order))  # combine observed + projected data and fix legend ordering

rm(obs_as, proj_as, anchor_obs_as, anchor_proj_as, offsets_as, obs_shifted_as, proj_series_as); gc()  # free intermediate objects and trigger garbage collection to control memory use

make_as_plot <- function(df_cell, title, y_label, out_path) {
  if (nrow(df_cell) == 0) return(invisible(NULL))  # skip plotting if there's no data for this country/age/sex cell
  proj_pts <- df_cell %>% filter(series != "Observed")  # subset of points to mark (projections only)
  p <- ggplot(df_cell, aes(x = year, y = value, color = series)) +
    geom_line(linewidth = 0.8) +  # draw a line for every series
    geom_point(data = proj_pts, size = 1) +  # add point markers on the projected portion only
    scale_color_manual(values = ssp_colors_obs, breaks = legend_order) +  # apply the fixed color scheme and legend ordering
    labs(title = title, x = "Year", y = y_label, color = NULL) +  # set title/axis labels and hide the legend title
    theme_bw(base_size = 12) +  # clean white-background theme, smaller font (many small plots)
    theme(legend.position = "bottom", plot.title = element_text(size = 11))  # legend below plot, smaller title text
  ggsave(out_path, p, width = 6, height = 4, dpi = 120)  # write the plot to disk as a PNG
}

n_saved <- 0  # counter for number of death plots saved in this section

for (cty in countries_plot) {
  cty_dir <- file.path(age_sex_dir, sanitize_filename(cty))  # per-country subfolder for this country's age/sex plots
  dir.create(cty_dir, showWarnings = FALSE)  # create that subfolder if it doesn't already exist

  df_cty       <- plot_df_as %>% filter(iso3 == cty)  # subset of data for this country
  ages_present <- age_order[age_order %in% unique(as.character(df_cty$age))]  # age bands present for this country, in canonical order

  for (ag in ages_present) {
    for (sx in c("Male", "Female")) {
      df_cell  <- df_cty %>% filter(age == ag, sex == sx)  # data for this specific country x age x sex combination
      out_path <- file.path(cty_dir, paste0(age_to_token(ag), "_", sx, ".png"))  # filename encoding the age band and sex
      make_as_plot(df_cell, paste0(sanitize_filename(cty), ": age ", ag, ", ", sx),
                  "CVD deaths", out_path)  # build and save the plot for this cell
      n_saved <- n_saved + 1  # increment the saved-plot counter
    }
  }
}

rm(plot_df_as); gc()  # free the large combined data frame and garbage collect before the next section

# =============================================================================
# 4. Country x age x sex: share
# =============================================================================

obs_as_share <- obs_raw %>%
  filter(allcause_death > 0) %>%  # avoid divide-by-zero when computing a share
  mutate(cardio_share_pct = 100 * death / allcause_death) %>%  # observed CVD share of all-cause deaths, per age/sex/year
  select(iso3, age, sex, year, cardio_share_pct)  # keep only the columns needed for plotting

proj_as_share <- proj_raw %>%
  filter(pred_allcause_deaths > 0) %>%  # avoid divide-by-zero when computing a share
  mutate(cardio_share_pct = 100 * pred_cardio_deaths / pred_allcause_deaths) %>%  # predicted CVD share of all-cause deaths, per age/sex/year/scenario
  select(scenario, iso3, age, sex, year, cardio_share_pct)  # keep only the columns needed for plotting

anchor_obs_share <- obs_as_share %>%
  filter(year == last_obs_year) %>%  # observed share in the last historical year, per group
  select(all_of(group_keys_as), anchor_obs = cardio_share_pct)  # keep the group keys plus the anchor value

anchor_proj_share <- proj_as_share %>%
  filter(year == first_proj_year) %>%  # projected shares in the first projection year, per group
  group_by(across(all_of(group_keys_as))) %>%
  summarise(anchor_proj = mean(cardio_share_pct, na.rm = TRUE), .groups = "drop")  # average across scenarios to get one anchor value per group

offsets_share <- anchor_obs_share %>%
  inner_join(anchor_proj_share, by = group_keys_as) %>%  # match observed and projected anchors by group
  mutate(offset = anchor_proj - anchor_obs) %>%  # compute the shift needed to connect observed to projected
  select(all_of(group_keys_as), offset)  # keep just the group keys and the offset

obs_shifted_share <- obs_as_share %>%
  inner_join(offsets_share, by = group_keys_as) %>%  # attach the group-specific offset to every observed row
  mutate(series = "Observed", value = cardio_share_pct + offset) %>%  # shift observed values and label the series
  select(all_of(group_keys_as), year, series, value)  # keep only the columns needed for plotting

proj_series_share <- proj_as_share %>%
  mutate(series = scenario, value = cardio_share_pct) %>%  # label projected rows by their scenario name
  select(all_of(group_keys_as), year, series, value)  # keep only the columns needed for plotting

plot_df_share <- bind_rows(obs_shifted_share, proj_series_share) %>%
  mutate(series = factor(series, levels = legend_order))  # combine observed + projected data and fix legend ordering

rm(obs_raw, proj_raw, obs_as_share, proj_as_share, anchor_obs_share, anchor_proj_share,
  offsets_share, obs_shifted_share, proj_series_share, obs_country, proj_country); gc()  # free everything no longer needed and garbage collect

n_saved_share <- 0  # counter for number of share plots saved in this section

for (cty in countries_plot) {
  cty_dir <- file.path(age_sex_share_dir, sanitize_filename(cty))  # per-country subfolder for this country's share plots
  dir.create(cty_dir, showWarnings = FALSE)  # create that subfolder if it doesn't already exist

  df_cty       <- plot_df_share %>% filter(iso3 == cty)  # subset of data for this country
  ages_present <- age_order[age_order %in% unique(as.character(df_cty$age))]  # age bands present for this country, in canonical order

  for (ag in ages_present) {
    for (sx in c("Male", "Female")) {
      df_cell  <- df_cty %>% filter(age == ag, sex == sx)  # data for this specific country x age x sex combination
      out_path <- file.path(cty_dir, paste0(age_to_token(ag), "_", sx, ".png"))  # filename encoding the age band and sex
      make_as_plot(df_cell, paste0(sanitize_filename(cty), ": age ", ag, ", ", sx),
                  "CVD share of all-cause deaths (%)", out_path)  # build and save the plot for this cell
      n_saved_share <- n_saved_share + 1  # increment the saved-plot counter
    }
  }
}

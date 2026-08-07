# =============================================================================
# AGE-STANDARDIZED RATE -- COUNTRY PLOTS
# Updated: 08/07/2026
# Author: Martina Otavova
# =============================================================================
#
#   country_plots/           -- age-standardized CVD mortality rate (per
#                                100,000, WHO World Standard Population,
#                                sexes pooled), one PNG per country
#   country_plots_allcause/  -- age-standardized all-cause mortality rate
#                                (per 100,000, sexes pooled), one PNG per
#                                country
#
# =============================================================================

.libPaths(c("/home/otavova/R/library"))  # point R at the custom library location for package lookups

library(dplyr)        # data wrangling (filter/mutate/group_by/summarise)
library(ggplot2)      # plotting
library(countrycode)  # convert ISO3 country codes to human-readable country names

out_dir <- "Results/AgeStandardizedRate"                            # root folder for age-standardized-rate outputs
country_plots_dir          <- file.path(out_dir, "country_plots")            # path for the per-country CVD ASMR plots
country_plots_allcause_dir <- file.path(out_dir, "country_plots_allcause")    # path for the per-country all-cause ASMR plots
for (d in c(country_plots_dir, country_plots_allcause_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)  # create each output directory if missing
}

ssp_colors_obs <- c(
  "Observed" = "black",
  "SSP1" = "#1B9E77", "SSP2" = "#D95F02", "SSP3" = "#7570B3",
  "SSP4" = "#E7298A", "SSP5" = "#66A61E"
)  # fixed color mapping so each scenario/series always renders in the same color across plots
legend_order <- c("Observed", "SSP1", "SSP2", "SSP3", "SSP4", "SSP5")  # controls the order series appear in plot legends

sanitize_filename <- function(iso) {
  name <- countrycode(iso, "iso3c", "country.name", warn = FALSE)  # look up the full country name from its ISO3 code
  if (is.na(name) || iso == "TWN") name <- if (iso == "TWN") "Taiwan" else iso  # fall back to "Taiwan" for TWN or to the raw code if lookup failed
  gsub("[^A-Za-z0-9]+", "_", name)  # replace any non-alphanumeric characters with underscores to make a safe filename
}

# =============================================================================
# 1. Load age-standardized rates (AgeStandardizedRate.R output, sexes pooled)
# =============================================================================

asmr <- readRDS(file.path(out_dir, "asmr_both_sexes.rds")) %>%  # load the pooled-sex age-standardized rate table (CVD + All-cause)
  mutate(iso3 = as.character(iso3))  # ensure iso3 is a plain character vector, not a factor

causes <- list(
  CVD        = list(dir = country_plots_dir,          y_label = "CVD ASMR (per 100,000)"),
  `All-cause` = list(dir = country_plots_allcause_dir, y_label = "ASMR (per 100,000)")
)  # per-cause output directory and y-axis label

# =============================================================================
# 2. Country-level age-standardized rate plot
# =============================================================================

make_asmr_plot <- function(iso, obs_cause, proj_cause, y_label, last_obs_year, first_proj_year, out_path) {
  obs_i  <- obs_cause  %>% filter(iso3 == iso) %>% arrange(year)              # observed series for this country
  proj_i <- proj_cause %>% filter(iso3 == iso) %>% arrange(scenario, year)    # projected series for this country
  if (nrow(obs_i) == 0 || nrow(proj_i) == 0) return(invisible(NULL))         # skip plotting if either series is missing

  anchor_obs  <- obs_i$asmr_per_100k[obs_i$year == last_obs_year]                              # observed rate in the last historical year
  anchor_proj <- mean(proj_i$asmr_per_100k[proj_i$year == first_proj_year], na.rm = TRUE)       # average projected rate across scenarios in the first projection year
  offset <- anchor_proj - anchor_obs                                                             # gap between observed and projected at the join point

  obs_plot  <- obs_i  %>% transmute(year, series = "Observed", value = asmr_per_100k + offset)  # shift observed series so it visually connects to the projections
  proj_plot <- proj_i %>% transmute(year, series = scenario,   value = asmr_per_100k)             # projected series, labeled by scenario name

  plot_df <- bind_rows(obs_plot, proj_plot) %>% mutate(series = factor(series, levels = legend_order))  # combine and fix factor levels for consistent legend order

  country_name <- countrycode(iso, "iso3c", "country.name", warn = FALSE)  # get full country name for the plot title
  if (is.na(country_name)) country_name <- iso  # fall back to the ISO code if the name lookup failed

  proj_plot <- proj_plot %>% mutate(series = factor(series, levels = legend_order))  # apply the same factor levels to the projection-only subset (used for point markers)

  p <- ggplot(plot_df, aes(x = year, y = value, color = series)) +
    geom_line(linewidth = 0.9) +  # draw a line for every series (observed + each scenario)
    geom_point(data = proj_plot, size = 1.3) +  # add point markers only on the projected portion
    scale_color_manual(values = ssp_colors_obs, breaks = legend_order) +  # apply the fixed color scheme and legend ordering
    labs(title = country_name, x = "Year", y = y_label, color = NULL) +  # set title/axis labels and hide the legend title
    theme_bw(base_size = 14) +  # clean white-background theme with larger base font
    theme(legend.position = "bottom", plot.title = element_text(face = "plain"))  # legend below plot, non-bold title

  ggsave(out_path, p, width = 7, height = 5, dpi = 150)  # write the plot to disk as a PNG
}

for (cause_label in names(causes)) {
  cause_cfg <- causes[[cause_label]]  # this cause's output directory + y-axis label

  asmr_cause <- asmr %>% filter(cause == cause_label)  # restrict to this cause (CVD or All-cause)
  obs_cause  <- asmr_cause %>% filter(scenario == "historical") %>% arrange(iso3, year)  # observed portion of the series
  proj_cause <- asmr_cause %>% filter(scenario != "historical")  %>% arrange(iso3, scenario, year)  # projected portion of the series

  countries_plot <- sort(unique(proj_cause$iso3))  # alphabetically sorted list of all countries that have projections
  last_obs_year   <- max(obs_cause$year)   # most recent year with observed data (used to anchor projections visually)
  first_proj_year <- min(proj_cause$year)  # first year with projected data (used to anchor projections visually)

  message("\n=== ", cause_label, ": ", length(countries_plot), " countries ===")

  for (iso in countries_plot) {
    make_asmr_plot(iso, obs_cause, proj_cause, cause_cfg$y_label, last_obs_year, first_proj_year,
                  file.path(cause_cfg$dir, paste0(sanitize_filename(iso), ".png")))  # save the ASMR plot for this country
  }

  cat("Saved", length(countries_plot), cause_label, "age-standardized rate country plots to", cause_cfg$dir, "\n")
}

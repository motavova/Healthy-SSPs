# =============================================================================
# Regional Trajectory Plots
# Updated: 17/08/2026
# Author: Martina Otavova
#
# Adds a time-series view (fig5) to the regional CVD-share figures: cardio
# deaths and cardio share by region across 1990-2100, historical + all SSPs
# overlaid, one panel per region. Complements the cross-sectional fig3/fig4
# bar charts already produced by RegionalSummaryPlots.R.
# =============================================================================

setwd("/home/otavova/Healthy-SSPs/SSPs_CVD_mortality")

.libPaths(c("/home/otavova/R/library"))

library(dplyr)
library(ggplot2)
library(scales)

summary_dir <- "Results/RegionalSummaryMeasures"
plot_dir    <- "Results/Plots/RegionalSummaries"

ssp_colours <- c("historical" = "grey35", "SSP1" = "#1b9e77", "SSP2" = "#d95f02",
                 "SSP3" = "#7570b3", "SSP4" = "#e7298a", "SSP5" = "#66a61e")

# WorldBankUpdate = World Bank's current country classification (with a
# High-Income bucket pulled out of geography). This used to also loop over a
# second "WorldBank" scheme (countrycode's region field, no High-Income
# carve-out), but that scheme's country membership was stale relative to the
# World Bank's current classification and has been retired; its last output
# is archived at Results/History/RegionalSummaries_old_worldbank/.
schemes <- list(
  WorldBankUpdate = list(
    summary_rds        = file.path(summary_dir, "regional_summary_worldbankupdate.rds"),
    out_dir             = file.path(plot_dir, "WorldBankUpdate"),
    label               = "World Bank regions",
    show_label_in_title = FALSE,        # titles stay generic, no "(World Bank regions)" suffix
    proj_start_year     = 2030,           # SSP lines are only drawn from this year onward
    show_historical     = TRUE,            # include the observed/historical line
    show_points         = FALSE,
    file_suffix         = "",
    recode        = c("High-Income" = "High-Income Countries"),
    region_levels = c(
      "East Asia & Pacific", "Europe & Central Asia", "High-Income Countries",
      "Latin America & Caribbean", "Middle East, North Africa, Afghanistan & Pakistan",
      "South Asia", "Sub-Saharan Africa"
    )
  ),
  WorldBankUpdate_full = list(
    summary_rds        = file.path(summary_dir, "regional_summary_worldbankupdate.rds"),
    out_dir             = file.path(plot_dir, "WorldBankUpdate"),
    label               = "World Bank regions",
    show_label_in_title = FALSE,
    proj_start_year     = 2025,          # earliest SSP year in the data (there is no 2024 row)
    show_historical     = TRUE,
    show_points         = FALSE,
    file_suffix         = "_full",         # keeps these separate from the 2030-onward fig5a/fig5b
    recode        = c("High-Income" = "High-Income Countries"),
    region_levels = c(
      "East Asia & Pacific", "Europe & Central Asia", "High-Income Countries",
      "Latin America & Caribbean", "Middle East, North Africa, Afghanistan & Pakistan",
      "South Asia", "Sub-Saharan Africa"
    )
  ),
  WorldBankUpdate_proj_dots = list(
    summary_rds        = file.path(summary_dir, "regional_summary_worldbankupdate.rds"),
    out_dir             = file.path(plot_dir, "WorldBankUpdate"),
    label               = "World Bank regions",
    show_label_in_title = FALSE,
    proj_start_year     = 2030,          # projections only, from 2030 onward
    show_historical     = FALSE,           # no observed/historical line
    show_points         = TRUE,            # dot marker on every 5-year data point
    file_suffix         = "_dots",
    recode        = c("High-Income" = "High-Income Countries"),
    region_levels = c(
      "East Asia & Pacific", "Europe & Central Asia", "High-Income Countries",
      "Latin America & Caribbean", "Middle East, North Africa, Afghanistan & Pakistan",
      "South Asia", "Sub-Saharan Africa"
    )
  )
)

for (scheme in schemes) {

  dir.create(scheme$out_dir, showWarnings = FALSE, recursive = TRUE)

  regional_summary <- readRDS(scheme$summary_rds)
  n_regions <- n_distinct(regional_summary$region)

  recode_region <- function(region) {
    if (length(scheme$recode) == 0) return(region)
    ifelse(region %in% names(scheme$recode), scheme$recode[region], region)
  }

  traj_data <- regional_summary %>%
    filter(scenario != "historical" | scheme$show_historical) %>%           # drop the observed line entirely when requested
    filter(scenario == "historical" | year >= scheme$proj_start_year) %>%   # drop SSP years before the requested projection start
    mutate(
      region = factor(recode_region(region), levels = scheme$region_levels),
      ssp    = factor(scenario, levels = c("historical", "SSP1", "SSP2", "SSP3", "SSP4", "SSP5"))
    )

  title_deaths <- if (scheme$show_label_in_title) {
    paste0("Cardiovascular deaths over time by region (", scheme$label, ")")
  } else {
    "Cardiovascular deaths over time by region"
  }
  title_share <- if (scheme$show_label_in_title) {
    paste0("Cardiovascular share of all-cause deaths over time by region (", scheme$label, ")")
  } else {
    "Cardiovascular share of all-cause deaths over time by region"
  }
  subtitle_text <- if (scheme$show_historical) {
    paste0("Historical (1990-2023) and SSP projections (", scheme$proj_start_year,
          "-2100); dashed line marks the projection start")
  } else {
    paste0("SSP projections, ", scheme$proj_start_year, "-2100")
  }
  x_breaks <- if (scheme$show_historical) seq(1990, 2100, by = 20) else seq(scheme$proj_start_year, 2100, by = 20)

  # ===========================================================================
  # fig5a: cardio deaths (millions) over time, by region -- historical + SSPs
  # ===========================================================================

  p_deaths <- ggplot(traj_data, aes(x = year, y = cardio_deaths_millions, colour = ssp)) +
    { if (scheme$show_historical) geom_vline(xintercept = scheme$proj_start_year, linetype = "dashed", colour = "grey70", linewidth = 0.4) } +
    geom_line(linewidth = 0.7) +
    { if (scheme$show_points) geom_point(size = 1.8) } +
    facet_wrap(~ region, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = ssp_colours,
                        breaks = c("historical", "SSP1", "SSP2", "SSP3", "SSP4", "SSP5"),
                        labels = c("Historical", "SSP1", "SSP2", "SSP3", "SSP4", "SSP5"),
                        name = NULL) +
    scale_x_continuous(breaks = x_breaks) +
    labs(title = title_deaths,
        subtitle = subtitle_text,
        x = NULL, y = "Cardiovascular deaths (millions)") +
    theme_bw(base_size = 16) +
    theme(
      plot.title         = element_text(face = "bold", size = 20),
      plot.subtitle      = element_text(size = 15),
      strip.background   = element_rect(fill = "grey85"),
      strip.text         = element_text(face = "bold", size = 15),
      axis.title         = element_text(size = 16),
      axis.text          = element_text(size = 13),
      legend.title       = element_text(size = 15),
      legend.text        = element_text(size = 14),
      panel.grid.minor   = element_blank(),
      legend.position    = "bottom"
    )

  deaths_stem   <- paste0("fig5a_trajectory_deaths_by_region", scheme$file_suffix)
  deaths_width  <- 12
  deaths_height <- ceiling(n_regions / 3) * 3 + 2

  ggsave(file.path(scheme$out_dir, paste0(deaths_stem, ".png")), p_deaths,
        width = deaths_width, height = deaths_height, dpi = 150, limitsize = FALSE)
  ggsave(file.path(scheme$out_dir, paste0(deaths_stem, ".svg")), p_deaths,
        width = deaths_width, height = deaths_height, device = grDevices::svg, limitsize = FALSE)

  # ===========================================================================
  # fig5b: cardio share of all-cause deaths over time, by region
  # ===========================================================================

  p_share <- ggplot(traj_data, aes(x = year, y = cardio_share, colour = ssp)) +
    { if (scheme$show_historical) geom_vline(xintercept = scheme$proj_start_year, linetype = "dashed", colour = "grey70", linewidth = 0.4) } +
    geom_line(linewidth = 0.7) +
    { if (scheme$show_points) geom_point(size = 1.8) } +
    facet_wrap(~ region, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = ssp_colours,
                        breaks = c("historical", "SSP1", "SSP2", "SSP3", "SSP4", "SSP5"),
                        labels = c("Historical", "SSP1", "SSP2", "SSP3", "SSP4", "SSP5"),
                        name = NULL) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(title = title_share,
        subtitle = subtitle_text,
        x = NULL, y = "Cardiovascular share of all-cause deaths") +
    theme_bw(base_size = 16) +
    theme(
      plot.title         = element_text(face = "bold", size = 20),
      plot.subtitle      = element_text(size = 15),
      strip.background   = element_rect(fill = "grey85"),
      strip.text         = element_text(face = "bold", size = 15),
      axis.title         = element_text(size = 16),
      axis.text          = element_text(size = 13),
      legend.title       = element_text(size = 15),
      legend.text        = element_text(size = 14),
      panel.grid.minor   = element_blank(),
      legend.position    = "bottom"
    )

  share_stem   <- paste0("fig5b_trajectory_share_by_region", scheme$file_suffix)
  share_width  <- 12
  share_height <- ceiling(n_regions / 3) * 3 + 2

  ggsave(file.path(scheme$out_dir, paste0(share_stem, ".png")), p_share,
        width = share_width, height = share_height, dpi = 150, limitsize = FALSE)
  ggsave(file.path(scheme$out_dir, paste0(share_stem, ".svg")), p_share,
        width = share_width, height = share_height, device = grDevices::svg, limitsize = FALSE)
}

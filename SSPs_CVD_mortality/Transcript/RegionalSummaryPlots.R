# =============================================================================
# Regional Summary Plots
# Updated: 14/08/2026
# Author: Martina Otavova
#
# Reproduces the regional CVD-share figures (sex x region x SSP, and
# region x SSP by year) for both regional divisions produced by
# RegionalSummaryMeasures.R -- Sellers (2020) and World Bank -- saving each
# set of figures into its own subfolder.
# =============================================================================

setwd("/home/otavova/Healthy-SSPs/SSPs_CVD_mortality")

.libPaths(c("/home/otavova/R/library"))

library(dplyr)
library(ggplot2)
library(scales)

summary_dir <- "Results/RegionalSummaryMeasures"
plot_dir    <- "Results/Plots/RegionalSummaries"

ssp_colours <- c("SSP1" = "#1b9e77", "SSP2" = "#d95f02", "SSP3" = "#7570b3",
                "SSP4" = "#e7298a", "SSP5" = "#66a61e")

schemes <- list(
  Sellers = list(
    summary_rds     = file.path(summary_dir, "regional_summary_sellers.rds"),
    summary_sex_rds = file.path(summary_dir, "regional_summary_sellers_sex.rds"),
    out_dir         = file.path(plot_dir, "Sellers"),
    label           = "Sellers regions",
    recode          = c("High-Income" = "High-Income Countries"),
    region_levels   = c(
      "High-Income Countries", "Sub-Saharan Africa", "South Asia",
      "Middle East & North Africa", "Latin America & Caribbean",
      "Europe & Central Asia", "East Asia & Pacific"
    )
  ),
  WorldBank = list(
    summary_rds     = file.path(summary_dir, "regional_summary_worldbank.rds"),
    summary_sex_rds = file.path(summary_dir, "regional_summary_worldbank_sex.rds"),
    out_dir         = file.path(plot_dir, "WorldBank"),
    label           = "World Bank regions",
    recode          = c(),
    region_levels   = c(
      "North America", "Sub-Saharan Africa", "South Asia",
      "Middle East & North Africa", "Latin America & Caribbean",
      "Europe & Central Asia", "East Asia & Pacific"
    )
  )
)

for (scheme in schemes) {

  dir.create(scheme$out_dir, showWarnings = FALSE, recursive = TRUE)

  regional_summary     <- readRDS(scheme$summary_rds)
  regional_summary_sex <- readRDS(scheme$summary_sex_rds)
  n_regions <- n_distinct(regional_summary$region)

  recode_region <- function(region) {
    if (length(scheme$recode) == 0) return(region)
    ifelse(region %in% names(scheme$recode), scheme$recode[region], region)
  }

  # ===========================================================================
  # fig3: CVD share by sex, region and SSP -- 2060 and 2100 side by side
  # ===========================================================================

  make_sex_share_data <- function(plot_year) {
    regional_summary_sex %>%
      filter(year == plot_year) %>%
      mutate(
        region   = factor(recode_region(region), levels = scheme$region_levels),
        ssp      = factor(scenario, levels = c("SSP5", "SSP4", "SSP3", "SSP2", "SSP1")),
        sex_year = paste0(sex, " ", plot_year)
      )
  }

  sex_share_data <- bind_rows(make_sex_share_data(2060), make_sex_share_data(2100)) %>%
    mutate(sex_year = factor(sex_year, levels = c("Female 2060", "Male 2060", "Female 2100", "Male 2100")))

  p_sex_share <- ggplot(sex_share_data, aes(x = region, y = cardio_share, fill = ssp)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    coord_flip() +
    facet_wrap(~ sex_year, nrow = 2) +
    scale_fill_manual(values = ssp_colours, breaks = c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5"),
                      labels = c("1", "2", "3", "4", "5"), name = "SSP") +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Cardiovascular share of all-cause deaths by sex, region and SSP",
        x = NULL, y = "Proportion of deaths") +
    theme_bw(base_size = 11) +
    theme(
      plot.title         = element_text(face = "bold", size = 12),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "right"
    )

  ggsave(file.path(scheme$out_dir, "fig3_sex_cvd_share_2060_2100.png"), p_sex_share,
        width = 13, height = ceiling(n_regions * 0.4) + 4, dpi = 150, limitsize = FALSE)

  # ===========================================================================
  # fig4: CVD share by region and SSP, both sexes combined
  #   -- 2030, 2050, 2080, 2100 individually, and 2060|2100 combined
  # ===========================================================================

  years <- c(2030, 2050, 2080, 2100)

  shared_max <- regional_summary %>%
    filter(year %in% years) %>%
    pull(cardio_share) %>%
    max(na.rm = TRUE)

  make_share_data <- function(plot_year) {
    regional_summary %>%
      filter(year == plot_year) %>%
      mutate(
        region = factor(recode_region(region), levels = scheme$region_levels),
        ssp    = factor(scenario, levels = c("SSP5", "SSP4", "SSP3", "SSP2", "SSP1"))
      )
  }

  for (plot_year in years) {
    year_data <- make_share_data(plot_year)

    p <- ggplot(year_data, aes(x = region, y = cardio_share, fill = ssp)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      coord_flip() +
      scale_fill_manual(values = ssp_colours, breaks = c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5"),
                        labels = c("1", "2", "3", "4", "5"), name = "SSP") +
      scale_y_continuous(labels = percent_format(accuracy = 1),
                        limits = c(0, shared_max * 1.05),
                        expand = expansion(mult = c(0, 0.02))) +
      labs(title = paste0("Cardiovascular share of all-cause deaths by region and SSP — ", plot_year),
          x = NULL, y = "Proportion of deaths") +
      theme_bw(base_size = 12) +
      theme(
        plot.title         = element_text(face = "bold", size = 12),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "right"
      )

    ggsave(file.path(scheme$out_dir, paste0("fig4_", plot_year, ".png")), p,
          width = 10, height = ceiling(n_regions * 0.4) + 2, dpi = 150, limitsize = FALSE)
  }

  combined_share_data <- bind_rows(make_share_data(2060), make_share_data(2100)) %>%
    mutate(year_label = as.character(year))

  p_combined_share <- ggplot(combined_share_data, aes(x = region, y = cardio_share, fill = ssp)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    coord_flip() +
    facet_wrap(~ year_label, ncol = 2) +
    scale_fill_manual(values = ssp_colours, breaks = c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5"),
                      labels = c("1", "2", "3", "4", "5"), name = "SSP") +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Cardiovascular share of all-cause deaths by region and SSP",
        x = NULL, y = "Proportion of deaths") +
    theme_bw(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 12),
      strip.background   = element_rect(fill = "grey85"),
      strip.text         = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "right"
    )

  ggsave(file.path(scheme$out_dir, "fig4_combined_2060_2100.png"), p_combined_share,
        width = 13, height = ceiling(n_regions * 0.4) + 2, dpi = 150, limitsize = FALSE)
}

# =============================================================================
# Sex Deaths Trajectory Plot
# Adds a time-series view of cardiovascular deaths by sex, region, and SSP
# (2025-2100), so that crossovers between male and female absolute death
# counts -- invisible in a fixed-year snapshot -- become visible directly.
# =============================================================================

setwd("/home/otavova/Healthy-SSPs/SSPs_CVD_mortality")

.libPaths(c("/home/otavova/R/library"))

library(dplyr)
library(ggplot2)

summary_dir <- "Results/RegionalSummaryMeasures"
plot_dir    <- "Results/Plots/RegionalSummaries/WorldBankUpdate"

ssp_colours <- c("SSP1" = "#1b9e77", "SSP2" = "#d95f02", "SSP3" = "#7570b3",
                 "SSP4" = "#e7298a", "SSP5" = "#66a61e")

region_levels <- c(
  "East Asia & Pacific", "Europe & Central Asia", "High-Income Countries",
  "Latin America & Caribbean", "Middle East, North Africa, Afghanistan & Pakistan",
  "South Asia", "Sub-Saharan Africa"
)

rs_sex <- readRDS(file.path(summary_dir, "regional_summary_worldbankupdate_sex.rds")) %>%
  filter(scenario != "historical", year >= 2030) %>%
  mutate(
    region = recode(region, "High-Income" = "High-Income Countries"),
    region = factor(region, levels = region_levels),
    ssp    = factor(scenario, levels = c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5")),
    cardio_deaths_millions = cardio_deaths / 1e6
  )

n_regions <- n_distinct(rs_sex$region)

p <- ggplot(rs_sex, aes(x = year, y = cardio_deaths_millions, colour = ssp, linetype = sex)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ region, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = ssp_colours, name = "SSP") +
  scale_linetype_manual(values = c("Male" = "solid", "Female" = "dashed"), name = "Sex") +
  labs(title = "Cardiovascular deaths by sex over time, by region and SSP",
      subtitle = "Solid = Male, dashed = Female. SSP projections, 2030-2100.",
      x = NULL, y = "Cardiovascular deaths (millions)") +
  theme_bw(base_size = 16) +
  theme(
    plot.title         = element_text(face = "bold", size = 20),
    plot.subtitle      = element_text(size = 15),
    strip.background   = element_rect(fill = "grey85"),
    strip.text         = element_text(face = "bold", size = 15),
    axis.title         = element_text(size = 16),
    axis.text          = element_text(size = 13),
    legend.title        = element_text(size = 15),
    legend.text        = element_text(size = 14),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom"
  )

ggsave(file.path(plot_dir, "fig6_sex_deaths_trajectory_by_region.png"), p,
      width = 14, height = ceiling(n_regions / 3) * 3.4 + 2.5, dpi = 150, limitsize = FALSE)

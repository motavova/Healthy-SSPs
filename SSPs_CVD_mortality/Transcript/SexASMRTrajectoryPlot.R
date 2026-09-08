# =============================================================================
# Sex Age-Standardized Mortality Rate (ASMR) Trajectory Plot
#
# Companion to fig6 (SexDeathsTrajectoryPlot.R), but using age-standardized
# cardiovascular mortality rates by sex instead of raw death counts. This
# removes the population-structure confound (e.g. differing numbers of
# elderly men vs women) so that any remaining male-female gap reflects a
# genuine difference in age-specific risk, not who happens to be alive.
#
# Every region and every year/SSP is standardized to the SAME reference
# age structure: the WHO World Standard Population (2000-2025), the
# convention used throughout WHO and Global Burden of Disease
# cardiovascular epidemiology, so these ASMRs are comparable not only
# across the regions/SSPs in this report but against externally published
# age-standardized CVD mortality rates.
#
# Source: Ahmad OB, Boschi-Pinto C, Lopez AD, Murray CJ, Lozano R, Inoue M.
# "Age standardization of rates: a new WHO standard." GPE Discussion Paper
# Series No. 31, World Health Organization, 2001. Weights below are the
# standard's published 5-year age-band values (per 100,000); the source
# table's terminal 95-99 and 100+ groups are summed into a single "95+"
# bucket to match this project's data (which tops out at an open 95+
# interval). These are widely reproduced values, but were transcribed from
# memory rather than the source PDF -- spot-check against the primary
# citation above before using this for a publication-facing figure.
# =============================================================================
who_std <- tibble::tribble(
  ~age,    ~who_weight,
  "0-4",       8860,
  "5-9",       8690,
  "10-14",     8600,
  "15-19",     8470,
  "20-24",     8220,
  "25-29",     7930,
  "30-34",     7610,
  "35-39",     7150,
  "40-44",     6590,
  "45-49",     6040,
  "50-54",     5370,
  "55-59",     4550,
  "60-64",     3720,
  "65-69",     2960,
  "70-74",     2210,
  "75-79",     1520,
  "80-84",      910,
  "85-89",      440,
  "90-94",      150,
  "95+",         45   # = 40 (95-99) + 5 (100+) in the original WHO table
)

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

rates <- readRDS(file.path(summary_dir, "age_sex_region_mortality_rates.rds")) %>%
  mutate(region = recode(region, "High-Income" = "High-Income Countries"))

# Standard age structure: ONE global reference, pooled (male + female)
# population across all regions combined, year 2030, SSP2 -- held fixed so
# it never reflects a particular scenario's future demographic assumptions,
# and applied identically to every region so regional magnitudes are
# directly comparable, not just each region's own male-female gap.
std_weights <- who_std %>%
  mutate(w = who_weight / sum(who_weight)) %>%   # normalize to sum to 1, self-correcting for any rounding in the source table
  select(age, w)

asmr <- rates %>%
  filter(scenario != "historical", year >= 2030) %>%
  left_join(std_weights, by = "age") %>%
  group_by(region, scenario, year, sex) %>%
  summarise(asmr_100k = sum(cardio_mortality_rate * w, na.rm = TRUE) * 1e5, .groups = "drop") %>%
  mutate(
    region = factor(region, levels = region_levels),
    ssp    = factor(scenario, levels = c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5"))
  )

n_regions <- n_distinct(asmr$region)

p <- ggplot(asmr, aes(x = year, y = asmr_100k, colour = ssp, linetype = sex)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ region, ncol = 3) +   # shared y-axis, so regional magnitudes are directly comparable
  scale_colour_manual(values = ssp_colours, name = "SSP") +
  scale_linetype_manual(values = c("Male" = "solid", "Female" = "dashed"), name = "Sex") +
  labs(title = "Age-standardized cardiovascular mortality rate by sex, region, and SSP",
      subtitle = "Solid = Male, dashed = Female. Standardized to the WHO World Standard Population (2000-2025). SSP projections, 2030-2100.",
      x = NULL, y = "Age-standardized cardiovascular deaths per 100,000") +
  theme_bw(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    strip.background   = element_rect(fill = "grey85"),
    strip.text         = element_text(face = "bold"),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom"
  )

ggsave(file.path(plot_dir, "fig7_sex_asmr_trajectory_by_region.png"), p,
      width = 14, height = ceiling(n_regions / 3) * 3.4 + 2.5, dpi = 150, limitsize = FALSE)

cat("Saved fig7_sex_asmr_trajectory_by_region.png\n")

# =============================================================================
# Global Figures
# Updated: 14/08/2026
# Author: Martina Otavova
# =============================================================================

.libPaths(c("/home/otavova/R/library"))

library(dplyr)
library(ggplot2)
library(scales)

out_dir <- "Results/GlobalFigures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ssp_colours <- c("SSP1" = "#1b9e77", "SSP2" = "#d95f02", "SSP3" = "#7570b3",
                "SSP4" = "#e7298a", "SSP5" = "#66a61e")

source("Transcript/Sellers_regions.R")

region_lookup <- bind_rows(
  sellers_regions,
  tibble(iso3 = c("CMR", "KIR"), region = c("Sub-Saharan Africa", "East Asia & Pacific"))
)

region_levels <- c(
  "High-Income Countries", "Sub-Saharan Africa", "South Asia",
  "Middle East & North Africa", "Latin America & Caribbean",
  "Europe & Central Asia", "East Asia & Pacific"
)

# =============================================================================
# 1. Load Samir 5-year recalculation, attach region
# =============================================================================

proj_samir <- readRDS("Results/53_Recalc_CVD_Deaths_SamirMx/proj_deaths_samir_5yr.rds") %>%
  mutate(iso3 = as.character(iso3), sex = as.character(sex)) %>%
  select(scenario, iso3, age, sex, year, pop,
        pred_cardio_deaths = cardio_deaths_samir,
        pred_allcause_deaths = allcause_deaths_samir) %>%
  left_join(region_lookup, by = "iso3")

# =============================================================================
# 2. Global and regional summaries (pooled-sex, and by sex)
# =============================================================================

global_samir5yr <- proj_samir %>%
  group_by(scenario, year) %>%
  summarise(
    cardio_deaths   = sum(pred_cardio_deaths,   na.rm = TRUE),
    allcause_deaths = sum(pred_allcause_deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cardio_share           = cardio_deaths / allcause_deaths,
    cardio_deaths_millions = cardio_deaths / 1e6
  )

saveRDS(global_samir5yr, file.path(out_dir, "global_summary_samir5yr.rds"))

regional_summary <- proj_samir %>%
  filter(!is.na(region)) %>%
  group_by(scenario, region, year) %>%
  summarise(
    cardio_deaths   = sum(pred_cardio_deaths,   na.rm = TRUE),
    allcause_deaths = sum(pred_allcause_deaths, na.rm = TRUE),
    n_countries     = n_distinct(iso3),
    .groups = "drop"
  ) %>%
  mutate(
    cardio_share           = cardio_deaths / allcause_deaths,
    cardio_deaths_millions = cardio_deaths / 1e6
  )

regional_summary_sex <- proj_samir %>%
  filter(!is.na(region)) %>%
  group_by(scenario, region, sex, year) %>%
  summarise(
    cardio_deaths   = sum(pred_cardio_deaths,   na.rm = TRUE),
    allcause_deaths = sum(pred_allcause_deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cardio_share = cardio_deaths / allcause_deaths)

saveRDS(regional_summary,     file.path(out_dir, "regional_summary_samir5yr.rds"))
saveRDS(regional_summary_sex, file.path(out_dir, "regional_summary_sex_samir5yr.rds"))

n_regions <- n_distinct(regional_summary$region)

# =============================================================================
# 3. fig_global: global CVD deaths by SSP (single panel, script 38's style)
# =============================================================================

p_deaths <- ggplot(global_samir5yr, aes(x = year, y = cardio_deaths_millions, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = ssp_colours, name = "SSP") +
  scale_x_continuous(breaks = c(2030, 2050, 2070, 2100)) +
  labs(title = "Global CVD deaths: Samir/WIC mx recalculation (native 5-year steps)",
      x = NULL, y = "CVD deaths (millions)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "fig_global_ssp_deaths_samir5yr.png"), p_deaths, width = 10, height = 6, dpi = 300)

# =============================================================================
# 4. fig3: CVD share by sex, region and SSP -- 2060 and 2100 side by side
# =============================================================================

make_sex_share_data <- function(plot_year) {
  regional_summary_sex %>%
    filter(year == plot_year) %>%
    mutate(
      region = ifelse(region == "High-Income", "High-Income Countries", region),
      region = factor(region, levels = region_levels),
      ssp    = factor(scenario, levels = c("SSP5", "SSP4", "SSP3", "SSP2", "SSP1")),
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
  labs(title = "Cardiovascular share of all-cause deaths by sex, region and SSP (Samir/WIC mx, 5yr)",
      x = NULL, y = "Proportion of deaths") +
  theme_bw(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right"
  )

ggsave(file.path(out_dir, "fig3_sex_cvd_share_2060_2100.png"), p_sex_share,
      width = 13, height = ceiling(n_regions * 0.4) + 4, dpi = 150, limitsize = FALSE)

# =============================================================================
# 5. fig4: CVD share by region and SSP, both sexes combined
#    -- 2030, 2050, 2080, 2100 individually, and 2060|2100 combined
# =============================================================================

years <- c(2030, 2050, 2080, 2100)

shared_max <- regional_summary %>%
  filter(year %in% years) %>%
  pull(cardio_share) %>%
  max(na.rm = TRUE)

make_share_data <- function(plot_year) {
  regional_summary %>%
    filter(year == plot_year) %>%
    mutate(
      region = ifelse(region == "High-Income", "High-Income Countries", region),
      region = factor(region, levels = region_levels),
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
    labs(title = paste0("Cardiovascular share of all-cause deaths by region and SSP (Samir/WIC mx, 5yr) — ", plot_year),
        x = NULL, y = "Proportion of deaths") +
    theme_bw(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 12),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "right"
    )

  ggsave(file.path(out_dir, paste0("fig4_", plot_year, ".png")), p,
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
  labs(title = "Cardiovascular share of all-cause deaths by region and SSP (Samir/WIC mx, 5yr)",
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

ggsave(file.path(out_dir, "fig4_combined_2060_2100.png"), p_combined_share,
      width = 13, height = ceiling(n_regions * 0.4) + 2, dpi = 150, limitsize = FALSE)

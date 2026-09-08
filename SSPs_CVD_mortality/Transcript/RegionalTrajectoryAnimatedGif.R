# =============================================================================
# Regional CVD Deaths Trajectory -- Animated GIF (SSPs reveal one at a time,
# across all region panels simultaneously)
#
# Companion animation to fig5a_trajectory_deaths_by_region_dots.png: the
# same 7-region faceted panel layout, projections only (2030-2100), dots at
# each 5-year point -- but revealed one SSP at a time. At each animation
# step, every region panel shows the same set of completed SSPs (fully
# drawn) plus the current SSP drawing in, in parallel across all panels.
#
# Panel order matches the "three trajectory shapes" narrative: already
# peaked (Europe & Central Asia); peaks mid-century (East Asia & Pacific,
# High-Income Countries, Latin America & Caribbean, South Asia); still
# rising at 2100 under the less favorable SSPs (Middle East, North Africa,
# Afghanistan & Pakistan; Sub-Saharan Africa).
# =============================================================================

setwd("/home/otavova/Healthy-SSPs/SSPs_CVD_mortality")

.libPaths(c("/home/otavova/R/library"))

library(dplyr)
library(ggplot2)

summary_dir <- "Results/RegionalSummaryMeasures"
out_dir     <- "Results/Plots/RegionalSummaries/WorldBankUpdate"
frame_dir   <- file.path(out_dir, "gif_frames_trajectory")
dir.create(frame_dir, showWarnings = FALSE, recursive = TRUE)

ssp_colours <- c("SSP1" = "#1b9e77", "SSP2" = "#d95f02", "SSP3" = "#7570b3",
                 "SSP4" = "#e7298a", "SSP5" = "#66a61e")
ssp_breaks  <- c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5")

# Ordered to match the three-trajectory-shape narrative: already peaked;
# peaks mid-century; still rising at 2100 under the less favorable SSPs.
region_levels <- c(
  "Europe & Central Asia",
  "East Asia & Pacific", "High-Income Countries", "Latin America & Caribbean", "South Asia",
  "Middle East, North Africa, Afghanistan & Pakistan", "Sub-Saharan Africa"
)

proj_start_year <- 2030

regional_summary <- readRDS(file.path(summary_dir, "regional_summary_worldbankupdate.rds")) %>%
  filter(scenario != "historical", year >= proj_start_year) %>%
  mutate(
    region = recode(region, "High-Income" = "High-Income Countries"),
    region = factor(region, levels = region_levels),
    scenario = factor(scenario, levels = ssp_breaks)
  )

proj_years <- sort(unique(regional_summary$year))
n_regions  <- n_distinct(regional_summary$region)
x_range    <- range(regional_summary$year)

# One dummy (NA-valued) row per region, so facet_wrap always has all 7 panels
# to draw even on frames where no SSP data has been revealed yet -- geom_line/
# geom_point silently skip NA y-values, so this renders as empty panels.
region_skeleton <- regional_summary %>%
  distinct(region) %>%
  mutate(year = x_range[1], cardio_deaths_millions = NA_real_,
        scenario = factor(NA, levels = ssp_breaks))

make_frame <- function(completed_ssps, current_ssp, reveal_up_to_year) {
  parts <- list(region_skeleton)
  for (s in completed_ssps) {
    parts[[length(parts) + 1]] <- regional_summary %>% filter(scenario == s)
  }
  if (!is.null(current_ssp)) {
    parts[[length(parts) + 1]] <- regional_summary %>%
      filter(scenario == current_ssp, year <= reveal_up_to_year)
  }
  plot_data <- bind_rows(parts) %>%
    mutate(scenario = factor(scenario, levels = ssp_breaks))

  ggplot(plot_data, aes(x = year, y = cardio_deaths_millions, colour = scenario)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8) +
    facet_wrap(~ region, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = ssp_colours, breaks = ssp_breaks, name = "SSP", drop = FALSE) +
    scale_x_continuous(limits = x_range, breaks = seq(proj_start_year, 2100, by = 20)) +
    labs(title = "Cardiovascular deaths over time by region",
        subtitle = paste0("SSP projections, ", proj_start_year, "-2100"),
        x = NULL, y = "Cardiovascular deaths (millions)") +
    theme_bw(base_size = 15) +
    theme(
      plot.title         = element_text(face = "bold", size = 18),
      plot.subtitle      = element_text(size = 15),
      strip.background   = element_rect(fill = "grey85"),
      strip.text         = element_text(face = "bold", size = 14),
      axis.text          = element_text(size = 13),
      axis.title         = element_text(size = 15),
      legend.text        = element_text(size = 14),
      legend.title       = element_text(size = 15),
      panel.grid.minor   = element_blank(),
      legend.position    = "bottom"
    )
}

frame_i <- 0
save_frame <- function(p) {
  frame_i <<- frame_i + 1
  ggsave(file.path(frame_dir, sprintf("frame_%03d.png", frame_i)), p,
        width = 14, height = ceiling(n_regions / 3) * 3.6 + 2.2, dpi = 110, limitsize = FALSE)
}

HOLD_START   <- 6
HOLD_BETWEEN <- 6
HOLD_END     <- 14

for (i in seq_len(HOLD_START)) save_frame(make_frame(character(0), NULL, NA))

completed <- character(0)
for (s in ssp_breaks) {
  for (yr in proj_years) {
    save_frame(make_frame(completed, s, yr))
  }
  completed <- c(completed, s)
  for (i in seq_len(HOLD_BETWEEN)) save_frame(make_frame(completed, NULL, NA))
}

for (i in seq_len(HOLD_END)) save_frame(make_frame(completed, NULL, NA))

cat("Saved", frame_i, "frames to", frame_dir, "\n")
cat("n proj_years per SSP:", length(proj_years), "\n")

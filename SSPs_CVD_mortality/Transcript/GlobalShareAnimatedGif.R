# =============================================================================
# Global CVD Share -- Animated GIF (SSPs reveal one at a time)
#
# Renders the fig_global_ssp_deaths plot as a frame-by-frame animation:
# historical draws in first, then each SSP line draws in, in turn (SSP1,
# then SSP2, ... SSP5), holding briefly once complete before the next SSP
# starts. Previously completed SSPs stay visible while the next one draws.
#
# Frames are rendered here as individual PNGs; a companion Python script
# (using Pillow, since gifski/ImageMagick/ffmpeg are unavailable in this
# environment -- no rustc/apt access) stitches them into the final GIF.
# =============================================================================

setwd("/home/otavova/Healthy-SSPs/SSPs_CVD_mortality")

.libPaths(c("/home/otavova/R/library"))

library(dplyr)
library(ggplot2)
library(scales)

out_dir   <- "Results/Plots/Global"
frame_dir <- file.path(out_dir, "gif_frames")
dir.create(frame_dir, showWarnings = FALSE, recursive = TRUE)

ssp_colours <- c("historical" = "grey35", "SSP1" = "#1B9E77", "SSP2" = "#D95F02",
                 "SSP3" = "#7570B3", "SSP4" = "#E7298A", "SSP5" = "#66A61E")
ssp_breaks  <- c("historical", "SSP1", "SSP2", "SSP3", "SSP4", "SSP5")
ssp_labels  <- c("Historical", "SSP1 - Sustainability", "SSP2 - Middle of the Road",
                 "SSP3 - Regional Rivalry", "SSP4 - Inequality", "SSP5 - Fossil-fueled Development")

last_hist_year  <- 2023
proj_start_year <- 2025

obs_and_proj <- readRDS("Results/Projections/obs_and_proj.rds") %>%
  mutate(iso3 = as.character(iso3), sex = as.character(sex)) %>%
  filter(scenario == "historical" | year >= proj_start_year)

global_summary <- obs_and_proj %>%
  group_by(scenario, year) %>%
  summarise(cardio_deaths = sum(cardio_deaths, na.rm = TRUE),
            allcause_deaths = sum(allcause_deaths, na.rm = TRUE), .groups = "drop") %>%
  mutate(cardio_share = cardio_deaths / allcause_deaths)

# connect each SSP to the historical line's endpoint, same technique as GlobalPlots.R
hist_anchor <- global_summary %>% filter(scenario == "historical", year == last_hist_year)
ssps <- setdiff(unique(global_summary$scenario), "historical")
anchors <- lapply(ssps, function(s) hist_anchor %>% mutate(scenario = s)) %>% bind_rows()
global_summary <- bind_rows(global_summary, anchors) %>% arrange(scenario, year)

hist_data <- global_summary %>% filter(scenario == "historical")
proj_years <- sort(unique(global_summary$year[global_summary$year >= last_hist_year &
                                              global_summary$scenario != "historical"]))

y_range <- range(global_summary$cardio_share)
x_range <- range(global_summary$year)

make_frame <- function(completed_ssps, current_ssp, reveal_up_to_year) {
  parts <- list(hist_data)
  for (s in completed_ssps) {
    parts[[length(parts) + 1]] <- global_summary %>% filter(scenario == s)
  }
  if (!is.null(current_ssp)) {
    parts[[length(parts) + 1]] <- global_summary %>%
      filter(scenario == current_ssp, year <= reveal_up_to_year)
  }
  plot_data <- bind_rows(parts) %>%
    mutate(scenario = factor(scenario, levels = ssp_breaks))

  ggplot(plot_data, aes(x = year, y = cardio_share, colour = scenario)) +
    geom_line(linewidth = 1.1) +
    scale_colour_manual(values = ssp_colours, breaks = ssp_breaks, labels = ssp_labels,
                        name = NULL, drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = y_range) +
    scale_x_continuous(limits = x_range, breaks = c(1990, 2010, 2025, 2050, 2075, 2100)) +
    labs(title = "Cardiovascular share of all-cause deaths",
        subtitle = "Historical (1990-2023) and SSP projections (2024-2100)",
        x = NULL, y = "Proportion of deaths") +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          axis.text = element_text(size = 13), axis.title = element_text(size = 14),
          legend.text = element_text(size = 13))
}

frame_i <- 0
save_frame <- function(p) {
  frame_i <<- frame_i + 1
  ggsave(file.path(frame_dir, sprintf("frame_%03d.png", frame_i)), p,
        width = 10, height = 6.5, dpi = 120)
}

HOLD_START  <- 6    # frames pausing on historical-only, before SSP1 starts
HOLD_BETWEEN <- 6    # frames pausing after each SSP completes, before the next starts
HOLD_END    <- 14   # frames pausing on the fully-complete plot, before the GIF loops

# 1. Hold on historical only
for (i in seq_len(HOLD_START)) save_frame(make_frame(character(0), NULL, NA))

# 2. Reveal each SSP in turn, holding on completion
completed <- character(0)
for (s in ssps) {
  for (yr in proj_years) {
    save_frame(make_frame(completed, s, yr))
  }
  completed <- c(completed, s)
  for (i in seq_len(HOLD_BETWEEN)) save_frame(make_frame(completed, NULL, NA))
}

# 3. Final hold on the complete plot
for (i in seq_len(HOLD_END)) save_frame(make_frame(completed, NULL, NA))

cat("Saved", frame_i, "frames to", frame_dir, "\n")

# =============================================================================
# Global CVD Deaths -- Animated GIF (SSPs reveal one at a time)
#
# Companion to GlobalShareAnimatedGif.R, but for the CVD deaths (millions)
# metric instead of the CVD share: projections only (no historical line),
# starting at 2025, with a dot marker appearing at each 5-year data point as
# each SSP line draws in, in turn (SSP1, then SSP2, ... SSP5), holding
# briefly once complete before the next SSP starts.
# =============================================================================

setwd("/home/otavova/Healthy-SSPs/SSPs_CVD_mortality")

.libPaths(c("/home/otavova/R/library"))

library(dplyr)
library(ggplot2)

out_dir   <- "Results/Plots/Global"
frame_dir <- file.path(out_dir, "gif_frames_deaths")
dir.create(frame_dir, showWarnings = FALSE, recursive = TRUE)

ssp_colours <- c("SSP1" = "#1B9E77", "SSP2" = "#D95F02", "SSP3" = "#7570B3",
                 "SSP4" = "#E7298A", "SSP5" = "#66A61E")
ssp_breaks  <- c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5")
ssp_labels  <- c("SSP1 - Sustainability", "SSP2 - Middle of the Road",
                 "SSP3 - Regional Rivalry", "SSP4 - Inequality", "SSP5 - Fossil-fueled Development")

proj_start_year <- 2025

global_summary <- readRDS("Results/Projections/obs_and_proj.rds") %>%
  mutate(iso3 = as.character(iso3), sex = as.character(sex)) %>%
  filter(scenario != "historical", year >= proj_start_year) %>%
  group_by(scenario, year) %>%
  summarise(cardio_deaths = sum(cardio_deaths, na.rm = TRUE), .groups = "drop") %>%
  mutate(cardio_deaths_millions = cardio_deaths / 1e6,
        scenario = factor(scenario, levels = ssp_breaks))

proj_years <- sort(unique(global_summary$year))
y_range <- range(global_summary$cardio_deaths_millions)
x_range <- range(global_summary$year)

make_frame <- function(completed_ssps, current_ssp, reveal_up_to_year) {
  parts <- list()
  for (s in completed_ssps) {
    parts[[length(parts) + 1]] <- global_summary %>% filter(scenario == s)
  }
  if (!is.null(current_ssp)) {
    parts[[length(parts) + 1]] <- global_summary %>%
      filter(scenario == current_ssp, year <= reveal_up_to_year)
  }
  plot_data <- if (length(parts) == 0) global_summary[0, ] else bind_rows(parts) %>%
    mutate(scenario = factor(scenario, levels = ssp_breaks))

  ggplot(plot_data, aes(x = year, y = cardio_deaths_millions, colour = scenario)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.4) +
    scale_colour_manual(values = ssp_colours, breaks = ssp_breaks, labels = ssp_labels,
                        name = NULL, drop = FALSE) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_y_continuous(limits = y_range) +
    scale_x_continuous(limits = x_range, breaks = seq(proj_start_year, 2100, by = 25)) +
    labs(title = "Global cardiovascular deaths",
        subtitle = paste0("SSP projections, ", proj_start_year, "-2100"),
        x = NULL, y = "CVD deaths (millions)") +
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

HOLD_START   <- 6    # frames pausing on an empty axes, before SSP1 starts
HOLD_BETWEEN <- 6    # frames pausing after each SSP completes, before the next starts
HOLD_END     <- 14   # frames pausing on the fully-complete plot, before the GIF loops

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

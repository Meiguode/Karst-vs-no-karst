library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(patchwork)
library(stringr)

data_files <- list(
  "Combined"  = "03_combined_with_koppen.csv",
  "Non-karst" = "03_nonkarst_with_koppen.csv",
  "Karst"     = "03_karst_with_koppen.csv"
)

climate_levels <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")

raw_data_list <- lapply(names(data_files), function(nm) {
  df <- read.csv(data_files[[nm]], stringsAsFactors = FALSE)
  df$Dataset <- nm
  df
})

all_raw_data <- bind_rows(raw_data_list) %>%
  mutate(
    Year = as.numeric(Study_midyear),
    kg_description = trimws(kg_description),
    ClimateGroup = case_when(
      grepl("^Tropical",   kg_description, ignore.case = TRUE) ~ "Tropical",
      grepl("^Arid",       kg_description, ignore.case = TRUE) ~ "Arid",
      grepl("^Temperate",  kg_description, ignore.case = TRUE) ~ "Temperate",
      grepl("^Cold",       kg_description, ignore.case = TRUE) ~ "Cold",
      grepl("^Polar",      kg_description, ignore.case = TRUE) ~ "Polar",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(ClimateGroup), !is.na(Rs_annual), !is.na(Year)) %>%
  mutate(
    ClimateGroup = factor(ClimateGroup, levels = climate_levels),
    Dataset_Label = factor(
      Dataset,
      levels = c("Combined", "Non-karst", "Karst"),
      labels = c("(a) Combined", "(b) Non-karst", "(c) Karst")
    )
  )

all_data_with_periods <- all_raw_data %>%
  mutate(
    Period = case_when(
      Dataset %in% c("Combined", "Non-karst") & Year < 2013 ~ "Pre",
      Dataset %in% c("Combined", "Non-karst") & Year <= 2016 ~ "Transition",
      Dataset %in% c("Combined", "Non-karst")               ~ "Post",
      Dataset == "Karst" & Year < 2009                      ~ "Pre",
      Dataset == "Karst" & Year <= 2014                     ~ "Transition",
      Dataset == "Karst"                                    ~ "Post",
      TRUE ~ NA_character_
    )
  )

annual_means <- all_data_with_periods %>%
  group_by(ClimateGroup, Dataset, Dataset_Label, Year, Period) %>%
  summarise(
    Rs_mean = mean(Rs_annual, na.rm = TRUE),
    Rs_sd   = sd(Rs_annual, na.rm = TRUE),
    N       = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Rs_se    = Rs_sd / sqrt(N),
    CI_lower = Rs_mean - 1.96 * Rs_se,
    CI_upper = Rs_mean + 1.96 * Rs_se
  )

period_boundaries <- annual_means %>%
  group_by(ClimateGroup, Dataset, Period) %>%
  summarise(
    xmin = min(Year),
    xmax = max(Year),
    .groups = "drop"
  )

trend_lines <- all_data_with_periods %>%
  group_by(ClimateGroup, Dataset, Dataset_Label, Period) %>%
  filter(n() >= 10) %>%
  do({
    fit <- lm(Rs_annual ~ Year, data = .)
    sm  <- summary(fit)

    p_value <- sm$coefficients["Year", "Pr(>|t|)"]

    # Generate smooth line across observed year range
    yrs <- seq(min(.$Year), max(.$Year), length.out = 30)

    data.frame(
      Year = yrs,
      Rs_trend = predict(fit, newdata = data.frame(Year = yrs)),
      Significant = ifelse(
        !is.na(p_value) & p_value <= 0.05,
        "Significant (p ≤ 0.05)",
        "Not significant (p > 0.05)"
      )
    )
  }) %>%
  ungroup()

period_colors <- c(
  "Pre" = "#E8F4F8",
  "Transition" = "#FFF2CC",
  "Post" = "#F8E8E8"
)

trend_colors <- c(
  "Significant (p ≤ 0.05)" = "#FFD700",
  "Not significant (p > 0.05)" = "#A9A9A9"
)

climate_zone_plots <- list()

for (cz in climate_levels) {

  plot_data   <- filter(annual_means, ClimateGroup == cz)
  period_data <- filter(period_boundaries, ClimateGroup == cz)
  trend_data  <- filter(trend_lines, ClimateGroup == cz)

  if (nrow(plot_data) == 0) next

  p <- ggplot() +
    geom_rect(
      data = period_data,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Period),
      alpha = 0.1
    ) +
    geom_ribbon(
      data = plot_data,
      aes(x = Year, ymin = Rs_mean - Rs_sd, ymax = Rs_mean + Rs_sd,
          group = Dataset_Label),
      fill = "gray70", alpha = 0.15
    ) +
    geom_line(
      data = plot_data,
      aes(x = Year, y = Rs_mean, color = Dataset_Label),
      linewidth = 1.1
    ) +
    geom_point(
      data = plot_data,
      aes(x = Year, y = Rs_mean, color = Dataset_Label),
      size = 2
    ) +
    geom_errorbar(
      data = plot_data,
      aes(x = Year, ymin = CI_lower, ymax = CI_upper, color = Dataset_Label),
      width = 0.4, alpha = 0.4
    ) +
    geom_line(
      data = trend_data,
      aes(
        x = Year, y = Rs_trend,
        group = interaction(Dataset_Label, Period),
        color = Significant
      ),
      linetype = "longdash",
      linewidth = 1
    ) +
    facet_wrap(~Dataset_Label, ncol = 3, scales = "free_y") +
    scale_color_manual(
      values = c(
        "(a) Combined" = "#50ad9f",
        "(b) Non-karst" = "#bc272d",
        "(c) Karst" = "#0000a2",
        trend_colors
      ),
      guide = "none"
    ) +
    scale_fill_manual(values = period_colors, guide = "none") +
    coord_cartesian(ylim = c(0, 3000)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = paste(cz, "Climate Zone"),
      x = "Year",
      y = expression(paste(R[s], " (g C m"^{-2}, " yr"^{-1}, ")"))
    )

  climate_zone_plots[[cz]] <- p
}

combined_figure <- wrap_plots(climate_zone_plots, ncol = 2)

ggsave(
  "Rs_ClimateZones_Breakpoints_Optimal.tiff",
  combined_figure,
  width = 16, height = 14, dpi = 600,
  compression = "lzw", bg = "white"
)

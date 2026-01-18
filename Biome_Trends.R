library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(cowplot)

source("statistical_prep.R")

fig_path <- "figures"
dir.create(fig_path, showWarnings = FALSE)

datasets_plot <- list(
  "Combined"   = combined,
  "Non-karst"  = nonkarst,
  "Karst"      = karst
)

biome_cols <- c(
  "Alpine"                 = "#1f78b4",
  "Boreal / Arctic"        = "#33a02c",
  "Temperate"              = "#e31a1c",
  "Tropical / Subtropical" = "#ff7f00"
)

all_data_list <- list()

for (i in seq_along(datasets_plot)) {

  panel_label <- letters[i]

  dat_clean <- datasets_plot[[i]] %>%
    mutate(
      Biome_clean = tolower(trimws(Biome)),
      Biome = case_when(
        grepl("tropical|subtropical", Biome_clean) ~ "Tropical / Subtropical",
        grepl("temperate|mediterranean|grassland", Biome_clean) ~ "Temperate",
        grepl("boreal|tundra|arctic", Biome_clean) ~ "Boreal / Arctic",
        grepl("alpine", Biome_clean) ~ "Alpine",
        TRUE ~ NA_character_
      ),
      Panel = panel_label
    ) %>%
    filter(!is.na(Biome))

  annual_stats <- dat_clean %>%
    group_by(Study_midyear, Biome, Panel) %>%
    summarise(
      Rs_mean = mean(Rs_annual, na.rm = TRUE),
      Rs_sd   = sd(Rs_annual, na.rm = TRUE),
      n       = n(),
      .groups = "drop"
    ) %>%
    mutate(
      Rs_se    = Rs_sd / sqrt(n),
      ci_lower = Rs_mean - 1.96 * Rs_se,
      ci_upper = Rs_mean + 1.96 * Rs_se
    )

  trend_lines <- dat_clean %>%
    group_by(Period, Biome, Panel) %>%
    group_modify(~{
      d <- .x
      if (nrow(d) < 3 || length(unique(d$Study_midyear)) < 2) return(tibble())
      m <- lm(Rs_annual ~ Study_midyear, d)
      s <- summary(m)$coefficients
      if (!("Study_midyear" %in% rownames(s))) return(tibble())
      if (s["Study_midyear", "Pr(>|t|)"] >= 0.05) return(tibble())
      yrs <- seq(min(d$Study_midyear), max(d$Study_midyear), length.out = 50)
      tibble(
        Study_midyear = yrs,
        Rs_trend = s["(Intercept)", "Estimate"] +
                   s["Study_midyear", "Estimate"] * yrs
      )
    }) %>% ungroup()

  all_data_list[[panel_label]] <- list(
    annual = annual_stats,
    trends = trend_lines
  )
}

all_annual <- bind_rows(lapply(all_data_list, `[[`, "annual"))
all_trends <- bind_rows(lapply(all_data_list, `[[`, "trends"))

bivariate_stats <- list()

for (panel in unique(all_data_list %>% names())) {

  dat_panel <- datasets_plot[[which(letters == panel)]] %>%
    mutate(
      Biome_clean = tolower(trimws(Biome)),
      Biome = case_when(
        grepl("tropical|subtropical", Biome_clean) ~ "Tropical / Subtropical",
        grepl("temperate|mediterranean|grassland", Biome_clean) ~ "Temperate",
        grepl("boreal|tundra|arctic", Biome_clean) ~ "Boreal / Arctic",
        grepl("alpine", Biome_clean) ~ "Alpine",
        TRUE ~ NA_character_
      ),
      Panel = panel
    ) %>%
    filter(!is.na(Biome), !is.na(Period))

  for (b in unique(dat_panel$Biome)) {
    for (p in unique(dat_panel$Period)) {

      dsub <- dat_panel %>%
        filter(Biome == b, Period == p)

      n_obs <- nrow(dsub)

      if (n_obs < 3 || length(unique(dsub$Study_midyear)) < 2) {
        bivariate_stats[[length(bivariate_stats) + 1]] <-
          data.frame(
            Panel = panel,
            Biome = b,
            Period = p,
            N = n_obs,
            Slope = NA,
            Slope_SE = NA,
            P_value = NA,
            R_squared = NA,
            Significant = "Insufficient data"
          )
        next
      }

      fit <- lm(Rs_annual ~ Study_midyear, data = dsub)
      sm  <- summary(fit)

      slope     <- sm$coefficients["Study_midyear", "Estimate"]
      slope_se  <- sm$coefficients["Study_midyear", "Std. Error"]
      pval      <- sm$coefficients["Study_midyear", "Pr(>|t|)"]
      r2        <- sm$r.squared

      bivariate_stats[[length(bivariate_stats) + 1]] <-
        data.frame(
          Panel = panel,
          Biome = b,
          Period = p,
          N = n_obs,
          Slope = slope,
          Slope_SE = slope_se,
          P_value = pval,
          R_squared = r2,
          Significant = ifelse(pval < 0.05, "Yes", "No")
        )
    }
  }
}

bivariate_stats_df <- bind_rows(bivariate_stats)

write.csv(
  bivariate_stats_df,
  file.path(fig_path, "Figure3_Bivariate_Regression_Stats.csv"),
  row.names = FALSE
)

base_panel <- function(panel_letter, show_xlab = FALSE, show_legend = FALSE) {

  annual_panel <- filter(all_annual, Panel == panel_letter)
  trends_panel <- filter(all_trends, Panel == panel_letter)

  y_max <- max(annual_panel$ci_upper, na.rm = TRUE) * 1.05
  y_min <- max(0, min(annual_panel$ci_lower, na.rm = TRUE) * 0.95)

  p <- ggplot() +
    geom_point(
      data = annual_panel,
      aes(x = Study_midyear, y = Rs_mean, color = Biome),
      size = 2.5
    ) +
    geom_errorbar(
      data = annual_panel,
      aes(x = Study_midyear, ymin = ci_lower, ymax = ci_upper, color = Biome),
      width = 0.5, linewidth = 0.6
    ) +
    geom_line(
      data = trends_panel,
      aes(x = Study_midyear, y = Rs_trend, color = Biome),
      linewidth = 1.3
    ) +
    scale_color_manual(
      values = biome_cols,
      guide = guide_legend(
        title = "Biome",
        nrow = 1,
        byrow = TRUE
      )
    ) +
    scale_x_continuous(breaks = seq(1990, 2020, 5),
                       limits = c(1990, 2022)) +
    scale_y_continuous(limits = c(y_min, y_max)) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.position = if (show_legend) "top" else "none",
      axis.title.y = element_text(face = "bold")
    ) +
    labs(y = expression(R[s]~"(g C m"^{-2}~"yr"^{-1}*")"))

  if (show_xlab) p <- p + labs(x = "Year")
  else p <- p + theme(axis.title.x = element_blank())
  p
}

legend_plot <- base_panel("a", show_legend = TRUE)
legend_grob <- get_legend(legend_plot)

panel_a <- base_panel("a")
panel_b <- base_panel("b")
panel_c <- base_panel("c", show_xlab = TRUE)

panels <- panel_a / panel_b / panel_c

final_plot <- plot_grid(
  legend_grob,
  panels,
  ncol = 1,
  rel_heights = c(0.08, 1)
)

ggsave(
  file.path(fig_path, "Figure3_Biome_Trends.tiff"),
  final_plot,
  width = 8.5,
  height = 10.5,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

legend_only_plot <- ggplot(
  data.frame(
    Biome = factor(names(biome_cols), levels = names(biome_cols)),
    x = seq_along(biome_cols),
    y = 1
  ),
  aes(x = x, y = y, color = Biome)
) +
  geom_point(size = 3) +
  geom_line(aes(group = Biome), linewidth = 1.2) +
  scale_color_manual(
    values = biome_cols,
    guide = guide_legend(
      title = "Biome",
      nrow = 1,
      byrow = TRUE
    )
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.4, "cm"),
    legend.margin = margin(6, 6, 6, 6)
  )

legend_grob_only <- get_legend(legend_only_plot)

legend_figure <- plot_grid(
  legend_grob_only,
  ncol = 1
)

ggsave(
  file.path(fig_path, "Figure3_Legend.tiff"),
  legend_figure,
  width = 8.5,
  height = 1.0,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

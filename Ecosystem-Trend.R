library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(forcats)
library(cowplot)

source("statistical_prep.R")

fig_path <- "figures"
dir.create(fig_path, showWarnings = FALSE)

primary_cols <- c(
  "#0000a2",
  "#bc272d",
  "#e9c716",
  "#50ad9f",
  "#9a5ba1",
  "#f4a300"
)

other_col <- "#bdbdbd"

period_cols <- c(
  "Pre-breakpoint" = "#0000a2",
  "Breakpoint" = "#50ad9f",
  "Post-breakpoint" = "#bc272d"
)

datasets_plot <- list(
  "Combined" = combined,
  "Non-karst" = nonkarst,
  "Karst" = karst
)

top_ecosystems_list <- lapply(datasets_plot, function(dat) {
  dat %>%
    filter(!is.na(Ecosystem_type), Outlier == 0) %>%
    count(Ecosystem_type, sort = TRUE) %>%
    slice_head(n = 6) %>%
    pull(Ecosystem_type)
})

prep_ecosystem_data <- function(dat, top_ecosystems) {
  dat %>%
    filter(
      !is.na(Ecosystem_type),
      !is.na(Rs_annual),
      Outlier == 0
    ) %>%
    mutate(
      Ecosystem_type = trimws(as.character(Ecosystem_type)),
      Ecosystem_grp = if_else(
        Ecosystem_type %in% top_ecosystems,
        Ecosystem_type,
        "Others"
      ),
      Ecosystem_grp = factor(
        Ecosystem_grp,
        levels = c(top_ecosystems, "Others")
      ),
      Period = factor(Period, 
                      levels = c("Pre-breakpoint", "Breakpoint", "Post-breakpoint"))
    )
}

prep_slope_data <- function(dat, top_ecosystems) {
  dat %>%
    filter(
      !is.na(Ecosystem_type),
      !is.na(Rs_annual),
      Outlier == 0,
      Ecosystem_type %in% top_ecosystems
    ) %>%
    mutate(
      Ecosystem_type = trimws(as.character(Ecosystem_type)),
      Ecosystem_grp = factor(Ecosystem_type, levels = top_ecosystems),
      Period = factor(Period, 
                      levels = c("Pre-breakpoint", "Breakpoint", "Post-breakpoint"))
    )
}

cat("\nTop 6 ecosystems by dataset:\n")
for (nm in names(datasets_plot)) {
  top_ecos <- top_ecosystems_list[[nm]]
  cat("\n", nm, ":\n")
  cat(paste0("  ", 1:6, ". ", top_ecos, "\n"), sep = "")
}


make_main_plot <- function(dat, panel_label, dataset_name) {
  
  annual_stats <- dat %>%
    group_by(Study_midyear, Ecosystem_grp, Period) %>%
    summarise(
      Rs_mean = mean(Rs_annual, na.rm = TRUE),
      Rs_sd = sd(Rs_annual, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  trend_models <- dat %>%
    group_by(Ecosystem_grp, Period) %>%
    filter(n() >= 3) %>%  # Need at least 3 points to fit a line
    group_modify(~ {
      m <- lm(Rs_annual ~ Study_midyear, .x)
      yrs <- seq(min(.x$Study_midyear), max(.x$Study_midyear), length.out = 50)
      pred <- predict(m, newdata = data.frame(Study_midyear = yrs), se.fit = TRUE)
      tibble(
        Study_midyear = yrs,
        fit = pred$fit,
        se = pred$se.fit
      )
    }) %>%
    ungroup()
  
  n_ecos <- length(levels(dat$Ecosystem_grp)) - 1
  cols <- c(primary_cols[1:n_ecos], other_col)
  names(cols) <- levels(dat$Ecosystem_grp)
  
  p <- ggplot() +
    geom_errorbar(
      data = annual_stats,
      aes(
        x = Study_midyear,
        ymin = Rs_mean - Rs_sd,
        ymax = Rs_mean + Rs_sd,
        color = Ecosystem_grp
      ),
      alpha = 0.3,
      width = 0.5
    ) +
    geom_point(
      data = annual_stats,
      aes(x = Study_midyear, y = Rs_mean, color = Ecosystem_grp),
      size = 2,
      alpha = 0.7
    ) +
    geom_ribbon(
      data = trend_models,
      aes(
        x = Study_midyear,
        ymin = fit - 1.96 * se,
        ymax = fit + 1.96 * se,
        fill = Ecosystem_grp,
        group = interaction(Ecosystem_grp, Period)
      ),
      alpha = 0.2
    ) +
    geom_line(
      data = trend_models,
      aes(
        x = Study_midyear, 
        y = fit, 
        color = Ecosystem_grp,
        linetype = Period,
        group = interaction(Ecosystem_grp, Period)
      ),
      linewidth = 1
    ) +
    scale_color_manual(values = cols, name = "Ecosystem") +
    scale_fill_manual(values = cols, name = "Ecosystem") +
    scale_linetype_manual(values = c("Pre-breakpoint" = "solid", 
                                     "Breakpoint" = "dashed", 
                                     "Post-breakpoint" = "dotted"),
                          name = "Period") +
    geom_vline(
      data = data.frame(
        dataset = dataset_name,
        breakpoint = if(dataset_name %in% c("Combined", "Non-karst")) {
          c(2012.5, 2016.5)
        } else {
          c(2008.5, 2014.5)
        }
      ),
      aes(xintercept = breakpoint),
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.5,
      alpha = 0.7
    ) +
    theme_classic(base_size = 13) +
    labs(
      x = "Year",
      y = expression(paste("Annual soil respiration (", R[s], " g C m"^{-2}, " yr"^{-1}, ")")),
      title = panel_label
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      legend.box = "vertical",
      legend.margin = margin(),
      legend.key.width = unit(1.5, "cm"),
      axis.title.y = element_text(size = 11)
    ) +
    guides(
      color = guide_legend(order = 1, ncol = 1),
      fill = guide_legend(order = 1, ncol = 1),
      linetype = guide_legend(order = 2, ncol = 1)
    )
  
  return(p)
}

make_bar_plot <- function(dat, dataset_name) {
  
  slope_data <- dat %>%
    group_by(Ecosystem_grp, Period) %>%
    filter(n() >= 3) %>%  # Need at least 3 points per period
    group_modify(~ {
      m <- lm(Rs_annual ~ Study_midyear, .x)
      sm <- summary(m)
      tibble(
        Slope = coef(m)[2],
        SE = sm$coefficients[2, 2],
        P = sm$coefficients[2, 4]
      )
    }) %>%
    ungroup() %>%
    mutate(
      CI_lo = Slope - 1.96 * SE,
      CI_hi = Slope + 1.96 * SE,
      Significance = case_when(
        P < 0.001 ~ "***",
        P < 0.01 ~ "**",
        P < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  p <- ggplot(slope_data, aes(x = Slope, y = Ecosystem_grp, fill = Period)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
    geom_text(
      data = slope_data %>% filter(Significance != ""),
      aes(label = Significance, 
          x = Slope + ifelse(Slope >= 0, max(abs(slope_data$Slope)) * 0.05, 
                            -max(abs(slope_data$Slope)) * 0.05)),
      position = position_dodge(width = 0.8),
      size = 4,
      vjust = 0.5,
      hjust = ifelse(slope_data %>% filter(Significance != "") %>% pull(Slope) >= 0, 0, 1)
    ) +
    scale_fill_manual(
      values = c(
        "Pre-breakpoint" = "#0000a2",
        "Breakpoint" = "#50ad9f",
        "Post-breakpoint" = "#bc272d"
      ),
      name = "Period"
    ) +
    scale_x_continuous(
      limits = function(x) {
        x_range <- range(x, na.rm = TRUE)
        padding <- diff(x_range) * 0.1
        c(x_range[1] - padding, x_range[2] + padding)
      },
      expand = expansion(mult = 0.05)
    ) +
    labs(
      x = expression(paste("R"[s], " Slope (g C m"^{-2}, " yr"^{-1}, ")")),
      y = "",
      title = paste0("(", letters[which(names(datasets_plot) == dataset_name) + 3], ")")
      #title = paste0("(", letters[which(names(datasets_plot) == dataset_name) + 3], ") ", dataset_name, " - Slope Analysis")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.title.x = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.direction = "horizontal",
      axis.text.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

cat("Creating Main Only plots (with top 6 ecosystems)...\n")

main_plots_only <- lapply(names(datasets_plot), function(nm) {
  top_ecos <- top_ecosystems_list[[nm]]
  dat <- prep_ecosystem_data(datasets_plot[[nm]], top_ecos)
  make_main_plot(dat, 
                paste0("(", letters[which(names(datasets_plot) == nm)], ") ", nm),
                dataset_name = nm)
})

cat("\nTop 6 ecosystems by dataset:\n")
for (nm in names(datasets_plot)) {
  top_ecos <- top_ecosystems_list[[nm]]
  cat("\n", nm, ":\n")
  cat(paste0("  ", 1:6, ". ", top_ecos, "\n"), sep = "")
}

final_main_plot <- wrap_plots(main_plots_only, nrow = 1, guides = "collect") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11, face = "bold"))

ggsave(
  file.path(fig_path, "Ecosystem_Trends_Top6_MainOnly.png"),
  final_main_plot,
  width = 18,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(fig_path, "Ecosystem_Trends_Top6_MainOnly.tiff"),
  final_main_plot,
  width = 18,
  height = 6.5,
  dpi = 300,
  bg = "white",
  compression = "lzw"
)

cat("\nCreating Combined figure\n")

all_main_plots <- list()
all_bar_plots <- list()

for (i in seq_along(names(datasets_plot))) {
  nm <- names(datasets_plot)[i]
  top_ecos <- top_ecosystems_list[[nm]]
  
  dat_main <- prep_ecosystem_data(datasets_plot[[nm]], top_ecos)
  main_p <- make_main_plot(dat_main, 
                          paste0("(", letters[i], ") ", nm),
                          dataset_name = nm) +
    theme(legend.position = "none",
          plot.title = element_text(size = 12, face = "bold"))
  
  dat_slope <- prep_slope_data(datasets_plot[[nm]], top_ecos)
  bar_p <- make_bar_plot(dat_slope, nm) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          legend.position = "none")
  
  all_main_plots[[i]] <- main_p
  all_bar_plots[[i]] <- bar_p
  
  cat("\nSlope analysis for", nm, ":\n")
  print(table(dat_slope$Ecosystem_grp, dat_slope$Period))
}

first_nm <- names(datasets_plot)[1]
top_ecos <- top_ecosystems_list[[first_nm]]
legend_dat <- prep_ecosystem_data(datasets_plot[[first_nm]], top_ecos)

ecosystem_legend_plot <- ggplot(data.frame(x = 1:length(levels(legend_dat$Ecosystem_grp)), 
                                           y = 1,
                                           Ecosystem = factor(levels(legend_dat$Ecosystem_grp), 
                                                            levels = levels(legend_dat$Ecosystem_grp)))) +
  geom_point(aes(x = x, y = y, color = Ecosystem), size = 4) +
  scale_color_manual(values = c(primary_cols[1:6], other_col),
                     name = "Ecosystem Type:",
                     labels = function(x) {
                       ifelse(x == "Others", "Others", x)
                     }) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9))

period_line_legend_plot <- ggplot(data.frame(x = 1:3, 
                                             y = 1,
                                             Period = factor(c("Pre-breakpoint", "Breakpoint", "Post-breakpoint"),
                                                           levels = c("Pre-breakpoint", "Breakpoint", "Post-breakpoint")))) +
  geom_line(aes(x = x, y = y, linetype = Period), linewidth = 1.5) +
  scale_linetype_manual(values = c("Pre-breakpoint" = "solid", 
                                   "Breakpoint" = "dashed", 
                                   "Post-breakpoint" = "dotted"),
                        name = "Trend Period:") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9),
        legend.key.width = unit(2, "cm"))

period_fill_legend_plot <- ggplot(data.frame(x = 1:3, 
                                             y = 1,
                                             Period = factor(c("Pre-breakpoint", "Breakpoint", "Post-breakpoint"),
                                                           levels = c("Pre-breakpoint", "Breakpoint", "Post-breakpoint")))) +
  geom_col(aes(x = x, y = y, fill = Period), width = 0.7) +
  scale_fill_manual(values = period_cols,
                    name = "Period (Bar Plots):") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9))

ecosystem_legend <- get_legend(ecosystem_legend_plot)
period_line_legend <- get_legend(period_line_legend_plot)
period_fill_legend <- get_legend(period_fill_legend_plot)

combined_legend <- plot_grid(
  ecosystem_legend,
  period_line_legend,
  period_fill_legend,
  nrow = 1,
  rel_widths = c(1.8, 1, 1)  # Wider for ecosystems (6 items + Others)
)

top_row <- plot_grid(plotlist = all_main_plots, nrow = 1, 
                     labels = NULL, align = "h", axis = "tb")

bottom_row <- plot_grid(plotlist = all_bar_plots, nrow = 1,
                       labels = NULL, align = "h", axis = "tb")

combined_plot <- plot_grid(
  top_row,
  bottom_row,
  combined_legend,
  ncol = 1,
  rel_heights = c(1, 1, 0.15),
  labels = c("", "", ""),
  label_size = 12
)

ggsave(
  file.path(fig_path, "Ecosystem_Trends_Top6_Combined.tiff"),
  combined_plot,
  width = 18,
  height = 12,
  dpi = 300,
  bg = "white",
  compression = "lzw"
)

cat("✓ Saved Combined plot\n")

top_ecosystems_info <- lapply(names(datasets_plot), function(nm) {
  top_ecos <- top_ecosystems_list[[nm]]
  dat <- prep_ecosystem_data(datasets_plot[[nm]], top_ecos)
  top_6 <- top_ecos
  others_count <- sum(dat$Ecosystem_grp == "Others")
  total_count <- nrow(dat)
  
  data.frame(
    Dataset = nm,
    Ecosystem_1 = top_6[1],
    Ecosystem_2 = top_6[2],
    Ecosystem_3 = top_6[3],
    Ecosystem_4 = top_6[4],
    Ecosystem_5 = top_6[5],
    Ecosystem_6 = top_6[6],
    Others_Count = others_count,
    Total_Count = total_count,
    Others_Percent = round(others_count/total_count * 100, 1)
  )
})

top_ecosystems_df <- do.call(rbind, top_ecosystems_info)
write.csv(top_ecosystems_df, 
          file.path(fig_path, "Top6_Ecosystems_Summary.csv"),
          row.names = FALSE)

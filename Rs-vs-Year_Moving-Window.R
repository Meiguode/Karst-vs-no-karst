library(dplyr)
library(ggplot2)
library(ggthemes)
library(mblm)
library(Kendall)
library(patchwork)

source("statistical_prep.R")

datasets <- list(
  "Non-karst" = nonkarst,
  "Karst"     = karst,
  "Combined"  = combined
)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

WINDOW <- 10
STEP   <- 1

mycolors <- colorRampPalette(
  c("#0000a2", "#bc272d", "#e9c716", "#50ad9f")
)

all_winresults <- list()

for (dname in names(datasets)) {

  df <- datasets[[dname]] %>%
    select(Study_midyear, Rs_annual) %>%
    filter(!is.na(Study_midyear), !is.na(Rs_annual))

  start_year <- min(df$Study_midyear)
  end_year   <- max(df$Study_midyear)

  out <- list()
  t <- start_year + WINDOW

  while (t <= end_year) {

    win <- df %>%
      filter(
        Study_midyear >= (t - WINDOW),
        Study_midyear < t
      )

    if (nrow(win) >= 5) {
      ts_fit <- mblm(Rs_annual ~ Study_midyear, win)

      slope     <- coef(ts_fit)[2]
      intercept <- coef(ts_fit)[1]

      kt <- Kendall(win$Study_midyear, win$Rs_annual)
      pval <- kt$sl

      r2 <- cor(win$Study_midyear, win$Rs_annual,
                method = "spearman", use = "complete.obs")^2

      out[[length(out) + 1]] <- data.frame(
        Dataset    = dname,
        Year.range = paste0(t - WINDOW, "-", t - 1),
        Year       = t - 1,
        Slope      = slope,
        Intercept  = intercept,
        R.square   = r2,
        P.value    = pval,
        Num        = nrow(win),
        Mean_Rs    = mean(win$Rs_annual, na.rm = TRUE),
        Median_Rs  = median(win$Rs_annual, na.rm = TRUE),
        SD_Rs      = sd(win$Rs_annual, na.rm = TRUE),
        SE_Rs      = sd(win$Rs_annual, na.rm = TRUE) / sqrt(nrow(win)),
        stringsAsFactors = FALSE
      )
    }

    t <- t + STEP
  }

  winresult <- bind_rows(out)

  winresult <- winresult %>%
    mutate(
      sig = case_when(
        P.value < 0.001 ~ "***",
        P.value < 0.01  ~ "**",
        P.value < 0.05  ~ "*",
        P.value < 0.1   ~ "·",
        TRUE            ~ ""
      ),
      Significance_Level = case_when(
        P.value < 0.001 ~ "p < 0.001",
        P.value < 0.01  ~ "p < 0.01",
        P.value < 0.05  ~ "p < 0.05",
        P.value < 0.1   ~ "p < 0.1",
        TRUE            ~ "Not significant"
      ),
      Num.italic = paste0("italic(", Num, ")")
    )

  write.csv(
    winresult,
    paste0("results/moving_window_comprehensive_stats_", dname, ".csv"),
    row.names = FALSE
  )

  significant_periods <- winresult %>% filter(P.value < 0.05)

  write.csv(
    significant_periods,
    paste0("results/significant_periods_stats_", dname, ".csv"),
    row.names = FALSE
  )

  all_winresults[[dname]] <- winresult

  cat("Finished moving-window analysis for:", dname, "\n")
}

all_slopes_vector <- unlist(lapply(all_winresults, function(x) x$Slope))
global_ymax <- max(abs(all_slopes_vector), na.rm = TRUE) * 1.3
y_N_global  <- global_ymax * 0.95

dataset_order <- c("Combined", "Non-karst", "Karst")
labels <- c("(a)", "(b)", "(c)")
plot_list <- list()

for (i in seq_along(dataset_order)) {
  dname <- dataset_order[i]
  winresult <- all_winresults[[dname]]
  
  winresult <- winresult %>%
    mutate(
      star_y = Slope + ifelse(Slope >= 0, 0.05 * global_ymax, -0.05 * global_ymax)
    )
  
  TS <- theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, size = 8, hjust = 1, vjust = 1),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "grey3", fill = NA, linewidth = 1),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),  # Left-aligned
    plot.title.position = "plot",  # Aligns title to the plot margin
    plot.margin = margin(5, 10, 5, 10, "pt")  # Adjusted margins
  )

  p <- ggplot(winresult, aes(x = Year, y = Slope)) +
    geom_bar(
      aes(fill = as.factor(Year)),
      stat = "identity",
      color = "navy",
      linewidth = 0.2,
      width = 0.7
    ) +
    geom_text(
      aes(y = star_y, label = sig),
      size = 3,
      vjust = ifelse(winresult$Slope >= 0, -0.5, 1.5)
    ) +
    geom_text(
      aes(y = y_N_global, label = Num.italic),
      size = 2.5,
      parse = TRUE
    ) +
    scale_fill_manual(values = mycolors(nrow(winresult))) +
    scale_x_continuous(
      breaks = winresult$Year,
      labels = winresult$Year
    ) +
    scale_y_continuous(
      limits = c(-global_ymax, global_ymax),
      breaks = pretty(c(-global_ymax, global_ymax), n = 6)
    ) +
    # THEME MUST BE ADDED WITH +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, size = 8, hjust = 1, vjust = 1),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "grey3", fill = NA, linewidth = 1),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 10),
      plot.title = element_text(size = 11, face = "bold", hjust = 0),
      plot.title.position = "plot",
      plot.margin = margin(5, 10, 5, 10, "pt")
    ) +
    labs(
      x = "Year",
      y = expression(Delta~R[s]~"(g C m"^{-2}~"yr"^{-2}*")"),
      title = labels[i]
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray40", linewidth = 0.3)

  plot_list[[dname]] <- p
  
  ggsave(
    paste0("figures/fig_moving_window_", dname, ".tiff"),
    p,
    width = 6.5,
    height = 2,
    units = "in",
    dpi = 900,
    compression = "lzw"
  )
}

combined_plot <- plot_list[["Combined"]] / plot_list[["Non-karst"]] / plot_list[["Karst"]] +
  plot_layout(heights = c(1, 1, 1))

ggsave(
  "figures/fig_moving_window_COMBINED_vertical.tiff",
  combined_plot,
  width = 7.5,
  height = 9,
  units = "in",
  dpi = 900,
  compression = "lzw"
)

ggsave(
  "figures/fig_moving_window_COMBINED_vertical.png",
  combined_plot,
  width = 7.5,
  height = 9,
  dpi = 300
)

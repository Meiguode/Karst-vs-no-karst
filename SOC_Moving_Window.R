library(dplyr)
library(ggplot2)
library(ggthemes)
library(mblm)
library(readr)
library(gridExtra)

datasets <- list(
  "Combined"  = "03_combined_with_koppen.csv",
  "Non-karst" = "03_nonkarst_with_koppen.csv",
  "Karst"     = "03_karst_with_koppen.csv"
)

prep_filter <- function(path, label) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(
      Study_midyear >= 1990,
      Study_midyear <= 2022,
      Meas_method %in% c("IRGA", "Gas Chromatography"),
      Ecosystem_type != "Cropland",
      Outlier == 0,
      !is.na(Rs_annual),
      !is.na(SOC_stock),
      Rs_annual > 0,
      Rs_annual <= 5000
    ) %>%
    mutate(
      Source = label,
      Year   = Study_midyear
    )
}

data_list <- lapply(names(datasets), function(name) {
  prep_filter(datasets[[name]], name)
})
names(data_list) <- names(datasets)

data_all      <- data_list[["Combined"]]
data_nonkarst <- data_list[["Non-karst"]]
data_karst    <- data_list[["Karst"]]

cat("\nDEBUG: rows after preprocessing =",
    nrow(bind_rows(data_list)), "\n")

t_start <- 60
t_end   <- 380
q       <- 60
p       <- 10

mblm_eqn <- function(df) {
  if (nrow(df) < 5) return(c(NA, NA))
  model <- mblm::mblm(Rs_annual ~ Study_midyear, df)
  c(coef(model)[2], summary(model)$coefficients[2, 4])
}

moving_window_soc <- function(data, label) {
  SOC.Area <- Slope <- P.value <- Num <- Mid.SOC <- SD <- SE <- c()
  t <- t_start
  while (t <= t_end) {
    win <- data %>%
      filter(SOC_stock >= (t - q) & SOC_stock < t)
    eq <- mblm_eqn(win)
    SOC.Area <- c(SOC.Area, paste0(t - q, "-", t))
    Mid.SOC  <- c(Mid.SOC, t - q / 2)
    Num      <- c(Num, nrow(win))
    Slope    <- c(Slope, eq[1])
    P.value  <- c(P.value, eq[2])
    sd_rs <- sd(win$Rs_annual, na.rm = TRUE)
    se_rs <- sd_rs / sqrt(nrow(win))
    SD <- c(SD, sd_rs)
    SE <- c(SE, se_rs)
    t <- t + p
  }
  out <- data.frame(SOC.Area, Slope, P.value, Num, Mid.SOC, SD, SE) %>%
    filter(!is.na(Slope), Num >= 40) %>%   # <-- apply cutoff here
    mutate(
      Dataset = label,
      SOC.Area = factor(SOC.Area, levels = SOC.Area),
      sig = case_when(
        P.value < 0.001 ~ "***",
        P.value < 0.01  ~ "**",
        P.value < 0.05  ~ "*",
        TRUE ~ ""
      )
    )
  return(out)
}

mw_all      <- moving_window_soc(data_all, "Combined")
mw_nonkarst <- moving_window_soc(data_nonkarst, "Non-karst")
mw_karst    <- moving_window_soc(data_karst, "Karst")

mw_table <- bind_rows(mw_all, mw_nonkarst, mw_karst)

supplementary_table <- mw_table %>%
  select(
    Dataset,
    `SOC Range (Mg ha⁻¹)` = SOC.Area,
    `Mid-SOC (Mg ha⁻¹)` = Mid.SOC,
    `Slope (g C m⁻² yr⁻²)` = Slope,
    `P-value` = P.value,
    `SD (g C m⁻² yr⁻¹)` = SD,
    `SE (g C m⁻² yr⁻¹)` = SE,
    `N` = Num,
    `Significance` = sig
  ) %>%
  mutate(
    `Slope (g C m⁻² yr⁻²)` = round(`Slope (g C m⁻² yr⁻²)`, 3),
    `P-value` = round(`P-value`, 4),
    `SD (g C m⁻² yr⁻¹)` = round(`SD (g C m⁻² yr⁻¹)`, 2),
    `SE (g C m⁻² yr⁻¹)` = round(`SE (g C m⁻² yr⁻¹)`, 2)
  ) %>%
  arrange(Dataset, `Mid-SOC (Mg ha⁻¹)`)

if (!dir.exists("figures")) dir.create("figures")
write.csv(supplementary_table, 
          "figures/Supplementary_Table_SOC_Moving_Window.csv", 
          row.names = FALSE)

cols <- c("#0000a2", "#bc272d", "#e9c716", "#50ad9f")
mycolors <- colorRampPalette(cols)

TS <- theme(
  text = element_text(size = 8),
  legend.position = "none",
  axis.text.x = element_text(angle = 45, size = 7, hjust = 1),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
  axis.text.y = element_text(size = 7),
  panel.grid = element_blank(),
  plot.title = element_blank(),  # No titles
  axis.title = element_text(size = 9),
  plot.margin = margin(5, 5, 5, 5, "pt")
)

plot_soc_consistent <- function(df, y_limits = NULL) {
  # Determine y-axis limits from all datasets if not provided
  if (is.null(y_limits)) {
    all_slopes <- c(mw_all$Slope, mw_nonkarst$Slope, mw_karst$Slope)
    y_max <- max(all_slopes, na.rm = TRUE) + 15
    y_min <- min(all_slopes, na.rm = TRUE) - 15
    y_limits <- c(y_min, y_max)
  }
  
  star_y <- df$Slope + 3
  count_y <- max(df$Slope, na.rm = TRUE) + 10
  
  p <- ggplot(df, aes(x = SOC.Area, y = Slope, fill = SOC.Area)) +
    geom_col(width = 0.8, alpha = 0.85, color = "black", linewidth = 0.1) +
    geom_text(aes(label = sig, y = star_y), size = 2.5) +
    geom_text(aes(label = Num, y = count_y), size = 2.0) +
    scale_fill_manual(values = mycolors(nrow(df))) +
    scale_y_continuous(limits = y_limits) +
    TS +
    labs(
      x = expression(paste("SOC Stock (Mg ha"^{-1}, ")")),
      y = expression(paste("Change in R"["s"], " (g C m"^{-2}, " yr"^{-2}, ")"))
    )
  
  return(p)
}

all_slopes <- c(mw_all$Slope, mw_nonkarst$Slope, mw_karst$Slope)
y_max <- max(all_slopes, na.rm = TRUE) + 15
y_min <- min(all_slopes, na.rm = TRUE) - 15
y_limits <- c(y_min, y_max)

p_all <- plot_soc_consistent(mw_all, y_limits)
p_nonkarst <- plot_soc_consistent(mw_nonkarst, y_limits)
p_karst <- plot_soc_consistent(mw_karst, y_limits)

combined_figure <- gridExtra::arrangeGrob(
  p_all,
  p_nonkarst,
  p_karst,
  nrow = 3,
  heights = c(1, 1, 1)
)

ggsave("figures/SOC_Moving_Window_All_Datasets.png", 
       plot = combined_figure, 
       width = 6.5, 
       height = 8, 
       dpi = 900,
       bg = "white")

# Save as TIFF
tiff("figures/SOC_Moving_Window_All_Datasets.tiff",
     width = 6.5, 
     height = 8, 
     units = "in", 
     res = 900,
     compression = "lzw",
     bg = "white")
grid::grid.draw(combined_figure)
dev.off()








#===================================SON BINS===========================================
library(dplyr)
library(ggplot2)
library(ggthemes)
library(readr)
library(gridExtra)
library(grid)

datasets <- list(
  "Combined"  = "03_combined_with_koppen.csv",
  "Non-karst" = "03_nonkarst_with_koppen.csv",
  "Karst"     = "03_karst_with_koppen.csv"
)

prep_filter <- function(path, label) {

  read_csv(path, show_col_types = FALSE) %>%
    filter(
      Study_midyear >= 1990,
      Study_midyear <= 2022,
      Meas_method %in% c("IRGA", "Gas Chromatography"),
      Ecosystem_type != "Cropland",
      Outlier == 0,
      !is.na(Rs_annual),
      !is.na(SOC_stock),

      Rs_annual > 0,
      Rs_annual <= 5000
    ) %>%
    mutate(
      Source = label,
      Year   = Study_midyear
    )
}

data_all <- bind_rows(
  lapply(names(datasets), function(name) {
    prep_filter(datasets[[name]], name)
  })
)

data_all <- data_all %>%
  mutate(
    SOC_bin = case_when(
      Source == "Karst" & SOC_stock < 100  ~ "0–100",
      Source == "Karst" & SOC_stock >= 100 ~ "100+",

      Source != "Karst" & SOC_stock < 100 ~ "0–100",
      Source != "Karst" & SOC_stock >= 100 & SOC_stock < 180 ~ "100–180",
      Source != "Karst" & SOC_stock >= 180 & SOC_stock < 270 ~ "180–270",
      Source != "Karst" & SOC_stock >= 270 ~ "270+",

      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(SOC_bin))

cat("\nDEBUG: SOC bin counts\n")
print(table(data_all$Source, data_all$SOC_bin))

results <- data.frame()

fit_soc_bin <- function(df, dataset_label, soc_label) {

  if (nrow(df) < 5) return(NULL)

  fit <- lm(Rs_annual ~ Year, data = df)
  sm  <- summary(fit)

  pval <- sm$coefficients["Year", "Pr(>|t|)"]

  results <<- rbind(
    results,
    data.frame(
      Dataset   = dataset_label,
      SOC_bin   = soc_label,
      Intercept = coef(fit)[1],
      Slope     = coef(fit)[2],
      P_value   = pval,
      SD        = sd(df$Rs_annual, na.rm = TRUE),
      SE        = sd(df$Rs_annual, na.rm = TRUE) / sqrt(nrow(df)),
      N         = nrow(df),
      Signif_code = ifelse(pval < 0.001, "***",
                           ifelse(pval < 0.01, "**",
                                  ifelse(pval < 0.05, "*", ""))),
      stringsAsFactors = FALSE
    )
  )
}

for (src in unique(data_all$Source)) {
  for (soc in unique(data_all$SOC_bin[data_all$Source == src])) {
    fit_soc_bin(
      data_all %>% filter(Source == src, SOC_bin == soc),
      src, soc
    )
  }
}

write.csv(results, "SOC_Rs_SD_SE_FixedBins_Karst2bins.csv", row.names = FALSE)

nonkarst_bins <- c("0–100", "100–180", "180–270", "270+")
karst_bins <- c("0–100", "100+")

nonkarst_bins_display <- paste0(nonkarst_bins, " (Mg C ha⁻¹)")
karst_bins_display <- paste0(karst_bins, " (Mg C ha⁻¹)")

soc_colors <- c(
  "0–100"   = "#0000a2",
  "100–180" = "#bc272d",
  "180–270" = "#e9c716",
  "270+"    = "#50ad9f",
  "100+"    = "#bc272d"
)

get_line_color <- function(dataset_label, soc_label) {
  row <- results %>% 
    filter(Dataset == dataset_label, as.character(SOC_bin) == soc_label)
  
  if (nrow(row) == 0) return("grey70")
  if (row$P_value[1] >= 0.05) return("grey70")
  
  return(soc_colors[[soc_label]])
}

plot_list_combined_nonkarst <- list()
datasets_plot1 <- c("Combined", "Non-karst")

for (dataset_idx in seq_along(datasets_plot1)) {
  dataset <- datasets_plot1[dataset_idx]
  
  for (soc_idx in seq_along(nonkarst_bins)) {
    soc <- nonkarst_bins[soc_idx]
    
    # Filter data for this combination
    plot_data <- data_all %>%
      filter(Source == dataset, as.character(SOC_bin) == soc)
    
    if (nrow(plot_data) >= 5) {
      # Create the plot
      p <- ggplot(plot_data, aes(x = Year, y = Rs_annual)) +
        geom_point(alpha = 0.25, size = 1.5, color = "grey40") +
        geom_smooth(
          method = "lm", 
          se = FALSE,
          color = get_line_color(dataset, soc),
          linewidth = 1.2
        ) +
        theme_few(base_size = 10) +
        labs(
          x = ifelse(dataset_idx == length(datasets_plot1), "Year", ""),
          y = ifelse(soc_idx == 1, 
                     expression(R[s]~"(g C m"^{-2}~"yr"^{-1}*")"), 
                     "")
        ) +
        scale_y_continuous(
          limits = c(0, 2500),
          breaks = seq(0, 2500, by = 500)
        ) +
        scale_x_continuous(
          limits = c(1990, 2022),
          breaks = seq(1990, 2020, by = 10)
        ) +
        theme(
          plot.margin = margin(5, 5, 5, 5, "pt"),
          plot.title = element_blank(),
          axis.title = element_text(size = 10)
        )
      
      plot_name <- paste0(dataset, "_", gsub("[^0-9A-Za-z]+", "_", soc))
      plot_list_combined_nonkarst[[plot_name]] <- p
    } else {
      # Create empty plot
      p <- ggplot() +
        theme_void() +
        theme(plot.background = element_rect(fill = "white", color = NA))
      
      plot_name <- paste0(dataset, "_", gsub("[^0-9A-Za-z]+", "_", soc))
      plot_list_combined_nonkarst[[plot_name]] <- p
    }
  }
}

layout_matrix_plot1 <- matrix(1:8, nrow = 2, ncol = 4, byrow = TRUE)

top_titles_plot1 <- lapply(nonkarst_bins_display, function(title) {
  grid::textGrob(
    title,
    gp = grid::gpar(fontsize = 11, fontface = "bold")
  )
})

top_row_plot1 <- gridExtra::arrangeGrob(
  grobs = top_titles_plot1,
  nrow = 1
)

main_grid_plot1 <- gridExtra::arrangeGrob(
  grobs = plot_list_combined_nonkarst,
  layout_matrix = layout_matrix_plot1,
  widths = rep(1, 4),
  heights = rep(1, 2)
)

figure1 <- gridExtra::arrangeGrob(
  # SOC bin titles
  top_row_plot1,
  
  # Main plots grid
  main_grid_plot1,
  
  nrow = 2,
  heights = c(0.5, 6)
)

ggsave(
  "SOC_Rs_Trends_Combined_Nonkarst.png",
  plot = figure1,
  width = 14,
  height = 7,
  dpi = 900,
  bg = "white"
)

ggsave(
  "SOC_Rs_Trends_Combined_Nonkarst.pdf",
  plot = figure1,
  width = 14,
  height = 7,
  bg = "white"
)

tiff(
  "SOC_Rs_Trends_Combined_Nonkarst.tiff",
  width = 14, 
  height = 7, 
  units = "in", 
  res = 900,
  compression = "lzw",
  bg = "white"
)
grid::grid.draw(figure1)
dev.off()

plot_list_karst <- list()

for (soc_idx in seq_along(karst_bins)) {
  soc <- karst_bins[soc_idx]
  
  plot_data <- data_all %>%
    filter(Source == "Karst", as.character(SOC_bin) == soc)
  
  if (nrow(plot_data) >= 5) {
    p <- ggplot(plot_data, aes(x = Year, y = Rs_annual)) +
      geom_point(alpha = 0.25, size = 1.5, color = "grey40") +
      geom_smooth(
        method = "lm", 
        se = FALSE,
        color = get_line_color("Karst", soc),
        linewidth = 1.2
      ) +
      theme_few(base_size = 10) +
      labs(
        x = "Year",
        y = expression(R[s]~"(g C m"^{-2}~"yr"^{-1}*")")
      ) +
      scale_y_continuous(
        limits = c(0, 2500),
        breaks = seq(0, 2500, by = 500)
      ) +
      scale_x_continuous(
        limits = c(1990, 2022),
        breaks = seq(1990, 2020, by = 10)
      ) +
      theme(
        plot.margin = margin(5, 5, 5, 5, "pt"),
        plot.title = element_blank(),
        axis.title = element_text(size = 10)
      )
    
    plot_name <- paste0("Karst_", gsub("[^0-9A-Za-z]+", "_", soc))
    plot_list_karst[[plot_name]] <- p
  } else {
    # Create empty plot
    p <- ggplot() +
      theme_void() +
      theme(plot.background = element_rect(fill = "white", color = NA))
    
    plot_name <- paste0("Karst_", gsub("[^0-9A-Za-z]+", "_", soc))
    plot_list_karst[[plot_name]] <- p
  }
}

for (empty_idx in (length(karst_bins) + 1):4) {
  plot_name <- paste0("Karst_empty_", empty_idx)
  plot_list_karst[[plot_name]] <- ggplot() + theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))
}

layout_matrix_plot2 <- matrix(1:4, nrow = 1, ncol = 4)

top_titles_plot2 <- list()

for (col in 1:4) {
  if (col <= length(karst_bins_display)) {
    title_text <- karst_bins_display[col]
  } else {
    title_text <- ""
  }
  
  top_titles_plot2[[col]] <- grid::textGrob(
    title_text,
    gp = grid::gpar(fontsize = 11, fontface = "bold")
  )
}

top_row_plot2 <- gridExtra::arrangeGrob(
  grobs = top_titles_plot2,
  nrow = 1
)

main_grid_plot2 <- gridExtra::arrangeGrob(
  grobs = plot_list_karst,
  layout_matrix = layout_matrix_plot2,
  widths = rep(1, 4),
  heights = 1
)

figure2 <- gridExtra::arrangeGrob(
  top_row_plot2,
  
  main_grid_plot2,
  
  nrow = 2,
  heights = c(0.5, 6)
)

ggsave(
  "SOC_Rs_Trends_Karst.png",
  plot = figure2,
  width = 14,
  height = 5,
  dpi = 900,
  bg = "white"
)

ggsave(
  "SOC_Rs_Trends_Karst.pdf",
  plot = figure2,
  width = 14,
  height = 5,
  bg = "white"
)

tiff(
  "SOC_Rs_Trends_Karst.tiff",
  width = 14, 
  height = 5, 
  units = "in", 
  res = 900,
  compression = "lzw",
  bg = "white"
)
grid::grid.draw(figure2)
dev.off()





#===================================SOC and Climate==============================================
library(dplyr)
library(ggplot2)
library(ggthemes)
library(readr)

datasets <- list(
  "Combined"  = "03_combined_with_koppen.csv",
  "Non-karst" = "03_nonkarst_with_koppen.csv",
  "Karst"     = "03_karst_with_koppen.csv"
)

climate_levels <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")

prep_filter <- function(path, label) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(
      Study_midyear >= 1990,
      Study_midyear <= 2022,
      Meas_method %in% c("IRGA", "Gas Chromatography"),
      Ecosystem_type != "Cropland",
      Outlier == 0,
      !is.na(Rs_annual),
      !is.na(SOC_stock),
      !is.na(kg_description),
      Rs_annual > 0,
      Rs_annual <= 5000
    ) %>%
    mutate(
      Source = label,
      Year   = Study_midyear,
      ClimateGroup = case_when(
        grepl("^Tropical",  kg_description, ignore.case = TRUE) ~ "Tropical",
        grepl("^Arid",      kg_description, ignore.case = TRUE) ~ "Arid",
        grepl("^Temperate", kg_description, ignore.case = TRUE) ~ "Temperate",
        grepl("^Cold",      kg_description, ignore.case = TRUE) ~ "Cold",
        grepl("^Polar",     kg_description, ignore.case = TRUE) ~ "Polar",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(ClimateGroup)) %>%
    mutate(ClimateGroup = factor(ClimateGroup, levels = climate_levels))
}

add_breakpoint_periods <- function(df, dataset_name) {
  if (dataset_name %in% c("Combined", "Non-karst")) {
    df %>% mutate(
      Period = case_when(
        Year < 2013 ~ "Pre-breakpoint",
        Year <= 2016 ~ "Breakpoint",
        TRUE ~ "Post-breakpoint"
      )
    )
  } else {
    df %>% mutate(
      Period = case_when(
        Year < 2009 ~ "Pre-breakpoint",
        Year <= 2014 ~ "Breakpoint",
        TRUE ~ "Post-breakpoint"
      )
    )
  }
}

data_all <- bind_rows(
  lapply(names(datasets), function(name) {
    df <- prep_filter(datasets[[name]], name)
    add_breakpoint_periods(df, name)
  })
)

data_all <- data_all %>%
  mutate(
    SOC_bin = case_when(
      Source == "Karst" & SOC_stock < 100  ~ "SOC 0–100",
      Source == "Karst" & SOC_stock >= 100 ~ "SOC 100+",
      Source != "Karst" & SOC_stock < 100 ~ "SOC 0–100",
      Source != "Karst" & SOC_stock >= 100 & SOC_stock < 180 ~ "SOC 100–180",
      Source != "Karst" & SOC_stock >= 180 & SOC_stock < 270 ~ "SOC 180–270",
      Source != "Karst" & SOC_stock >= 270 ~ "SOC 270+",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(SOC_bin))

results <- data.frame()

fit_model <- function(df, src, soc, clim, per) {
  if (nrow(df) < 5) return(NULL)
  fit <- lm(Rs_annual ~ Year, data = df)
  sm  <- summary(fit)
  results <<- rbind(
    results,
    data.frame(
      Dataset = src,
      ClimateZone = clim,
      SOC_bin = soc,
      Period = per,
      Intercept = coef(fit)[1],
      Slope = coef(fit)[2],
      P_value = sm$coefficients["Year","Pr(>|t|)"],
      SD = sd(df$Rs_annual, na.rm = TRUE),
      SE = sd(df$Rs_annual, na.rm = TRUE) / sqrt(nrow(df)),
      N = nrow(df),
      stringsAsFactors = FALSE
    )
  )
}

for (src in unique(data_all$Source)) {
  for (clim in climate_levels) {
    for (soc in unique(data_all$SOC_bin)) {
      for (per in unique(data_all$Period)) {
        fit_model(
          data_all %>% filter(Source==src, ClimateGroup==clim,
                              SOC_bin==soc, Period==per),
          src, soc, clim, per
        )
      }
    }
  }
}

write.csv(results,
          "SOC_Rs_Trends_ByClimate_SOC_Breakpoints.csv",
          row.names = FALSE)

results_clean <- results %>%
  select(-Intercept, -SD) %>%
  mutate(
    Slope = round(Slope, 2),
    SE = round(SE, 2),
    P_value = round(P_value, 4),
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01 ~ "**",
      P_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(Dataset, ClimateZone, SOC_bin, Period)

write.csv(results_clean, 
          "Supplementary_Table_SOC_Rs_Trends_Clean.csv", 
          row.names = FALSE)

period_cols <- c(
  "Pre-breakpoint"  = "blue",
  "Breakpoint"      = "red",
  "Post-breakpoint" = "darkgreen"
)

plot_fun <- function(df, src, clim, soc) {

  sig_periods <- results %>%
    filter(Dataset==src, ClimateZone==clim,
           SOC_bin==soc, P_value <= 0.05) %>%
    pull(Period)

  if (length(sig_periods)==0) return(NULL)

  ggplot(df, aes(Year, Rs_annual)) +
    geom_point(alpha=0.25, size=1, color="grey70") +
    geom_smooth(
      data = df %>% filter(Period %in% sig_periods),
      aes(color = Period),
      method="lm", se=FALSE, linewidth=1
    ) +
    scale_color_manual(values = period_cols) +
    scale_y_continuous(limits=c(0,2500)) +
    theme_few() +
    theme(legend.position="none") +
    labs(
      title=paste(src, clim, soc),
      x="Year",
      y=expression(R[s]~"(g C m"^{-2}~"yr"^{-1}*")")
    )
}

dir.create("figures_SOC_climate_breakpoints", showWarnings = FALSE)

for (src in unique(data_all$Source)) {
  for (clim in climate_levels) {
    for (soc in unique(data_all$SOC_bin)) {

      dfp <- data_all %>%
        filter(Source==src, ClimateGroup==clim, SOC_bin==soc)

      p <- plot_fun(dfp, src, clim, soc)

      if (!is.null(p)) {
        ggsave(
          paste0("figures_SOC_climate_breakpoints/",
                 src,"_",clim,"_",gsub("[^A-Za-z0-9]+","_",soc),".png"),
          plot=p, width=2.8, height=2.2, dpi=900
        )
      }
    }
  }
}

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


source("functions_plot.R")
source("statistical_prep.R")

cols <- c("#0000a2", "#bc272d", "#e9c716", "#50ad9f")
mycolors <- colorRampPalette(cols)

rcdata = dplyr::select(work.data, SOC_stock, Rs_annual, Study_midyear) 
t = q = 60
p = 10
i = 1
ord = c()
SOC.Area = c()
Slope = c()
R.square = c()
P.value = c()
Num = c()
Mid.SOC = c()

while (t <= 380) {
  ord = c(ord, i)
  win = filter(rcdata, SOC_stock >= (t - q) & SOC_stock < (t))
  eq = mblm_eqn(df = win)
  SOC.Area = c(SOC.Area, paste0((t - q), '-', (t)))
  Slope = c(Slope, eq[1])
  P.value = c(P.value, eq[2])
  Num = c(Num, nrow(win))
  Mid.SOC = c(Mid.SOC, t - (0.5 * q))
  t = t + p
  i = i + 1
}

winresult = data.frame(matrix(nrow = length(ord), ncol = 0))
winresult$Order = ord
winresult$SOC.Area = SOC.Area
winresult$Slope = as.numeric(Slope)
winresult$P.value = as.numeric(P.value)
winresult$Num = Num
winresult$Mid.SOC = Mid.SOC

winresult$SOC.Area <- factor(winresult$SOC.Area, levels = c(winresult$SOC.Area))
winresult$sig = ''
winresult[winresult$P.value < 0.1, which(colnames(winresult) == 'sig')] = ''
winresult[winresult$P.value < 0.05, which(colnames(winresult) == 'sig')] = '*'
winresult[winresult$P.value < 0.001, which(colnames(winresult) == 'sig')] = '**'
winresult[winresult$P.value < 0.001, which(colnames(winresult) == 'sig')] = '***'

TS <- theme(text = element_text(size = 6),
            legend.position = "none",
            axis.text.x = element_text(angle = 45, size = 5, hjust = 1, vjust = 1),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
            axis.text.y = element_text(size = 6))

n_rows <- nrow(winresult)
star_p <- as.numeric(Slope)

if (n_rows >= 11) {
  star_p[1:11] <- star_p[1:11] - 4.5
}
if (n_rows >= 33) {
  star_p[12:33] <- star_p[12:33] + 1
}
if (n_rows >= 27) {
  star_p[20:27] <- star_p[20:27] - 3.5
}

star <- geom_text(aes(x = SOC.Area, y = star_p, label = sig), size = 2.1)
winresult$Num.italic = paste0('italic(', winresult$Num, ')')
count <- geom_text(aes(x = SOC.Area, y = 100, label = Num.italic), size = 1.7, parse = TRUE)

SFM = scale_fill_manual(values = mycolors(nrow(winresult)))

barplot = ggplot(winresult, aes(x = SOC.Area, y = Slope)) + 
  geom_bar(aes(fill = SOC.Area), stat = "identity", alpha = 0.8, color = "navy", linewidth = 0.1) +
  TS + 
  SFM + 
  star + 
  count + 
  labs(
    x = expression(paste("SOC Stock (g C m"^{-2}, ")")),
    y = expression(paste("Change in R"["s"], " (g C m"^{-2}, " yr"^{-2}, ")"))
  ) +
  ggtitle("(c)") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 6, margin = margin(t = 2)),
    axis.title.y = element_text(size = 6, margin = margin(r = 2))
  )

if (!dir.exists("figures")) {
  dir.create("figures")
}

ggsave("figures/fig3_c.pdf", barplot, width = 6.5, height = 2, units = 'in', dpi = 900)
ggsave("figures/fig3_c.png", barplot, width = 6.5, height = 2, units = 'in', dpi = 900)
ggsave("figures/fig3_c.tiff", barplot, width = 6.5, height = 2, units = 'in', dpi = 900, compression = "lzw")

cat("SOC window analysis plot created successfully!\n")
cat("Number of SOC bins:", n_rows, "\n")
cat("- X-axis: SOC Stock (g C m⁻²)\n")
cat("- Y-axis: Change in Rₛ (g C m⁻² yr⁻²)\n")
cat("Files saved:\n")
cat("- figures/fig3_c.pdf\n")
cat("- figures/fig3_c.png\n")
cat("- figures/fig3_c.tiff\n")

#blue_colors <- c('#E6F3FF', '#80BFFF', '#4DA6FF', '#0073E6')

blue_colors <- c('#1B4F72', '#2874A6', '#5DADE2', '#85C1E9')

year_min <- floor(min(work.data$Study_midyear, na.rm = TRUE))
year_max <- ceiling(max(work.data$Study_midyear, na.rm = TRUE))
cat("Actual year range in your data:", year_min, "to", year_max, "\n")

year_breaks <- seq(1990, 2020, by = 10)  # 1990, 2000, 2010, 2020

# 0-100
SOC1 = filter(work.data, SOC_stock < 100)
summary(lm(Rs_annual ~ Study_midyear, data = SOC1))
(p <- ggplot(aes_string(x = "Study_midyear", y = "Rs_annual"), data = SOC1) + 
    geom_point(color = blue_colors[1], size = 1.6, alpha = 0.3, shape = 16, stroke = 0.05) + 
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = blue_colors[1], size = 1.3) +
    theme_few() +
    labs(
      x = "Year",
      y = expression(paste("R"["s"], " (g C m"^{-2}, " yr"^{-1}, ")"))
    ) +
    ggtitle(expression(paste("(a) SOC 0-100 g C m"^{-2}))) +
    scale_x_continuous(breaks = year_breaks, limits = c(year_min, year_max)) +
    scale_y_continuous(breaks = c(0, 800, 1600, 2400, 3200), limits = c(0, 3400)) +
    theme(text = element_text(size = 8), 
          panel.border = element_rect(color = "grey3", linewidth = 0.5)) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in")))
ggsave("figures/SOC_0-100.pdf", p, width = 2.2, height = 2, units = 'in', dpi = 900)
ggsave("figures/SOC_0-100.png", p, width = 2.2, height = 2, units = 'in', dpi = 900)

# 100-180
SOC2 = filter(work.data, SOC_stock >= 100 & SOC_stock < 180)
summary(lm(Rs_annual ~ Study_midyear, data = SOC2))
(p <- ggplot(aes_string(x = "Study_midyear", y = "Rs_annual"), data = SOC2) + 
    geom_point(color = blue_colors[2], size = 1.6, alpha = 0.3, shape = 16, stroke = 0.05) + 
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = blue_colors[2], size = 1.3) +
    theme_few() +
    labs(
      x = "Year",
      y = expression(paste("R"["s"], " (g C m"^{-2}, " yr"^{-1}, ")"))
    ) +
    ggtitle(expression(paste("(b) SOC 100-180 g C m"^{-2}))) +
    scale_x_continuous(breaks = year_breaks, limits = c(year_min, year_max)) +
    scale_y_continuous(breaks = c(0, 700, 1400, 2100, 2800), limits = c(0, 2900)) +
    theme(text = element_text(size = 8), 
          panel.border = element_rect(color = "grey3", linewidth = 0.5)) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in")))
ggsave("figures/SOC_100-180.pdf", p, width = 2.2, height = 2, units = 'in', dpi = 900)
ggsave("figures/SOC_100-180.png", p, width = 2.2, height = 2, units = 'in', dpi = 900)

# 180-270
SOC3 = filter(work.data, SOC_stock >= 180 & SOC_stock < 270)
summary(lm(Rs_annual ~ Study_midyear, data = SOC3))
(p <- ggplot(aes_string(x = "Study_midyear", y = "Rs_annual"), data = SOC3) + 
    geom_point(color = blue_colors[3], size = 1.6, alpha = 0.4, shape = 16, stroke = 0.05) + 
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = blue_colors[3], size = 1.3, linetype = "dashed") +
    theme_few() +
    labs(
      x = "Year",
      y = expression(paste("R"["s"], " (g C m"^{-2}, " yr"^{-1}, ")"))
    ) +
    ggtitle(expression(paste("(c) SOC 180-270 g C m"^{-2}))) +
    scale_x_continuous(breaks = year_breaks, limits = c(year_min, year_max)) +
    scale_y_continuous(breaks = c(0, 700, 1400, 2100, 2800), limits = c(0, 2900)) +
    theme(text = element_text(size = 8), 
          panel.border = element_rect(color = "grey3", linewidth = 0.5)) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in")))
ggsave("figures/SOC_180-270.pdf", p, width = 2.2, height = 2, units = 'in', dpi = 900)
ggsave("figures/SOC_180-270.png", p, width = 2.2, height = 2, units = 'in', dpi = 900)

# 270+
SOC4 = filter(work.data, SOC_stock >= 270)
summary(lm(Rs_annual ~ Study_midyear, data = SOC4))
(p <- ggplot(aes_string(x = "Study_midyear", y = "Rs_annual"), data = SOC4) + 
    geom_point(color = blue_colors[4], size = 1.6, alpha = 0.6, shape = 16, stroke = 0.05) + 
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = blue_colors[4], size = 1.3) +
    theme_few() +
    labs(
      x = "Year",
      y = expression(paste("R"["s"], " (g C m"^{-2}, " yr"^{-1}, ")"))
    ) +
    ggtitle(expression(paste("(d) SOC 270+ g C m"^{-2}))) +
    scale_x_continuous(breaks = year_breaks, limits = c(year_min, year_max)) +
    scale_y_continuous(breaks = c(0, 700, 1400, 2100, 2800), limits = c(0, 2900)) +
    theme(text = element_text(size = 8), 
          panel.border = element_rect(color = "grey3", linewidth = 0.5)) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in")))
ggsave("figures/SOC_270+.pdf", p, width = 2.2, height = 2, units = 'in', dpi = 900)
ggsave("figures/SOC_270+.png", p, width = 2.2, height = 2, units = 'in', dpi = 900)

cat("SOC interval plots created with proper time range (", year_min, "-", year_max, "):\n")
cat("- X-axis: Year (with breaks at:", paste(year_breaks, collapse = ", "), ")\n")
cat("- Y-axis: Rₛ (g C m⁻² yr⁻¹)\n")
cat("- Titles include SOC ranges with proper units\n")
cat("Files saved in 'figures' folder:\n")
cat("- SOC_0-100.pdf/png\n")
cat("- SOC_100-180.pdf/png\n")
cat("- SOC_180-270.pdf/png\n")
cat("- SOC_270+.pdf/png\n")

# SOC interval analysis with time period breakdown
# Define time periods including full period
time_periods <- list(
  list(name = "1990-2000", start = 1990, end = 2000),
  list(name = "2000-2010", start = 2000, end = 2010),
  list(name = "2010-2022", start = 2010, end = 2022),
  list(name = "1990-2022", start = 1990, end = 2022)
)

soc_intervals <- list(
  list(name = "SOC 0-100 g C m⁻²", data = SOC1),
  list(name = "SOC 100-180 g C m⁻²", data = SOC2),
  list(name = "SOC 180-270 g C m⁻²", data = SOC3),
  list(name = "SOC 270+ g C m⁻²", data = SOC4)
)

extract_regression_results <- function(soc_data, interval_name) {
  if(nrow(soc_data) > 2) {
    model <- lm(Rs_annual ~ Study_midyear, data = soc_data)
    summary_model <- summary(model)
    
    return(data.frame(
      SOC_Interval = interval_name,
      Slope = coef(model)[2],
      P_Value = summary_model$coefficients[2, 4],
      N_Observations = nrow(soc_data),
      R_Squared = summary_model$r.squared,
      Significant = ifelse(summary_model$coefficients[2, 4] < 0.05, "YES", "NO"),
      stringsAsFactors = FALSE
    ))
  }
  return(NULL)
}

# Function to analyze by time period
analyze_time_period <- function(soc_data, interval_name, period) {
  period_data <- soc_data %>% 
    filter(Study_midyear >= period$start & Study_midyear <= period$end)
  
  if(nrow(period_data) > 2) {
    model <- lm(Rs_annual ~ Study_midyear, data = period_data)
    summary_model <- summary(model)
    
    return(data.frame(
      SOC_Interval = interval_name,
      Time_Period = period$name,
      Slope = coef(model)[2],
      P_Value = summary_model$coefficients[2, 4],
      N_Observations = nrow(period_data),
      R_Squared = summary_model$r.squared,
      Significant = ifelse(summary_model$coefficients[2, 4] < 0.05, "YES", "NO"),
      stringsAsFactors = FALSE
    ))
  }
  return(NULL)
}

results_full <- do.call(rbind, lapply(1:4, function(i) {
  extract_regression_results(soc_intervals[[i]]$data, soc_intervals[[i]]$name)
}))

write.csv(results_full, "SOC_interval_full_period_results.csv", row.names = FALSE)

time_period_results <- do.call(rbind, lapply(1:4, function(i) {
  soc_interval <- soc_intervals[[i]]
  period_results <- do.call(rbind, lapply(time_periods, function(period) {
    analyze_time_period(soc_interval$data, soc_interval$name, period)
  }))
  period_results
}))

for(period in time_periods) {
  period_data <- time_period_results %>% filter(Time_Period == period$name)
  filename <- paste0("SOC_interval_", gsub("-", "_", period$name), "_results.csv")
  write.csv(period_data, filename, row.names = FALSE)
}

write.csv(time_period_results, "SOC_interval_regression_results_comprehensive.csv", row.names = FALSE)

detailed_summary <- do.call(rbind, lapply(1:4, function(i) {
  soc_data <- soc_intervals[[i]]$data
  model <- lm(Rs_annual ~ Study_midyear, data = soc_data)
  summary_model <- summary(model)
  
  data.frame(
    SOC_Interval = soc_intervals[[i]]$name,
    Slope = coef(model)[2],
    P_Value = summary_model$coefficients[2, 4],
    N_Observations = nrow(soc_data),
    R_Squared = summary_model$r.squared,
    Significant = ifelse(summary_model$coefficients[2, 4] < 0.05, "YES", "NO"),
    Mean_Rs = mean(soc_data$Rs_annual, na.rm = TRUE),
    SD_Rs = sd(soc_data$Rs_annual, na.rm = TRUE),
    Years_Covered = paste(round(min(soc_data$Study_midyear, na.rm = TRUE), 1), 
                         "-", round(max(soc_data$Study_midyear, na.rm = TRUE), 1)),
    stringsAsFactors = FALSE
  )
}))

write.csv(detailed_summary, "SOC_interval_detailed_summary.csv", row.names = FALSE)

summary_table <- time_period_results %>%
  select(SOC_Interval, Time_Period, Slope, P_Value, Significant, N_Observations, R_Squared)

write.csv(summary_table, "SOC_interval_summary_table.csv", row.names = FALSE)

soc_summary <- results_full %>%
  mutate(SOC_Interval = gsub("SOC ", "", SOC_Interval)) %>%
  group_by(SOC_Interval) %>%
  summarise(
    Total_Observations = sum(N_Observations),
    Total_Rs_Percentage = NA, # Placeholder for additional calculations
    Significant_Relationships = sum(Significant == "YES"),
    Mean_Slope = mean(Slope)
  )

write.csv(soc_summary, "SOC_interval_summary.csv", row.names = FALSE)

time_summary <- time_period_results %>%
  group_by(Time_Period) %>%
  summarise(
    Total_Observations = sum(N_Observations),
    Total_Rs_Percentage = NA, # Placeholder for additional calculations
    Significant_Relationships = sum(Significant == "YES"),
    Mean_Slope = mean(Slope, na.rm = TRUE)
  )

write.csv(time_summary, "SOC_time_period_summary.csv", row.names = FALSE)

comprehensive_time_analysis <- do.call(rbind, lapply(1:4, function(i) {
  soc_interval <- soc_intervals[[i]]
  do.call(rbind, lapply(time_periods, function(period) {
    period_data <- soc_interval$data %>% 
      filter(Study_midyear >= period$start & Study_midyear <= period$end)
    
    if(nrow(period_data) > 2) {
      model <- lm(Rs_annual ~ Study_midyear, data = period_data)
      summary_model <- summary(model)
      
      # Calculate additional metrics
      mean_rs <- mean(period_data$Rs_annual, na.rm = TRUE)
      sd_rs <- sd(period_data$Rs_annual, na.rm = TRUE)
      
      return(data.frame(
        Time_Period = period$name,
        SOC_Interval = soc_interval$name,
        Slope = coef(model)[2],
        P_Value = summary_model$coefficients[2, 4],
        R_Squared = summary_model$r.squared,
        N_Observations = nrow(period_data),
        Significant = ifelse(summary_model$coefficients[2, 4] < 0.05, "YES", "NO"),
        Mean_Rs = mean_rs,
        SD_Rs = sd_rs,
        Rs_Percentage = NA, # Placeholder
        Top_Biomes = NA,    # Placeholder - would require biome data
        Year_Range = paste(min(period_data$Study_midyear, na.rm = TRUE), 
                          "-", max(period_data$Study_midyear, na.rm = TRUE)),
        stringsAsFactors = FALSE
      ))
    }
    return(NULL)
  }))
}))

write.csv(comprehensive_time_analysis, "SOC_time_period_comprehensive_results.csv", row.names = FALSE)

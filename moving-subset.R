#Actually use this one
#fig1
source("0_functions_plot.R")
source("4_statistical_prep.R")

# OLS estimator===============================================================================
### Moving subset analysis for year windows

# Define blue color palette to match karst aquifer map
mycolors = colorRampPalette(c("#0000a2", "#bc272d", "#e9c716", "#50ad9f"))

# Use dplyr::select to avoid conflicts
rcdata = dplyr::select(work.data, Latitude, Rs_annual, Study_midyear)

# Initialize variables for the window analysis
t = q = 10
p = 1
i = 1

# Initialize empty vectors to store the results
ord = c()
Lat.Area = c()
Slope = c()
R.square = c()
P.value = c()
Num = c()
Mid.Lat = c()
Year.range = c()

# Additional vectors for comprehensive statistics
Mean_Rs = c()
Median_Rs = c()
SD_Rs = c()
SE_Rs = c()
Intercept = c()
Adj_R.square = c()

# Adjust the year calculation to focus on the range from 1990 to 2022
while (t <= max(rcdata$Study_midyear - 1990 + p + 0.1)) {
  i = i + 1
  win = filter(rcdata, Study_midyear - 1990 >= (t - q) & (Study_midyear - 1990) < t)
  
  # Calculate comprehensive statistics for each window
  if (nrow(win) > 0) {
    eq = lm_eqn(df = win)
    Year.range = c(Year.range, paste0((1990 + t - q), '-', (1990 + t - 1)))
    Slope = c(Slope, eq[1])
    R.square = c(R.square, eq[2])
    P.value = c(P.value, eq[3])
    Num = c(Num, nrow(win))
    
    # Calculate additional statistics
    Mean_Rs = c(Mean_Rs, mean(win$Rs_annual, na.rm = TRUE))
    Median_Rs = c(Median_Rs, median(win$Rs_annual, na.rm = TRUE))
    SD_Rs = c(SD_Rs, sd(win$Rs_annual, na.rm = TRUE))
    SE_Rs = c(SE_Rs, sd(win$Rs_annual, na.rm = TRUE) / sqrt(nrow(win)))
    
    # Fit linear model to get intercept and adjusted R-squared
    if (nrow(win) > 1) {
      lm_model <- lm(Rs_annual ~ Study_midyear, data = win)
      Intercept = c(Intercept, coef(lm_model)[1])
      Adj_R.square = c(Adj_R.square, summary(lm_model)$adj.r.squared)
    } else {
      Intercept = c(Intercept, NA)
      Adj_R.square = c(Adj_R.square, NA)
    }
  }
  t = t + p
}

# Create the comprehensive result data frame
winresult = data.frame(
  Year.range, 
  Slope = as.numeric(Slope),
  R.square = as.numeric(R.square),
  Adj_R.square = as.numeric(Adj_R.square),
  P.value = as.numeric(P.value),
  Num = as.numeric(Num),
  Mean_Rs = as.numeric(Mean_Rs),
  Median_Rs = as.numeric(Median_Rs),
  SD_Rs = as.numeric(SD_Rs),
  SE_Rs = as.numeric(SE_Rs),
  Intercept = as.numeric(Intercept)
)

# Add significance stars based on p-value
winresult$sig = ''
winresult[winresult$P.value < 0.1, which(colnames(winresult) == 'sig')] = '·'
winresult[winresult$P.value < 0.05, which(colnames(winresult) == 'sig')] = '*'
winresult[winresult$P.value < 0.01, which(colnames(winresult) == 'sig')] = '**'
winresult[winresult$P.value < 0.001, which(colnames(winresult) == 'sig')] = '***'

# Add significance level description
winresult$Significance_Level = case_when(
  winresult$P.value < 0.001 ~ "p < 0.001",
  winresult$P.value < 0.01 ~ "p < 0.01", 
  winresult$P.value < 0.05 ~ "p < 0.05",
  winresult$P.value < 0.1 ~ "p < 0.1",
  TRUE ~ "Not significant"
)

# Filter significant periods (p < 0.05)
significant_periods <- winresult %>% 
  filter(P.value < 0.05)

# Print significant periods
cat("=== SIGNIFICANT PERIODS (p < 0.05) ===\n")
if (nrow(significant_periods) > 0) {
  print(significant_periods)
} else {
  cat("No periods with p < 0.05 significance\n")
}

# Print all periods with comprehensive statistics
cat("\n=== ALL PERIODS - COMPREHENSIVE STATISTICS ===\n")
print(winresult)

# Save comprehensive results to CSV
write.csv(winresult, "moving_window_comprehensive_stats.csv", row.names = FALSE)
cat("\nComprehensive statistics saved to: moving_window_comprehensive_stats.csv\n")

# Save significant periods to separate CSV
write.csv(significant_periods, "significant_periods_stats.csv", row.names = FALSE)
cat("Significant periods statistics saved to: significant_periods_stats.csv\n")

# Create summary table of significant periods
if (nrow(significant_periods) > 0) {
  significant_summary <- significant_periods %>%
    dplyr::select(Year.range, Slope, P.value, R.square, Adj_R.square, Num, Mean_Rs, SD_Rs, SE_Rs, Significance_Level)
  
  cat("\n=== SUMMARY OF SIGNIFICANT PERIODS ===\n")
  print(significant_summary)
  
  # Save summary
  write.csv(significant_summary, "significant_periods_summary.csv", row.names = FALSE)
  cat("Significant periods summary saved to: significant_periods_summary.csv\n")
}

# Continue with the original plotting code...
# Create a custom theme for the plot
TS <- theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, size = 6, hjust = 1, vjust = 1),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "grey3", fill = NA, size = 1),
            axis.text.y = element_text(size = 6))

# Create 'star_p' values based on the slope, ensuring the correct number of elements
star_p = c(as.numeric(Slope)[1:7] + 1, 
           as.numeric(Slope)[8:12] - 3, 
           as.numeric(Slope)[13:18] + 1, 
           as.numeric(Slope)[19:22] + 1)

# Define the geom_text layers for the plot
star <- geom_text(aes(x = Year.range, y = star_p, label = sig), size = 5)
winresult$Num.italic = paste0('italic(', winresult$Num, ')')
count <- geom_text(aes(x = Year.range, y = 54, label = Num.italic), size = 2, parse = TRUE)

# Set color scale and y-axis limits
SFM = scale_fill_manual(values = mycolors(nrow(winresult)))
SYC = scale_y_continuous(breaks = c(-15, 0, 15, 30, 45), limits = c(-15, 55))

# Create the barplot with proper axis labels
barplot = ggplot(winresult, aes(x = Year.range, y = Slope)) + 
  geom_bar(aes(fill = Year.range), stat = "identity", color = "navy", linewidth = 0.2) +  # Blue outline
  TS +  # Apply custom theme
  SFM +  # Apply blue color scale
  star +  # Add significance stars
  count +  # Add counts in italic
  SYC +  # Adjust y-axis
  labs(
    x = "Year",
    y = expression(paste(Delta, "R"["s"], " (g C m"^{-2}, " yr"^{-2}, ")"))
  ) +
  #ggtitle("(d)") +  # Add title
  # Additional styling to match your karst map aesthetic
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 8, margin = margin(t = 5)),
    axis.title.y = element_text(size = 8, margin = margin(r = 5))
  )

# Create the 'figures' folder if it doesn't exist
if (!dir.exists("figures")) {
  dir.create("figures")
}

# Save the plot as a PDF
ggsave("figures/fig1_d.pdf", plot = barplot, width = 6.5, height = 2, units = 'in', dpi = 900)

# Save the plot as a PNG (additional format requested)
ggsave("figures/fig1_d.png", plot = barplot, width = 6.5, height = 2, units = 'in', dpi = 900)

# Also save as high-resolution TIFF for publications
ggsave("figures/fig1_d.tiff", plot = barplot, width = 6.5, height = 2, units = 'in', dpi = 900, 
       compression = "lzw")

cat("\nPlot created with proper axis labels:\n")
cat("- X-axis: Year\n")
cat("- Y-axis: Change in Rₛ (g C m⁻² yr⁻²)\n")
cat("Files saved:\n")
cat("- figures/fig1_d.pdf\n")
cat("- figures/fig1_d.png\n")
cat("- figures/fig1_d.tiff\n")
cat("- moving_window_comprehensive_stats.csv\n")
cat("- significant_periods_stats.csv\n")
if (nrow(significant_periods) > 0) {
  cat("- significant_periods_summary.csv\n")
}


## Moving subset analysis for SOC windows

# Define blue color palette to match karst aquifer map
cols <- c("#0000a2", "#bc272d", "#e9c716", "#50ad9f")  # Light to dark blue gradient
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

winresult$SOC.Area <- factor(winresult$SOC.Area, levels = c(winresult$SOC.Area))  # prevent the reorder of levels 
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

# FIXED: Dynamic star_p calculation based on actual data length
n_rows <- nrow(winresult)
star_p <- as.numeric(Slope)  # Start with actual slope values

# Apply adjustments based on position in the sequence
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

# Create the barplot with proper axis labels
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

# Create the 'figures' folder if it doesn't exist
if (!dir.exists("figures")) {
  dir.create("figures")
}

# Save the plot
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















# Rs~year in SOC intervals with blue color scheme and proper labels

# Define blue color palette
#blue_colors <- c('#E6F3FF', '#80BFFF', '#4DA6FF', '#0073E6')  # Light to dark blue

blue_colors <- c('#1B4F72', '#2874A6', '#5DADE2', '#85C1E9')

# Get the actual year range from your data
year_min <- floor(min(work.data$Study_midyear, na.rm = TRUE))
year_max <- ceiling(max(work.data$Study_midyear, na.rm = TRUE))
cat("Actual year range in your data:", year_min, "to", year_max, "\n")

# Calculate appropriate breaks for your time range
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


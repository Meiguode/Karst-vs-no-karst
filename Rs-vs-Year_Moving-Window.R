source("functions_plot.R")
source("statistical_prep.R")

# Moving subset analysis for year windows
mycolors = colorRampPalette(c("#0000a2", "#bc272d", "#e9c716", "#50ad9f"))

rcdata = dplyr::select(work.data, Latitude, Rs_annual, Study_midyear)

t = q = 10
p = 1
i = 1

ord = c()
Lat.Area = c()
Slope = c()
R.square = c()
P.value = c()
Num = c()
Mid.Lat = c()
Year.range = c()

Mean_Rs = c()
Median_Rs = c()
SD_Rs = c()
SE_Rs = c()
Intercept = c()
Adj_R.square = c()

while (t <= max(rcdata$Study_midyear - 1990 + p + 0.1)) {
  i = i + 1
  win = filter(rcdata, Study_midyear - 1990 >= (t - q) & (Study_midyear - 1990) < t)
  
  if (nrow(win) > 0) {
    eq = lm_eqn(df = win)
    Year.range = c(Year.range, paste0((1990 + t - q), '-', (1990 + t - 1)))
    Slope = c(Slope, eq[1])
    R.square = c(R.square, eq[2])
    P.value = c(P.value, eq[3])
    Num = c(Num, nrow(win))
    
    Mean_Rs = c(Mean_Rs, mean(win$Rs_annual, na.rm = TRUE))
    Median_Rs = c(Median_Rs, median(win$Rs_annual, na.rm = TRUE))
    SD_Rs = c(SD_Rs, sd(win$Rs_annual, na.rm = TRUE))
    SE_Rs = c(SE_Rs, sd(win$Rs_annual, na.rm = TRUE) / sqrt(nrow(win)))
    
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

winresult$sig = ''
winresult[winresult$P.value < 0.1, which(colnames(winresult) == 'sig')] = '·'
winresult[winresult$P.value < 0.05, which(colnames(winresult) == 'sig')] = '*'
winresult[winresult$P.value < 0.01, which(colnames(winresult) == 'sig')] = '**'
winresult[winresult$P.value < 0.001, which(colnames(winresult) == 'sig')] = '***'

winresult$Significance_Level = case_when(
  winresult$P.value < 0.001 ~ "p < 0.001",
  winresult$P.value < 0.01 ~ "p < 0.01", 
  winresult$P.value < 0.05 ~ "p < 0.05",
  winresult$P.value < 0.1 ~ "p < 0.1",
  TRUE ~ "Not significant"
)

significant_periods <- winresult %>% 
  filter(P.value < 0.05)

cat("=== SIGNIFICANT PERIODS (p < 0.05) ===\n")
if (nrow(significant_periods) > 0) {
  print(significant_periods)
} else {
  cat("No periods with p < 0.05 significance\n")
}

cat("\n=== ALL PERIODS - COMPREHENSIVE STATISTICS ===\n")
print(winresult)

write.csv(winresult, "moving_window_comprehensive_stats.csv", row.names = FALSE)
cat("\nComprehensive statistics saved to: moving_window_comprehensive_stats.csv\n")

write.csv(significant_periods, "significant_periods_stats.csv", row.names = FALSE)
cat("Significant periods statistics saved to: significant_periods_stats.csv\n")

if (nrow(significant_periods) > 0) {
  significant_summary <- significant_periods %>%
    dplyr::select(Year.range, Slope, P.value, R.square, Adj_R.square, Num, Mean_Rs, SD_Rs, SE_Rs, Significance_Level)
  
  cat("\n=== SUMMARY OF SIGNIFICANT PERIODS ===\n")
  print(significant_summary)
  
  write.csv(significant_summary, "significant_periods_summary.csv", row.names = FALSE)
  cat("Significant periods summary saved to: significant_periods_summary.csv\n")
}

TS <- theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, size = 6, hjust = 1, vjust = 1),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "grey3", fill = NA, size = 1),
            axis.text.y = element_text(size = 6))

star_p = c(as.numeric(Slope)[1:7] + 1, 
           as.numeric(Slope)[8:12] - 3, 
           as.numeric(Slope)[13:18] + 1, 
           as.numeric(Slope)[19:22] + 1)

star <- geom_text(aes(x = Year.range, y = star_p, label = sig), size = 5)
winresult$Num.italic = paste0('italic(', winresult$Num, ')')
count <- geom_text(aes(x = Year.range, y = 54, label = Num.italic), size = 2, parse = TRUE)

SFM = scale_fill_manual(values = mycolors(nrow(winresult)))
SYC = scale_y_continuous(breaks = c(-15, 0, 15, 30, 45), limits = c(-15, 55))

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
  #ggtitle("(d)") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 8, margin = margin(t = 5)),
    axis.title.y = element_text(size = 8, margin = margin(r = 5))
  )

if (!dir.exists("figures")) {
  dir.create("figures")
}

ggsave("figures/fig1_d.pdf", plot = barplot, width = 6.5, height = 2, units = 'in', dpi = 900)
ggsave("figures/fig1_d.png", plot = barplot, width = 6.5, height = 2, units = 'in', dpi = 900)
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


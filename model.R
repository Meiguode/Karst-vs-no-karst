source("0_functions_plot.R")
source("4_statistical_prep.R")
library(dplyr)

workingpath <- getwd()
dir.create(paste0(workingpath, "/ModelTable"), showWarnings = FALSE)
res.data <- read.csv("03_processed_data_complete_final_non_karst.csv", check.names = FALSE)
res.data <- filter(
  res.data,
  Study_midyear >= 1987,
  Meas_method %in% c("IRGA", "Gas Chromatography"),
  Ecosystem_type != "Agriculture",
  Outlier == 0
)

dv <- c("Rs_annual")
idv <- c(
  "Study_midyear", "Latitude", "Biome",
  "MAT", "MAP", "delta.Mean.Temp.", "delta.Mean.Precip.", "SOC_stock",
  "Meas_method", "Ecosystem_type", "Partition_method", "Stage", "Elevation"
)

work.data <- na.omit(res.data[, c(idv, dv)])
work.data <- droplevels(work.data)

work.data$Year <- work.data$Study_midyear
work.data$Year2 <- work.data$Study_midyear^2

lmformula <- paste(
  "Rs_annual ~ Year + I(Year^2) +",
  "Biome + Meas_method + Latitude + Elevation + Stage + Ecosystem_type + SOC_stock +",
  "MAT + MAP + delta.Mean.Temp. + delta.Mean.Precip. +",
  "Year*Biome + Year*SOC_stock + Year*Stage + Year*Ecosystem_type + Year*Elevation +",
  "Year*Meas_method + Year*Latitude +",
  "Biome*MAT + Biome*MAP + Biome*MAT*MAP + Biome*delta.Mean.Temp. +",
  "Biome*delta.Mean.Precip. + Biome*delta.Mean.Temp.*delta.Mean.Precip. +",
  "MAT*MAP*Biome + delta.Mean.Temp.*delta.Mean.Precip.*Biome"
)

multi.lm <- lm(lmformula, data = work.data, na.action = "na.omit")

full.model.1 <- model.cof(multi.lm)
full.model.anova <- anova(multi.lm)

full.model.1$effect <- rownames(full.model.1)
full.model.anova$effect <- rownames(full.model.anova)
full.model.1 <- merge(full.model.anova, full.model.1, by = "effect", all = TRUE)

write.csv(
  full.model.1,
  paste0(workingpath, "/ModelTable/Rs.fullmodel.year.nonkarst.csv"),
  row.names = FALSE
)

library(car)

results_path <- paste0(workingpath, "/ModelTable/")
dir.create(results_path, showWarnings = FALSE)

sink(paste0(results_path, "Model_Summary_nonKarst.txt"))
cat("===== LINEAR MODEL SUMMARY =====\n\n")
print(summary(multi.lm))
sink()

write.csv(full.model.anova,
          paste0(results_path, "ANOVA_Table_nonKarst.csv"),
          row.names = TRUE)

write.csv(coef(summary(multi.lm)),
          paste0(results_path, "Coefficients_Table_nonKarst.csv"),
          row.names = TRUE)

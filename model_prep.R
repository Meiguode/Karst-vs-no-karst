library(dplyr)

datasets <- list(
  "Non-karst" = "03_nonkarst_with_koppen.csv",
  "Karst"     = "03_karst_with_koppen.csv",
  "Combined"  = "03_combined_with_koppen.csv"
)

model_vars <- c(
  "Rs_annual", "Study_midyear", "YearsOfData",
  "Biome", "kg_code",
  "Meas_method", "Stage", "Ecosystem_type",
  "Latitude", "Elevation", "SOC_stock",
  "MAT", "MAP",
  "delta.Mean.Temp.", "delta.Mean.Precip."
)

prep_filter <- function(file) {

  df_raw <- read.csv(file, stringsAsFactors = FALSE)

  cat("\nFILE:", file, "\n")
  cat("Rows on read:", nrow(df_raw), "\n")

  df <- df_raw[
    df_raw$Study_midyear >= 1990 &
    df_raw$Meas_method %in% c("IRGA", "Gas Chromatography") &
    df_raw$Ecosystem_type != "Cropland" &
    df_raw$Outlier == 0,
  ]

  df$Study_midyear <- as.numeric(df$Study_midyear)
  cat("Study_midyear range:", range(df$Study_midyear, na.rm = TRUE), "\n")

  cat("Rows after filters:", nrow(df), "\n")

  na_counts <- colSums(is.na(df[, model_vars]))
  cat("NA counts (model vars only):\n")
  print(na_counts[na_counts > 0])

  df <- df[, model_vars]
  df <- df[complete.cases(df), ]

  cat("Rows after NA removal (model vars only):", nrow(df), "\n")

  return(df)
}

prepped <- lapply(datasets, prep_filter)

add_breakpoint_periods <- function(df, dataset_name) {

  if (dataset_name %in% c("Combined", "Non-karst")) {
    df$Period <- ifelse(
      df$Study_midyear < 2013, "Pre-breakpoint (1990–2012)",
      ifelse(df$Study_midyear <= 2016,
             "Breakpoint (2013–2016)",
             "Post-breakpoint (2017–2020)")
    )
  } else {
    df$Period <- ifelse(
      df$Study_midyear < 2009, "Pre-breakpoint (1990–2008)",
      ifelse(df$Study_midyear <= 2014,
             "Breakpoint (2009–2014)",
             "Post-breakpoint (2015–2019)")
    )
  }

  df
}

prepped <- mapply(add_breakpoint_periods, prepped, names(prepped),
                  SIMPLIFY = FALSE)

cat("\nFINAL ROW COUNTS:\n")
print(lapply(prepped, nrow))

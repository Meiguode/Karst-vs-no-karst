library(dplyr)

datasets <- list(
  "Non-karst" = "03_processed_data_complete_final_nonkarst.csv",
  "Karst"     = "03_processed_data_complete_final_karst.csv",
  "Combined"  = "03_processed_data_complete_final.csv"
)

prep_filter <- function(file) {

  df <- read.csv(file, stringsAsFactors = FALSE)

  df %>%
    filter(
      Study_midyear >= 1990,
      Meas_method %in% c("IRGA", "Gas Chromatography"),
      Ecosystem_type != "Cropland",
      Outlier == 0,
      !is.na(Rs_annual),
      !is.na(Study_midyear),

      Rs_annual > 0,
      Rs_annual <= 5000
    )
}
prepped <- lapply(datasets, prep_filter)

add_breakpoint_periods <- function(df, dataset_name) {

  if (dataset_name %in% c("Combined", "Non-karst")) {

    df %>%
      mutate(
        Period = case_when(
          Study_midyear < 2013 ~ "Pre-breakpoint",
          Study_midyear <= 2016 ~ "Breakpoint",
          TRUE ~ "Post-breakpoint"
        )
      )

  } else if (dataset_name == "Karst") {

    df %>%
      mutate(
        Period = case_when(
          Study_midyear < 2009 ~ "Pre-breakpoint",
          Study_midyear <= 2014 ~ "Breakpoint",
          TRUE ~ "Post-breakpoint"
        )
      )
  }
}

prepped <- mapply(
  add_breakpoint_periods,
  prepped,
  names(prepped),
  SIMPLIFY = FALSE
)

nonkarst <- prepped[["Non-karst"]]
karst    <- prepped[["Karst"]]
combined <- prepped[["Combined"]]

cat("\nRows after prep:\n")
cat("Non-karst:", nrow(nonkarst), "\n")
cat("Karst:", nrow(karst), "\n")
cat("Combined:", nrow(combined), "\n")

cat("\nPeriod distribution:\n")
print(lapply(list(nonkarst = nonkarst,
                  karst = karst,
                  combined = combined),
             function(df) table(df$Period)))

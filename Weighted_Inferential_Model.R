source("functions_plot.R")
source("model_prep.R")

library(dplyr)
library(car)
library(broom)
library(tibble)
library(caret)

set.seed(123)

workingpath <- getwd()
results_path <- file.path(workingpath, "ModelTable")

dir.create(results_path, showWarnings = FALSE, recursive = TRUE)

lm_base <- Rs_annual ~
  Study_midyear + I(Study_midyear^2) +
  Biome +
  Meas_method + Stage + Ecosystem_type +
  Latitude + Elevation + SOC_stock +
  MAT + MAP +
  delta.Mean.Temp. + delta.Mean.Precip. +
  Study_midyear * Biome +
  Study_midyear * SOC_stock +
  Study_midyear * Stage +
  Study_midyear * Ecosystem_type

lm_koppen <- update(lm_base, . ~ . + kg_code)

run_kfold_cv <- function(dat_model, formula, label, k = 10) {
  set.seed(123)
  fold_ids <- sample(rep(1:k, length.out = nrow(dat_model)))
  
  cv_results <- list()
  
  for (i in 1:k) {
    train_idx <- which(fold_ids != i)
    test_idx <- which(fold_ids == i)
    
    train_data <- dat_model[train_idx, ]
    test_data <- dat_model[test_idx, ]
    
    train_actual <- train_data$Rs_annual
    test_actual <- test_data$Rs_annual
    
    factor_vars <- all.vars(formula)
    factor_vars <- factor_vars[factor_vars %in% colnames(dat_model)]
    
    rows_to_keep <- rep(TRUE, nrow(test_data))
    for (var in factor_vars) {
      if (is.factor(train_data[[var]])) {
        train_levels <- levels(train_data[[var]])
        test_levels <- as.character(test_data[[var]])
        new_level_rows <- !test_levels %in% train_levels
        if (any(new_level_rows)) {
          rows_to_keep[new_level_rows] <- FALSE
        }
      }
    }
    
    test_data <- test_data[rows_to_keep, ]
    test_actual <- test_actual[rows_to_keep]
    
    if (nrow(test_data) == 0) {
      cv_results[[i]] <- data.frame(
        fold = i,
        train_n = nrow(train_data),
        test_n = 0,
        train_r2 = NA,
        test_r2 = NA,
        r2_reduction = NA,
        reduction_pct = NA
      )
      next
    }
    
    for (var in factor_vars) {
      if (is.factor(train_data[[var]])) {
        train_levels <- levels(train_data[[var]])
        test_data[[var]] <- factor(test_data[[var]], levels = train_levels)
      }
    }
    
    if ("weights" %in% colnames(train_data)) {
      model_cv <- lm(formula, data = train_data, weights = weights)
    } else {
      model_cv <- lm(formula, data = train_data)
    }
    
    train_pred <- predict(model_cv, newdata = train_data)
    if (sd(train_pred) > 0 && sd(train_actual) > 0) {
      train_r2 <- cor(train_pred, train_actual)^2
    } else {
      train_r2 <- NA
    }
    
    test_pred <- tryCatch({
      predict(model_cv, newdata = test_data)
    }, error = function(e) {
      return(rep(NA, nrow(test_data)))
    })
    
    if (any(is.na(test_pred)) || length(test_pred) != length(test_actual)) {
      test_r2 <- NA
    } else if (sd(test_pred) > 0 && sd(test_actual) > 0) {
      test_r2 <- cor(test_pred, test_actual)^2
    } else {
      test_r2 <- NA
    }
    
    if (!is.na(train_r2) && !is.na(test_r2)) {
      r2_reduction <- train_r2 - test_r2
      reduction_pct <- (r2_reduction / train_r2) * 100
    } else {
      r2_reduction <- NA
      reduction_pct <- NA
    }
    
    cv_results[[i]] <- data.frame(
      fold = i,
      train_n = nrow(train_data),
      test_n = nrow(test_data),
      train_r2 = train_r2,
      test_r2 = test_r2,
      r2_reduction = r2_reduction,
      reduction_pct = reduction_pct
    )
  }
  
  cv_df <- do.call(rbind, cv_results)
  
  cv_df_valid <- cv_df[!is.na(cv_df$train_r2) & !is.na(cv_df$test_r2), ]
  
  if (nrow(cv_df_valid) == 0) {
    cv_summary <- data.frame(
      mean_train_r2 = NA,
      mean_test_r2 = NA,
      mean_reduction_pct = NA,
      sd_train_r2 = NA,
      sd_test_r2 = NA
    )
  } else {
    cv_summary <- cv_df_valid %>%
      summarise(
        mean_train_r2 = mean(train_r2),
        mean_test_r2 = mean(test_r2),
        mean_reduction_pct = mean(reduction_pct),
        sd_train_r2 = sd(train_r2),
        sd_test_r2 = sd(test_r2)
      )
  }
  
  write.csv(
    cv_df,
    file.path(results_path, paste0("CV_fold_results_", label, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    cv_summary,
    file.path(results_path, paste0("CV_summary_", label, ".csv")),
    row.names = FALSE
  )
  
  return(cv_summary)
}

run_rs_model <- function(dat, label, k = 10) {

  cat("RUNNING DATASET:", label, "\n")

  cat("Rows received from statistical_prep:", nrow(dat), "\n")

  
  # CREATE WEIGHTS
  if ("YearsOfData" %in% colnames(dat)) {
    dat$weights <- dat$YearsOfData
    cat("✓ Using YearsOfData for weights (range:", 
        range(dat$YearsOfData, na.rm = TRUE), ")\n")
  } else {
    cat("⚠️ No 'YearsOfData' column found. Using equal weights.\n")
    dat$weights <- 1
  }

  fac_vars <- c("Biome","Ecosystem_type","Meas_method","Stage","kg_code")
  dat[fac_vars] <- lapply(dat[fac_vars], factor)


  kg_freq <- sort(table(dat$kg_code), decreasing = TRUE)
  print(kg_freq)

  min_kg_n <- max(10, floor(0.05 * nrow(dat)))
  rare_kg <- names(kg_freq[kg_freq < min_kg_n])

  dat_k <- dat
  dat_k$kg_code <- as.character(dat_k$kg_code)
  dat_k$kg_code[dat_k$kg_code %in% rare_kg] <- "Rare"
  dat_k$kg_code <- factor(dat_k$kg_code)

  use_koppen <- nlevels(dat_k$kg_code) >= 2

  if (!use_koppen) {
    cat("⚠️ Köppen removed (insufficient levels)\n")
    dat_model <- dat[, setdiff(names(dat), "kg_code")]
  } else {
    cat("✅ Köppen retained with", nlevels(dat_k$kg_code), "levels\n")
    dat_model <- dat_k
  }
  
  m_base <- lm(lm_base, dat_model, weights = weights)
  s_base <- summary(m_base)

  if (use_koppen) {
    m_k <- lm(lm_koppen, dat_model, weights = weights)
    s_k <- summary(m_k)
    print(anova(m_base, m_k))
  }
  
  cat("\n--- Running k-fold cross-validation (k =", k, ") ---\n")
    
    cv_base <- run_kfold_cv(dat_model, lm_base, paste0(label, "_base"), k = k)
    cat("Base model - Training R²:", round(cv_base$mean_train_r2, 3), 
        "Validation R²:", round(cv_base$mean_test_r2, 3),
        "Reduction:", round(cv_base$mean_reduction_pct, 1), "%\n")
    
    if (use_koppen) {
      cv_k <- run_kfold_cv(dat_model, lm_koppen, paste0(label, "_koppen"), k = k)
      cat("Köppen model - Training R²:", round(cv_k$mean_train_r2, 3),
          "Validation R²:", round(cv_k$mean_test_r2, 3),
          "Reduction:", round(cv_k$mean_reduction_pct, 1), "%\n")
    }

  # ANOVA

  if (use_koppen) {
    anova_mod <- anova(m_k)
    anova_label <- "with_Koppen"
  } else {
    anova_mod <- anova(m_base)
    anova_label <- "base"
  }

  anova_tbl <- broom::tidy(anova_mod)

  write.csv(
    anova_tbl,
    file.path(
      results_path,
      paste0("ANOVA_", label, "_", anova_label, ".csv")
    ),
    row.names = FALSE
  )

  cat("\nANOVA table saved for", label, "(", anova_label, "model )\n")
  
  tibble(
      Dataset = label,
      N = nobs(m_base),
      Koppen_Used = use_koppen,
      R2_base = s_base$r.squared,
      Adj_R2_base = s_base$adj.r.squared,
      R2_koppen = if (use_koppen) s_k$r.squared else NA,
      Adj_R2_koppen = if (use_koppen) s_k$adj.r.squared else NA,
      # NEW: Cross-validation metrics
      CV_Train_R2 = cv_base$mean_train_r2,
      CV_Test_R2 = cv_base$mean_test_r2,
      CV_Reduction_Pct = cv_base$mean_reduction_pct,
      CV_Train_R2_k = if (use_koppen) cv_k$mean_train_r2 else NA,
      CV_Test_R2_k = if (use_koppen) cv_k$mean_test_r2 else NA,
      CV_Reduction_Pct_k = if (use_koppen) cv_k$mean_reduction_pct else NA
    )
  }

model_table <- bind_rows(
  run_rs_model(prepped[["Non-karst"]], "Non-karst"),
  run_rs_model(prepped[["Karst"]],     "Karst"),
  run_rs_model(prepped[["Combined"]],  "Combined")
)

write.csv(
  model_table,
  file.path(results_path, "Table_Model_Performance_Karst_vs_NonKarst.csv"),
  row.names = FALSE
)

print(model_table)

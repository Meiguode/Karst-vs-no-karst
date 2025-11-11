source("0_functions_plot.R")

complete_data <- read.csv("srdb-data.csv", check.names = FALSE)

cat("Total records in srdb-data.csv:", nrow(complete_data), "\n")
cat("Study_midyear range:", range(complete_data$Study_midyear, na.rm = TRUE), "\n")
cat("Unique Meas_method values:", unique(complete_data$Meas_method), "\n")
cat("Unique Ecosystem_type values:", unique(complete_data$Ecosystem_type), "\n")
cat("Unique Quality_flag values:", unique(complete_data$Quality_flag), "\n")

cat("\n--- Filtering steps ---\n")
cat("Records with Study_midyear >= 1990:", sum(complete_data$Study_midyear >= 1990, na.rm = TRUE), "\n")
cat("Records with Study_midyear <= 2022:", sum(complete_data$Study_midyear <= 2022, na.rm = TRUE), "\n")
cat("Records with IRGA or Gas Chromatography:", sum(complete_data$Meas_method %in% c('IRGA','Gas Chromatography'), na.rm = TRUE), "\n")
cat("Records not Agriculture:", sum(complete_data$Ecosystem_type != 'Agriculture', na.rm = TRUE), "\n")
cat("Records with Quality_flag == 0:", sum(complete_data$Quality_flag == 0, na.rm = TRUE), "\n")

res.data <- complete_data %>%
  filter(Study_midyear >= 1990 & Study_midyear <= 2022) %>%
  dplyr::select(Record_number, Study_midyear, YearsOfData, Latitude, Longitude, 
         Site_name, Study_number, Elevation, Biome, Ecosystem_state, Ecosystem_type,
         Leaf_habit, MAT, MAP, Meas_method, Partition_method, Stage,
         Rs_annual, Rh_annual, Ra_annual,
         Manipulation, Soil_BD, Soil_CN, Soil_sand, Soil_silt, Soil_clay,
         Quality_flag)

cat("\nRecords after time filter:", nrow(res.data), "\n")

res.data <- res.data %>%
  filter(Meas_method %in% c('IRGA','Gas Chromatography') | is.na(Meas_method))

cat("Records after measurement method filter:", nrow(res.data), "\n")

res.data <- res.data %>%
  filter(Ecosystem_type != 'Agriculture' | is.na(Ecosystem_type))

cat("Records after ecosystem filter:", nrow(res.data), "\n")

res.data <- res.data %>%
  filter(Quality_flag == 0 | is.na(Quality_flag))

cat("Records after quality filter:", nrow(res.data), "\n")

res.data <- res.data %>%
  mutate(
    Ecosystem_type2 = case_when(
      Ecosystem_type %in% c('Evergreen Forest', 'Deciduous Forest', 'Mixed Forest') ~ 'Forest',
      Ecosystem_type %in% c('Grassland', 'Savanna') ~ 'Grassland',
      Ecosystem_type %in% c('Shrubland', 'Tundra') ~ 'Shrubland',
      TRUE ~ as.character(Ecosystem_type)
    ),
    Outlier_rh = 0,
    Outlier_ra = 0
  )

res.data <- res.data %>%
  filter(!is.na(Rs_annual) | !is.na(Rh_annual) | !is.na(Ra_annual))

cat("Final records after removing NA respiration data:", nrow(res.data), "\n")

if (nrow(res.data) == 0) {
  cat("\nNo data after all filters. Let's check what we have with just time range:\n")
  minimal_data <- complete_data %>%
    filter(Study_midyear >= 1990 & Study_midyear <= 2022) %>%
    filter(!is.na(Rs_annual) | !is.na(Rh_annual) | !is.na(Ra_annual))
  
  cat("Records with just time range and respiration data:", nrow(minimal_data), "\n")
  
  res.data <- minimal_data %>%
    dplyr::select(Record_number, Study_midyear, YearsOfData, Latitude, Longitude, 
           Site_name, Study_number, Elevation, Biome, Ecosystem_state, Ecosystem_type,
           Leaf_habit, MAT, MAP, Meas_method, Partition_method, Stage,
           Rs_annual, Rh_annual, Ra_annual,
           Manipulation, Soil_BD, Soil_CN, Soil_sand, Soil_silt, Soil_clay) %>%
    mutate(
      Ecosystem_type2 = case_when(
        Ecosystem_type %in% c('Evergreen Forest', 'Deciduous Forest', 'Mixed Forest') ~ 'Forest',
        Ecosystem_type %in% c('Grassland', 'Savanna') ~ 'Grassland',
        Ecosystem_type %in% c('Shrubland', 'Tundra') ~ 'Shrubland',
        TRUE ~ as.character(Ecosystem_type)
      ),
      Outlier_rh = 0,
      Outlier_ra = 0
    )
}

write.csv(res.data, "complete_data_rhra.csv", row.names = FALSE)

cat("\nPreprocessing complete!\n")
cat("Final number of records:", nrow(res.data), "\n")
cat("Number of records with Rh_annual:", sum(!is.na(res.data$Rh_annual)), "\n")
cat("Number of records with Ra_annual:", sum(!is.na(res.data$Ra_annual)), "\n")
cat("Year range:", range(res.data$Study_midyear, na.rm = TRUE), "\n")
cat("File saved as: complete_data_rhra.csv\n")

#-
source("0_functions_plot.R")
source("4_statistical_prep.R")
library(dplyr)
library(ggplot2)
library(ggthemes)

workingpath <- getwd()
dir.create(paste0(workingpath, "/ModelTable"), showWarnings = FALSE)

cat("=== STEP 1: LOADING AND MERGING SOC_STOCK DATA ===\n")

res.data <- read.csv("record-number_rhra_with_climate_soc.csv", check.names = FALSE)
cat("Main data dimensions:", dim(res.data), "\n")

sample.data <- read.csv("02_sampledataset.csv", check.names = FALSE)
cat("Sample data dimensions:", dim(sample.data), "\n")

soc_data <- sample.data %>%
  dplyr::select(Record_number, SOC_stock) %>%
  filter(!is.na(SOC_stock))

cat("SOC data to merge:", nrow(soc_data), "records with SOC_stock values\n")

res.data <- res.data %>%
  left_join(soc_data, by = "Record_number")

cat("After merging SOC_stock:\n")
cat("Main data dimensions:", dim(res.data), "\n")
cat("Non-NA SOC_stock values in main data:", sum(!is.na(res.data$SOC_stock)), "\n")

cat("\n=== STEP 2: DATA PREPARATION ===\n")

desired_vars <- c(
  "Study_midyear", "Latitude", "Biome",
  "MAT", "MAP", "delta.Mean.Temp.", "delta.Mean.Precip.", "SOC_stock",
  "Meas_method", "Ecosystem_type", "Partition_method", "Stage", "Elevation"
)

available_vars <- desired_vars[desired_vars %in% names(res.data)]
cat("✅ Available variables:", paste(available_vars, collapse = ", "), "\n")

idv <- available_vars

res.data.rh <- res.data %>%
  dplyr::filter(
    Study_midyear >= 1990,
    Study_midyear <= 2022,
    Meas_method %in% c("IRGA", "Gas Chromatography"),
    Ecosystem_type != "Agriculture",
    Outlier_rh == 0
  )

cat("Rh data after filtering:", nrow(res.data.rh), "rows\n")

res.data.ra <- res.data %>%
  dplyr::filter(
    Study_midyear >= 1990,
    Study_midyear <= 2022,
    Meas_method %in% c("IRGA", "Gas Chromatography"),
    Ecosystem_type != "Agriculture",
    Outlier_ra == 0
  )

cat("Ra data after filtering:", nrow(res.data.ra), "rows\n")

create_complete_dataset <- function(data, response_var, dataset_name) {
  cat("Creating complete dataset for", dataset_name, "...\n")
  
  selected_vars <- c(available_vars, response_var)
  
  subset_data <- data %>% dplyr::select(all_of(selected_vars))
  
  cat("  Dataset dimensions before removing NAs:", dim(subset_data), "\n")
  
  complete_data <- na.omit(subset_data)
  
  cat("  Dataset dimensions after removing NAs:", dim(complete_data), "\n")
  
  factor_vars <- c("Biome", "Meas_method", "Ecosystem_type", "Partition_method", "Stage")
  vars_to_remove <- c()
  
  for (var in factor_vars) {
    if (var %in% names(complete_data)) {
      unique_levels <- length(unique(complete_data[[var]]))
      cat("  ", var, "has", unique_levels, "levels\n")
      if (unique_levels < 2) {
        cat("  ⚠️ Removing", var, "- only one level\n")
        vars_to_remove <- c(vars_to_remove, var)
        complete_data[[var]] <- NULL
      } else {
        # Convert to factor if not already
        complete_data[[var]] <- as.factor(complete_data[[var]])
      }
    }
  }
  
  current_available_vars <- available_vars[!available_vars %in% vars_to_remove]
  
  complete_data$Year <- complete_data$Study_midyear
  complete_data$Year2 <- complete_data$Study_midyear^2
  
  return(list(data = complete_data, available_vars = current_available_vars))
}

cat("\n--- Rh Data ---\n")
rh_result <- create_complete_dataset(res.data.rh, "Rh_annual", "Rh")
work.data.rh <- rh_result$data
rh_available_vars <- rh_result$available_vars

cat("\n--- Ra Data ---\n")
ra_result <- create_complete_dataset(res.data.ra, "Ra_annual", "Ra")
work.data.ra <- ra_result$data
ra_available_vars <- ra_result$available_vars

cat("Rh complete cases:", nrow(work.data.rh), "rows\n")
cat("Ra complete cases:", nrow(work.data.ra), "rows\n")

cat("\n=== STEP 3: BUILDING MODEL FORMULA ===\n")

build_formula <- function(available_vars, response_var, data) {
  # Base terms that should always be included
  base_terms <- c("Year", "I(Year^2)")
  
  main_effects <- available_vars[!available_vars %in% c("Study_midyear")]
  base_terms <- c(base_terms, main_effects)
  
  final_terms <- c()
  for (term in base_terms) {
    # For simple terms
    if (!grepl("[:*]", term)) {
      if (term %in% names(data)) {
        # Check if it's a factor with sufficient levels
        if (is.factor(data[[term]])) {
          if (length(levels(data[[term]])) >= 2) {
            final_terms <- c(final_terms, term)
          }
        } else {
          final_terms <- c(final_terms, term)
        }
      }
    } else {
      # For interaction terms, check if all components exist
      components <- unlist(strsplit(term, "[:*]"))
      if (all(components %in% names(data))) {
        # Check if any factor components have sufficient levels
        valid_term <- TRUE
        for (comp in components) {
          if (is.factor(data[[comp]]) && length(levels(data[[comp]])) < 2) {
            valid_term <- FALSE
            break
          }
        }
        if (valid_term) {
          final_terms <- c(final_terms, term)
        }
      }
    }
  }
  
  # Create formula
  if (length(final_terms) > 0) {
    formula_str <- paste(response_var, "~", paste(final_terms, collapse = " + "))
  } else {
    formula_str <- paste(response_var, "~ Year")  # Fallback to simple model
  }
  
  return(formula_str)
}

rh_formula <- build_formula(rh_available_vars, "Rh_annual", work.data.rh)
ra_formula <- build_formula(ra_available_vars, "Ra_annual", work.data.ra)

cat("Rh model formula:\n", rh_formula, "\n\n")
cat("Ra model formula:\n", ra_formula, "\n\n")

cat("=== STEP 4: FITTING MODELS ===\n")

# Function to fit model and handle errors
fit_model <- function(formula_str, data, model_name) {
  cat("Fitting", model_name, "...\n")
  
  tryCatch({
    model <- lm(formula_str, data = data, na.action = "na.omit")
    cat("✅", model_name, "fitted successfully!\n")
    cat("   Sample size:", nrow(data), "observations\n")
    cat("   Number of predictors:", length(coef(model)) - 1, "\n")
    cat("   R-squared:", summary(model)$r.squared, "\n")
    
    # Check ANOVA effects
    anova_result <- anova(model)
    cat("   Effects in ANOVA:", paste(rownames(anova_result), collapse = ", "), "\n\n")
    
    return(model)
  }, error = function(e) {
    cat("❌", model_name, "fitting failed:", e$message, "\n")
    
    # Try simpler model
    cat("   Trying simpler model...\n")
    simple_formula <- paste(all.vars(as.formula(formula_str))[1], "~ Year")
    tryCatch({
      model <- lm(simple_formula, data = data, na.action = "na.omit")
      cat("   ✅ Simple model fitted successfully\n")
      return(model)
    }, error = function(e2) {
      cat("   ❌ Even simple model failed:", e2$message, "\n\n")
      return(NULL)
    })
  })
}

rh_model <- fit_model(rh_formula, work.data.rh, "Rh model")

ra_model <- fit_model(ra_formula, work.data.ra, "Ra model")

cat("=== STEP 5: GENERATING OUTPUTS ===\n")

generate_model_outputs <- function(model, model_name, response_var, data) {
  if (is.null(model)) {
    cat("❌ No model for", model_name, "- skipping outputs\n")
    return(NULL)
  }
  
  cat("Generating outputs for", model_name, "...\n")
  
  full.model.1 <- model.cof(model)
  full.model.anova <- anova(model)
  
  full.model.1$effect <- rownames(full.model.1)
  full.model.anova$effect <- rownames(full.model.anova)
  full.model.combined <- merge(full.model.anova, full.model.1, by = "effect", all = TRUE)
  
  write.csv(
    full.model.combined,
    paste0(workingpath, "/ModelTable/", model_name, ".fullmodel.year.csv"),
    row.names = FALSE
  )
  
  results_path <- paste0(workingpath, "/ModelTable/")
  
  sink(paste0(results_path, model_name, "_Model_Summary.txt"))
  cat("===== LINEAR MODEL SUMMARY =====\n\n")
  cat("Response variable:", response_var, "\n")
  cat("Sample size:", nrow(data), "\n")
  cat("Model formula:", paste(deparse(formula(model)), collapse = "\n"), "\n\n")
  print(summary(model))
  sink()
  
  write.csv(full.model.anova,
            paste0(results_path, model_name, "_ANOVA_Table.csv"),
            row.names = TRUE)
  
  write.csv(coef(summary(model)),
            paste0(results_path, model_name, "_Coefficients_Table.csv"),
            row.names = TRUE)
  
  tryCatch({
    pdf(paste0(results_path, model_name, "_Diagnostic_Plots.pdf"), width = 8, height = 8)
    par(mfrow = c(2, 2))
    plot(model, main = paste("Diagnostic Plots for", response_var))
    dev.off()
  }, error = function(e) {
    cat("   ⚠️ Diagnostic plots failed:", e$message, "\n")
  })
  
  cat("✅ Outputs generated for", model_name, "\n")
  return(list(combined = full.model.combined, anova = full.model.anova))
}

rh_results <- generate_model_outputs(rh_model, "Rh", "Rh_annual", work.data.rh)
ra_results <- generate_model_outputs(ra_model, "Ra", "Ra_annual", work.data.ra)

cat("\n=== STEP 6: CREATING SUPPLEMENTARY TABLE 3 ===\n")

if (!is.null(rh_results) && !is.null(ra_results)) {
  # Extract ANOVA results
  rh_anova <- rh_results$anova
  ra_anova <- ra_results$anova
  
  clean_effect_names <- function(effects) {
    effects <- gsub("Study_midyear", "Year", effects)
    effects <- gsub("I(Study_midyear^2)", "Year2", effects, fixed = TRUE)
    effects <- gsub("Elevation", "Altitude", effects)
    effects <- gsub("Ecosystem_type", "Ecosystem type", effects)
    effects <- gsub("Meas_method", "Measurement method", effects)
    effects <- gsub("Partition_method", "Partitioning method", effects)
    effects <- gsub("delta.Mean.Temp.", "ΔMAT", effects)
    effects <- gsub("delta.Mean.Precip.", "ΔMAP", effects)
    effects <- gsub("SOC_stock", "SOC", effects)
    effects <- gsub(":", "×", effects)
    return(effects)
  }
  
  rh_effects <- clean_effect_names(rownames(rh_anova))
  ra_effects <- clean_effect_names(rownames(ra_anova))
  
  all_effects <- unique(c(rh_effects, ra_effects))
  
  comprehensive_table <- data.frame(
    Effect = all_effects,
    Degree_of_freedom = NA,
    Rh_F = NA,
    Rh_P = NA,
    Ra_F = NA,
    Ra_P = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in 1:length(rh_effects)) {
    effect <- rh_effects[i]
    idx <- which(comprehensive_table$Effect == effect)
    if (length(idx) > 0) {
      comprehensive_table$Degree_of_freedom[idx] <- rh_anova$Df[i]
      comprehensive_table$Rh_F[idx] <- round(rh_anova$`F value`[i], 3)
      p_value <- rh_anova$`Pr(>F)`[i]
      comprehensive_table$Rh_P[idx] <- ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
    }
  }
  
  for (i in 1:length(ra_effects)) {
    effect <- ra_effects[i]
    idx <- which(comprehensive_table$Effect == effect)
    if (length(idx) > 0) {
      comprehensive_table$Ra_F[idx] <- round(ra_anova$`F value`[i], 3)
      p_value <- ra_anova$`Pr(>F)`[i]
      comprehensive_table$Ra_P[idx] <- ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
    }
  }
  
  comprehensive_table <- comprehensive_table[!is.na(comprehensive_table$Rh_F) | !is.na(comprehensive_table$Ra_F), ]
  comprehensive_table <- comprehensive_table[order(-comprehensive_table$Rh_F, na.last = TRUE), ]
  
  write.csv(comprehensive_table,
            paste0(workingpath, "/ModelTable/Supplementary_Table_3_ANOVA_Summary.csv"),
            row.names = FALSE)
  
  sink(paste0(workingpath, "/ModelTable/Supplementary_Table_3_Formatted.txt"))
  cat("Supplementary Table 3. Summary of the effects in the linear models of Rh and Ra\n")
  cat("in 1990-2022.\n\n")
  cat("Effect\t\t\tDegree of\tRh\t\tRa\n")
  cat("\t\t\tfreedom\t\tF\tP\tF\tP\n")
  
  for (i in 1:nrow(comprehensive_table)) {
    effect <- comprehensive_table$Effect[i]
    if (nchar(effect) < 16) effect <- paste0(effect, "\t\t")
    else if (nchar(effect) < 24) effect <- paste0(effect, "\t")
    
    cat(effect, 
        comprehensive_table$Degree_of_freedom[i], "\t",
        comprehensive_table$Rh_F[i], "\t", comprehensive_table$Rh_P[i], "\t",
        comprehensive_table$Ra_F[i], "\t", comprehensive_table$Ra_P[i], "\n")
  }
  
  cat("\nA total of", nrow(work.data.rh), "observations of soil Rh and", 
      nrow(work.data.ra), "observations of soil Ra were analyzed\n")
  cat("with the linear model controlled by year-weighted variables.\n")
  sink()
  
  cat("✅ Supplementary Table 3 created successfully!\n")
} else {
  cat("❌ Cannot create Supplementary Table 3 - one or both models failed\n")
}

cat("\n=== STEP 7: COMPREHENSIVE DECADAL ANALYSIS ===\n")

decades <- list(
  "1990-2000" = c(1990, 2000),
  "2000-2010" = c(2000, 2010),
  "2010-2022" = c(2010, 2022),
  "1990-2022" = c(1990, 2022)
)

decadal_results <- data.frame()

for (decade_name in names(decades)) {
  start_year <- decades[[decade_name]][1]
  end_year <- decades[[decade_name]][2]
  
  cat("Analyzing decade:", decade_name, "\n")
  
  # Filter data for this decade
  Rh_decade <- res.data.rh %>% 
    dplyr::filter(Study_midyear >= start_year & Study_midyear <= end_year)
  
  Ra_decade <- res.data.ra %>% 
    dplyr::filter(Study_midyear >= start_year & Study_midyear <= end_year)
  
  cat("  Rh data points:", nrow(Rh_decade), "\n")
  cat("  Ra data points:", nrow(Ra_decade), "\n")
  
  if (nrow(Rh_decade) >= 3) {
    mean_rh <- mean(Rh_decade$Rh_annual, na.rm = TRUE)
    sd_rh <- sd(Rh_decade$Rh_annual, na.rm = TRUE)
    
    rh_trend <- lm(Rh_annual ~ Study_midyear, data = Rh_decade)
    rh_coef <- coef(summary(rh_trend))
    
    slope_rh <- ifelse("Study_midyear" %in% rownames(rh_coef), rh_coef["Study_midyear", "Estimate"], NA)
    pvalue_rh <- ifelse("Study_midyear" %in% rownames(rh_coef), rh_coef["Study_midyear", "Pr(>|t|)"], NA)
    
    if (!is.na(mean_rh) && mean_rh != 0 && !is.na(slope_rh)) {
      percent_change_rh <- (slope_rh * 10 / mean_rh) * 100
    } else {
      percent_change_rh <- NA
    }
    
    if (!is.na(slope_rh) && !is.na(pvalue_rh)) {
      conf_int <- confint(rh_trend)
      slope_lower <- ifelse("Study_midyear" %in% rownames(conf_int), conf_int["Study_midyear", 1], NA)
      slope_upper <- ifelse("Study_midyear" %in% rownames(conf_int), conf_int["Study_midyear", 2], NA)
    } else {
      slope_lower <- NA
      slope_upper <- NA
    }
    
    decadal_results <- rbind(decadal_results, data.frame(
      Period = decade_name,
      Response = "Rh_annual",
      Sample_Size = nrow(Rh_decade),
      Mean_Value = mean_rh,
      Std_Dev = sd_rh,
      Slope = slope_rh,
      Slope_Lower_CI = slope_lower,
      Slope_Upper_CI = slope_upper,
      Year_Pvalue = pvalue_rh,
      Percent_Change_Per_Decade = percent_change_rh,
      Significance = ifelse(pvalue_rh < 0.05, "SIGNIFICANT", 
                           ifelse(pvalue_rh < 0.1, "MARGINAL", "NOT SIGNIFICANT"))
    ))
  }
  
  if (nrow(Ra_decade) >= 3) {
    mean_ra <- mean(Ra_decade$Ra_annual, na.rm = TRUE)
    sd_ra <- sd(Ra_decade$Ra_annual, na.rm = TRUE)
    
    # Simple linear trend
    ra_trend <- lm(Ra_annual ~ Study_midyear, data = Ra_decade)
    ra_coef <- coef(summary(ra_trend))
    
    slope_ra <- ifelse("Study_midyear" %in% rownames(ra_coef), ra_coef["Study_midyear", "Estimate"], NA)
    pvalue_ra <- ifelse("Study_midyear" %in% rownames(ra_coef), ra_coef["Study_midyear", "Pr(>|t|)"], NA)
    
    # Calculate percent change per decade
    if (!is.na(mean_ra) && mean_ra != 0 && !is.na(slope_ra)) {
      percent_change_ra <- (slope_ra * 10 / mean_ra) * 100
    } else {
      percent_change_ra <- NA
    }
    
    # Calculate confidence intervals
    if (!is.na(slope_ra) && !is.na(pvalue_ra)) {
      conf_int <- confint(ra_trend)
      slope_lower <- ifelse("Study_midyear" %in% rownames(conf_int), conf_int["Study_midyear", 1], NA)
      slope_upper <- ifelse("Study_midyear" %in% rownames(conf_int), conf_int["Study_midyear", 2], NA)
    } else {
      slope_lower <- NA
      slope_upper <- NA
    }
    
    decadal_results <- rbind(decadal_results, data.frame(
      Period = decade_name,
      Response = "Ra_annual",
      Sample_Size = nrow(Ra_decade),
      Mean_Value = mean_ra,
      Std_Dev = sd_ra,
      Slope = slope_ra,
      Slope_Lower_CI = slope_lower,
      Slope_Upper_CI = slope_upper,
      Year_Pvalue = pvalue_ra,
      Percent_Change_Per_Decade = percent_change_ra,
      Significance = ifelse(pvalue_ra < 0.05, "SIGNIFICANT", 
                           ifelse(pvalue_ra < 0.1, "MARGINAL", "NOT SIGNIFICANT"))
    ))
  }
}

write.csv(decadal_results,
          paste0(workingpath, "/ModelTable/Comprehensive_Decadal_Analysis.csv"),
          row.names = FALSE)

summary_decadal <- decadal_results %>%
  dplyr::select(Period, Response, Sample_Size, Slope, Year_Pvalue, Significance) %>%
  mutate(
    Slope = round(Slope, 4),
    Year_Pvalue = round(Year_Pvalue, 4)
  )

write.csv(summary_decadal,
          paste0(workingpath, "/ModelTable/Decadal_Trends_Summary.csv"),
          row.names = FALSE)

cat("✅ Comprehensive decadal analysis saved!\n")

cat("\n=== STEP 8: CREATING DECADAL TREND PLOTS ===\n")

tryCatch({
  p_rh <- ggplot() +
    geom_point(data = res.data.rh, aes(x = Study_midyear, y = Rh_annual), 
               alpha = 0.3, color = "red") +
    geom_smooth(data = res.data.rh, aes(x = Study_midyear, y = Rh_annual), 
                method = "lm", color = "darkred", se = TRUE) +
    facet_wrap(~cut(Study_midyear, breaks = c(1990, 2000, 2010, 2022), 
                   labels = c("1990-2000", "2000-2010", "2010-2022"))) +
    labs(title = "Rh Annual Trends by Decade",
         x = "Year", y = "Rh Annual (g C m⁻²)") +
    theme_few()
  
  ggsave(paste0(workingpath, "/ModelTable/Rh_Decadal_Trends.png"), 
         p_rh, width = 12, height = 6, dpi = 300)
  
  p_ra <- ggplot() +
    geom_point(data = res.data.ra, aes(x = Study_midyear, y = Ra_annual), 
               alpha = 0.3, color = "blue") +
    geom_smooth(data = res.data.ra, aes(x = Study_midyear, y = Ra_annual), 
                method = "lm", color = "darkblue", se = TRUE) +
    facet_wrap(~cut(Study_midyear, breaks = c(1990, 2000, 2010, 2022), 
                   labels = c("1990-2000", "2000-2010", "2010-2022"))) +
    labs(title = "Ra Annual Trends by Decade",
         x = "Year", y = "Ra Annual (g C m⁻²)") +
    theme_few()
  
  ggsave(paste0(workingpath, "/ModelTable/Ra_Decadal_Trends.png"), 
         p_ra, width = 12, height = 6, dpi = 300)
  
  cat("✅ Decadal trend plots created successfully!\n")
}, error = function(e) {
  cat("⚠️ Decadal plots failed:", e$message, "\n")
})

cat("\n=== STEP 9: FINAL SUMMARY ===\n")

# List all generated files
output_files <- list.files(paste0(workingpath, "/ModelTable"))
cat("📊 ALL GENERATED FILES:\n")
for (file in output_files) {
  file_info <- file.info(paste0(workingpath, "/ModelTable/", file))
  cat(sprintf("• %-45s (%d bytes)\n", file, file_info$size))
}

cat("\n📈 DECADAL TRENDS SUMMARY:\n")
cat("Period,Response,Sample_Size,Slope,Year_Pvalue,Significance\n")
for (i in 1:nrow(summary_decadal)) {
  row <- summary_decadal[i, ]
  cat(paste0('"', row$Period, '","', row$Response, '",', 
             row$Sample_Size, ",", row$Slope, ",", 
             row$Year_Pvalue, ',"', row$Significance, '"\n'))
}

cat("\n🎯 ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("Key outputs include:\n")
cat("• Supplementary Table 3 with comprehensive ANOVA results\n")
cat("• Full model outputs for Rh and Ra\n")
cat("• Comprehensive decadal analysis with slopes and p-values\n")
cat("• Decadal trend visualizations\n")
cat("• All year combinations: 1990-2000, 2000-2010, 2010-2022, 1990-2022\n")

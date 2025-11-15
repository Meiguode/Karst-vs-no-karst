# Load required libraries
library(dplyr)
library(tidyr)
library(broom)
library(stringr)

# Read your data
data <- read.csv("complete_data.csv")

# Define the time periods: 1990-2000, 2000-2010, 2010-2022
data <- data %>%
  mutate(period = case_when(
    Study_midyear >= 1990 & Study_midyear <= 2000 ~ "1990-2000",
    Study_midyear >= 2000 & Study_midyear <= 2010 ~ "2000-2010", 
    Study_midyear >= 2010 & Study_midyear <= 2022 ~ "2010-2022",
    Study_midyear >= 1990 & Study_midyear <= 2022 ~ "1990-2022"
  ))

# Function to perform linear regression and extract slope, p-value, and count
get_regression_stats <- function(df) {
  if (nrow(df) < 2) {
    return(data.frame(Slope = NA, P = NA, count = 0))
  }
  tryCatch({
    model <- lm(Rs_annual ~ Study_midyear, data = df)
    tidy_model <- tidy(model)
    data.frame(
      Slope = round(tidy_model$estimate[2], 3),
      P = round(tidy_model$p.value[2], 3),
      count = nrow(df)
    )
  }, error = function(e) {
    return(data.frame(Slope = NA, P = NA, count = nrow(df)))
  })
}

# Create the main analysis for each ecosystem type and period
results <- data %>%
  filter(!is.na(Rs_annual), !is.na(Study_midyear), !is.na(period)) %>%
  group_by(Ecosystem_type2, period) %>%
  group_modify(~ get_regression_stats(.x)) %>%
  ungroup() %>%
  pivot_wider(
    names_from = period,
    values_from = c(Slope, P, count),
    names_glue = "{period}_{.value}"
  )

# Ensure all required columns exist and handle missing combinations
required_columns <- c(
  "1990-2000_Slope", "1990-2000_P", "1990-2000_count",
  "2000-2010_Slope", "2000-2010_P", "2000-2010_count", 
  "2010-2022_Slope", "2010-2022_P", "2010-2022_count",
  "1990-2022_Slope", "1990-2022_P", "1990-2022_count"
)

# Add missing columns with NA values
for (col in required_columns) {
  if (!col %in% colnames(results)) {
    results[[col]] <- NA
  }
}

# Reorder columns to match the desired format
results <- results %>%
  select(
    Ecosystem_type = Ecosystem_type2,
    `1990-2000_Slope`, `1990-2000_P`, `1990-2000_count`,
    `2000-2010_Slope`, `2000-2010_P`, `2000-2010_count`,
    `2010-2022_Slope`, `2010-2022_P`, `2010-2022_count`,
    `1990-2022_Slope`, `1990-2022_P`, `1990-2022_count`
  )

# Create the Forest row (combination of deciduous, evergreen, mixed forests)
forest_data <- data %>%
  filter(Ecosystem_type2 %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest"),
         !is.na(Rs_annual), !is.na(Study_midyear))

forest_results <- forest_data %>%
  group_by(period) %>%
  group_modify(~ get_regression_stats(.x)) %>%
  ungroup() %>%
  pivot_wider(
    names_from = period,
    values_from = c(Slope, P, count),
    names_glue = "{period}_{.value}"
  ) %>%
  mutate(Ecosystem_type = "Forest")

# Ensure all required columns exist for forest results
for (col in required_columns) {
  if (!col %in% colnames(forest_results)) {
    forest_results[[col]] <- NA
  }
}

forest_results <- forest_results %>%
  select(
    Ecosystem_type,
    `1990-2000_Slope`, `1990-2000_P`, `1990-2000_count`,
    `2000-2010_Slope`, `2000-2010_P`, `2000-2010_count`,
    `2010-2022_Slope`, `2010-2022_P`, `2010-2022_count`,
    `1990-2022_Slope`, `1990-2022_P`, `1990-2022_count`
  )

# Combine all results
final_table <- bind_rows(forest_results, results) %>%
  arrange(factor(Ecosystem_type, levels = c(
    "Forest", "Deciduous Forest", "Evergreen Forest", "Mixed Forest",
    "Grassland", "Shrubland", "Wetland", "Others"
  )))

# Convert all numeric columns to character first before replacing NAs
final_table <- final_table %>%
  mutate(across(everything(), as.character))

# Now replace NA values with empty strings for cleaner output
final_table[is.na(final_table)] <- ""

# Format p-values to match the example (with <0.001 notation)
final_table <- final_table %>%
  mutate(across(contains("_P"), ~ case_when(
    .x == "" ~ "",
    as.numeric(.x) < 0.001 ~ "<0.001",
    as.numeric(.x) >= 0.001 ~ as.character(round(as.numeric(.x), 3)),
    TRUE ~ ""
  )))

# Rename columns to match the exact format in the example
colnames(final_table) <- c(
  "Ecosystem type",
  "1990-2000 Slope", "1990-2000 P", "1990-2000 count",
  "2000-2010 Slope", "2000-2010 P", "2000-2010 count", 
  "2010-2022 Slope", "2010-2022 P", "2010-2022 count",
  "1990-2022 Slope", "1990-2022 P", "1990-2022 count"
)

# Save the table to CSV
write.csv(final_table, "Supplementary_Table_1_RS_temporal_changes.csv", 
          row.names = FALSE, na = "")

# ============================================================================
# EXTRACT INFORMATION FOR THE DESCRIPTIVE PARAGRAPH
# ============================================================================

# Calculate the predominant ecosystems in your data
ecosystem_counts <- data %>%
  filter(!is.na(Ecosystem_type2), !is.na(Rs_annual)) %>%
  count(Ecosystem_type2) %>%
  arrange(desc(n)) %>%
  mutate(percentage = round(n / sum(n) * 100, 1))

# Get the top 2 predominant ecosystems
predominant_ecosystems <- head(ecosystem_counts, 2)
predominant_text <- paste(
  predominant_ecosystems$Ecosystem_type2[1], "and", 
  predominant_ecosystems$Ecosystem_type2[2]
)

total_percentage <- sum(predominant_ecosystems$percentage)

# Extract specific results for the paragraph
get_value <- function(ecosystem, period, stat) {
  col_name <- paste(period, stat, sep = " ")
  row <- final_table[final_table$`Ecosystem type` == ecosystem, col_name]
  if (length(row) == 0 || row == "") return(NA)
  return(row)
}

# Grassland results for full period
grassland_slope_1990_2022 <- get_value("Grassland", "1990-2022", "Slope")
grassland_p_1990_2022 <- get_value("Grassland", "1990-2022", "P")

# Forest results for full period
forest_p_1990_2022 <- get_value("Forest", "1990-2022", "P")

# Individual forest types for different periods
evergreen_slope_2000_2010 <- get_value("Evergreen Forest", "2000-2010", "Slope")
evergreen_p_2000_2010 <- get_value("Evergreen Forest", "2000-2010", "P")

deciduous_p_2000_2010 <- get_value("Deciduous Forest", "2000-2010", "P")

mixed_slope_2000_2010 <- get_value("Mixed Forest", "2000-2010", "Slope")
mixed_p_2000_2010 <- get_value("Mixed Forest", "2000-2010", "P")

# Early period forest results
evergreen_p_1990_2000 <- get_value("Evergreen Forest", "1990-2000", "P")
deciduous_p_1990_2000 <- get_value("Deciduous Forest", "1990-2000", "P")
mixed_p_1990_2000 <- get_value("Mixed Forest", "1990-2000", "P")

# Helper function to format p-value statements
format_p_statement <- function(p_value) {
  if (is.na(p_value) || p_value == "") return("P value not available")
  p_num <- as.numeric(gsub("<", "", p_value))
  if (p_value == "<0.001" || p_num < 0.001) {
    return("P < 0.001")
  } else if (p_num > 0.050) {
    return("P > 0.050")
  } else {
    return(paste0("P = ", p_value))
  }
}

# Build the descriptive paragraph
paragraph <- paste(
  "The RS data in our dataset are mainly collected from ", predominant_ecosystems$Ecosystem_type2[1], 
  " and ", predominant_ecosystems$Ecosystem_type2[2], ", which cover ", total_percentage, 
  "% of the total observations. ",
  sep = ""
)

# Add grassland information if available
if (!is.na(grassland_slope_1990_2022) && !is.na(grassland_p_1990_2022)) {
  p_stmt <- format_p_statement(grassland_p_1990_2022)
  if (grepl("< 0.05", p_stmt) || (grepl("=", p_stmt) && as.numeric(gsub("P = ", "", p_stmt)) < 0.05)) {
    change_dir <- ifelse(as.numeric(grassland_slope_1990_2022) > 0, "increased", "decreased")
    paragraph <- paste(
      paragraph,
      "RS in grasslands worldwide ", change_dir, " at a rate of ", abs(as.numeric(grassland_slope_1990_2022)), 
      " g C m⁻² yr⁻² during 1990-2022 (", p_stmt, ", Supplementary Table 1). ",
      sep = ""
    )
  } else {
    paragraph <- paste(
      paragraph,
      "RS in grasslands worldwide remained unchanged during 1990-2022 (", p_stmt, "). ",
      sep = ""
    )
  }
}

# Add forest information
if (!is.na(forest_p_1990_2022)) {
  p_stmt <- format_p_statement(forest_p_1990_2022)
  if (grepl("> 0.050", p_stmt)) {
    paragraph <- paste(paragraph, "In contrast, RS in forests worldwide remained unchanged (", p_stmt, "). ", sep = "")
  } else {
    change_dir <- ifelse(as.numeric(get_value("Forest", "1990-2022", "Slope")) > 0, "increased", "decreased")
    paragraph <- paste(
      paragraph, "In contrast, RS in forests worldwide ", change_dir, " (", p_stmt, "). ", 
      sep = ""
    )
  }
}

# Add early period forest results
paragraph <- paste(
  paragraph,
  "When forests were divided into evergreen forests, deciduous forests, and mixed forests, we found that RS ",
  sep = ""
)

early_period_changes <- c()
if (!is.na(evergreen_p_1990_2000) && !grepl("> 0.050", format_p_statement(evergreen_p_1990_2000))) {
  early_period_changes <- c(early_period_changes, "evergreen forests")
}
if (!is.na(deciduous_p_1990_2000) && !grepl("> 0.050", format_p_statement(deciduous_p_1990_2000))) {
  early_period_changes <- c(early_period_changes, "deciduous forests")
}
if (!is.na(mixed_p_1990_2000) && !grepl("> 0.050", format_p_statement(mixed_p_1990_2000))) {
  early_period_changes <- c(early_period_changes, "mixed forests")
}

if (length(early_period_changes) == 0) {
  paragraph <- paste(paragraph, "remained unchanged in 1990-2000 in all three forest types (P > 0.050, Supplementary Table 1). ", sep = "")
} else {
  paragraph <- paste(
    paragraph,
    "remained unchanged in 1990-2000 in some forest types but showed changes in ", 
    paste(early_period_changes, collapse = " and "), 
    " (Supplementary Table 1). ",
    sep = ""
  )
}

# Add middle period forest results
middle_period_results <- c()
if (!is.na(evergreen_p_2000_2010) && !grepl("> 0.050", format_p_statement(evergreen_p_2000_2010))) {
  change_dir <- ifelse(as.numeric(evergreen_slope_2000_2010) > 0, "increased", "decreased")
  middle_period_results <- c(middle_period_results, 
    paste("RS ", change_dir, " (", format_p_statement(evergreen_p_2000_2010), ") at a rate of ", 
          abs(as.numeric(evergreen_slope_2000_2010)), " g C m⁻² yr⁻² in evergreen forests", sep = ""))
}

if (!is.na(deciduous_p_2000_2010) && grepl("> 0.050", format_p_statement(deciduous_p_2000_2010))) {
  middle_period_results <- c(middle_period_results, "RS remained unchanged in deciduous forests")
} else if (!is.na(deciduous_p_2000_2010)) {
  change_dir <- ifelse(as.numeric(get_value("Deciduous Forest", "2000-2010", "Slope")) > 0, "increased", "decreased")
  middle_period_results <- c(middle_period_results, 
    paste("RS ", change_dir, " in deciduous forests (", format_p_statement(deciduous_p_2000_2010), ")", sep = ""))
}

if (!is.na(mixed_p_2000_2010) && !grepl("> 0.050", format_p_statement(mixed_p_2000_2010))) {
  change_dir <- ifelse(as.numeric(mixed_slope_2000_2010) > 0, "increased", "decreased")
  middle_period_results <- c(middle_period_results, 
    paste("RS ", change_dir, " (", format_p_statement(mixed_p_2000_2010), ") at a rate of ", 
          abs(as.numeric(mixed_slope_2000_2010)), " g C m⁻² yr⁻² in mixed forests", sep = ""))
} else if (!is.na(mixed_p_2000_2010)) {
  middle_period_results <- c(middle_period_results, "RS remained unchanged in mixed forests")
}

if (length(middle_period_results) > 0) {
  paragraph <- paste(
    paragraph,
    "In 2000-2010, ",
    paste(middle_period_results, collapse = ", "),
    ". ",
    sep = ""
  )
}

# Final sentence about geographical patterns
paragraph <- paste(
  paragraph,
  "These results provide insights into the temporal changes of RS across different ecosystems during 1990-2022.",
  sep = ""
)

# Create a data frame for the paragraph and save it
paragraph_df <- data.frame(
  Description = "Results summary paragraph",
  Content = paragraph
)

# Save the paragraph to CSV
write.csv(paragraph_df, "RS_temporal_changes_summary_paragraph.csv", row.names = FALSE)

# Also save the ecosystem counts
write.csv(ecosystem_counts, "ecosystem_distribution_stats.csv", row.names = FALSE)

# Print everything
cat("=== SUPPLEMENTARY TABLE 1 ===\n")
print(final_table, width = Inf)

cat("\n=== ECOSYSTEM DISTRIBUTION ===\n")
print(ecosystem_counts)

cat("\n=== DESCRIPTIVE PARAGRAPH ===\n")
cat(paragraph)
cat("\n\n")

cat("Files saved:\n")
cat("1. Supplementary_Table_1_RS_temporal_changes.csv\n")
cat("2. RS_temporal_changes_summary_paragraph.csv\n")
cat("3. ecosystem_distribution_stats.csv\n")





# PLOT BIT
# FULL, SELF-CONTAINED, ROBUST SCRIPT - BAR PLOTS FINAL VERSION
# Copy-paste and run as-is (requires tidyverse)

library(tidyverse)

# ------------ USER CONTROLS ------------
csv_file <- "Supplementary_Table_1_RS_temporal_changes.csv"
TOP_N <- 5   # set to 3 or 5 depending on how many ecosystems you want to show
# ---------------------------------------

# 1) read all columns as character to avoid type-mixing problems
df_raw <- read_csv(csv_file, col_types = cols(.default = "c"), show_col_types = FALSE)

# 2) select only the 3 periods you said are clean
wanted_cols <- c(
  "Ecosystem type",
  "1990-2000 Slope", "1990-2000 P", "1990-2000 count",
  "2000-2010 Slope", "2000-2010 P", "2000-2010 count",
  "2010-2022 Slope", "2010-2022 P", "2010-2022 count"
)
present_cols <- intersect(wanted_cols, names(df_raw))
df <- df_raw %>% select(all_of(present_cols))

# 3) remove rows with empty ecosystem type
df <- df %>% filter(!is.na(`Ecosystem type`) & str_trim(`Ecosystem type`) != "")

# 4) pivot: safe reshape because all columns are character
long <- df %>%
  pivot_longer(
    cols = -`Ecosystem type`,
    names_to = c("Period", "Stat"),
    names_pattern = "(.*)\\s+(Slope|P|count)$",
    values_to = "value",
    values_drop_na = FALSE
  ) %>%
  pivot_wider(names_from = Stat, values_from = value)

# 5) clean numeric columns: handle "<0.001" and blanks
clean_num <- function(x) {
  x <- ifelse(x == "" | is.na(x), NA, x)
  x <- ifelse(x == "<0.001", "0.0009", x)
  suppressWarnings(as.numeric(x))
}

plot_df <- long %>%
  mutate(
    Slope = clean_num(Slope),
    count = clean_num(count),
    P = clean_num(P)
  )

# 6) create significance flag; drop rows lacking Slope or count
plot_df <- plot_df %>%
  mutate(
    Significance = case_when(
      P < 0.001 ~ "***",
      P < 0.01 ~ "**",
      P < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  filter(!is.na(Slope) & !is.na(count))

# 7) choose top N ecosystems by total sample count across the three periods
# BUT REMOVE THE GENERIC "FOREST" CATEGORY
ecosum <- plot_df %>%
  filter(`Ecosystem type` != "Forest") %>%  # Remove generic Forest
  group_by(`Ecosystem type`) %>%
  summarise(total_count = sum(count, na.rm = TRUE)) %>%
  arrange(desc(total_count))

top_ecos <- head(ecosum$`Ecosystem type`, TOP_N)

plot_sel <- plot_df %>% 
  filter(`Ecosystem type` %in% top_ecos) %>%
  filter(`Ecosystem type` != "Forest")  # Double-check removal

# 8) ensure Period ordering and ALPHABETICAL order for ecosystems
period_levels <- c("1990-2000", "2000-2010", "2010-2022")
plot_sel <- plot_sel %>% 
  mutate(
    Period = factor(Period, levels = intersect(period_levels, unique(Period))),
    `Ecosystem type` = factor(`Ecosystem type`, levels = sort(unique(`Ecosystem type`), decreasing = TRUE))
  )

# 9) HORIZONTAL BAR PLOT WITH YOUR COLOR SCHEME
p_bar <- ggplot(plot_sel, aes(x = Slope, y = `Ecosystem type`, fill = Period)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
  # Add asterisks only for significant values
  geom_text(
    data = plot_sel %>% filter(Significance != ""),
    aes(label = Significance, x = Slope + ifelse(Slope >= 0, max(Slope)*0.05, min(Slope)*0.05)),
    position = position_dodge(width = 0.8),
    size = 4,
    vjust = 0.5,
    hjust = ifelse(plot_sel %>% filter(Significance != "") %>% pull(Slope) >= 0, 0, 1)
  ) +
  # Your established color scheme with red as highest
  scale_fill_manual(
    values = c(
      "1990-2000" = "#0000a2",  # Your dark blue
      "2000-2010" = "#50ad9f",  # Your teal  
      "2010-2022" = "#bc272d"   # Your red for most recent/highest
    ),
    name = ""
  ) +
  # Fix x-axis scale to prevent cut-off zeros
  scale_x_continuous(
    limits = c(min(plot_sel$Slope) * 1.1, max(plot_sel$Slope) * 1.1),
    expand = expansion(mult = 0.05)
  ) +
  labs(
    x = expression(paste("R"[s], " Slope (g C m"^{-2}, " yr"^{-1}, ")")),
    y = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),  # Remove grid lines
    panel.grid.minor = element_blank(),  # Remove grid lines
    legend.position = "bottom",
    axis.text.y = element_text(face = "bold"),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),  # Add x-axis line
    axis.ticks.x = element_line(color = "black", linewidth = 0.5)   # Add x-axis ticks
  )

print(p_bar)

# Save the bar plot
ggsave("rs_trends_barplot_final.png", p_bar, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("rs_trends_barplot_final.pdf", p_bar, width = 8, height = 6, bg = "white")
ggsave("rs_trends_barplot_final.tiff", p_bar, width = 8, height = 6, dpi = 300, bg = "white")
# Show what ecosystems made the final cut
cat("Final ecosystems in the plot (alphabetical order):\n")
print(sort(unique(plot_sel$`Ecosystem type`)))

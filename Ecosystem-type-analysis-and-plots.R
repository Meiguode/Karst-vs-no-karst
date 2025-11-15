library(dplyr)
library(tidyr)
library(broom)
library(stringr)

data <- read.csv("complete_data.csv")
data <- data %>%
  mutate(period = case_when(
    Study_midyear >= 1990 & Study_midyear <= 2000 ~ "1990-2000",
    Study_midyear >= 2000 & Study_midyear <= 2010 ~ "2000-2010", 
    Study_midyear >= 2010 & Study_midyear <= 2022 ~ "2010-2022",
    Study_midyear >= 1990 & Study_midyear <= 2022 ~ "1990-2022"
  ))

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

required_columns <- c(
  "1990-2000_Slope", "1990-2000_P", "1990-2000_count",
  "2000-2010_Slope", "2000-2010_P", "2000-2010_count", 
  "2010-2022_Slope", "2010-2022_P", "2010-2022_count",
  "1990-2022_Slope", "1990-2022_P", "1990-2022_count"
)

for (col in required_columns) {
  if (!col %in% colnames(results)) {
    results[[col]] <- NA
  }
}

results <- results %>%
  select(
    Ecosystem_type = Ecosystem_type2,
    `1990-2000_Slope`, `1990-2000_P`, `1990-2000_count`,
    `2000-2010_Slope`, `2000-2010_P`, `2000-2010_count`,
    `2010-2022_Slope`, `2010-2022_P`, `2010-2022_count`,
    `1990-2022_Slope`, `1990-2022_P`, `1990-2022_count`
  )

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

final_table <- bind_rows(forest_results, results) %>%
  arrange(factor(Ecosystem_type, levels = c(
    "Forest", "Deciduous Forest", "Evergreen Forest", "Mixed Forest",
    "Grassland", "Shrubland", "Wetland", "Others"
  )))

final_table <- final_table %>%
  mutate(across(everything(), as.character))

final_table[is.na(final_table)] <- ""

final_table <- final_table %>%
  mutate(across(contains("_P"), ~ case_when(
    .x == "" ~ "",
    as.numeric(.x) < 0.001 ~ "<0.001",
    as.numeric(.x) >= 0.001 ~ as.character(round(as.numeric(.x), 3)),
    TRUE ~ ""
  )))

colnames(final_table) <- c(
  "Ecosystem type",
  "1990-2000 Slope", "1990-2000 P", "1990-2000 count",
  "2000-2010 Slope", "2000-2010 P", "2000-2010 count", 
  "2010-2022 Slope", "2010-2022 P", "2010-2022 count",
  "1990-2022 Slope", "1990-2022 P", "1990-2022 count"
)

write.csv(final_table, "Supplementary_Table_1_RS_temporal_changes.csv", 
          row.names = FALSE, na = "")



# PLOT BIT
library(tidyverse)

csv_file <- "Supplementary_Table_1_RS_temporal_changes.csv"
TOP_N <- 5

df_raw <- read_csv(csv_file, col_types = cols(.default = "c"), show_col_types = FALSE)

wanted_cols <- c(
  "Ecosystem type",
  "1990-2000 Slope", "1990-2000 P", "1990-2000 count",
  "2000-2010 Slope", "2000-2010 P", "2000-2010 count",
  "2010-2022 Slope", "2010-2022 P", "2010-2022 count"
)
present_cols <- intersect(wanted_cols, names(df_raw))
df <- df_raw %>% select(all_of(present_cols))

df <- df %>% filter(!is.na(`Ecosystem type`) & str_trim(`Ecosystem type`) != "")

long <- df %>%
  pivot_longer(
    cols = -`Ecosystem type`,
    names_to = c("Period", "Stat"),
    names_pattern = "(.*)\\s+(Slope|P|count)$",
    values_to = "value",
    values_drop_na = FALSE
  ) %>%
  pivot_wider(names_from = Stat, values_from = value)

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

# Remove generic forest category to reduce ambiguity
ecosum <- plot_df %>%
  filter(`Ecosystem type` != "Forest") %>%  # Remove generic Forest
  group_by(`Ecosystem type`) %>%
  summarise(total_count = sum(count, na.rm = TRUE)) %>%
  arrange(desc(total_count))

top_ecos <- head(ecosum$`Ecosystem type`, TOP_N)

plot_sel <- plot_df %>% 
  filter(`Ecosystem type` %in% top_ecos) %>%
  filter(`Ecosystem type` != "Forest")  # Double-check removal

period_levels <- c("1990-2000", "2000-2010", "2010-2022")
plot_sel <- plot_sel %>% 
  mutate(
    Period = factor(Period, levels = intersect(period_levels, unique(Period))),
    `Ecosystem type` = factor(`Ecosystem type`, levels = sort(unique(`Ecosystem type`), decreasing = TRUE))
  )

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

ggsave("rs_trends_barplot_final.png", p_bar, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("rs_trends_barplot_final.pdf", p_bar, width = 8, height = 6, bg = "white")
ggsave("rs_trends_barplot_final.tiff", p_bar, width = 8, height = 6, dpi = 300, bg = "white")

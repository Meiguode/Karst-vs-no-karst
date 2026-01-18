library(dplyr)
library(ggplot2)
library(ggthemes)
library(mblm)
library(patchwork)

source("statistical_prep.R")

datasets <- list(
  "Combined" = combined,
  "Non-karst" = nonkarst,
  "Karst" = karst
)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

MAP_LIMS <- c(-500, 1500)
MAT_LIMS <- c(-3, 8)
RS_LIMS  <- c(0, 3000)

cols <- list(
  MAP = "#0000a2",
  MAT = "#50ad9f"
)

all_plots <- list()

for (i in seq_along(names(datasets))) {
  dname <- names(datasets)[i]
  df <- datasets[[dname]]
  
  if ("delta.Mean.Precip." %in% names(df)) {
    dat <- df %>%
      filter(
        !is.na(Rs_annual),
        !is.na(delta.Mean.Precip.),
        is.finite(Rs_annual),
        is.finite(delta.Mean.Precip.),
        delta.Mean.Precip. >= MAP_LIMS[1],
        delta.Mean.Precip. <= MAP_LIMS[2],
        Rs_annual <= RS_LIMS[2]
      )
    
    if (nrow(dat) >= 20) {
      ts_fit <- tryCatch({
        mblm(Rs_annual ~ delta.Mean.Precip., dat)
      }, error = function(e) {
        cat("Warning: mblm failed for", dname, "MAP:", e$message, "\n")
        lm(Rs_annual ~ delta.Mean.Precip., dat)
      })

      ts_slope <- coef(ts_fit)[2]

      library(Kendall)
      tau_test <- MannKendall(dat$delta.Mean.Precip.)
      ts_pvalue <- tau_test$sl[1]

      if (!exists("theilsen_results")) {
        theilsen_results <- data.frame(
          Dataset = character(),
          Variable = character(),
          N = integer(),
          Slope = numeric(),
          P_value = numeric(),
          stringsAsFactors = FALSE
        )
      }

      theilsen_results <- rbind(theilsen_results, data.frame(
        Dataset = dname,
        Variable = "Delta MAP",  # Fixed: use "Delta MAP" for this section
        N = nrow(dat),
        Slope = ts_slope,
        P_value = ts_pvalue
      ))
      
      p <- ggplot(dat, aes(delta.Mean.Precip., Rs_annual)) +
        geom_point(color = cols$MAP, alpha = 0.25, size = 1.4) +
        geom_abline(
          slope = coef(ts_fit)[2],
          intercept = coef(ts_fit)[1],
          linetype = "dashed",
          linewidth = 1
        ) +
        theme_few() +
        labs(
          x = expression(Delta~MAP~"(mm)"),
          y = expression(R[s]~"(g C m"^{-2}~"yr"^{-1}*")")
        ) +
        scale_x_continuous(limits = MAP_LIMS) +
        scale_y_continuous(limits = RS_LIMS) +
        theme(
          text = element_text(size = 9),
          panel.border = element_rect(color = "grey20", fill = NA),
          plot.title = element_blank()
        )
      
      # Store with column label as the key
      all_plots[[paste0("MAP_", dname)]] <- p
    }
  }
  
  if ("delta.Mean.Temp." %in% names(df)) {
    dat <- df %>%
      filter(
        !is.na(Rs_annual),
        !is.na(delta.Mean.Temp.),
        is.finite(Rs_annual),
        is.finite(delta.Mean.Temp.),
        delta.Mean.Temp. >= MAT_LIMS[1],
        delta.Mean.Temp. <= MAT_LIMS[2],
        Rs_annual <= RS_LIMS[2]
      )
    
    if (nrow(dat) >= 20) {
      ts_fit <- tryCatch({
        mblm(Rs_annual ~ delta.Mean.Temp., dat)
      }, error = function(e) {
        cat("Warning: mblm failed for", dname, "MAT:", e$message, "\n")
        lm(Rs_annual ~ delta.Mean.Temp., dat)
      })

      ts_slope <- coef(ts_fit)[2]

      library(Kendall)
      tau_test <- MannKendall(dat$delta.Mean.Temp.)
      ts_pvalue <- tau_test$sl[1]

      if (!exists("theilsen_results")) {
        theilsen_results <- data.frame(
          Dataset = character(),
          Variable = character(),
          N = integer(),
          Slope = numeric(),
          P_value = numeric(),
          stringsAsFactors = FALSE
        )
      }

      theilsen_results <- rbind(theilsen_results, data.frame(
        Dataset = dname,
        Variable = "Delta MAT",
        N = nrow(dat),
        Slope = ts_slope,
        P_value = ts_pvalue
      ))
      
      p <- ggplot(dat, aes(delta.Mean.Temp., Rs_annual)) +
        geom_point(color = cols$MAT, alpha = 0.25, size = 1.4) +
        geom_abline(
          slope = coef(ts_fit)[2],
          intercept = coef(ts_fit)[1],
          linetype = "dashed",
          linewidth = 1
        ) +
        theme_few() +
        labs(
          x = expression(Delta~MAT~"(" * degree * "C)"),
          y = expression(R[s]~"(g C m"^{-2}~"yr"^{-1}*")")
        ) +
        scale_x_continuous(limits = MAT_LIMS) +
        scale_y_continuous(limits = RS_LIMS) +
        theme(
          text = element_text(size = 9),
          panel.border = element_rect(color = "grey20", fill = NA),
          plot.title = element_blank()
        )
      
      all_plots[[paste0("MAT_", dname)]] <- p
    }
  }
}

column_order <- c("Combined", "Non-karst", "Karst")
column_labels <- c("(a)", "(b)", "(c)")

add_column_label <- function(plot, label) {
  plot + 
    theme(plot.margin = margin(t = 20)) +  # Add space at top for label
    # Add column label at the top
    annotate("text", x = mean(MAP_LIMS), y = RS_LIMS[2] * 1.05,
             label = label, size = 4, fontface = "bold", vjust = 0)
}

if (length(all_plots) >= 6) {
  
  top_row_plots <- list()
  bottom_row_plots <- list()
  
  for (i in seq_along(column_order)) {
    dname <- column_order[i]
    label <- column_labels[i]
    
    map_key <- paste0("MAP_", dname)
    if (map_key %in% names(all_plots)) {
      top_row_plots[[i]] <- add_column_label(all_plots[[map_key]], label)
    }
    
    mat_key <- paste0("MAT_", dname)
    if (mat_key %in% names(all_plots)) {
      bottom_row_plots[[i]] <- all_plots[[mat_key]]
    }
  }
  
  top_row <- wrap_plots(top_row_plots, ncol = 3)
  
  bottom_row <- wrap_plots(bottom_row_plots, ncol = 3)
  
  combined_figure <- top_row / bottom_row
  
  ggsave(
    "figures/combined_delta_climate_all.tiff",
    combined_figure,
    width = 9,
    height = 5,
    units = "in",
    dpi = 900,
    compression = "lzw"
  )

results_df <- data.frame()

add_results <- function(dataset, variable, period, method, n, slope, p, r2) {
  results_df <<- rbind(
    results_df,
    data.frame(
      Dataset = dataset,
      Variable = variable,
      Period = period,
      Method = method,
      N = n,
      Slope = slope,
      P_value = p,
      R_squared = r2,
      stringsAsFactors = FALSE
    )
  )
}

for (dname in names(datasets)) {
  df <- datasets[[dname]]
  
  for (per in c("All", unique(df$Period))) {
    dfp <- if (per == "All") df else filter(df, Period == per)
    
    if ("delta.Mean.Precip." %in% names(dfp)) {
      dat <- dfp %>% filter(
        !is.na(Rs_annual), !is.na(delta.Mean.Precip.),
        delta.Mean.Precip. >= MAP_LIMS[1],
        delta.Mean.Precip. <= MAP_LIMS[2]
      )
      
      if (nrow(dat) >= 20) {
        fit <- lm(Rs_annual ~ delta.Mean.Precip., dat)
        add_results(dname, "Delta MAP", per, "OLS",
                    nrow(dat), coef(fit)[2],
                    summary(fit)$coefficients[2,4],
                    summary(fit)$r.squared)
      }
    }
    
    if ("delta.Mean.Temp." %in% names(dfp)) {
      dat <- dfp %>% filter(
        !is.na(Rs_annual), !is.na(delta.Mean.Temp.),
        delta.Mean.Temp. >= MAT_LIMS[1],
        delta.Mean.Temp. <= MAT_LIMS[2]
      )
      
      if (nrow(dat) >= 20) {
        fit <- lm(Rs_annual ~ delta.Mean.Temp., dat)
        add_results(dname, "Delta MAT", per, "OLS",
                    nrow(dat), coef(fit)[2],
                    summary(fit)$coefficients[2,4],
                    summary(fit)$r.squared)
      }
    }
  }
}

write.csv(
  results_df,
  "results/statistical_results_delta_climate_fullperiod_plus_breakpoints.csv",
  row.names = FALSE
)

library(sf)
library(dplyr)
library(readr)
library(ggplot2)
library(future)
library(furrr)
library(purrr)

cat("Reading karst aquifers shapefile...\n")
karst_aquifers <- st_read("/mnt/d/Karst aquifers/WHYMAP_WOKAM/shp/whymap_karst__v1_poly.shp", quiet = TRUE)

cat("Checking geometry validity...\n")
invalid_geometries <- st_is_valid(karst_aquifers)
invalid_geometries_count <- sum(!invalid_geometries, na.rm = TRUE)
cat("Number of invalid geometries before fixing:", invalid_geometries_count, "\n")

if (invalid_geometries_count > 0) {
  cat("Fixing invalid geometries...\n")
  
  karst_aquifers_valid <- tryCatch({
    st_make_valid(karst_aquifers)
  }, error = function(e) {
    cat("st_make_valid failed:", e$message, "\n")
    return(karst_aquifers)
  })
  
  invalid_after_make_valid <- sum(!st_is_valid(karst_aquifers_valid), na.rm = TRUE)
  cat("Invalid geometries after st_make_valid:", invalid_after_make_valid, "\n")
  
  if (invalid_after_make_valid > 0) {
    cat("Trying geometry simplification...\n")
    karst_aquifers_simple <- tryCatch({
      st_simplify(karst_aquifers_valid, preserveTopology = TRUE, dTolerance = 0.00001)
    }, error = function(e) {
      cat("Simplification failed:", e$message, "\n")
      return(karst_aquifers_valid)
    })
    
    karst_aquifers <- karst_aquifers_simple
  } else {
    karst_aquifers <- karst_aquifers_valid
  }
  
  final_invalid_count <- sum(!st_is_valid(karst_aquifers), na.rm = TRUE)
  cat("Invalid geometries after simplification:", final_invalid_count, "\n")
  
  if (final_invalid_count > 0) {
    cat("Removing remaining invalid geometries...\n")
    valid_indices <- which(st_is_valid(karst_aquifers))
    karst_aquifers <- karst_aquifers[valid_indices, ]
    cat("Retained", nrow(karst_aquifers), "valid geometries out of", length(invalid_geometries), "\n")
  }
}

final_check <- sum(!st_is_valid(karst_aquifers), na.rm = TRUE)
cat("Final number of invalid geometries:", final_check, "\n")

if (final_check > 0) {
  cat("Warning: Some geometries remain invalid. Proceeding with valid subset only.\n")
  karst_aquifers <- karst_aquifers[st_is_valid(karst_aquifers), ]
}

cat("Reading soil respiration data...\n")
soil_data <- read_csv("/mnt/d/Soil_Respiration_Analysis/complete_data.csv")

cat("Converting soil data to spatial object...\n")
soil_data_sf <- st_as_sf(soil_data, coords = c("Longitude", "Latitude"), crs = 4326)

if (st_crs(karst_aquifers) != st_crs(soil_data_sf)) {
  cat("Transforming CRS...\n")
  soil_data_sf <- st_transform(soil_data_sf, st_crs(karst_aquifers))
}

cat("Splitting data into chunks...\n")
chunk_size <- ceiling(nrow(soil_data_sf) / 5)
soil_data_chunks <- split(soil_data_sf, rep(1:5, each = chunk_size, length.out = nrow(soil_data_sf)))

check_intersection <- function(chunk, karst_aquifers) {
  cat("Processing chunk with", nrow(chunk), "points...\n")
  
  result <- tryCatch({
    intersects <- st_intersects(chunk, karst_aquifers, sparse = TRUE)
    intersecting_indices <- which(apply(intersects, 1, any))
    if (length(intersecting_indices) > 0) {
      chunk[intersecting_indices, ]
    } else {
      NULL
    }
  }, error = function(e) {
    cat("Standard intersection failed:", e$message, "\n")
    
    result2 <- tryCatch({
      cat("Trying sparse = FALSE method...\n")
      intersects_matrix <- st_intersects(chunk, karst_aquifers, sparse = FALSE)
      intersecting_indices <- which(apply(intersects_matrix, 1, any))
      if (length(intersecting_indices) > 0) {
        chunk[intersecting_indices, ]
      } else {
        NULL
      }
    }, error = function(e2) {
      cat("Sparse = FALSE also failed:", e2$message, "\n")
      
      cat("Trying distance-based approach...\n")
      result3 <- tryCatch({
        # Find points within a small distance (10 meters)
        nearest <- st_is_within_distance(chunk, karst_aquifers, dist = 10)
        intersecting_indices <- which(apply(nearest, 1, any))
        if (length(intersecting_indices) > 0) {
          chunk[intersecting_indices, ]
        } else {
          NULL
        }
      }, error = function(e3) {
        cat("All intersection methods failed for this chunk\n")
        return(NULL)
      })
      return(result3)
    })
    return(result2)
  })
  
  if (!is.null(result)) {
    cat("Found", nrow(result), "intersecting points in this chunk\n")
  } else {
    cat("No intersecting points found in this chunk\n")
  }
  
  return(result)
}

cat("Starting intersection analysis...\n")
karst_points_list <- list()

for (i in seq_along(soil_data_chunks)) {
  cat(paste0("Processing chunk ", i, "/", length(soil_data_chunks), "...\n"))
  chunk_result <- check_intersection(soil_data_chunks[[i]], karst_aquifers)
  if (!is.null(chunk_result)) {
    karst_points_list[[i]] <- chunk_result
  }
}

cat("Combining results from all chunks...\n")
if (length(karst_points_list) > 0) {
  karst_points <- do.call(rbind, karst_points_list)
  cat("Total karst points found:", nrow(karst_points), "\n")
} else {
  karst_points <- NULL
  cat("No karst points found in any chunk\n")
}

if (!is.null(karst_points) && nrow(karst_points) > 0) {
  cat("Adding karst aquifer attributes...\n")
  
  karst_points_with_attr <- tryCatch({
    st_join(karst_points, karst_aquifers)
  }, error = function(e) {
    cat("Could not add aquifer attributes with st_join:", e$message, "\n")
    # Fallback: manual attribute assignment
    cat("Trying manual attribute assignment...\n")
    return(karst_points)
  })
  
  cat("Columns in karst_points_with_attr:", paste(names(karst_points_with_attr), collapse = ", "), "\n")
  
  # Step 11: Save karst points with ALL attributes
  cat("Saving karst points...\n")
  write_csv(st_drop_geometry(karst_points_with_attr), "/mnt/d/Karst_aquifers_with_soil_data.csv")
} else {
  karst_points_with_attr <- NULL
}

cat("Finding non-karst points...\n")

karst_indices <- numeric(0)

if (!is.null(karst_points) && nrow(karst_points) > 0) {
  # Get the indices of karst points in the original soil_data_sf
  # We need to match by coordinates or unique identifier
  
  if ("ID" %in% names(soil_data_sf)) {
    karst_indices <- which(soil_data_sf$ID %in% karst_points$ID)
  } else {
    # Method 2: Match by coordinates (more reliable)
    cat("Matching points by coordinates...\n")
    soil_coords <- st_coordinates(soil_data_sf)
    karst_coords <- st_coordinates(karst_points)
    
    # Create a unique key for coordinates (rounded to 6 decimal places for matching)
    soil_keys <- apply(round(soil_coords, 6), 1, paste, collapse = "_")
    karst_keys <- apply(round(karst_coords, 6), 1, paste, collapse = "_")
    
    karst_indices <- which(soil_keys %in% karst_keys)
  }
  
  cat("Found", length(karst_indices), "karst point indices in original data\n")
}

if (length(karst_indices) > 0) {
  non_karst_points <- soil_data_sf[-karst_indices, ]
  cat("Non-karst points found:", nrow(non_karst_points), "\n")
} else {
  non_karst_points <- soil_data_sf
  cat("No karst points found, so all points are non-karst:", nrow(non_karst_points), "\n")
}

if (!is.null(non_karst_points) && nrow(non_karst_points) > 0) {
  cat("Saving non-karst points...\n")
  write_csv(st_drop_geometry(non_karst_points), "/mnt/d/Non_karst_aquifers_with_soil_data.csv")
} else {
  cat("No non-karst points to save\n")
}

cat("Creating visualization...\n")
tryCatch({
  p <- ggplot() +
    geom_sf(data = karst_aquifers, fill = "lightblue", color = "blue", alpha = 0.5) +
    geom_sf(data = soil_data_sf, color = "gray", size = 0.3, alpha = 0.6) +
    theme_minimal() +
    ggtitle("Karst Aquifers and Soil Data Points")
  
  if (!is.null(karst_points) && nrow(karst_points) > 0) {
    p <- p + geom_sf(data = karst_points, color = "red", size = 0.5)
  }
  
  print(p)
}, error = function(e) {
  cat("Visualization failed:", e$message, "\n")
})

cat("\n=== PROCESSING COMPLETE ===\n")
if (!is.null(karst_points_with_attr) && nrow(karst_points_with_attr) > 0) {
  cat("Karst points found:", nrow(karst_points_with_attr), "\n")
  cat("Karst data saved to: /mnt/d/Karst_aquifers_with_soil_data.csv\n")
  
  karst_cols <- names(karst_points_with_attr)
  cat("Available attribute columns in karst data:", paste(karst_cols, collapse = ", "), "\n")
} else {
  cat("No karst points found.\n")
}

if (!is.null(non_karst_points) && nrow(non_karst_points) > 0) {
  cat("Non-karst points found:", nrow(non_karst_points), "\n")
  cat("Non-karst data saved to: /mnt/d/Non_karst_aquifers_with_soil_data.csv\n")
} else {
  cat("No non-karst points found.\n")
}

cat("Total soil data points processed:", nrow(soil_data_sf), "\n")
cat("Breakdown:\n")
cat("- Karst points:", ifelse(!is.null(karst_points), nrow(karst_points), 0), "\n")
cat("- Non-karst points:", ifelse(!is.null(non_karst_points), nrow(non_karst_points), 0), "\n")
cat("- Sum check:", ifelse(!is.null(karst_points) && !is.null(non_karst_points), 
                           nrow(karst_points) + nrow(non_karst_points), "N/A"), 
    "(should equal", nrow(soil_data_sf), ")\n")

tryCatch({
  publication_map <- ggplot() +
    # Plot non-karst points first (in background)
    geom_sf(data = non_karst_points, color = "gray60", size = 0.8, alpha = 0.7, shape = 16) +
    # Plot karst aquifers
    geom_sf(data = karst_aquifers, fill = "lightblue", color = "blue", alpha = 0.3, linewidth = 0.1) +
    # Plot karst points on top (in foreground)
    geom_sf(data = karst_points, color = "red", size = 0.8, alpha = 0.8, shape = 16) +
    # Use a clean theme suitable for publications
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "none"  # No legend for clean look
    ) +
    coord_sf(datum = st_crs(karst_aquifers))
  
  ggsave(
    filename = "/mnt/d/Karst_Aquifers_Soil_Data_Map.tiff",
    plot = publication_map,
    device = "tiff",
    width = 8,           # inches
    height = 6,          # inches
    units = "in",
    dpi = 600,           # High resolution for publications
    compression = "lzw"   # Lossless compression
  )
  
  ggsave(
    filename = "/mnt/d/Karst_Aquifers_Soil_Data_Map.pdf",
    plot = publication_map,
    device = "pdf",
    width = 8,
    height = 6
  )
  ggsave(
    filename = "/mnt/d/Karst_Aquifers_Soil_Data_Map.png",
    plot = publication_map,
    device = "png",
    width = 8,
    height = 6,
    dpi = 300
  )
  
}, error = function(e) {
  cat("Error creating publication map:", e$message, "\n")
})

#-
library(dplyr)

workingpath <- getwd()

cat("=== STEP 1: LOADING ALL DATASETS ===\n")

main_data <- read.csv("complete_data_with_climate_soc.csv", check.names = FALSE)
cat("Main data dimensions:", dim(main_data), "\n")
cat("Main data columns:", paste(names(main_data)[1:10], collapse = ", "), "...\n\n")

karst_data <- read.csv("Karst_aquifers_with_soil_data.csv", check.names = FALSE)
cat("Karst data dimensions:", dim(karst_data), "\n")
cat("Karst data columns:", paste(names(karst_data)[1:10], collapse = ", "), "...\n\n")

nonkarst_data <- read.csv("Non_karst_aquifers_with_soil_data.csv", check.names = FALSE)
cat("Non-karst data dimensions:", dim(nonkarst_data), "\n")
cat("Non-karst data columns:", paste(names(nonkarst_data)[1:10], collapse = ", "), "...\n\n")

cat("=== STEP 2: EXTRACTING RECORD NUMBERS ===\n")

karst_record_numbers <- unique(karst_data$Record_number)
nonkarst_record_numbers <- unique(nonkarst_data$Record_number)

cat("Unique karst record numbers:", length(karst_record_numbers), "\n")
cat("Unique non-karst record numbers:", length(nonkarst_record_numbers), "\n")

overlap <- intersect(karst_record_numbers, nonkarst_record_numbers)
cat("Overlapping record numbers:", length(overlap), "\n")
if (length(overlap) > 0) {
  cat("Overlapping records:", paste(overlap, collapse = ", "), "\n")
}

cat("\n=== STEP 3: CREATING KARST AND NON-KARST SUBSETS ===\n")

karst_subset <- main_data %>%
  dplyr::filter(Record_number %in% karst_record_numbers)

cat("Karst subset dimensions:", dim(karst_subset), "\n")
cat("Karst records found in main data:", nrow(karst_subset), "/", length(karst_record_numbers), 
    sprintf("(%.1f%%)\n", nrow(karst_subset)/length(karst_record_numbers)*100))

nonkarst_subset <- main_data %>%
  dplyr::filter(Record_number %in% nonkarst_record_numbers)

cat("Non-karst subset dimensions:", dim(nonkarst_subset), "\n")
cat("Non-karst records found in main data:", nrow(nonkarst_subset), "/", length(nonkarst_record_numbers), 
    sprintf("(%.1f%%)\n", nrow(nonkarst_subset)/length(nonkarst_record_numbers)*100))

all_record_numbers <- unique(main_data$Record_number)
classified_records <- union(karst_record_numbers, nonkarst_record_numbers)
unclassified_records <- setdiff(all_record_numbers, classified_records)

cat("Unclassified records (in main data but not in karst/non-karst):", length(unclassified_records), "\n")

cat("\n=== STEP 4: ADDING IDENTIFIER COLUMNS ===\n")

main_data_with_karst <- main_data %>%
  mutate(
    Karst_type = case_when(
      Record_number %in% karst_record_numbers ~ "Karst",
      Record_number %in% nonkarst_record_numbers ~ "Non_karst",
      TRUE ~ "Unclassified"
    )
  )

karst_distribution <- main_data_with_karst %>%
  count(Karst_type)

cat("Karst type distribution in main data:\n")
print(karst_distribution)

cat("\n=== STEP 5: SAVING SUBSETS ===\n")

write.csv(karst_subset, "complete_data_rhra_karst.csv", row.names = FALSE)

write.csv(nonkarst_subset, "complete_data_rhra_nonkarst.csv", row.names = FALSE)

cat("\n=== STEP 6: SUMMARY STATISTICS ===\n")

print_summary_stats <- function(data, dataset_name) {
  cat("\n📊", dataset_name, "Summary:\n")
  cat("  Total records:", nrow(data), "\n")
  cat("  Records with Rh_annual:", sum(!is.na(data$Rh_annual)), "\n")
  cat("  Records with Ra_annual:", sum(!is.na(data$Ra_annual)), "\n")
  cat("  Year range:", min(data$Study_midyear, na.rm = TRUE), "-", max(data$Study_midyear, na.rm = TRUE), "\n")
  
  if ("Rh_annual" %in% names(data)) {
    cat("  Rh_annual mean:", mean(data$Rh_annual, na.rm = TRUE), "\n")
    cat("  Rh_annual sd:", sd(data$Rh_annual, na.rm = TRUE), "\n")
  }
  
  if ("Ra_annual" %in% names(data)) {
    cat("  Ra_annual mean:", mean(data$Ra_annual, na.rm = TRUE), "\n")
    cat("  Ra_annual sd:", sd(data$Ra_annual, na.rm = TRUE), "\n")
  }
}

print_summary_stats(karst_subset, "KARST DATA")
print_summary_stats(nonkarst_subset, "NON-KARST DATA")
print_summary_stats(main_data, "ORIGINAL DATA")

cat("\n=== STEP 7: DATA QUALITY CHECK ===\n")

common_vars <- c("Study_midyear", "Latitude", "Biome", "MAT", "MAP", 
                 "Meas_method", "Ecosystem_type", "Rh_annual", "Ra_annual")

cat("Variable completeness in karst data:\n")
for (var in common_vars) {
  if (var %in% names(karst_subset)) {
    non_na <- sum(!is.na(karst_subset[[var]]))
    pct <- round(non_na/nrow(karst_subset)*100, 1)
    cat(sprintf("  %-20s: %d/%d (%5.1f%%)\n", var, non_na, nrow(karst_subset), pct))
  }
}

cat("\nVariable completeness in non-karst data:\n")
for (var in common_vars) {
  if (var %in% names(nonkarst_subset)) {
    non_na <- sum(!is.na(nonkarst_subset[[var]]))
    pct <- round(non_na/nrow(nonkarst_subset)*100, 1)
    cat(sprintf("  %-20s: %d/%d (%5.1f%%)\n", var, non_na, nrow(nonkarst_subset), pct))
  }
}

cat("\n=== STEP 8: COMPARISON SUMMARY ===\n")
comparison_summary <- data.frame(
  Metric = c("Total Records", "Records with Rh", "Records with Ra", 
             "Mean Rh", "Mean Ra", "Year Range"),
  Karst = c(
    nrow(karst_subset),
    sum(!is.na(karst_subset$Rh_annual)),
    sum(!is.na(karst_subset$Ra_annual)),
    round(mean(karst_subset$Rh_annual, na.rm = TRUE), 2),
    round(mean(karst_subset$Ra_annual, na.rm = TRUE), 2),
    paste(round(min(karst_subset$Study_midyear, na.rm = TRUE), 1), 
          "-", round(max(karst_subset$Study_midyear, na.rm = TRUE), 1))
  ),
  NonKarst = c(
    nrow(nonkarst_subset),
    sum(!is.na(nonkarst_subset$Rh_annual)),
    sum(!is.na(nonkarst_subset$Ra_annual)),
    round(mean(nonkarst_subset$Rh_annual, na.rm = TRUE), 2),
    round(mean(nonkarst_subset$Ra_annual, na.rm = TRUE), 2),
    paste(round(min(nonkarst_subset$Study_midyear, na.rm = TRUE), 1), 
          "-", round(max(nonkarst_subset$Study_midyear, na.rm = TRUE), 1))
  )
)

print(comparison_summary)

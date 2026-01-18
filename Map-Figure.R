library(ggplot2)
library(dplyr)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

karst_data <- read.csv("03_karst_with_koppen.csv", stringsAsFactors = FALSE)
nonkarst_data <- read.csv("03_nonkarst_with_koppen.csv", stringsAsFactors = FALSE)

karst_data <- karst_data %>%
  mutate(
    geology = "Karst",
    Latitude = as.numeric(Latitude),
    Longitude = as.numeric(Longitude)
  )

nonkarst_data <- nonkarst_data %>%
  mutate(
    geology = "Non-Karst",
    Latitude = as.numeric(Latitude),
    Longitude = as.numeric(Longitude)
  )

combined_data <- bind_rows(karst_data, nonkarst_data) %>%
  filter(!is.na(Latitude), !is.na(Longitude), !is.na(kg_short))

combined_data <- combined_data %>%
  mutate(
    climate_zone = case_when(
      substr(kg_short, 1, 1) == "A" ~ "Tropical",
      substr(kg_short, 1, 1) == "B" ~ "Arid",
      substr(kg_short, 1, 1) == "C" ~ "Temperate",
      substr(kg_short, 1, 1) == "D" ~ "Cold",
      substr(kg_short, 1, 1) == "E" ~ "Polar",
      TRUE ~ NA_character_
    )
  )

combined_data <- combined_data %>%
  mutate(geology = factor(geology, levels = c("Karst", "Non-Karst")))

koppen_raster <- rast("/mnt/d/Review-Chk/1991_2020/koppen_geiger_0p1.tif")

koppen_df <- as.data.frame(koppen_raster, xy = TRUE, na.rm = TRUE)
names(koppen_df) <- c("x", "y", "kg_code")

koppen_df <- koppen_df %>%
  mutate(
    climate_zone = case_when(
      kg_code %in% 1:3   ~ "Tropical",
      kg_code %in% 4:7   ~ "Arid",
      kg_code %in% 8:16  ~ "Temperate",
      kg_code %in% 17:28 ~ "Cold",
      kg_code %in% 29:30 ~ "Polar",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(climate_zone))

climate_colors <- c(
  "Tropical"   = "#66C2A5",
  "Arid"       = "#FC8D62", 
  "Temperate"  = "#8DA0CB",
  "Cold"       = "#E78AC3",
  "Polar"      = "#A6D854"
)

geology_colors <- c(
  "Karst" = "#E41A1C",
  "Non-Karst" = "#377EB8"
)

geology_shapes <- c(
  "Karst" = 17,
  "Non-Karst" = 19
)

world <- ne_countries(scale = "medium", returnclass = "sf")

final_map <- ggplot() +
  
  geom_tile(
    data = koppen_df,
    aes(x = x, y = y, fill = climate_zone),
    alpha = 0.7
  ) +
  
  geom_sf(
    data = world,
    fill = NA,
    color = "gray30",
    linewidth = 0.25,
    alpha = 0.9
  ) +
  
  geom_point(
    data = filter(combined_data, geology == "Non-Karst"),
    aes(
      x = Longitude,
      y = Latitude,
      color = geology,
      shape = geology
    ),
    size = 2.2,
    alpha = 0.9,
    stroke = 0.5
  ) +
  
  geom_point(
    data = filter(combined_data, geology == "Karst"),
    aes(
      x = Longitude,
      y = Latitude,
      color = geology,
      shape = geology
    ),
    size = 2.2,
    alpha = 0.9,
    stroke = 0.5
  ) +
  
  scale_fill_manual(
    values = climate_colors,
    name = "Climate Zone",
    guide = guide_legend(
      order = 1,
      override.aes = list(
        alpha = 0.8, 
        size = 3, 
        shape = NA,
        color = NA
      )
    )
  ) +
  
  scale_color_manual(
    values = geology_colors,
    name = "Geology",
    guide = guide_legend(
      order = 2,
      override.aes = list(
        size = 3, 
        alpha = 1, 
        shape = c(17, 19)
      )
    )
  ) +
  
  scale_shape_manual(
    values = geology_shapes,
    guide = "none"
  ) +
  
  coord_sf(
    xlim = c(-180, 180),
    ylim = c(-60, 85),
    expand = FALSE,
    crs = st_crs(4326)
  ) +
  
  annotation_scale(
    location = "bl",
    width_hint = 0.2,
    bar_cols = c("black", "white"),
    text_cex = 0.8,
    line_width = 1,
    height = unit(0.15, "cm")
  ) +
  
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    pad_x = unit(0.2, "in"),
    pad_y = unit(0.2, "in"),
    style = north_arrow_fancy_orienteering(
      fill = c("black", "white"),
      line_col = "black",
      text_size = 8
    ),
    height = unit(0.8, "cm"),
    width = unit(0.8, "cm")
  ) +
  
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    
    legend.position = c(0.02, 0.08),
    legend.justification = c(0, 0),
    
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.y = unit(0.15, "cm"),
    
    legend.title = element_blank(),
    legend.text = element_text(
      size = 9,
      margin = margin(r = 8)
    )
  )

ggsave(
  "figure1_clean_map.tiff",
  plot = final_map,
  width = 10,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

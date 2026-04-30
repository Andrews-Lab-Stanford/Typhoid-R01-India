#### TYPHOID R01: PHAGE ABUNDANCE EXPLORATION (39 SITES) ####
# Purpose: Visualize phage abundance (plaques) by time and space for independent catchments.
rm(list = ls())

library(tidyverse)
library(readxl)
library(lubridate)
library(sf)
library(viridis)
library(ggspatial)

# 1. CONFIGURE PATHS
data_dir <- "Data"
manuscript_dir <- "R01-Typhoid-Phage-Wastewater-Manuscript/Analysis Data"

# 2. LOAD SHARED DATA (FROM SAP ANALYSIS)
data_rdata <- file.path("Data", "cleaned_analysis_data.RData")
if (!file.exists(data_rdata)) {
  stop("Cleaned data not found. Please run TyphoidR01_SAP_Analysis.R first.")
}
load(data_rdata) # Loads 'res_intermediate'

# 3. SPATIAL EXPLORATION: PHAGE ABUNDANCE BY SPACE
# Get coordinates for the 38 sites
sites_coords <- read.csv(file.path("Data", "Site_54.csv")) %>%
  mutate(Site = as.numeric(Site)) %>%
  filter(!is.na(long_), !is.na(lat))

spatial_df <- res_intermediate %>%
  group_by(Site) %>%
  summarise(
    avg_abundance = mean(total_plaques, na.rm = TRUE),
    max_abundance = max(total_plaques, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  inner_join(sites_coords, by = "Site") %>%
  st_as_sf(coords = c("long_", "lat"), crs = 4326)

# 4. TEMPORAL EXPLORATION: PHAGE ABUNDANCE BY TIME
temporal_df <- res_intermediate %>%
  group_by(month) %>%
  summarise(
    mean_plaques = mean(total_plaques, na.rm = TRUE),
    total_plaques_global = sum(total_plaques, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(month))

p_time <- ggplot(temporal_df, aes(x = month, y = mean_plaques)) +
  geom_col(fill = "darkred", alpha = 0.6) +
  theme_minimal() +
  labs(title = "Temporal Trend of Phage Abundance (38 Independent Sites)",
       subtitle = "Mean plaques per site-month",
       x = "Month-Year",
       y = "Mean Plaques (Direct Counts)") +
  theme(plot.title = element_text(face = "bold", size = 14))
p_time
ggsave("Phage_Abundance_Time.png", p_time, width = 10, height = 8)

# Configure map limits
bbox <- st_bbox(spatial_df)
xlim_val <- as.numeric(c(bbox[1] - 0.01, bbox[3] + 0.01))
ylim_val <- as.numeric(c(bbox[2] - 0.01, bbox[4] + 0.01))

p_space <- ggplot() +
  annotation_map_tile(type = "osm", zoom = 14) +
  annotate("rect", xmin = xlim_val[1], xmax = xlim_val[2], ymin = ylim_val[1], ymax = ylim_val[2], 
           fill = "white", alpha = 0.5) +
  geom_sf(data = spatial_df, aes(fill = avg_abundance), 
          shape = 21, color = "black", stroke = 0.5, alpha = 0.8, size = 5) +
  scale_fill_viridis_c(option = "magma", name = "Avg Plaques") +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tr", which_north = "true") +
  coord_sf(xlim = xlim_val, ylim = ylim_val, crs = 4326) +
  theme_bw() +
  labs(title = "Spatial Distribution of Phage Abundance (38 Independent Sites)",
       subtitle = "Average plaques per site across study period",
       x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(face = "bold", size = 14))
p_space 

# 5. SAVE ARTIFACTS
ggsave("Phage_Abundance_Space.png", p_space, width = 10, height = 8)
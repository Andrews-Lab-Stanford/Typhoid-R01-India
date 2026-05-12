#### TYPHOID R01: BIVARIATE SPATIAL CORRELATION ANALYSIS ####
# Purpose: Quantify geographic association between wastewater phage and clinical typhoid.
# Methodology: Moran's I (Univariate Clustering) and Lee's L (Bivariate Correlation).

library(tidyverse)
library(readxl)
library(sf)
library(gstat)
library(sp)
library(patchwork)
library(viridis)
library(ggspatial) 
library(spdep) # For Moran's I and Lee's L
library(prettymapr) 

# 2. LOAD DATA
# 2.1 Site Metadata (Coordinates)
data_dir <- "Data"
sites <- read.csv(file.path(data_dir, "Site_54.csv")) %>%
  mutate(Site = as.numeric(Site)) %>%
  filter(!is.na(long_), !is.na(lat))

# 2.2 LOAD SHARED DATA (FROM SAP ANALYSIS)
data_rdata <- file.path("Data", "cleaned_analysis_data.RData")
if (!file.exists(data_rdata)) {
  stop("Cleaned data not found. Please run TyphoidR01_SAP_Analysis.R first.")
}
load(data_rdata) # Loads 'res_intermediate'

# Define Periods accurately
res_p <- res_intermediate %>%
  mutate(
    Period = case_when(
      month >= as.Date("2024-01-01") & month <= as.Date("2025-01-31") ~ "Jan 2024 - Jan 2025",
      month >= as.Date("2025-02-01") & month <= as.Date("2026-01-31") ~ "Feb 2025 - Jan 2026",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Period))

# 3. AGGREGATE STATS BY PERIOD
stats_agg_period <- res_p %>%
  group_by(Period, Site) %>%
  summarise(
    positivity_rate = sum(positive_samples, na.rm=T) / sum(total_samples, na.rm=T),
    total_cases = sum(clinical_cases, na.rm=T),
    Pop = max(Pop, na.rm=T),
    .groups = "drop"
  ) %>%
  mutate(incidence_1000 = (total_cases / Pop) * 1000)

stats_agg_total <- res_p %>%
  group_by(Site) %>%
  summarise(
    positivity_rate = sum(positive_samples, na.rm=T) / sum(total_samples, na.rm=T),
    total_cases = sum(clinical_cases, na.rm=T),
    Pop = max(Pop, na.rm=T),
    Period = "Overall (2024-2026)",
    .groups = "drop"
  ) %>%
  mutate(incidence_1000 = (total_cases / Pop) * 1000)

spatial_combined <- bind_rows(stats_agg_period, stats_agg_total)

# 4. SPATIAL DATA MERGE
spatial_sf_combined <- sites %>%
  inner_join(spatial_combined, by = "Site") %>%
  mutate(Period = factor(Period, levels = c("Jan 2024 - Jan 2025", "Feb 2025 - Jan 2026", "Overall (2024-2026)"))) %>%
  st_as_sf(coords = c("long_", "lat"), crs = 4326) %>%
  st_transform(32644) 

# 4. UNIVARIATE & BIVARIATE STATS
periods <- sort(unique(spatial_sf_combined$Period))

spatial_stats <- map_df(periods, function(p) {
  df_p <- spatial_sf_combined %>% filter(Period == p)
  coords_p <- st_coordinates(df_p)
  
  # Handle identical points
  if (any(duplicated(coords_p))) coords_p <- jitter(coords_p, amount = 0.0001)
  
  # Spatial weights list
  nb <- knn2nb(knearneigh(coords_p, k = 1))
  listw <- nb2listw(nb, style = "W")
  
  # Moran's I (Phage clustering)
  m_phage <- moran.mc(as.numeric(df_p$positivity_rate), listw, nsim = 999)
  
  # Moran's I (Clinical clustering)
  m_clin <- moran.mc(as.numeric(df_p$incidence_1000), listw, nsim = 999)
  
  # Lee's L (Bivariate correlation)
  l_biv <- lee.mc(as.numeric(df_p$positivity_rate), as.numeric(df_p$incidence_1000), listw, nsim = 999)
  
  data.frame(
    Period = p,
    Moran_Phage_I = m_phage$statistic,
    Moran_Phage_p = m_phage$p.value,
    Moran_Clin_I = m_clin$statistic,
    Moran_Clin_p = m_clin$p.value,
    Lee_L = l_biv$statistic,
    Lee_p = l_biv$p.value
  )
})

cat("\n--- SPATIAL AUTO-CORRELATION & BIVARIATE RESULTS ---\n")
print(spatial_stats)

### 5. MAPPING
# Refined limits for cleaner publication look
sites_utm <- st_transform(spatial_sf_combined, 32644)
bbox_utm <- st_bbox(st_buffer(sites_utm, 1500))
bbox_wgs84 <- st_bbox(st_transform(st_as_sfc(bbox_utm), 4326))

# Use numeric indices for robustness
xlim_val <- as.numeric(c(bbox_wgs84[1], bbox_wgs84[3]))
ylim_val <- as.numeric(c(bbox_wgs84[2], bbox_wgs84[4]))
XMIN <- xlim_val[1]
XMAX <- xlim_val[2]
YMIN <- ylim_val[1]
YMAX <- ylim_val[2]

# Mapping Helper for Univariate (Publication Style - Muted Basemap)
create_univariate_map <- function(site_data, var_name, title, legend_name, stat_df, stat_prefix, palette = "mako", use_size = FALSE) {
  # Prepare stats label
  stat_col <- paste0("Moran_", stat_prefix, "_I")
  p_col    <- paste0("Moran_", stat_prefix, "_p")
  
  # Create label string in the stat_df and ENSURE factor levels match site_data
  facet_labs <- stat_df %>% 
    mutate(Period = factor(Period, levels = levels(site_data$Period))) %>%
    mutate(lab = sprintf("Moran's I = %.2f (p = %.3f)", !!sym(stat_col), !!sym(p_col)))
  
  p <- ggplot() +
    annotation_map_tile(type = "osm", zoom = 13) +
    # Wash out background with a white rectangle
    annotate("rect", xmin = XMIN, xmax = XMAX, ymin = YMIN, ymax = YMAX, 
             fill = "white", alpha = 0.5)
  
  if (use_size) {
    p <- p + geom_sf(data = st_transform(site_data, 4326), 
                     aes(size = !!sym(var_name)), 
                     shape = 21, fill = "#BD0026", color = "black", stroke = 0.5, alpha = 0.7) +
      scale_size_continuous(name = legend_name, breaks = seq(0, 12, 4), range = c(1, 10))
  } else {
    p <- p + geom_sf(data = st_transform(site_data, 4326), 
                     aes(fill = !!sym(var_name)), 
                     shape = 21, size = 4, color = "black", stroke = 0.5, alpha = 0.8) +
      scale_fill_viridis_c(option = palette, name = legend_name, labels = scales::percent_format(accuracy = 1))
  }
  
  p + 
    geom_text(data = facet_labs, aes(x = XMIN + 0.005, y = YMAX - 0.005, label = lab), 
              size = 5, fontface = "bold", inherit.aes = FALSE, hjust = 0, vjust = 1, color = "black") +
    facet_wrap(~Period) +
    annotation_scale(location = "bl", width_hint = 0.4) +
    annotation_north_arrow(location = "tr", which_north = "true", 
                           style = north_arrow_fancy_orienteering) +
    coord_sf(xlim = xlim_val, ylim = ylim_val, expand = FALSE, crs = 4326) +
    theme_bw() + 
    labs(title = title, x = "Longitude", y = "Latitude") + 
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(face = "bold", size = 12))
}

# Plan: 3 maps
p_phage <- create_univariate_map(spatial_sf_combined, "positivity_rate", "Wastewater Phage Positivity (%)", "Phage %", 
                                 spatial_stats, "Phage", "mako")
p_clin  <- create_univariate_map(spatial_sf_combined, "incidence_1000", "Clinical Typhoid Incidence (per 1,000)", "Cases/1,000", 
                                 spatial_stats, "Clin", "rocket", use_size = TRUE)

# Combined Bivariate Overlay (Publication Style - Muted)
facet_labs <- spatial_stats %>% 
  mutate(Period = factor(Period, levels = levels(spatial_sf_combined$Period))) %>%
  mutate(lab = sprintf("Lee's L = %.2f (p = %.3f)", Lee_L, Lee_p))
p_biv <- ggplot() +
  annotation_map_tile(type = "osm", zoom = 13) +
  annotate("rect", xmin = XMIN, xmax = XMAX, ymin = YMIN, ymax = YMAX, 
           fill = "white", alpha = 0.5) +
  geom_sf(data = st_transform(spatial_sf_combined, 4326), 
          aes(fill = positivity_rate, size = incidence_1000), 
          shape = 21, color = "black", stroke = 0.3, alpha = 0.8) +
  scale_fill_viridis_c(option = "turbo", name = "Phage Positivity (%)", labels = scales::percent_format(accuracy = 1)) +
  scale_size_continuous(name = "Clinical Incidence\n(per 1,000 population)", breaks = seq(0, 12, 4), range = c(1, 10)) +
  geom_text(data = facet_labs, aes(x = XMIN + 0.002, y = YMAX - 0.005, label = lab), 
            size = 4.5, fontface = "bold", inherit.aes = FALSE, hjust = 0, vjust = 1, color = "black") +
  facet_wrap(~Period) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = xlim_val, ylim = ylim_val, expand = FALSE, crs = 4326) +
  theme_bw() + 
  labs(title = "Bivariate Correlation: Wastewater Phage Positivity vs. Clinical Typhoid Incidence",
       subtitle = "Spatial overlap of environmental detection (% of samples that are positive) and disease incidence (cases per 1,000)",
       x = "Longitude", y = "Latitude") + 
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold", size = 12))

# SAVE
suppressWarnings({
  ggsave("Typhoid_Spatial_Univariate_Phage.png", p_phage, width = 12, height = 7, dpi = 300)
  ggsave("Typhoid_Spatial_Univariate_Clinical.png", p_clin, width = 12, height = 7, dpi = 300)
  ggsave("Typhoid_Spatial_Bivariate_Combined.png", p_biv, width = 12, height = 7, dpi = 300)
  write.csv(spatial_stats, "Spatial_Correlation_Stats.csv", row.names = FALSE)
})

# 6. FINAL PLOTTING (SHOW ON CONSOLE)
p_phage
p_clin
p_biv


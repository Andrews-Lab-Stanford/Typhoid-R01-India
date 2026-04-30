#### TYPHOID R01 PHAGE ASSAYS ####
# User: Esther Jung
# Last Updated: February 19, 2026

rm(list = ls())

#### CONFIGURE ####
# Core data manipulation and visualization
library(tidyverse) # Includes dplyr, ggplot2, readr, etc.
library(readxl) # Excel file reading
library(lubridate) # Date handling
library(tidyselect) # Column selection helpers

# Visualization and plotting
library(ggplot2)
library(ggpubr) # Publication-ready plots
library(patchwork) # Plot arrangement
library(gridExtra)
library(grid)

# Mapping
library(leaflet)
library(leafsync) # Synced leaflet maps
library(htmlwidgets)
library(htmltools)

# Statistical analysis
library(gmodels)
library(MASS) # Negative binomial models
library(pscl) # Zero-inflated models

# Other utilities
library(here) # Project paths
library(gt) # Tables
library(REDCapR) # REDCap data access
library(devtools)
library(anytime)
library(shiny)
library(DiagrammeR)
library(rsconnect)

setwd("~/Google Drive/My Drive/Typhoid R01/Data")
sites <- read.csv("Site_54.csv") # Maps the catchments out
phage <- read.csv("REDCapPhageData.csv") # REDCap phage assay data until January 2026
pcr <- read.csv("typhi_sitePos_20251231.csv") # PCR by site until 12/31/2025
monthly_cases <- read_excel("Site_monthlytyphoidcases_FINAL.xlsx") # Clinical cases by site and month until January 2026

#### DATA CLEANING ============================================================

# Create case_site from monthly_cases: aggregate total cases by site
case_site <- monthly_cases %>%
  mutate(
    # Sum all monthly case columns
    total_cases = rowSums(across(matches("^(2024|2025|2026)[_-]")), na.rm = TRUE),
    case_positivity_percent = round(100 * total_cases / Pop, 2)
  ) %>%
  dplyr::select(Site, Pop, total_cases, case_positivity_percent) %>%
  arrange(desc(case_positivity_percent))

# Create monthly_cases_long for time series analysis
monthly_cases_long <- monthly_cases %>%
  pivot_longer(
    cols = matches("^(2024|2025|2026)[_-]"), # All month columns
    names_to = "month",
    values_to = "cases"
  ) %>%
  mutate(
    month = ym(month), # Convert to date for proper ordering
    positivity_percent = round((cases / Pop) * 100, 2)
  )

# Sites with cases
case_y_site <- case_site %>% filter(total_cases > 0)

names(pcr)

pcr <- pcr %>%
  rename(Site = site, Positive = MS, Unknown = NA.) %>%
  dplyr::select(c("Site", "Positive", "Unknown")) # PCR positivity

phage <- phage %>%
  mutate(result_status = case_when(
    amplification_tw == 1 | direct_tw == 1 ~ "Positive",
    amplification_tw == 2 & direct_tw == 2 ~ "Negative",
    (amplification_tw == 2 & is.na(direct_tw)) | (direct_tw == 2 & is.na(amplification_tw)) ~ "Negative",
    TRUE ~ "Unknown"
  ))
table(phage$result_status)

freq_table_df <- phage %>% # Phage positivity
  rename(Site = site_id_tw) %>%
  dplyr::select(Site, result_status) %>%
  group_by(Site) %>%
  summarise(
    Positive = sum(result_status == "Positive"),
    Total = n(),
    `Positivity %` = round(100 * (Positive / Total), 2),
    .groups = "drop"
  )

#### TABLES AND FIGURES ####
# Amplification = Enrichment assay | 1 = +, 2 = -
# Direct = Direct plating assay | 1 = +, 2 = -

min(phage$date_received_tw) # First sample received 2024-01-07, 2 year collection
max(phage$date_results_tw)

# Table 1
length(phage$record_id)
table(phage$amplification_tw, useNA = "ifany") # enrichment = qualitative
table(phage$direct_tw, useNA = "ifany") # direct assay
table(phage$direct_counts_tw) # direct abundance counts

sum(phage$amplification_tw == 1 & phage$direct_tw == 1, na.rm = TRUE) # both positive
sum(phage$amplification_tw == 2 & phage$direct_tw == 2, na.rm = TRUE) # both negative
sum(is.na(phage$amplification_tw) & is.na(phage$direct_tw)) # both NA

sum(phage$amplification_tw == 1 | phage$direct_tw == 1, na.rm = TRUE) # positive on either
both_negative <- phage$amplification_tw == 2 & phage$direct_tw == 2
one_negative_one_missing <- (phage$amplification_tw == 2 & is.na(phage$direct_tw)) |
  (phage$direct_tw == 2 & is.na(phage$amplification_tw))
sum(both_negative | one_negative_one_missing, na.rm = TRUE) # total number of negatives

# Table 2
addmargins(table(phage$amplification_tw, phage$direct_tw, useNA = "always"))

# Phage positivity by month - clean date and status
phage_month <- phage %>%
  mutate(
    date_sample = as.Date(ifelse(
      substr(date_results_tw, 1, 2) %in% c("29", "22") | substr(date_results_tw, 1, 4) == "2004", # key to correct years
      paste0("2024", substr(date_results_tw, 5, 10)),
      date_results_tw
    ), format = "%Y-%m-%d"),
    month = floor_date(date_sample, "month"),
    status = case_when(
      amplification_tw == 1 | direct_tw == 1 ~ "Positive",
      amplification_tw == 2 & direct_tw == 2 ~ "Negative",
      (amplification_tw == 2 & is.na(direct_tw)) |
        (is.na(amplification_tw) & direct_tw == 2) ~ "Negative",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(status), !is.na(site_id_tw))

# Table 3: Monthly summary of phage positives
phage_summary <- phage_month %>%
  group_by(month_year = month) %>%
  summarise(
    `Direct assay` = sum(direct_tw == 1, na.rm = TRUE),
    `Enrichment assay` = sum(amplification_tw == 1, na.rm = TRUE),
    Total = n(),
    .groups = "drop"
  ) %>%
  complete(
    month_year = seq.Date(
      from = min(phage_month$month, na.rm = TRUE),
      to = max(phage_month$month, na.rm = TRUE),
      by = "month"
    ),
    fill = list(`Direct assay` = 0, `Enrichment assay` = 0, Total = 0)
  )
sum(phage_summary$Total)

# Figure 1: Total Typhi Phage Positives by Month
phage_plot_data <- phage_summary %>%
  mutate(
    month_year = as.Date(month_year),
    total_positive = `Direct assay` + `Enrichment assay`,
    positivity_pct = ifelse(Total > 0, round(100 * total_positive / Total, 1), 0)
  ) %>%
  filter(month_year >= as.Date("2024-03-01") & month_year <= as.Date("2026-01-31"))
# 1a. Absolute Numbers Plot
phage_count_plot <- ggplot(phage_plot_data, aes(x = month_year, y = total_positive)) +
  geom_col(fill = "#7260a7") +
  geom_text(aes(label = total_positive), vjust = -0.4, size = 3.5) +
  scale_x_date(breaks = sort(unique(phage_plot_data$month_year)), date_labels = "%b %Y") +
  theme_pubr(base_size = 14) +
  labs(
    x = "Month",
    y = expression("Number of samples positive for" ~ italic("S.") ~ "Typhi phage")
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Monthly_Phage_Positives.png", phage_count_plot, width = 10, height = 4, dpi = 300)
#### MAPS BY CATCHMENTS ========================================================

# Prepare case map data with positivity calculations
case_map_data <- sites %>%
  left_join(case_site %>% dplyr::select(Site, total_cases, Pop, Positivity = case_positivity_percent),
    by = "Site"
  )
case_pal <- colorNumeric("YlOrRd", case_map_data$Positivity)

case_map <- leaflet(case_map_data) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addProviderTiles("Esri.WorldStreetMap") %>%
  addCircleMarkers(
    ~long_, ~lat,
    color = "black",
    fillColor = ~ case_pal(Positivity),
    fillOpacity = 0.9,
    radius = ~ sqrt(Pop) * 0.09, # Adjust as needed for sensible radii
    stroke = TRUE, weight = 1, opacity = 1,
    label = ~ paste0(
      "Site: ", Site,
      " | Cases: ", total_cases,
      " | Population: ", Pop,
      " | Positivity (%): ", Positivity
    ),
    labelOptions(textsize = 15)
  ) %>%
  addLegend(
    "bottomleft",
    pal = case_pal,
    values = case_map_data$Positivity,
    title = "Blood Culture Positivity (# cases/pop*100)"
  )

pcr_map_data <- pcr %>%
  left_join(freq_table_df %>% dplyr::select(Site, Total), by = "Site") %>%
  mutate(Positivity = round(100 * (Positive / Total), 2)) %>%
  left_join(sites, by = "Site")

pcr_pal <- colorNumeric("GnBu", pcr_map_data$Positivity)

pcr_map <- leaflet(pcr_map_data) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addProviderTiles("Esri.WorldStreetMap") %>%
  addCircleMarkers(
    ~long_, ~lat,
    color = "black",
    fillColor = ~ pcr_pal(Positivity),
    fillOpacity = 0.9,
    radius = ~ sqrt(Total) * 2,
    stroke = TRUE, weight = 1, opacity = 1,
    label = ~ paste0(
      "Site: ", Site,
      " | Positive Samples: ", Positive,
      " | Total Samples: ", Total,
      " | Positivity (%): ", Positivity
    ),
    labelOptions(textsize = 15)
  ) %>%
  addLegend(
    "bottomleft",
    pal = pcr_pal,
    values = pcr_map_data$Positivity,
    title = "PCR Positivity (# positive/n*100)"
  )

phage_map_data <- merge(sites, freq_table_df, by = "Site", all.x = TRUE)
phage_pal <- colorNumeric("BuPu", phage_map_data$`Positivity %`)

phage_map <- leaflet(phage_map_data) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addProviderTiles("Esri.WorldStreetMap") %>%
  addCircleMarkers(
    ~long_, ~lat,
    color = "black",
    fillColor = ~ phage_pal(`Positivity %`),
    fillOpacity = 0.9,
    radius = ~ sqrt(Total) * 2,
    stroke = TRUE, weight = 1, opacity = 1,
    label = ~ paste0(
      "Site: ", Site,
      " | Positive Samples: ", Positive,
      " | Total Samples: ", Total,
      " | Positivity (%): ", `Positivity %`
    ),
    labelOptions(textsize = 15)
  ) %>%
  addLegend(
    "bottomleft",
    pal = phage_pal,
    values = phage_map_data$`Positivity %`,
    title = "Phage Positivity (# positive/n*100)"
  )

# Synced maps object
synced_maps <- sync(phage_map, case_map)

# Save using save_html instead of saveWidget
saveWidget(case_map, "case_map.html", selfcontained = TRUE)
saveWidget(pcr_map, "pcr_map.html", selfcontained = TRUE)
saveWidget(phage_map, "phage_map.html", selfcontained = TRUE)
save_html(synced_maps, file = "Typhoid_Positivity_Maps.html")

#### SCATTERPLOTS =====================================================

# Improved Scatterplot Function
create_scatterplot <- function(data, x_var, x_label, color, title) {
  # Sort Site by x_var (positivity) for better visual trend
  data_sorted <- data %>%
    mutate(Site = fct_reorder(factor(Site), {{ x_var }}))

  ggplot(data_sorted, aes(x = {{ x_var }}, y = Site)) +
    geom_point(color = color, size = 3, alpha = 0.8) +
    # Add a subtle light grey line to help eye follow
    geom_segment(aes(x = 0, xend = {{ x_var }}, y = Site, yend = Site),
      color = "grey90", linewidth = 0.2
    ) +
    theme_pubr(base_size = 12) +
    labs(title = title, x = x_label, y = "Catchment Site") +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 11),
      axis.text.y = element_text(size = 9),
      panel.grid.major.x = element_line(color = "grey95")
    )
}

# 1. Blood Culture Scatterplot
blood_culture_plot <- create_scatterplot(
  case_site, case_positivity_percent,
  "Incidence (cases per 100 population)", "#0077b6",
  "Blood Culture Positivity"
)

# 2. PCR Scatterplot
pcr_plot <- create_scatterplot(
  pcr_map_data, Positivity,
  "Triple Positive Samples (%)", "#386641",
  "PCR Positivity"
)

# 3. Phage Scatterplot
phage_plot <- create_scatterplot(
  phage_map_data, `Positivity %`,
  "Samples Positive by Enrichment or Direct (%)", "#7209b7",
  "Phage Positivity"
)

combined_plots <- blood_culture_plot + pcr_plot + phage_plot +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Comparison of Typhoid Positivity Proportions by Catchment Site",
    subtitle = "Sites ordered by positivity within each assay type",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5)
    )
  )

combined_plots

ggsave("Positivity_Scatterplots.png", combined_plots, width = 15, height = 10)

# Merge datasets into one by Site for correlation plots
plot_df <- case_site %>%
  dplyr::select(Site, blood = case_positivity_percent) %>%
  left_join(pcr_map_data %>% dplyr::select(Site, pcr = Positivity), by = "Site") %>%
  left_join(phage_map_data %>% dplyr::select(Site, phage = `Positivity %`), by = "Site")

corr_scatter <- function(data, x, y, xlab, ylab, color = "black") {
  ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(size = 3, alpha = 0.6, color = color) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1, color = "darkred", fill = "pink", alpha = 0.2) +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 5) +
    theme_pubr(base_size = 14) +
    labs(x = xlab, y = ylab)
}

p1 <- corr_scatter(
  plot_df, "blood", "pcr",
  "Blood Culture Incidence (%)",
  "PCR Site Positivity (%)",
  "#0077b6"
)

p2 <- corr_scatter(
  plot_df, "blood", "phage",
  "Blood Culture Incidence (%)",
  "Phage Site Positivity (%)",
  "#0077b6"
)

p3 <- corr_scatter(
  plot_df, "pcr", "phage",
  "PCR Site Positivity (%)",
  "Phage Site Positivity (%)",
  "#386641"
)

# Arrange correlation plots
combined_correlations <- (p1 | p2) /
  (p3 | plot_spacer()) +
  plot_annotation(
    title = "Correlations between Environmental and Clinical Positivity",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

combined_correlations
ggsave("Positivity_Correlations.png", combined_correlations, width = 12, height = 10)


#### BARPLOTS =====================================================
blood_culture_data <- case_site %>%
  mutate(Binned = cut(case_positivity_percent, breaks = c(-Inf, 0, 0.2, 0.5, 1, 10)))

pcr_map_data <- pcr_map_data %>%
  mutate(Binned = cut(Positivity, breaks = c(-Inf, 0, 5, 10, 15, 20, 40)))

phage_map_data <- phage_map_data %>%
  mutate(Binned = cut(`Positivity %`, breaks = c(-Inf, 0, 5, 10, 25, 50, 85)))

# Improved Binned Barplot function
create_binned_barplot <- function(data, binned_var, title, fill_color) {
  ggplot(data, aes(x = {{ binned_var }})) +
    geom_bar(fill = fill_color, alpha = 0.8) +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 4) +
    labs(title = title, x = "Positivity Category", y = "Number of Sites") +
    theme_pubr(base_size = 12) +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# 1. Blood Culture Barplot
blood_culture_binned_plot <- create_binned_barplot(
  blood_culture_data, Binned,
  "Blood Culture Catchment Incidence", "#0077b6"
)

# 2. PCR Barplot
pcr_binned_plot <- create_binned_barplot(
  pcr_map_data, Binned,
  "PCR Site Positivity", "#386641"
)

# 3. Phage Barplot
phage_binned_plot <- create_binned_barplot(
  phage_map_data, Binned,
  "Phage Site Positivity", "#7209b7"
)

binned_plots <- blood_culture_binned_plot + pcr_binned_plot + phage_binned_plot +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Distribution of Positivity Categories across Catchment Sites",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

binned_plots
ggsave("Positivity_Barcharts.png", binned_plots, width = 12, height = 6)

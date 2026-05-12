#### TYPHOID R01 PHAGE ASSAYS ####
# User: Esther Jung                                                         
# Last Updated: April 4, 2025                                                

rm(list = ls())

#### CONFIGURE ####
library(readxl)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(devtools)
library(lubridate)
library(tidyselect)
library(gmodels)
library(anytime)
library(shiny)
library(DiagrammeR)
library(rsconnect)
library(REDCapR)
library(gt)
library(here)
library(leaflet)
library(readr)
library(gridExtra)
library(grid)
library(leafsync)
library(htmlwidgets)
library(htmltools)
library(ggpubr)

setwd("~/Google Drive/My Drive/Typhoid R01/Data")
sites <- read.csv("Site_54.csv")
phage <- read.csv("NovelToolsTyphoidVac_DATA_2025-12-16_0120.csv") # Replace manually
pcr <- read.csv("typhi_positives_site.csv") # PCR by site until 10/08/2025
monthly_cases <- read_excel("Site_monthlytyphoidcases_Sep2025.xlsx") # Clinical cases by site and month until September 2025

#### DATA CLEANING ============================================================

names(pcr)

pcr <- pcr %>% rename(Site = site, Positive = MS, Unknown = NA.) %>% select(c("Site", "Positive", "Unknown")) # PCR positivity
phage <- phage %>%
  mutate(result_status = case_when(
    amplification_tw == 1 | direct_tw == 1 ~ "Positive",  
    amplification_tw == 2 & direct_tw == 2 ~ "Negative",  
    (amplification_tw == 2 & is.na(direct_tw)) | (direct_tw == 2 & is.na(amplification_tw)) ~ "Negative",  
    TRUE ~ "Unknown"
  ))
table(phage$result_status)

freq_table_df <- phage %>% # Phage positivity
  select(c(site_id_tw, result_status)) %>%
  group_by(site_id_tw) %>%
  summarise(
    Positive = sum(result_status == "Positive"),
    Total = n(),
    Positivity = round(100*(Positive / Total),2)
  )
colnames(freq_table_df) <- c("Site", "Positive", "Total", "Positivity %")
write.csv(freq_table_df, "freq_table_df.csv", row.names = FALSE)


#### TABLES AND FIGURES ####
# Amplification = Enrichment assay | 1 = +, 2 = -
# Direct = Direct plating assay | 1 = +, 2 = -

min(phage$date_received_tw) # First sample received 2024-01-07, 2 year collection 

# Table 1
length(phage$record_id)
table(phage$amplification_tw, useNA = "ifany") # enrichment = qualitative
table(phage$direct_tw, useNA = "ifany") # direct = quantitative
table(phage$direct_counts_tw) # direct abundance

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

# Positivity by month
phage_month <- phage %>%
  mutate(
    date_sample = as.Date(ifelse(
      substr(date_results_tw, 1, 2) %in% c("29", "22") | substr(date_results_tw, 1, 4) == "2004",
      paste0("2024", substr(date_results_tw, 5, 10)),
      date_results_tw
    ), format = "%Y-%m-%d"),
    status = case_when(
      amplification_tw == 1 | direct_tw == 1 ~ "Positive",
      (amplification_tw == 2 & direct_tw == 2) |
        (amplification_tw == 2 & is.na(direct_tw)) |
        (is.na(amplification_tw) & direct_tw == 2) ~ "Negative",
      TRUE ~ NA_character_
    )
  ) %>%
  select(c("record_id", "date_sample", "amplification_tw", "direct_tw", "status")) %>%
  filter(!is.na(status))

# Table 3
phage_summary <- phage_month %>%
  mutate(
    month_year = floor_date(date_sample, "month"),
    assay_type = case_when(
      amplification_tw == 1 ~ "Enrichment assay",
      direct_tw == 1 ~ "Direct assay",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(month_year) %>%
  summarise(
    `Direct assay` = sum(direct_tw == 1, na.rm = TRUE),
    `Enrichment assay` = sum(amplification_tw == 1, na.rm = TRUE),
    Total = n()
  ) %>%
  complete(
    month_year = seq.Date(
      from = min(phage_month$date_sample, na.rm = TRUE),
      to = max(phage_month$date_sample, na.rm = TRUE),
      by = "month"
    ),
    fill = list(`Direct assay` = 0, `Enrichment assay` = 0, Total = 0)
  )

print(phage_summary, n = 28)
sum(phage_summary$Total)

# Figure 1
phage_plot_data <- phage_summary %>%
  filter(month_year <= "2025-11-01") %>%
  pivot_longer(cols = c("Direct assay", "Enrichment assay"),
               names_to = "assay_type", values_to = "count") %>%
  mutate(month_year_format = format(month_year, "%b %Y"))

phage_plot <- ggplot(phage_plot_data, aes(x = month_year_format, y = count, fill = assay_type)) +
  geom_bar(stat = "identity", fill = "#484848") +
  facet_wrap(~assay_type, scales = "fixed") +
  labs(x = "Date", y = "Typhi positives") +
  scale_x_discrete(limits = unique(phage_plot_data$month_year_format)) +  # Corrected line
  theme_gray(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, hjust = 0))

phage_plot_data <- phage_summary %>%
  filter(month_year <= as.Date("2025-11-01")) %>%
  mutate(month_year = as.Date(month_year)) %>%
  # Assuming your columns are named exactly "Direct assay" and "Enrichment assay"
  mutate(total_positive = (`Direct assay` + `Enrichment assay`)) %>%
  mutate(month_year_format = format(month_year, "%b %Y"))

# Create the singular barplot
phage_plot <- ggplot(phage_plot_data, aes(x = month_year_format, y = total_positive)) +
  geom_col(fill = "#484848") +
  geom_text(aes(label = total_positive), vjust = -0.4, size = 3) +
  labs(
    title = "Total Typhi Phage Positives by Month (Direct or Enrichment)",
    x = "Month",
    y = "Number of Positive Samples"
  ) +
  scale_x_discrete(limits = unique(phage_plot_data$month_year_format)) +
  theme_gray(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

phage_plot

#### MAPS BY CATCHMENTS ========================================================

cases <- case_site %>%
  mutate(Positivity = round(100 * (total_cases/ Pop), 2))

case_map_data <- merge(sites, cases, by = "Site", all.x = TRUE)
case_pal <- colorNumeric("YlOrRd", case_map_data$Positivity)

case_map <- leaflet(case_map_data) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addProviderTiles("Esri.WorldStreetMap") %>%
  addCircleMarkers(
    ~long_, ~lat,
    color = "black",
    fillColor = ~case_pal(Positivity),
    fillOpacity = 0.9, 
    radius = ~sqrt(Pop)*0.09,  # Adjust as needed for sensible radii
    stroke = TRUE, weight = 1, opacity = 1,
    label = ~paste0("Site: ", Site, 
                    " | Cases: ", total_cases, 
                    " | Population: ", Pop, 
                    " | Positivity (%): ", Positivity),
    labelOptions(textsize = 15)
  ) %>%
  addLegend(
    "bottomleft",
    pal = case_pal,
    values = case_map_data$Positivity,
    title = "Blood Culture Positivity (# cases/pop*100)"
  )

pcr_map_data <- pcr %>%
  left_join(freq_table_df %>% select(Site, Total), by = "Site") %>%
  mutate(Positivity = round(100 * (Positive / Total), 2)) %>%
  left_join(sites, by = "Site")

pcr_pal <- colorNumeric("GnBu", pcr_map_data$Positivity)

pcr_map <- leaflet(pcr_map_data) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addProviderTiles("Esri.WorldStreetMap") %>%
  addCircleMarkers(
    ~long_, ~lat,
    color = "black",
    fillColor = ~pcr_pal(Positivity),
    fillOpacity = 0.9, 
    radius = ~sqrt(Total)*2,
    stroke = TRUE, weight = 1, opacity = 1,
    label = ~paste0("Site: ", Site, 
                    " | Positive Samples: ", Positive, 
                    " | Total Samples: ", Total, 
                    " | Positivity (%): ", Positivity),
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
    fillColor = ~phage_pal(`Positivity %`),  
    fillOpacity = 0.9, 
    radius = ~sqrt(Total)*2,
    stroke = TRUE, weight = 1, opacity = 1,  
    label = ~paste0("Site: ", Site, 
                    " | Positive Samples: ", Positive, 
                    " | Total Samples: ", Total, 
                    " | Positivity (%): ", `Positivity %`),
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

create_scatterplot <- function(data, x_var, x_label, color, title) {
  ggplot(data, aes(x = .data[[deparse(substitute(x_var))]], y = factor(Site))) +
    geom_point(color = color, size = 3) +
    theme_minimal() +
    labs(title = title, x = x_label, y = "Site") +
    scale_y_discrete(labels = levels(factor(data$Site))) +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8))
}

# Blood Culture Scatterplot
blood_culture_plot <- create_scatterplot(case_site, case_positivity_percent, 
                                         "Blood Culture Positivity (cases per 100 population)", "blue", 
                                         "Blood Culture Positivity")

# PCR Scatterplot
pcr_plot <- create_scatterplot(pcr_map_data, Positivity, "PCR Positivity (% of samples triple positive by Moore Swab)", "darkgreen",
                               "PCR Positivity (%)")

# Phage Scatterplot
phage_plot <- create_scatterplot(phage_map_data, `Positivity %`, "Phage Positivity (% of samples positive by enrichment or direct)", "purple3", 
                                 "Phage Positivity (%)")

combined_plots <- blood_culture_plot + pcr_plot + phage_plot +
  plot_layout(ncol = 3) +  # three columns (side by side)
  plot_annotation(title = "Comparison of Positivity by Site")
combined_plots

ggsave("Positivity_Scatterplots.png", combined_plots, width = 15, height = 10)

library(ggpubr)

# Merge datasets into one by Site
plot_df <- case_site %>%
  select(Site, blood = case_positivity_percent) %>%
  left_join(pcr_map_data %>% select(Site, pcr = Positivity), by = "Site") %>%
  left_join(phage_map_data %>% select(Site, phage = `Positivity %`), by = "Site")

corr_scatter <- function(data, x, y, xlab, ylab) {
  ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
    stat_cor(method = "pearson", label.x = 0.05 * max(data[[x]]), label.y = 0.95 * max(data[[y]])) +
    theme_minimal() +
    labs(x = xlab, y = ylab)
}

p1 <- corr_scatter(plot_df, "blood", "pcr",
                   "Blood Culture Positivity (% of population)",
                   "PCR Positivity (% of samples)")

p2 <- corr_scatter(plot_df, "blood", "phage",
                   "Blood Culture Positivity (% of population)",
                   "Phage Positivity (% of samples)")

p3 <- corr_scatter(plot_df, "pcr", "phage",
                   "PCR Positivity (% of samples)",
                   "Phage Positivity (% of samples)")

library(patchwork)
combined <- (p1 | p2) /
  (p3 | plot_spacer())

combined

#### BARPLOTS =====================================================
blood_culture_data <- case_site %>%
  mutate(Binned = cut(case_positivity_percent, breaks = c(-Inf, 0, 0.2, 0.5, 1, 10)))

pcr_map_data <- pcr_map_data %>%
  mutate(Binned = cut(Positivity, breaks = c(-Inf, 0, 5, 10, 15, 20, 40)))

phage_map_data <- phage_map_data %>%
  mutate(Binned = cut(`Positivity %`, breaks = c(-Inf, 0, 5, 10, 25, 50, 85)))

create_binned_barplot <- function(data, binned_var, title, fill_color) {
  ggplot(data, aes(x = {{binned_var}})) +
    geom_bar(fill = fill_color) +
    labs(title = title, x = "Positivity Category", y = "Number of Sites") +
    theme_minimal(base_size = 10)
}
blood_culture_binned_plot <- create_binned_barplot(blood_culture_data, Binned, "Blood Culture Positivity (cases per 100 people)", "blue")
pcr_binned_plot <- create_binned_barplot(pcr_map_data, Binned, "PCR Positivity (% of samples)", "darkgreen")
phage_binned_plot <- create_binned_barplot(phage_map_data, Binned, "Phage Positivity (% of samples)", "purple3")
binned_plots <- blood_culture_binned_plot + pcr_binned_plot + phage_binned_plot + 
  plot_layout(ncol = 3)

ggsave("Positivity_Barcharts.png", binned_plots, width = 12, height = 5)

#### CASE BY CATCHMENT AND PHAGE BY CATCHMENT OVER TIME ========================

# Clean and reshape
monthly_cases_long <- monthly_cases %>%
  pivot_longer(cols = -c(Site, Pop), 
               names_to = "month", 
               values_to = "cases") %>%
  mutate(
    month = ym(month),
    positivity_percent = round((cases / Pop) * 100,2)
  )

case_site <- monthly_cases %>%
  group_by(Site) %>%
  summarise(
    Pop = Pop,
    total_cases = rowSums(across(matches("^(2024|2025)_")), na.rm = TRUE),
    case_positivity_percent = round(100 * total_cases / Pop, 2)) %>%
  arrange(desc(case_positivity_percent))

case_y_site <- case_site %>% filter(total_cases > 0)

monthly_cases_long <- monthly_cases %>%
  pivot_longer(
    cols = matches("^(2024|2025)_"),  # all month columns
    names_to = "month",
    values_to = "cases"
  ) %>%
  mutate(
    month = ym(month)  # convert to date for proper ordering
  )

# 2. Aggregate across sites for each month
cases_by_month <- monthly_cases_long %>%
  group_by(month) %>%
  summarise(
    total_cases = sum(cases, na.rm = TRUE)
  )

  # Plot total cases by month
  ggplot(cases_by_month, aes(x = month, y = total_cases)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = total_cases), vjust = -0.4, size = 3) +
    scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = "Total Culture-Confirmed Typhoid Cases by Month",
      x = "Month",
      y = "Total Cases"
    )

# Clean & aggregate phage sample data with plaque counts
phage_month <- phage %>%
  mutate(
    date_sample = as.Date(ifelse(
      substr(date_results_tw, 1, 2) %in% c("29", "22") | substr(date_results_tw, 1, 4) == "2004",  
      paste0("2024", substr(date_results_tw, 5, 10)),
      date_results_tw
    ), format = "%Y-%m-%d"),
    month = floor_date(date_sample, "month"),
    status = case_when(
      amplification_tw == 1 | direct_tw == 1 ~ "Positive",
      amplification_tw == 2 & direct_tw == 2 ~ "Negative",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(status), !is.na(site_id_tw))
 
# Aggregate to site-month level
phage_by_site_month <- phage_month %>%
  group_by(Site = site_id_tw, month) %>%
  summarise(
    total_samples = n(),
    positive_samples = sum(status == "Positive", na.rm = TRUE),
    total_plaques = sum(direct_counts_tw, na.rm = TRUE), 
    phage_positivity = round(100 * positive_samples / total_samples, 2),
    .groups = "drop"
  ) 

# Aggregate to site level
phage_site <- phage_month %>%
  group_by(Site = site_id_tw) %>%
  summarise(
    total_samples = n(),
    positive_samples = sum(status == "Positive", na.rm = TRUE),
    total_plaques = sum(direct_counts_tw, na.rm = TRUE), 
    phage_positivity = round(100 * positive_samples / total_samples, 2),
    .groups = "drop"
  ) 

# Final dataset
merged_df <- phage_by_site_month %>%
  left_join(monthly_cases_long, by = c("Site", "month")) %>%
  filter(!is.na(cases), !is.na(total_samples)) 
length(unique(merged_df$Site)) # should be 54 sites

#### Median positivity comparison ==============================================
data <- left_join(case_site, phage_site, by = c("Site"))

# Create a grouping variable based on whether there are any cases
data <- data %>%
  mutate(case_group = ifelse(total_cases > 0, "Typhoid cases", "No typhoid cases"))

data <- data %>%
  mutate(case_category = case_when(
    total_cases == 0 ~ "0 cases",
    total_cases >= 1 & total_cases <= 9 ~ "1-9 cases",
    total_cases >= 10 ~ "≥10 cases"
  ),
  # Make it a factor with ordered levels
  case_category = factor(case_category, 
                         levels = c("0 cases", "1-9 cases", "≥10 cases")))

# Calculate positivity proportions for Phage
data <- data %>%
  mutate(phage_positivity = positive_samples / total_samples)

# Plot 1: Phage Positivity
p1 <- ggplot(data, aes(x = case_group, y = phage_positivity, fill = case_group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 0.8) +
theme_minimal() +
  labs(x = "", y = "Phage positivity", fill = "")
p1

data %>% group_by(case_group) %>%
  summarise(
    n_sites = n(), 
    median_phage_positivity = median(phage_positivity, na.rm = TRUE)
  )
# Plot with categorical breakdown
p2 <- ggplot(data, aes(x = case_category, y = phage_positivity, fill = case_category)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  stat_compare_means(method = "kruskal.test", label.y = max(data$phage_positivity, na.rm = TRUE) * 1.05) +
  theme_minimal() +
  labs(x = "Typhoid case category", 
       y = "Phage positivity", 
       fill = "Case category") +
  scale_fill_brewer(palette = "Set2")

p2

data %>%
  group_by(case_category) %>%
  summarise(
    n_sites = n(),
    median_phage_positivity = median(phage_positivity, na.rm = TRUE)
  )

# Optional: Pairwise comparisons between groups
library(ggpubr)
my_comparisons <- list(c("0 cases", "1-9 cases"), 
                       c("1-9 cases", "≥10 cases"), 
                       c("0 cases", "≥10 cases"))

p3 <- ggplot(data, aes(x = case_category, y = phage_positivity, fill = case_category)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  stat_compare_means(label.y = max(data$phage_positivity, na.rm = TRUE) * 1.15) +
  theme_minimal() +
  labs(x = "Typhoid case category", 
       y = "Phage positivity", 
       fill = "Case category") +
  scale_fill_brewer(palette = "Set2")

p3

# Summary statistics by category
data %>%
  group_by(case_category) %>%
  summarise(
    n_sites = n(),
    median_positivity = median(phage_positivity, na.rm = TRUE),
    IQR_positivity = IQR(phage_positivity, na.rm = TRUE),
    mean_positivity = mean(phage_positivity, na.rm = TRUE),
    sd_positivity = sd(phage_positivity, na.rm = TRUE)
  )

#Plot 2: Moore Swab PCR Positive
p2 <- ggplot(data, aes(x = case_group, y = MS..PCR., fill = case_group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 8) +
  theme_minimal() +
  labs(x = "", y = "MS PCR positive", fill = "")

#Plot 3: Grab PCR Positive (NA PCR)
p3 <- ggplot(data, aes(x = case_group, y = NA..PCR., fill = case_group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 3) +
  theme_minimal() +
  labs(x = "", y = "Grab PCR positive", fill = "")

# Arrange plots
env_surv <- ggarrange(p1, ncol = 1, nrow = 1)
ggsave("env_surv20251017.png", env_surv, width = 12, height = 4, dpi = 300)

#### MAIN MODEL ================================================================

library(dplyr)
library(lubridate)

merged_df <- merged_df %>%
  arrange(Site, month) %>%
  group_by(Site) %>%
  mutate(
    cases_lag1 = lag(cases, 1),
    cases_lag2 = lag(cases, 2)  
    ) %>%
  ungroup()

# Basic model: Does previous-month case count predict phage positivity?
model_lag1 <- glm(
  cbind(positive_samples, total_samples - positive_samples) ~ 
    cases_lag1 + as.factor(Site),
  data = merged_df,
  family = binomial
)
exp(coef(model_lag1))         # Odds ratios
exp(confint(model_lag1))      # Confidence intervals

# 2-month lag model
model_lag2 <- glm(
  cbind(positive_samples, total_samples - positive_samples) ~ 
    cases_lag2 + as.factor(Site),
  data = merged_df,
  family = binomial
)
exp(coef(model_lag2))
exp(confint(model_lag2))

# Plaques
mean(merged_df$total_plaques, na.rm = TRUE)
var(merged_df$total_plaques, na.rm = TRUE)
mean(merged_df$total_plaques == 0, na.rm = TRUE)  # majority are 0's

library(pscl)
library(MASS)

# Standard neg binom model
model_nb <- glm.nb(
  total_plaques ~ cases_lag1 + as.factor(Site),
  data = merged_df
)

summary(model_nb)

# Zero-Inflated neg binom model
model_zinb <- zeroinfl(
  total_plaques ~ cases_lag2 + as.factor(Site) | cases_lag1,
  dist = "negbin",
  data = merged_df
)

summary(model_zinb)
exp(coef(model_zinb))
exp(confint(model_zinb))

AIC(model_nb, model_zinb)

# REVERSE
merged_df <- merged_df %>%
  arrange(Site, month) %>%
  group_by(Site) %>%
  mutate(
    phage_lag1 = lag(phage_pos_rate, 1),
    phage_lag2 = lag(phage_pos_rate, 2),
    plaques_lag1 = lag(total_plaques, 1),
    plaques_lag2 = lag(total_plaques, 2)
  ) %>%
  ungroup()

model_cases_lag1 <- glm(
  I(cases > 0) ~ phage_lag1 + as.factor(Site),
  data = merged_df,
  family = binomial
)

exp(coef(model_cases_lag1))  # odds ratios
exp(confint(model_cases_lag1)) # CIs

model_cases_lag2 <- glm(
  I(cases > 0) ~ phage_lag2 + as.factor(Site),
  data = merged_df,
  family = binomial
)

exp(coef(model_cases_lag2))  # odds ratios
exp(confint(model_cases_lag2)) # CIs

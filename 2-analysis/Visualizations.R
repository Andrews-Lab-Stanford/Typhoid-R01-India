#### TYPHOID R01 HIGH-LEVEL VISUALIZATIONS ####
# Purpose: Generate Monthly Trends, Correlations, and Spatial Hotspots
# Derived from: TyphoidR01_SAP_Analysis.R and shared slide examples

library(tidyverse)
library(readxl)
library(lubridate)
library(gridExtra)
library(patchwork) # for combining plots

# 1. CONFIGURE PATHS
data_dir <- "Data"
manuscript_dir <- "R01-Typhoid-Phage-Wastewater-Manuscript/Analysis Data"

# 2. LOAD SHARED DATA (FROM SAP ANALYSIS)
data_rdata <- file.path("Data", "cleaned_analysis_data.RData")
if (!file.exists(data_rdata)) {
  stop("Cleaned data not found. Please run TyphoidR01_SAP_Analysis.R first.")
}
load(data_rdata) # Loads 'res_intermediate'

# 3. GENERATE VISUALIZATIONS

# 3.1 SLIDE 1 LEFT: GLOBAL MONTHLY TRENDS
global_trends <- res_intermediate %>%
  group_by(month) %>%
  summarise(
    total_cases = sum(clinical_cases, na.rm = TRUE),
    n = sum(total_samples, na.rm = TRUE),
    pos = sum(positive_samples, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(positivity_pct = (pos / n) * 100)

# Dual-axis scaling factor
scale_fac <- max(global_trends$total_cases, na.rm=T) / 100

trend_plot <- ggplot(global_trends, aes(x = month)) +
  geom_col(aes(y = total_cases), fill = "steelblue", alpha = 0.7) +
  geom_line(aes(y = positivity_pct * scale_fac), color = "darkred", linewidth = 1.2) +
  geom_point(aes(y = positivity_pct * scale_fac), color = "darkred", size = 2) +
  scale_y_continuous(
    name = "Total Monthly Typhoid Cases",
    sec.axis = sec_axis(~./scale_fac, name = "S. Typhi Detection in WW (%)")
  ) +
  theme_minimal() +
  labs(title = "Monthly Trends of Typhoid Fever and Wastewater Detection",
       x = "Month-Year") +
  theme(axis.title.y.right = element_text(color = "darkred"),
        axis.text.y.right = element_text(color = "darkred"),
        plot.title = element_text(face = "bold", size = 14))

trend_plot

# 3.2 SLIDE 1 RIGHT: CORRELATION (Independent Catchments)
corr_df <- res_intermediate %>%
  group_by(Site) %>%
  summarise(
    total_cases = sum(clinical_cases, na.rm = TRUE),
    Pop = max(Pop, na.rm = TRUE),
    positivity_rate = sum(positive_samples, na.rm = TRUE) / sum(total_samples, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(positivity_rate))

library(ggpubr)

corr_plot <- ggplot(corr_df, aes(x = total_cases, y = positivity_rate)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.7) +
  geom_smooth(method = "loess", span = 0.7, color = "darkred", fill = "pink", alpha = 0.2) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4) +
  theme_minimal() +
  labs(title = "Correlation: Wastewater Detection vs. Disease Burden",
       subtitle = "Unified Analysis across Non-overlapping Catchments (38 Sites)",
       x = "Total Number of Clinical Typhoid Cases",
       y = "Wastewater Detection (Overall Positivity Ratio)") +
  theme(plot.title = element_text(face = "bold", size = 14))

# 3.3 SUPPLEMENTARY: DATA DISTRIBUTIONS (Site-Month Comparison)
dist_df <- res_intermediate %>%
  mutate(clinical_cases = replace_na(clinical_cases, 0),
         total_plaques = replace_na(total_plaques, 0))

p1 <- ggplot(dist_df, aes(x = clinical_cases)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  theme_minimal() + labs(title = "A) Clinical Cases: All Site-Months", x = "Cases", y = "Freq")

p2 <- ggplot(dist_df %>% filter(clinical_cases > 0), aes(x = clinical_cases)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white", boundary = 0.5) +
  theme_minimal() + labs(title = "B) Clinical Cases: Positives Only", x = "Cases", y = "Freq") +
  xlim(0.5, NA) 

p3 <- ggplot(dist_df, aes(x = total_plaques)) +
  geom_histogram(bins = 30, fill = "darkred", color = "white") +
  theme_minimal() + labs(title = "C) Phage Abundance: All Site-Months", x = "Plaques", y = "Freq")

p4 <- ggplot(dist_df %>% filter(total_plaques > 0), aes(x = total_plaques)) +
  geom_histogram(binwidth = 10, fill = "darkred", color = "white", boundary = 0) +
  theme_minimal() + labs(title = "D) Phage Abundance: Positives Only", x = "Plaques", y = "Freq") +
  xlim(0.5, NA)

dist_plot_clinical <- p1 / p2 + plot_annotation(title = "Clinical Case Distributions (38 Sites)")
dist_plot_phage <- p3 / p4 + plot_annotation(title = "Phage Abundance Distributions (38 Sites)")

# 3.4 SUPPLEMENTARY: COVARIATE TRENDS OVER TIME
cov_trends <- res_intermediate %>%
  group_by(month) %>%
  summarise(
    mean_temp = mean(avg_temp, na.rm = TRUE),
    mean_rain = mean(total_rain, na.rm = TRUE),
    mean_flow = mean(avg_flow, na.rm = TRUE),
    .groups = "drop"
  )

ct1 <- ggplot(res_intermediate, aes(x = month, y = avg_temp)) +
  geom_line(aes(group = Site), color = "gray80", alpha = 0.3) +
  geom_line(data = cov_trends, aes(y = mean_temp), color = "darkorange", linewidth = 1.2) +
  theme_minimal() + labs(title = "Average Temperature (°C)", x = NULL, y = "Temp")

ct2 <- ggplot(res_intermediate, aes(x = month, y = total_rain)) +
  geom_line(aes(group = Site), color = "gray80", alpha = 0.3) +
  geom_line(data = cov_trends, aes(y = mean_rain), color = "blue", linewidth = 1.2) +
  theme_minimal() + labs(title = "Total Monthly Rainfall (mm)", x = NULL, y = "Rain")

ct3 <- ggplot(res_intermediate, aes(x = month, y = avg_flow)) +
  geom_line(aes(group = Site), color = "gray80", alpha = 0.3) +
  geom_line(data = cov_trends, aes(y = mean_flow), color = "darkgreen", linewidth = 1.2) +
  theme_minimal() + labs(title = "Average Monthly Flow (m3/h)", x = "Month", y = "Flow")

cov_plot <- (ct1 / ct2 / ct3) + plot_annotation(title = "Environmental Covariate Trends (38 Sites)",
                                               subtitle = "Gray lines = individual sites | Heavy lines = global monthly mean")
# 4. SAVE ARTIFACTS
ggsave("Typhoid_Global_Trends.png", trend_plot, width = 10, height = 6)
ggsave("Typhoid_Correlations.png", corr_plot, width = 10, height = 8)
ggsave("Typhoid_Distributions_Clinical.png", dist_plot_clinical, width = 10, height = 10)
ggsave("Typhoid_Distributions_Phage.png", dist_plot_phage, width = 10, height = 10)
ggsave("Typhoid_Covariate_Trends.png", cov_plot, width = 10, height = 12)
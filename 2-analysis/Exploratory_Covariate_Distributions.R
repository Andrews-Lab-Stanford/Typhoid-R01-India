#### TYPHOID R01 — EXPLORATORY COVARIATE DISTRIBUTIONS ####
# Purpose: Review the distribution of site-level demographic covariates (Age, WASH, SES).
# This helps determine if our current binning/summarization captures the variance across sites.

library(tidyverse)
library(patchwork)

# ── 1. Load Data ──────────────────────────────────────────────────────────────
base_dir   <- "/Users/esthjung/Library/CloudStorage/GoogleDrive-esthjung@stanford.edu/My Drive/Typhoid R01"
rdata_path <- file.path(base_dir, "Data/analysis_ready_weekly.RData")
load(rdata_path) # loads: census_site, site_status

# Join status to census metadata
site_meta <- census_site %>%
  inner_join(site_status, by = "Site") %>%
  # Rename columns to match plot logic if necessary, or just use them
  rename(
    pct_water = pct_water_improved,
    pct_toilet = pct_toilet_improved
  )

# ── 2. Visualizations (Grouped by Site Status) ─────────────────────────────
# Set colors
group_colors <- c("FALSE" = "#95a5a6", "TRUE" = "#e74c3c")

plot_compare <- function(df, var, title, ylab) {
  ggplot(df, aes(x = has_cases, y = !!sym(var), fill = has_cases)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_fill_manual(values = group_colors) +
    labs(title = title, x = "Has Clinical Cases", y = ylab) +
    theme_minimal() +
    theme(legend.position = "none")
}

p1 <- plot_compare(site_meta, "pct_under_15", "Age: % Under 15", "% of Population")
p2 <- plot_compare(site_meta, "pct_high_ses", "SES: % High (Class I+II+III)", "% of Population")
p3 <- plot_compare(site_meta, "pct_water",    "WASH: % Improved Water", "% of Population")
p4 <- plot_compare(site_meta, "pct_toilet",   "WASH: % Improved Toilet", "% of Population")

combined_plot <- (p1 | p2) / (p3 | p4) + 
  plot_annotation(title = "Comparing Covariate Distributions: Sites With vs. Without Clinical Cases (N=38)")
 
combined_plot

# Print Group Means
cat("\n--- Mean Values by Site Group ---\n")
site_meta %>%
  group_by(has_cases) %>%
  summarise(
    mean_under_15 = mean(pct_under_15),
    mean_high_ses = mean(pct_high_ses),
    mean_water    = mean(pct_water),
    mean_toilet   = mean(pct_toilet),
    .groups = "drop"
  ) %>%
  print()

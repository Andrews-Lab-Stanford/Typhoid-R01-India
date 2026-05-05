#### TYPHOID R01 — STEP 1: DATA CLEANING & PREPARATION ####
# Last Updated: May 2026
# Purpose: Load pre-processed RData, prepare analysis-ready datasets, and save.
#
# INPUT:  Data/cleaned_analysis_data.RData  (contains 'res_intermediate')
# OUTPUT: Data/analysis_ready.RData
#         — merged_data   (full site-month panel with all covariates)
#         — res_table1    (descriptive dataset; burn-in excluded)
#         — analysis_df   (analytic dataset; scaled, complete cases)
#
# Run this FIRST before 2_Main_Models.R or 3_Table_1.R


library(tidyverse)
library(lubridate)

# ── Directories ────────────────────────────────────────────────────────────────
base_dir  <- "/Users/esthjung/Library/CloudStorage/GoogleDrive-esthjung@stanford.edu/My Drive/Typhoid R01"
data_dir  <- file.path(base_dir, "Data")

cat("=== TYPHOID R01: DATA PREPARATION ===\n")

# ── 1. Load Pre-Cleaned Snapshot ───────────────────────────────────────────────
rdata_path <- file.path(data_dir, "cleaned_analysis_data.RData")
if (!file.exists(rdata_path)) {
  stop("cleaned_analysis_data.RData not found at: ", rdata_path,
       "\nPlease ensure the raw data pipeline has been run to generate this file.")
}
load(rdata_path)  # loads: res_intermediate
cat("Loaded res_intermediate:", nrow(res_intermediate), "rows,",
    n_distinct(res_intermediate$Site), "sites\n")

# ── 2. Column Mapping ──────────────────────────────────────────────────────────
# Standardise column names used throughout the analysis
merged_data <- res_intermediate %>%
  mutate(
    temp     = avg_temp,
    rainfall = total_rain,
    flow     = avg_flow,
    hf183    = hf183_abs
  )

# ── 3. Site-Level Case Status ──────────────────────────────────────────────────
# has_cases = TRUE if the site ever recorded ≥1 clinical typhoid case
site_status <- merged_data %>%
  group_by(Site) %>%
  summarise(has_cases = any(clinical_cases > 0, na.rm = TRUE), .groups = "drop")

merged_data <- merged_data %>%
  dplyr::select(-any_of("has_cases")) %>%
  inner_join(site_status, by = "Site")

cat("  Sites with cases:   ", sum(site_status$has_cases), "\n")
cat("  Sites without cases:", sum(!site_status$has_cases), "\n")

# ── 4. Descriptive Dataset (Table 1) ──────────────────────────────────────────
# Exclude burn-in months (first 2 months per site needed to compute 2-month lag)
res_table1 <- merged_data %>% filter(!is.na(cases_lag2))
cat("Table 1 dataset:", nrow(res_table1), "site-months\n")

# ── 5. Analytic Dataset (Models) ──────────────────────────────────────────────
# Further filter to complete cases for weather & HF183, then z-score scale covariates
analysis_df <- res_table1 %>%
  filter(!is.na(temp), !is.na(hf183)) %>%
  mutate(across(
    c(flow, rainfall, temp, hf183,
      pct_under_15, pct_water_improved, pct_toilet_improved, pct_high_ses),
    scale
  ))
cat("Analytic dataset:  ", nrow(analysis_df), "site-months\n")

# ── 6. Save Analysis-Ready Data ────────────────────────────────────────────────
output_path <- file.path(data_dir, "analysis_ready.RData")
save(merged_data, res_table1, analysis_df, file = output_path)
cat("\nSUCCESS: Saved to", output_path, "\n")
cat("Objects saved: merged_data, res_table1, analysis_df\n")

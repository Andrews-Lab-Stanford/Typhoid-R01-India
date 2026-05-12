#### TYPHOID R01 — STEP 1: WEEKLY DATA CLEANING & PREPARATION ####
# Last Updated: May 2026 by Esther Jung
# Purpose: Load raw data, aggregate to weekly levels, and prepare analysis-ready weekly datasets.
#
# INPUT:  data/REDCapPhageData.csv (phage dataset)
#         data/Typhoidcases_site_monthly weekly.xlsx (sheet 2: weekly cases)
#         R01-Typhoid-Phage-Wastewater-Manuscript/Analysis Data/weather_flow_hf183.xlsx (weather, ww flow, and HF183 covariates)
#         data/census_data.xlsx (for site-level census data)
# OUTPUT: Data/analysis_ready_weekly.RData
#
# Run this FIRST before 2_Main_Models_Weekly.R and 3_Table_1_Weekly.R

library(tidyverse)
library(readxl)
library(lubridate)

# ── Directories ────────────────────────────────────────────────────────────────
base_dir       <- "/Users/esthjung/Library/CloudStorage/GoogleDrive-esthjung@stanford.edu/My Drive/Typhoid R01"
root_data_dir  <- file.path(base_dir, "Data")
repo_data_dir  <- file.path(base_dir, "Typhoid-R01-India/data")
manuscript_dir <- repo_data_dir

# ── 1. Load & Process Weekly Clinical Cases ────────────────────────────────────
cases_file <- file.path(repo_data_dir, "Typhoidcases_site_monthly weekly.xlsx")
cat("Loading weekly cases from:", cases_file, "\n")

# Load overlap status from the 'monthly' sheet
overlap_status <- read_excel(cases_file, sheet = "monthly") %>%
  dplyr::select(Site = ESsite, overlap)

cases_raw <- read_excel(cases_file, sheet = "weekly")

# Pivot long, parse weeks, and filter for overlap == 1
cases_weekly <- cases_raw %>%
  dplyr::select(-cumulative_cases) %>%
  pivot_longer(cols = -ESsite, names_to = "week_str", values_to = "clinical_cases") %>%
  rename(Site = ESsite) %>%
  mutate(Site = as.character(as.numeric(Site))) %>%
  inner_join(overlap_status %>% mutate(Site = as.character(as.numeric(Site))), by = "Site") %>%
  filter(overlap == 1) %>%
  mutate(
    # Parse ISO 8601 week string (YYYY-WXX) to the Monday of that week
    # ISO week 1 always contains Jan 4th.
    week_date = map_vec(week_str, function(s) {
      parts <- strsplit(s, "-W")[[1]]
      y <- as.numeric(parts[1])
      w <- as.numeric(parts[2])
      # Jan 4 is always in ISO Week 1
      jan4 <- as.Date(paste0(y, "-01-04"))
      # Find the Monday of that week
      w1_monday <- floor_date(jan4, unit = "week", week_start = 1)
      # Add the number of weeks
      return(w1_monday + weeks(w - 1))
    })
  ) %>%
  filter(!is.na(week_date)) %>%
  # Ensure chronological completeness for each site before lagging
  group_by(Site) %>%
  complete(week_date = seq.Date(min(week_date), max(week_date), by = "week"), 
           fill = list(clinical_cases = 0)) %>%
  arrange(Site, week_date) %>%
  mutate(
    cases_lag0 = clinical_cases,
    cases_lag1 = lag(clinical_cases, 1),
    cases_lag2 = lag(clinical_cases, 2),
    cases_lag3 = lag(clinical_cases, 3),
    cases_lag4 = lag(clinical_cases, 4)
  ) %>%
  ungroup()

# ── 3. Load & Process Phage Data (Weekly) ──────────────────────────────────────
phage_file <- file.path(repo_data_dir, "REDCapPhageData.csv")
cat("Loading phage data from:", phage_file, "\n")

phage_raw <- read.csv(phage_file)
phage_weekly <- phage_raw %>%
  mutate(
    date_sample = as.Date(parse_date_time(date_results_tw, c("mdy", "ymd"))),
    # Month-year correction for 2004/pre-2024 dates
    date_sample = as.Date(ifelse(
      !is.na(date_sample) & (year(date_sample) == 2004 | year(date_sample) < 2024),
      paste0("2024", substr(as.character(date_sample), 5, 10)),
      as.character(date_sample)
    )),
    # Align to Monday of the week
    week_date = floor_date(date_sample, unit = "week", week_start = 1),
    status = case_when(
      amplification_tw == 1 | direct_tw == 1 ~ "Positive",
      amplification_tw == 2 & direct_tw == 2 ~ "Negative",
      TRUE ~ "Negative"
    )
  ) %>%
  group_by(Site = as.character(as.numeric(site_id_tw)), week_date) %>%
  summarise(
    total_samples    = n(),
    positive_samples = sum(status == "Positive", na.rm = TRUE),
    total_plaques    = sum(direct_counts_tw, na.rm = TRUE),
    .groups = "drop"
  )

# ── 4. Load & Process Weather/HF183 Data (Weekly) ─────────────────────────────
weather_file <- file.path(repo_data_dir, "weather_flow_hf183.xlsx")
cat("Loading weather/HF183 from:", weather_file, "\n")

weather_raw <- read_excel(weather_file)
weather_weekly <- weather_raw %>%
  mutate(
    Site = as.character(as.numeric(substr(as.character(Sample_ID), 1, 3))),
    # Align to Monday of the week
    week_date = floor_date(as.Date(date_collection_new), unit = "week", week_start = 1),
    rainfall  = ifelse(tolower(`x24_hours_rainfall_mm_recorded_from_0830_hrs_ist_of_yesterday_to_0830_hrs_ist_of_today`) %in% c("nil", "trace"), 0, 
                       as.numeric(`x24_hours_rainfall_mm_recorded_from_0830_hrs_ist_of_yesterday_to_0830_hrs_ist_of_today`)),
    temp      = (as.numeric(maximum_temp_o_c_recorded) + as.numeric(minimum_temp_o_c_recorded)) / 2,
    flow      = as.numeric(g_per_hour_in_m3),
    hf183     = as.numeric(hf_183_tr)
  ) %>%
  group_by(Site, week_date) %>%
  summarise(
    temp     = mean(temp, na.rm = TRUE),
    rainfall = sum(rainfall, na.rm = TRUE),
    flow     = mean(flow, na.rm = TRUE),
    hf183    = mean(hf183, na.rm = TRUE),
    .groups = "drop"
  )

# ── 5. Load Site-Level Census Data ───────────────────────────────────────────
census_file_raw <- file.path(repo_data_dir, "census_data.xlsx")
cat("Loading raw census data from:", census_file_raw, "\n")

# Read the four domains from separate sheets
age_df    <- read_excel(census_file_raw, sheet = "Age category")
water_df  <- read_excel(census_file_raw, sheet = "Source of drinking water")
toilet_df <- read_excel(census_file_raw, sheet = "Toilet_source")
ses_df    <- read_excel(census_file_raw, sheet = "SES")

# Also need the 'Pop' variable (total catchment population)
# This usually comes from a separate sheet or the Total_n of the Age sheet
# In the RData it was named 'Pop'. Let's use Total_n from Age sheet as Pop.
census_combined <- age_df %>%
  dplyr::select(Site = CA, Pop = Total_n, `0-4_n`, `5-14_n`, `15-49_n`, `≥50_n`, Total_n) %>%
  left_join(water_df %>% dplyr::select(Site = CA, Improved_n_water = Improved_n, Total_n_water = Total_n), by = "Site") %>%
  left_join(toilet_df %>% dplyr::select(Site = CA, Improved_n_toilet = Improved_n, Total_n_toilet = Total_n), by = "Site") %>%
  left_join(ses_df %>% dplyr::select(Site = CA, `Class I_n`, `Class II_n`, `Class III_n`, `Class IV_n`, `Class V_n`, Total_n_ses = Total_n), by = "Site")

census_site <- census_combined %>%
  mutate(Site = as.character(as.numeric(Site))) %>%
  group_by(Site) %>%
  summarise(across(everything(), first), .groups = "drop") %>%
  mutate(
    pct_under_15       = 100 * (as.numeric(`0-4_n`) + as.numeric(`5-14_n`)) / as.numeric(Pop),
    pct_water_improved = 100 * as.numeric(Improved_n_water) / as.numeric(Total_n_water),
    pct_toilet_improved= 100 * as.numeric(Improved_n_toilet) / as.numeric(Total_n_toilet),
    # High SES = Class 1, 2, and 3
    pct_high_ses       = 100 * (as.numeric(`Class I_n`) + as.numeric(`Class II_n`) + as.numeric(`Class III_n`)) / as.numeric(Total_n_ses)
  )

# ── 6. Merge Weekly Data ──────────────────────────────────────────────────────
merged_weekly <- cases_weekly %>%
  inner_join(phage_weekly, by = c("Site", "week_date")) %>%
  left_join(weather_weekly, by = c("Site", "week_date")) %>%
  left_join(census_site,    by = "Site")

cat("Merged Weekly Dataset:", nrow(merged_weekly), "site-weeks\n")

# ── 7. Scaling & Site Status ──────────────────────────────────────────────────
# IMPORTANT: 'has_cases' status should be based on the FULL clinical surveillance panel,
# not just the weeks where we happened to take a wastewater sample.
site_status <- cases_weekly %>%
  group_by(Site) %>%
  summarise(has_cases = any(clinical_cases > 0, na.rm = TRUE), .groups = "drop")

merged_weekly <- merged_weekly %>%
  inner_join(site_status, by = "Site")

# Final analytic set: filter complete cases for weather/HF183 and scale
analysis_df_weekly <- merged_weekly %>%
  filter(!is.na(temp), !is.na(hf183), !is.na(cases_lag4)) %>%
  mutate(across(
    c(flow, rainfall, temp, hf183,
      pct_under_15, pct_water_improved, pct_toilet_improved, pct_high_ses),
    scale
  ))

cat("Final Weekly Analytic dataset:", nrow(analysis_df_weekly), "site-weeks\n")

# ── 8. Save ───────────────────────────────────────────────────────────────────
output_path <- file.path(root_data_dir, "analysis_ready_weekly.RData")
save(cases_weekly, merged_weekly, analysis_df_weekly, census_site, site_status, file = output_path)

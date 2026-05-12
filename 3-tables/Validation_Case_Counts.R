#### TYPHOID R01 — ANALYTIC DATA VALIDATION ####
# Purpose: Briefly verify the final analytic dataset for Table 1 (Population, Sites, and Cases).

library(tidyverse)

# ── 1. Load Cleaned Analytic Data ─────────────────────────────────────────────
base_dir   <- "/Users/esthjung/Library/CloudStorage/GoogleDrive-esthjung@stanford.edu/My Drive/Typhoid R01"
rdata_path <- file.path(base_dir, "Data/analysis_ready_weekly.RData")
load(rdata_path) # loads: cases_weekly, merged_weekly, analysis_df_weekly

# ── 2. Site & Case Validation (Full Panel) ────────────────────────────────────
cat("\n--- SITE CLASSIFICATION (Full Clinical Panel) ---\n")
site_summary <- cases_weekly %>%
  group_by(Site) %>%
  summarise(
    total_cases = sum(clinical_cases, na.rm = TRUE),
    has_any_cases = total_cases > 0,
    .groups = "drop"
  )

print(table(has_any_cases = site_summary$has_any_cases))

cat("\n--- TOTAL CLINICAL CASES (Full Panel) ---\n")
cat("Total Cases:", sum(site_summary$total_cases), "\n")

# ── 3. Population Validation ──────────────────────────────────────────────────
cat("\n--- POPULATION BY SITE GROUP ---\n")
# Extract site-level population from merged data
site_pop <- merged_weekly %>%
  group_by(Site) %>%
  summarise(
    Pop = first(Pop),
    has_cases = first(has_cases),
    .groups = "drop"
  )

pop_audit <- site_pop %>%
  group_by(has_cases) %>%
  summarise(
    n_sites   = n(),
    total_pop = sum(Pop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(pct_pop = round(100 * total_pop / sum(total_pop), 1))

print(pop_audit)

cat("\n--- GRAND TOTAL POPULATION (38 Sites) ---\n")
cat("Total Population:", sum(site_pop$Pop), "\n")
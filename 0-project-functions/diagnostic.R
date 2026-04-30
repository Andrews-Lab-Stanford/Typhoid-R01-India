library(readxl)
library(tidyverse)
library(lubridate)

base_dir <- "/Users/esthjung/Library/CloudStorage/GoogleDrive-esthjung@stanford.edu/My Drive/Typhoid R01"
manuscript_dir <- file.path(base_dir, "R01-Typhoid-Phage-Wastewater-Manuscript/Analysis Data")
data_dir <- file.path(base_dir, "Data")

# 1. Phage
phage_raw <- read.csv(file.path(data_dir, "REDCapPhageData.csv"))
cat("Phage records:", nrow(phage_raw), "\n")
print(head(phage_raw$date_results_tw))

# 2. Cases
cases_file <- file.path(manuscript_dir, "Site_monthlytyphoidcases_Jan 2026.xlsx")
case_sheets <- excel_sheets(cases_file)
cat("Case sheets:", length(case_sheets), "\n")

# 3. Weather
weather_raw <- read_excel(file.path(manuscript_dir, "weather_flow_hf183.xlsx"))
cat("Weather records:", nrow(weather_raw), "\n")

# 4. Census
census_file <- file.path(manuscript_dir, "sample_census_data.xlsx")
age_raw <- read_excel(census_file, sheet = "Age category")
cat("Census Age records:", nrow(age_raw), "\n")

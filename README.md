# Typhoid R01 — Vellore, India Analysis Pipeline

This repository contains the analytical pipeline for the Typhoid R01 study, a collaboration between Stanford, CMC Vellore, MGH, and UC Davis. This guide outlines the steps required to replicate the manuscript results, from raw data cleaning to the final regression models.

---

## Getting Started (Core Pipeline)

Follow these steps in order to generate the primary manuscript outputs. All scripts should be run from the repository root.

### **Step 1: Data Preparation & Cleaning**
Run the cleaning script to ingest raw Excel/REDCap data, synchronize weekly site dates, and generate the analytic environment.
- **Script**: `1-data-cleaning/1_Data_Cleaning_Weekly.R`
- **Input**: Raw files in `data/` (Clinical cases, Census, and Weather covariates).
- **Output**: `Data/analysis_ready_weekly.RData` (Saves objects to your local drive for subsequent steps).

### **Step 2: Table 1 (Descriptive Statistics)**
Generate the manuscript Table 1 comparing sites with vs. without typhoid cases.
- **Script**: `3-tables/3_Table_1_Weekly.R`
- **Output**: `results/Table1_Weekly.csv`
- **Note**: This script pulls clinical totals from the full surveillance panel (70 cases) and environmental medians from sampled weeks.

### **Step 3: Main Regression Models**
Run the primary negative binomial and binary outcome models for the manuscript.
- **Script**: `2-analysis/2_Main_Models_Weekly.R`
- **Output**: Model coefficients and statistics printed to console/results.

---

## Data Verification
To independently verify the sites with and without cases, the total population denominator, and clinical case counts, run the brief validation tool:
- **Script**: `3-tables/Validation_Case_Counts.R`

---

## Directory Structure

### **1-data-cleaning/**
- `1_Data_Cleaning_Weekly.R`: Primary cleaning script for the weekly manuscript pipeline.
- `1_Data_Cleaning.R`: (Legacy) Original monthly granularity cleaning script.

### **2-analysis/**
- `2_Main_Models_Weekly.R`: **Primary manuscript models.**
- `R01_Team_Meeting_Analyses.R`: Scripts used for collaborative status updates and team meetings.
- `Exploratory_...`: Various scripts for mapping, model testing, and trend analysis.

### **3-tables/**
- `3_Table_1_Weekly.R`: **Primary manuscript Table 1 generator.**
- `Validation_Case_Counts.R`: Audit tool for data verification.

### **data/**
Contains raw datasets including:
- `census_data.xlsx`: Raw demographic and SES data.
- `Typhoidcases_site_monthly weekly.xlsx`: Raw clinical surveillance data.
- `weather_flow_hf183.xlsx`: Environmental and phage density covariates.

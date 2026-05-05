#### TYPHOID R01 SAP ANALYSIS ####
# Last Updated: April 24, 2026
# Reverted to optimized fast-processing version with updated Table 1 variables.

library(tidyverse)
library(readxl)
library(lubridate)
library(glmmTMB)
library(broom.mixed)

# Set correct directories
base_dir <- "/Users/esthjung/Library/CloudStorage/GoogleDrive-esthjung@stanford.edu/My Drive/Typhoid R01"
manuscript_dir <- file.path(base_dir, "R01-Typhoid-Phage-Wastewater-Manuscript/Analysis Data")
data_dir <- file.path(base_dir, "Data")

#### 1. DATA LOADING & CLEANING (FAST MODE) ####
cat("Loading pre-cleaned data from RData...\n")
rdata_path <- file.path(data_dir, "cleaned_analysis_data.RData")

if (file.exists(rdata_path)) {
  load(rdata_path) # Loads 'res_intermediate'
  
  # Ensure the data structure matches the expected names for the rest of the script
  # Map column names if they differ in the RData (e.g., hf183_abs vs hf183)
  merged_data <- res_intermediate %>%
    mutate(
      temp = avg_temp,
      rainfall = total_rain,
      flow = avg_flow,
      hf183 = hf183_abs
    )
  
  # Define Site-level case status (ever positive)
  site_status <- merged_data %>% 
    group_by(Site) %>% 
    summarise(has_cases = any(clinical_cases > 0, na.rm = TRUE), .groups = "drop")
  
  merged_data <- merged_data %>% 
    dplyr::select(-any_of("has_cases")) %>% 
    inner_join(site_status, by = "Site")
  
  # Descriptive Dataset (Exclude burn-in)
  res_table1 <- merged_data %>% filter(!is.na(cases_lag2))
  
  # Analytic Dataset (Exclude missing weather/HF183)
  analysis_df <- res_table1 %>%
    filter(!is.na(temp), !is.na(hf183)) %>%
    mutate(across(c(flow, rainfall, temp, hf183, pct_under_15, pct_water_improved, pct_toilet_improved, pct_high_ses), scale))
  
  cat("SUCCESS: Data loaded and prepared.\n")
} else {
  stop("Pre-cleaned data file not found at: ", rdata_path)
}

#### 3. TABLE 1: DESCRIPTIVE SUMMARY ####

get_val <- function(df, var, type = "n_percent", overall_val = NULL) {
  if (nrow(df) == 0) return("0")
  if (type == "sum") return(format(sum(df[[var]], na.rm=T), big.mark=","))
  if (type == "median_iqr") {
    return(paste0(round(median(df[[var]], na.rm=T), 1), " (", round(quantile(df[[var]], 0.25, na.rm=T), 1), "-", round(quantile(df[[var]], 0.75, na.rm=T), 1), ")"))
  }
  if (type == "n_percent_pop") {
    total_var <- case_when(grepl("_n$|years$", var) ~ "Total_n", grepl("water", var) ~ "Total_n_water", grepl("toilet", var) ~ "Total_n_toilet", grepl("Class", var) ~ "Total_n_ses", TRUE ~ "Total_n")
    site_vals <- df %>% group_by(Site) %>% summarise(v = first(!!sym(var)), t = first(!!sym(total_var)), .groups="drop")
    v_sum <- sum(site_vals$v, na.rm=T); t_sum <- sum(site_vals$t, na.rm=T)
    return(paste0(format(v_sum, big.mark=","), " (", round(100*v_sum/t_sum, 1), "%)"))
  }
  if (type == "sum_n_percent") {
    n <- if(var == "Pop") sum((df %>% group_by(Site) %>% summarise(v = first(Pop)))$v) else sum(df[[var]], na.rm=T)
    pct <- round(100 * n / as.numeric(overall_val), 1)
    return(paste0(format(n, big.mark=","), " (", pct, "%)"))
  }
  if (type == "n_percent_samples") {
    pos <- sum(df$positive_samples, na.rm=T); tot <- sum(df$total_samples, na.rm=T)
    return(paste0(pos, " (", round(100*pos/tot, 1), "%)"))
  }
  return("")
}

get_p_value <- function(var_name, data_type = "continuous") {
  if (data_type == "continuous") return(format.pval(wilcox.test(res_table1[[var_name]] ~ res_table1$has_cases)$p.value, digits=3))
  if (data_type == "continuous_site") {
    site_data <- res_table1 %>% group_by(Site) %>% summarise(v = first(!!sym(var_name)), hc = first(has_cases), .groups="drop")
    return(format.pval(wilcox.test(site_data$v ~ site_data$hc)$p.value, digits=3))
  }
  if (data_type == "chisq_pop") {
    total_var <- case_when(grepl("_n$|years$", var_name) ~ "Total_n", grepl("water", var_name) ~ "Total_n_water", grepl("toilet", var_name) ~ "Total_n_toilet", grepl("Class", var_name) ~ "Total_n_ses", TRUE ~ "Total_n")
    site_data <- res_table1 %>% group_by(Site) %>% summarise(v = first(!!sym(var_name)), t = first(!!sym(total_var)), hc = first(has_cases), .groups="drop")
    tab <- matrix(c(sum(site_data$v[site_data$hc], na.rm=T), sum(site_data$t[site_data$hc], na.rm=T) - sum(site_data$v[site_data$hc], na.rm=T), sum(site_data$v[!site_data$hc], na.rm=T), sum(site_data$t[!site_data$hc], na.rm=T) - sum(site_data$v[!site_data$hc], na.rm=T)), nrow=2)
    return(format.pval(fisher.test(tab)$p.value, digits=3))
  }
  if (data_type == "fisher_samples") {
    pos_w <- sum(res_table1$positive_samples[res_table1$has_cases], na.rm=T); tot_w <- sum(res_table1$total_samples[res_table1$has_cases], na.rm=T)
    pos_wo <- sum(res_table1$positive_samples[!res_table1$has_cases], na.rm=T); tot_wo <- sum(res_table1$total_samples[!res_table1$has_cases], na.rm=T)
    return(format.pval(fisher.test(matrix(c(pos_w, tot_w-pos_w, pos_wo, tot_wo-pos_wo), nrow=2))$p.value, digits=3))
  }
  return("-")
}

overall_pop <- sum((res_table1 %>% group_by(Site) %>% summarise(v = first(Pop)))$v)
overall_samples <- sum(res_table1$total_samples)

overall_summary <- data.frame(
  Characteristic = c("Total population, n (%)", "0-4 years", "5-14 years", "15-49 years", ">=50 years", "Class 1 (>=Rs.9130)", "Class 2 (Rs.4565-9129)", "Class 3 (Rs.2739-4564)", "Class 4 (Rs.1369-2738)", "Class 5 (<Rs.1369)", "Improved drinking water", "Improved toilet", "Median population density (IQR)", "Total clinical typhoid cases, n", "Total wastewater samples, n (%)", "Total positive samples, n (%)", "Median HF183 ACN (IQR)", "Median cumulative rainfall, mm (IQR)", "Median temperature, C (IQR)", "Median flow rate, m3/h (IQR)"),
  Overall = c(format(overall_pop, big.mark=","), get_val(res_table1, "0-4_n", "n_percent_pop"), get_val(res_table1, "5-14_n", "n_percent_pop"), get_val(res_table1, "15-49_n", "n_percent_pop"), get_val(res_table1, "≥50_n", "n_percent_pop"), get_val(res_table1, "Class I_n", "n_percent_pop"), get_val(res_table1, "Class II_n", "n_percent_pop"), get_val(res_table1, "Class III_n", "n_percent_pop"), get_val(res_table1, "Class IV_n", "n_percent_pop"), get_val(res_table1, "Class V_n", "n_percent_pop"), get_val(res_table1, "Improved_n_water", "n_percent_pop"), get_val(res_table1, "Improved_n_toilet", "n_percent_pop"), "-", get_val(res_table1, "clinical_cases", "sum"), format(overall_samples, big.mark=","), get_val(res_table1, "positive_samples", "n_percent_samples"), get_val(res_table1, "hf183", "median_iqr"), get_val(res_table1, "rainfall", "median_iqr"), get_val(res_table1, "temp", "median_iqr"), get_val(res_table1, "flow", "median_iqr")),
  With_Cases = c(get_val(res_table1 %>% filter(has_cases), "Pop", "sum_n_percent", overall_pop), get_val(res_table1 %>% filter(has_cases), "0-4_n", "n_percent_pop"), get_val(res_table1 %>% filter(has_cases), "5-14_n", "n_percent_pop"), get_val(res_table1 %>% filter(has_cases), "15-49_n", "n_percent_pop"), get_val(res_table1 %>% filter(has_cases), "≥50_n", "n_percent_pop"), get_val(res_table1 %>% filter(has_cases), "Class I_n", "n_percent_pop"), get_val(res_table1 %>% filter(has_cases), "Class II_n", "n_percent_pop"), get_val(res_table1 %>% filter(has_cases), "Class III_n", "n_percent_pop"), get_val(res_table1 %>% filter(has_cases), "Class IV_n", "n_percent_pop"), get_val(res_table1 %>% filter(has_cases), "Class V_n", "n_percent_pop"), get_val(res_table1 %>% filter(has_cases), "Improved_n_water", "n_percent_pop"), get_val(res_table1 %>% filter(has_cases), "Improved_n_toilet", "n_percent_pop"), "-", get_val(res_table1 %>% filter(has_cases), "clinical_cases", "sum"), get_val(res_table1 %>% filter(has_cases), "total_samples", "sum_n_percent", overall_samples), get_val(res_table1 %>% filter(has_cases), "positive_samples", "n_percent_samples"), get_val(res_table1 %>% filter(has_cases), "hf183", "median_iqr"), get_val(res_table1 %>% filter(has_cases), "rainfall", "median_iqr"), get_val(res_table1 %>% filter(has_cases), "temp", "median_iqr"), get_val(res_table1 %>% filter(has_cases), "flow", "median_iqr")),
  Without_Cases = c(get_val(res_table1 %>% filter(!has_cases), "Pop", "sum_n_percent", overall_pop), get_val(res_table1 %>% filter(!has_cases), "0-4_n", "n_percent_pop"), get_val(res_table1 %>% filter(!has_cases), "5-14_n", "n_percent_pop"), get_val(res_table1 %>% filter(!has_cases), "15-49_n", "n_percent_pop"), get_val(res_table1 %>% filter(!has_cases), "≥50_n", "n_percent_pop"), get_val(res_table1 %>% filter(!has_cases), "Class I_n", "n_percent_pop"), get_val(res_table1 %>% filter(!has_cases), "Class II_n", "n_percent_pop"), get_val(res_table1 %>% filter(!has_cases), "Class III_n", "n_percent_pop"), get_val(res_table1 %>% filter(!has_cases), "Class IV_n", "n_percent_pop"), get_val(res_table1 %>% filter(!has_cases), "Class V_n", "n_percent_pop"), get_val(res_table1 %>% filter(!has_cases), "Improved_n_water", "n_percent_pop"), get_val(res_table1 %>% filter(!has_cases), "Improved_n_toilet", "n_percent_pop"), "-", "0", get_val(res_table1 %>% filter(!has_cases), "total_samples", "sum_n_percent", overall_samples), get_val(res_table1 %>% filter(!has_cases), "positive_samples", "n_percent_samples"), get_val(res_table1 %>% filter(!has_cases), "hf183", "median_iqr"), get_val(res_table1 %>% filter(!has_cases), "rainfall", "median_iqr"), get_val(res_table1 %>% filter(!has_cases), "temp", "median_iqr"), get_val(res_table1 %>% filter(!has_cases), "flow", "median_iqr")),
  P_Value = c(get_p_value("Pop", "continuous_site"), get_p_value("0-4_n", "chisq_pop"), get_p_value("5-14_n", "chisq_pop"), get_p_value("15-49_n", "chisq_pop"), get_p_value("≥50_n", "chisq_pop"), get_p_value("Class I_n", "chisq_pop"), get_p_value("Class II_n", "chisq_pop"), get_p_value("Class III_n", "chisq_pop"), get_p_value("Class IV_n", "chisq_pop"), get_p_value("Class V_n", "chisq_pop"), get_p_value("Improved_n_water", "chisq_pop"), get_p_value("Improved_n_toilet", "chisq_pop"), "-", "-", "-", get_p_value("positive_samples", "fisher_samples"), get_p_value("hf183", "continuous"), get_p_value("rainfall", "continuous"), get_p_value("temp", "continuous"), get_p_value("flow", "continuous"))
)

print(overall_summary)

#### 4. STATISTICAL MODELS ####
# Main Forward Model: Phage Positivity ~ Lagged Cases
model_adj <- glmmTMB(cbind(positive_samples, total_samples - positive_samples) ~ cases_lag2 + flow + rainfall + temp + hf183 + pct_under_15 + pct_water_improved + pct_toilet_improved + pct_high_ses + (1 | Site), family = betabinomial, data = analysis_df)
summary(model_adj)

# IRR Model: Phage Counts ~ Lagged Cases
model_zip <- glmmTMB(total_plaques ~ cases_lag2 + flow + rainfall + temp + hf183 + pct_under_15 + pct_water_improved + pct_toilet_improved + pct_high_ses + (1 | Site), ziformula = ~1, family = poisson, data = analysis_df)
summary(model_zip)

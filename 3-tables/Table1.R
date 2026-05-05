#### TYPHOID R01 — TABLE 1: DESCRIPTIVE STATISTICS ####
# Last Updated: May 2026
# Purpose: Generate Table 1 comparing site-months with vs. without clinical typhoid cases
# Requires: cleaned_analysis_data.RData (run 1-data-cleaning/DataCleaning.R first)
# Output:   results/Table1.csv

library(tidyverse)
library(lubridate)

# ── Directories ────────────────────────────────────────────────────────────────
base_dir   <- "/Users/esthjung/Library/CloudStorage/GoogleDrive-esthjung@stanford.edu/My Drive/Typhoid R01"
data_dir   <- file.path(base_dir, "Data")
results_dir <- file.path(base_dir, "Typhoid-R01-India/results")
dir.create(results_dir, showWarnings = FALSE)

#### 1. LOAD DATA ####
cat("Loading cleaned data...\n")
rdata_path <- file.path(data_dir, "cleaned_analysis_data.RData")
if (!file.exists(rdata_path)) stop("Run 1-data-cleaning/DataCleaning.R first to generate: ", rdata_path)
load(rdata_path)  # loads 'res_intermediate'

# Map column names
merged_data <- res_intermediate %>%
  mutate(
    temp     = avg_temp,
    rainfall = total_rain,
    flow     = avg_flow,
    hf183    = hf183_abs
  )

# Define site-level case status (ever positive site)
site_status <- merged_data %>%
  group_by(Site) %>%
  summarise(has_cases = any(clinical_cases > 0, na.rm = TRUE), .groups = "drop")

merged_data <- merged_data %>%
  dplyr::select(-any_of("has_cases")) %>%
  inner_join(site_status, by = "Site")

# Exclude burn-in (need 2-month lag)
res_table1 <- merged_data %>% filter(!is.na(cases_lag2))

cat("Table 1 dataset: ", nrow(res_table1), "rows,", n_distinct(res_table1$Site), "sites\n")
cat("  Sites with cases:    ", sum(site_status$has_cases), "\n")
cat("  Sites without cases: ", sum(!site_status$has_cases), "\n\n")

#### 2. HELPER FUNCTIONS ####

get_val <- function(df, var, type = "n_percent", overall_val = NULL) {
  if (nrow(df) == 0) return("0")
  if (type == "sum") return(format(sum(df[[var]], na.rm = TRUE), big.mark = ","))
  if (type == "median_iqr") {
    return(paste0(
      round(median(df[[var]], na.rm = TRUE), 1), " (",
      round(quantile(df[[var]], 0.25, na.rm = TRUE), 1), "-",
      round(quantile(df[[var]], 0.75, na.rm = TRUE), 1), ")"
    ))
  }
  if (type == "n_percent_pop") {
    total_var <- case_when(
      grepl("_n$|years$", var)  ~ "Total_n",
      grepl("water", var)       ~ "Total_n_water",
      grepl("toilet", var)      ~ "Total_n_toilet",
      grepl("Class", var)       ~ "Total_n_ses",
      TRUE                      ~ "Total_n"
    )
    site_vals <- df %>%
      group_by(Site) %>%
      summarise(v = first(!!sym(var)), t = first(!!sym(total_var)), .groups = "drop")
    v_sum <- sum(site_vals$v, na.rm = TRUE)
    t_sum <- sum(site_vals$t, na.rm = TRUE)
    return(paste0(format(v_sum, big.mark = ","), " (", round(100 * v_sum / t_sum, 1), "%)"))
  }
  if (type == "sum_n_percent") {
    n <- if (var == "Pop") {
      sum((df %>% group_by(Site) %>% summarise(v = first(Pop)))$v)
    } else {
      sum(df[[var]], na.rm = TRUE)
    }
    pct <- round(100 * n / as.numeric(overall_val), 1)
    return(paste0(format(n, big.mark = ","), " (", pct, "%)"))
  }
  if (type == "n_percent_samples") {
    pos <- sum(df$positive_samples, na.rm = TRUE)
    tot <- sum(df$total_samples, na.rm = TRUE)
    return(paste0(pos, " (", round(100 * pos / tot, 1), "%)"))
  }
  return("")
}

get_p_value <- function(var_name, data_type = "continuous") {
  tryCatch({
    if (data_type == "continuous") {
      return(format.pval(wilcox.test(res_table1[[var_name]] ~ res_table1$has_cases)$p.value, digits = 3))
    }
    if (data_type == "continuous_site") {
      site_data <- res_table1 %>%
        group_by(Site) %>%
        summarise(v = first(!!sym(var_name)), hc = first(has_cases), .groups = "drop")
      return(format.pval(wilcox.test(site_data$v ~ site_data$hc)$p.value, digits = 3))
    }
    if (data_type == "chisq_pop") {
      total_var <- case_when(
        grepl("_n$|years$", var_name) ~ "Total_n",
        grepl("water", var_name)      ~ "Total_n_water",
        grepl("toilet", var_name)     ~ "Total_n_toilet",
        grepl("Class", var_name)      ~ "Total_n_ses",
        TRUE                          ~ "Total_n"
      )
      site_data <- res_table1 %>%
        group_by(Site) %>%
        summarise(v = first(!!sym(var_name)), t = first(!!sym(total_var)), hc = first(has_cases), .groups = "drop")
      tab <- matrix(c(
        sum(site_data$v[site_data$hc],  na.rm = TRUE),
        sum(site_data$t[site_data$hc],  na.rm = TRUE) - sum(site_data$v[site_data$hc],  na.rm = TRUE),
        sum(site_data$v[!site_data$hc], na.rm = TRUE),
        sum(site_data$t[!site_data$hc], na.rm = TRUE) - sum(site_data$v[!site_data$hc], na.rm = TRUE)
      ), nrow = 2)
      return(format.pval(fisher.test(tab)$p.value, digits = 3))
    }
    if (data_type == "fisher_samples") {
      pos_w  <- sum(res_table1$positive_samples[res_table1$has_cases],  na.rm = TRUE)
      tot_w  <- sum(res_table1$total_samples[res_table1$has_cases],     na.rm = TRUE)
      pos_wo <- sum(res_table1$positive_samples[!res_table1$has_cases], na.rm = TRUE)
      tot_wo <- sum(res_table1$total_samples[!res_table1$has_cases],    na.rm = TRUE)
      return(format.pval(fisher.test(matrix(c(pos_w, tot_w - pos_w, pos_wo, tot_wo - pos_wo), nrow = 2))$p.value, digits = 3))
    }
    return("-")
  }, error = function(e) "-")
}

#### 3. ASSEMBLE TABLE 1 ####
cat("Building Table 1...\n")

overall_pop     <- sum((res_table1 %>% group_by(Site) %>% summarise(v = first(Pop)))$v)
overall_samples <- sum(res_table1$total_samples)

overall_summary <- data.frame(
  Characteristic = c(
    "Total population, n (%)",
    "  0-4 years", "  5-14 years", "  15-49 years", "  ≥50 years",
    "SES Class 1 (≥Rs.9130)", "SES Class 2 (Rs.4565-9129)", "SES Class 3 (Rs.2739-4564)",
    "SES Class 4 (Rs.1369-2738)", "SES Class 5 (<Rs.1369)",
    "Improved drinking water", "Improved toilet",
    "Median population density (IQR)",
    "Total clinical typhoid cases, n",
    "Total wastewater samples, n (%)",
    "Positive wastewater samples, n (%)",
    "Median HF183 ACN (IQR)",
    "Median cumulative rainfall, mm (IQR)",
    "Median temperature, °C (IQR)",
    "Median flow rate, m³/h (IQR)"
  ),
  Overall = c(
    format(overall_pop, big.mark = ","),
    get_val(res_table1, "0-4_n",          "n_percent_pop"),
    get_val(res_table1, "5-14_n",         "n_percent_pop"),
    get_val(res_table1, "15-49_n",        "n_percent_pop"),
    get_val(res_table1, "≥50_n",          "n_percent_pop"),
    get_val(res_table1, "Class I_n",      "n_percent_pop"),
    get_val(res_table1, "Class II_n",     "n_percent_pop"),
    get_val(res_table1, "Class III_n",    "n_percent_pop"),
    get_val(res_table1, "Class IV_n",     "n_percent_pop"),
    get_val(res_table1, "Class V_n",      "n_percent_pop"),
    get_val(res_table1, "Improved_n_water",  "n_percent_pop"),
    get_val(res_table1, "Improved_n_toilet", "n_percent_pop"),
    "-",
    get_val(res_table1, "clinical_cases", "sum"),
    format(overall_samples, big.mark = ","),
    get_val(res_table1, "positive_samples", "n_percent_samples"),
    get_val(res_table1, "hf183",    "median_iqr"),
    get_val(res_table1, "rainfall", "median_iqr"),
    get_val(res_table1, "temp",     "median_iqr"),
    get_val(res_table1, "flow",     "median_iqr")
  ),
  With_Cases = c(
    get_val(res_table1 %>% filter(has_cases), "Pop",              "sum_n_percent", overall_pop),
    get_val(res_table1 %>% filter(has_cases), "0-4_n",           "n_percent_pop"),
    get_val(res_table1 %>% filter(has_cases), "5-14_n",          "n_percent_pop"),
    get_val(res_table1 %>% filter(has_cases), "15-49_n",         "n_percent_pop"),
    get_val(res_table1 %>% filter(has_cases), "≥50_n",           "n_percent_pop"),
    get_val(res_table1 %>% filter(has_cases), "Class I_n",       "n_percent_pop"),
    get_val(res_table1 %>% filter(has_cases), "Class II_n",      "n_percent_pop"),
    get_val(res_table1 %>% filter(has_cases), "Class III_n",     "n_percent_pop"),
    get_val(res_table1 %>% filter(has_cases), "Class IV_n",      "n_percent_pop"),
    get_val(res_table1 %>% filter(has_cases), "Class V_n",       "n_percent_pop"),
    get_val(res_table1 %>% filter(has_cases), "Improved_n_water",  "n_percent_pop"),
    get_val(res_table1 %>% filter(has_cases), "Improved_n_toilet", "n_percent_pop"),
    "-",
    get_val(res_table1 %>% filter(has_cases), "clinical_cases",    "sum"),
    get_val(res_table1 %>% filter(has_cases), "total_samples",     "sum_n_percent", overall_samples),
    get_val(res_table1 %>% filter(has_cases), "positive_samples",  "n_percent_samples"),
    get_val(res_table1 %>% filter(has_cases), "hf183",    "median_iqr"),
    get_val(res_table1 %>% filter(has_cases), "rainfall", "median_iqr"),
    get_val(res_table1 %>% filter(has_cases), "temp",     "median_iqr"),
    get_val(res_table1 %>% filter(has_cases), "flow",     "median_iqr")
  ),
  Without_Cases = c(
    get_val(res_table1 %>% filter(!has_cases), "Pop",             "sum_n_percent", overall_pop),
    get_val(res_table1 %>% filter(!has_cases), "0-4_n",          "n_percent_pop"),
    get_val(res_table1 %>% filter(!has_cases), "5-14_n",         "n_percent_pop"),
    get_val(res_table1 %>% filter(!has_cases), "15-49_n",        "n_percent_pop"),
    get_val(res_table1 %>% filter(!has_cases), "≥50_n",          "n_percent_pop"),
    get_val(res_table1 %>% filter(!has_cases), "Class I_n",      "n_percent_pop"),
    get_val(res_table1 %>% filter(!has_cases), "Class II_n",     "n_percent_pop"),
    get_val(res_table1 %>% filter(!has_cases), "Class III_n",    "n_percent_pop"),
    get_val(res_table1 %>% filter(!has_cases), "Class IV_n",     "n_percent_pop"),
    get_val(res_table1 %>% filter(!has_cases), "Class V_n",      "n_percent_pop"),
    get_val(res_table1 %>% filter(!has_cases), "Improved_n_water",  "n_percent_pop"),
    get_val(res_table1 %>% filter(!has_cases), "Improved_n_toilet", "n_percent_pop"),
    "-",
    "0",
    get_val(res_table1 %>% filter(!has_cases), "total_samples",    "sum_n_percent", overall_samples),
    get_val(res_table1 %>% filter(!has_cases), "positive_samples", "n_percent_samples"),
    get_val(res_table1 %>% filter(!has_cases), "hf183",    "median_iqr"),
    get_val(res_table1 %>% filter(!has_cases), "rainfall", "median_iqr"),
    get_val(res_table1 %>% filter(!has_cases), "temp",     "median_iqr"),
    get_val(res_table1 %>% filter(!has_cases), "flow",     "median_iqr")
  ),
  P_Value = c(
    get_p_value("Pop",              "continuous_site"),
    get_p_value("0-4_n",           "chisq_pop"),
    get_p_value("5-14_n",          "chisq_pop"),
    get_p_value("15-49_n",         "chisq_pop"),
    get_p_value("≥50_n",           "chisq_pop"),
    get_p_value("Class I_n",       "chisq_pop"),
    get_p_value("Class II_n",      "chisq_pop"),
    get_p_value("Class III_n",     "chisq_pop"),
    get_p_value("Class IV_n",      "chisq_pop"),
    get_p_value("Class V_n",       "chisq_pop"),
    get_p_value("Improved_n_water",  "chisq_pop"),
    get_p_value("Improved_n_toilet", "chisq_pop"),
    "-", "-", "-",
    get_p_value("positive_samples", "fisher_samples"),
    get_p_value("hf183",    "continuous"),
    get_p_value("rainfall", "continuous"),
    get_p_value("temp",     "continuous"),
    get_p_value("flow",     "continuous")
  ),
  stringsAsFactors = FALSE
)

#### 4. PRINT & SAVE ####
print(overall_summary)

output_path <- file.path(results_dir, "Table1.csv")
write.csv(overall_summary, output_path, row.names = FALSE)
cat("\nSUCCESS: Table 1 saved to:", output_path, "\n")

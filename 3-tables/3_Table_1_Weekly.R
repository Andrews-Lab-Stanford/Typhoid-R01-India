#### TYPHOID R01 — STEP 3: WEEKLY TABLE 1 (DESCRIPTIVE STATISTICS) ####
# Last Updated: May 2026
# Purpose: Generate Table 1 comparing sites with vs. without typhoid cases based on the WEEKLY analytic set.
#
# INPUT:  Data/analysis_ready_weekly.RData  (run 1-data-cleaning/1_Data_Cleaning_Weekly.R first)
# OUTPUT: results/Table1_Weekly.csv

library(tidyverse)
library(readxl)

# ── Directories ────────────────────────────────────────────────────────────────
base_dir    <- "/Users/esthjung/Library/CloudStorage/GoogleDrive-esthjung@stanford.edu/My Drive/Typhoid R01"
data_dir    <- file.path(base_dir, "Data")
repo_dir    <- file.path(base_dir, "Typhoid-R01-India")
results_dir <- file.path(repo_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load Data ──────────────────────────────────────────────────────────────────
rdata_path <- file.path(data_dir, "analysis_ready_weekly.RData")
if (!file.exists(rdata_path)) {
  stop("Run 1-data-cleaning/1_Data_Cleaning_Weekly.R first to generate: ", rdata_path)
}
load(rdata_path)   # loads: cases_weekly, merged_weekly, analysis_df_weekly
cat("Full clinical panel:", nrow(cases_weekly), "site-weeks\n")
cat("Merged (sampled weeks):", nrow(merged_weekly), "site-weeks\n")

# ── Build res_table: Full clinical panel + census vars (site-level) ────────────
# Pull census + has_cases from merged_weekly (one row per site)
site_meta <- merged_weekly %>%
  group_by(Site) %>%
  summarise(
    across(c(Pop, `0-4_n`, `5-14_n`, `15-49_n`, `≥50_n`,
             `Class I_n`, `Class II_n`, `Class III_n`, `Class IV_n`, `Class V_n`,
             Improved_n_water, Total_n_water,
             Improved_n_toilet, Total_n_toilet,
             Total_n_ses, has_cases),
           first),
    .groups = "drop"
  )

# res_table = every week × site from the full clinical panel, with census joined
res_table <- cases_weekly %>%
  inner_join(site_meta, by = "Site")

# env_table = only sampled weeks (unscaled) for environmental/wastewater rows
env_table <- merged_weekly

cat("res_table rows:", nrow(res_table), "\n")
cat("Total clinical cases in res_table:", sum(res_table$clinical_cases, na.rm = TRUE), "\n")

# ── Helper: get_val() ──────────────────────────────────────────────────────────
get_val <- function(df, var, type = "n_percent", overall_val = NULL) {
  if (nrow(df) == 0) return("0")

  if (type == "sum")
    return(format(sum(df[[var]], na.rm = TRUE), big.mark = ","))

  if (type == "median_iqr")
    return(paste0(
      round(median(df[[var]], na.rm = TRUE), 1), " (",
      round(quantile(df[[var]], 0.25, na.rm = TRUE), 1), "-",
      round(quantile(df[[var]], 0.75, na.rm = TRUE), 1), ")"
    ))

  if (type == "n_percent_pop") {
    total_var <- case_when(
      grepl("_n$|years$|^.0-4_n", var) ~ "Pop",
      grepl("water",  var)     ~ "Total_n_water",
      grepl("toilet", var)     ~ "Total_n_toilet",
      grepl("Class",  var)     ~ "Total_n_ses",
      TRUE                     ~ "Pop"
    )
    sv  <- df %>% group_by(Site) %>%
           summarise(v = first(!!sym(var)), t = first(!!sym(total_var)), .groups = "drop")
    vs  <- sum(sv$v, na.rm = TRUE); ts <- sum(sv$t, na.rm = TRUE)
    return(paste0(format(vs, big.mark = ","), " (", round(100 * vs / ts, 1), "%)"))
  }

  if (type == "sum_n_percent") {
    n <- if (var == "Pop")
      sum((df %>% group_by(Site) %>% summarise(v = first(Pop)))$v)
    else
      sum(df[[var]], na.rm = TRUE)
    pct <- round(100 * n / as.numeric(overall_val), 1)
    return(paste0(format(n, big.mark = ","), " (", pct, "%)"))
  }

  if (type == "n_percent_samples") {
    pos <- sum(df$positive_samples, na.rm = TRUE)
    tot <- sum(df$total_samples,    na.rm = TRUE)
    return(paste0(pos, " (", round(100 * pos / tot, 1), "%)"))
  }

  return("")
}

# ── Helper: get_p_value() ──────────────────────────────────────────────────────
get_p_value <- function(var_name, data_type = "continuous") {
  # Use env_table for environmental/wastewater vars, res_table for the rest
  df_use <- if (var_name %in% c("total_samples", "positive_samples", "hf183", "rainfall", "temp", "flow")) env_table else res_table
  
  tryCatch({
    if (data_type == "continuous")
      return(format.pval(
        wilcox.test(df_use[[var_name]] ~ df_use$has_cases)$p.value, digits = 3))

    if (data_type == "continuous_site") {
      sd <- df_use %>% group_by(Site) %>%
            summarise(v = first(!!sym(var_name)), hc = first(has_cases), .groups = "drop")
      return(format.pval(wilcox.test(sd$v ~ sd$hc)$p.value, digits = 3))
    }

    if (data_type == "chisq_pop") {
      total_var <- case_when(
        grepl("_n$|years$|^.0-4_n", var_name) ~ "Pop",
        grepl("water",  var_name)     ~ "Total_n_water",
        grepl("toilet", var_name)     ~ "Total_n_toilet",
        grepl("Class",  var_name)     ~ "Total_n_ses",
        TRUE                          ~ "Pop"
      )
      sd  <- df_use %>% group_by(Site) %>%
             summarise(v = first(!!sym(var_name)), t = first(!!sym(total_var)),
                       hc = first(has_cases), .groups = "drop")
      tab <- matrix(c(
        sum(sd$v[ sd$hc], na.rm = TRUE), sum(sd$t[ sd$hc], na.rm = TRUE) - sum(sd$v[ sd$hc], na.rm = TRUE),
        sum(sd$v[!sd$hc], na.rm = TRUE), sum(sd$t[!sd$hc], na.rm = TRUE) - sum(sd$v[!sd$hc], na.rm = TRUE)
      ), nrow = 2)
      return(format.pval(fisher.test(tab)$p.value, digits = 3))
    }

    if (data_type == "fisher_samples") {
      # This always uses env_table
      pw  <- sum(env_table$positive_samples[ env_table$has_cases], na.rm = TRUE)
      tw  <- sum(env_table$total_samples[    env_table$has_cases], na.rm = TRUE)
      pwo <- sum(env_table$positive_samples[!env_table$has_cases], na.rm = TRUE)
      two <- sum(env_table$total_samples[   !env_table$has_cases], na.rm = TRUE)
      return(format.pval(
        fisher.test(matrix(c(pw, tw - pw, pwo, two - pwo), nrow = 2))$p.value, digits = 3))
    }
    return("-")
  }, error = function(e) "-")
}

# ── Totals ─────────────────────────────────────────────────────────────────────
overall_pop     <- sum((res_table %>% group_by(Site) %>% summarise(v = first(Pop)))$v)
overall_samples <- sum(env_table$total_samples)

# ── Build Table 1 ──────────────────────────────────────────────────────────────
vars <- list(
  list("Number of sites, n",                 "Site",              "n_distinct",       "-"),
  list("Total population, n (%)",            "Pop",               "sum_n_percent",    "continuous_site"),
  list("  0-4 years, n (%)",                 "0-4_n",             "n_percent_pop",    "chisq_pop"),
  list("  5-14 years, n (%)",                "5-14_n",            "n_percent_pop",    "chisq_pop"),
  list("  15-49 years, n (%)",               "15-49_n",           "n_percent_pop",    "chisq_pop"),
  list("  \u226550 years, n (%)",            "≥50_n",             "n_percent_pop",    "chisq_pop"),
  list("SES Class 1 (\u2265Rs.9130), n (%)", "Class I_n",         "n_percent_pop",    "chisq_pop"),
  list("SES Class 2 (Rs.4565-9129), n (%)",  "Class II_n",        "n_percent_pop",    "chisq_pop"),
  list("SES Class 3 (Rs.2739-4564), n (%)",  "Class III_n",       "n_percent_pop",    "chisq_pop"),
  list("SES Class 4 (Rs.1369-2738), n (%)",  "Class IV_n",        "n_percent_pop",    "chisq_pop"),
  list("SES Class 5 (<Rs.1369), n (%)",      "Class V_n",         "n_percent_pop",    "chisq_pop"),
  list("Improved drinking water, n (%)",      "Improved_n_water",  "n_percent_pop",    "chisq_pop"),
  list("Improved toilet, n (%)",             "Improved_n_toilet", "n_percent_pop",    "chisq_pop"),
  list("Total clinical typhoid cases, n",    "clinical_cases",    "sum",              "-"),
  list("Total wastewater samples, n (%)",    "total_samples",     "sum_n_percent",    "-"),
  list("Positive wastewater samples, n (%)", "positive_samples",  "n_percent_samples","fisher_samples"),
  list("Median HF183 ACN (IQR)",             "hf183",             "median_iqr",       "continuous"),
  list("Median rainfall, mm (IQR)",          "rainfall",          "median_iqr",       "continuous"),
  list("Median temperature, \u00b0C (IQR)",  "temp",              "median_iqr",       "continuous"),
  list("Median flow rate, m\u00b3/h (IQR)",  "flow",              "median_iqr",       "continuous")
)

build_col <- function(df_full, df_env, overall_val_pop, overall_val_samp) {
  sapply(vars, function(v) {
    var   <- v[[2]]; type <- v[[3]]
    if (type == "n_distinct") return(as.character(n_distinct(df_full$Site)))
    
    # Decide which dataset to use
    df_use <- if (var %in% c("total_samples", "positive_samples", "hf183", "rainfall", "temp", "flow")) df_env else df_full
    
    ov <- if (var == "Pop") overall_val_pop else if (var == "total_samples") overall_val_samp else NULL
    tryCatch(get_val(df_use, var, type, ov), error = function(e) "-")
  })
}

overall_summary <- data.frame(
  Characteristic = sapply(vars, `[[`, 1),
  Overall        = build_col(res_table, env_table,
                             overall_pop, overall_samples),
  With_Cases     = build_col(res_table %>% filter( has_cases), env_table %>% filter( has_cases),
                             overall_pop, overall_samples),
  Without_Cases  = build_col(res_table %>% filter(!has_cases), env_table %>% filter(!has_cases),
                             overall_pop, overall_samples),
  P_Value        = sapply(vars, function(v) {
    if (v[[4]] == "-") return("-")
    get_p_value(v[[2]], v[[4]])
  }),
  stringsAsFactors = FALSE
)

# Sites without cases have 0 clinical cases by definition
# (This is already handled by the logic above, but for total clarity:)
# overall_summary$Without_Cases[overall_summary$Characteristic == "Total clinical typhoid cases, n"] <- "0"

print(overall_summary, right = FALSE)

# ── Save ───────────────────────────────────────────────────────────────────────
out <- file.path(results_dir, "Table1_Weekly.csv")
write.csv(overall_summary, out, row.names = FALSE)
cat("\nSUCCESS: Weekly Table 1 saved to", out, "\n")

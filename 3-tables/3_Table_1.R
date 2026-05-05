#### TYPHOID R01 — STEP 3: TABLE 1 (DESCRIPTIVE STATISTICS) ####
# Last Updated: May 2026
# Purpose: Generate Table 1 comparing sites with vs. without typhoid cases.
#
# INPUT:  Data/analysis_ready.RData  (run 1-data-cleaning/1_Data_Cleaning.R first)
# OUTPUT: results/Table1.csv

library(tidyverse)

# ── Directories ────────────────────────────────────────────────────────────────
base_dir    <- "/Users/esthjung/Library/CloudStorage/GoogleDrive-esthjung@stanford.edu/My Drive/Typhoid R01"
data_dir    <- file.path(base_dir, "Data")
results_dir <- file.path(base_dir, "Typhoid-R01-India/results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load Data ──────────────────────────────────────────────────────────────────
rdata_path <- file.path(data_dir, "analysis_ready.RData")
if (!file.exists(rdata_path)) {
  stop("Run 1-data-cleaning/1_Data_Cleaning.R first to generate: ", rdata_path)
}
load(rdata_path)   # loads: merged_data, res_table1, analysis_df
cat("Loaded res_table1:", nrow(res_table1), "site-months\n")

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
      grepl("_n$|years$", var) ~ "Total_n",
      grepl("water",  var)     ~ "Total_n_water",
      grepl("toilet", var)     ~ "Total_n_toilet",
      grepl("Class",  var)     ~ "Total_n_ses",
      TRUE                     ~ "Total_n"
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
  tryCatch({
    if (data_type == "continuous")
      return(format.pval(
        wilcox.test(res_table1[[var_name]] ~ res_table1$has_cases)$p.value, digits = 3))

    if (data_type == "continuous_site") {
      sd <- res_table1 %>% group_by(Site) %>%
            summarise(v = first(!!sym(var_name)), hc = first(has_cases), .groups = "drop")
      return(format.pval(wilcox.test(sd$v ~ sd$hc)$p.value, digits = 3))
    }

    if (data_type == "chisq_pop") {
      total_var <- case_when(
        grepl("_n$|years$", var_name) ~ "Total_n",
        grepl("water",  var_name)     ~ "Total_n_water",
        grepl("toilet", var_name)     ~ "Total_n_toilet",
        grepl("Class",  var_name)     ~ "Total_n_ses",
        TRUE                          ~ "Total_n"
      )
      sd  <- res_table1 %>% group_by(Site) %>%
             summarise(v = first(!!sym(var_name)), t = first(!!sym(total_var)),
                       hc = first(has_cases), .groups = "drop")
      tab <- matrix(c(
        sum(sd$v[ sd$hc], na.rm = TRUE), sum(sd$t[ sd$hc], na.rm = TRUE) - sum(sd$v[ sd$hc], na.rm = TRUE),
        sum(sd$v[!sd$hc], na.rm = TRUE), sum(sd$t[!sd$hc], na.rm = TRUE) - sum(sd$v[!sd$hc], na.rm = TRUE)
      ), nrow = 2)
      return(format.pval(fisher.test(tab)$p.value, digits = 3))
    }

    if (data_type == "fisher_samples") {
      pw  <- sum(res_table1$positive_samples[ res_table1$has_cases], na.rm = TRUE)
      tw  <- sum(res_table1$total_samples[    res_table1$has_cases], na.rm = TRUE)
      pwo <- sum(res_table1$positive_samples[!res_table1$has_cases], na.rm = TRUE)
      two <- sum(res_table1$total_samples[   !res_table1$has_cases], na.rm = TRUE)
      return(format.pval(
        fisher.test(matrix(c(pw, tw - pw, pwo, two - pwo), nrow = 2))$p.value, digits = 3))
    }
    return("-")
  }, error = function(e) "-")
}

# ── Totals ─────────────────────────────────────────────────────────────────────
overall_pop     <- sum((res_table1 %>% group_by(Site) %>% summarise(v = first(Pop)))$v)
overall_samples <- sum(res_table1$total_samples)

# ── Build Table 1 ──────────────────────────────────────────────────────────────
vars <- list(
  list("Total population, n (%)",            "Pop",               "sum_n_percent",    "continuous_site"),
  list("  0-4 years, n (%)",                 "0-4_n",             "n_percent_pop",    "chisq_pop"),
  list("  5-14 years, n (%)",                "5-14_n",            "n_percent_pop",    "chisq_pop"),
  list("  15-49 years, n (%)",               "15-49_n",           "n_percent_pop",    "chisq_pop"),
  list("  \u226550 years, n (%)",            "\u226550_n",        "n_percent_pop",    "chisq_pop"),
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

build_col <- function(df, overall_val_pop, overall_val_samp) {
  sapply(vars, function(v) {
    var   <- v[[2]]; type <- v[[3]]
    ov    <- if (var == "Pop") overall_val_pop else if (var == "total_samples") overall_val_samp else NULL
    tryCatch(get_val(df, var, type, ov), error = function(e) "-")
  })
}

overall_summary <- data.frame(
  Characteristic = sapply(vars, `[[`, 1),
  Overall        = build_col(res_table1,
                             overall_pop, overall_samples),
  With_Cases     = build_col(res_table1 %>% filter( has_cases),
                             overall_pop, overall_samples),
  Without_Cases  = build_col(res_table1 %>% filter(!has_cases),
                             overall_pop, overall_samples),
  P_Value        = sapply(vars, function(v) {
    if (v[[4]] == "-") return("-")
    get_p_value(v[[2]], v[[4]])
  }),
  stringsAsFactors = FALSE
)

# Fix "0" for clinical cases in Without_Cases column (sites without cases = 0 by definition)
overall_summary$Without_Cases[overall_summary$Characteristic == "Total clinical typhoid cases, n"] <- "0"

print(overall_summary, right = FALSE)

# ── Save ───────────────────────────────────────────────────────────────────────
out <- file.path(results_dir, "Table1.csv")
write.csv(overall_summary, out, row.names = FALSE)
cat("\nSUCCESS: Table 1 saved to", out, "\n")

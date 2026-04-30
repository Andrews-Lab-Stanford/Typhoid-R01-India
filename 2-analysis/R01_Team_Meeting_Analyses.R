#### TYPHOID R01 PHAGE ASSAYS ####
# User: Esther Jung
# Last Updated: February 19, 2026

rm(list = ls())

#### CONFIGURE ####
# Core data manipulation and visualization
library(tidyverse) # Includes dplyr, ggplot2, readr, etc.
library(readxl) # Excel file reading
library(lubridate) # Date handling
library(tidyselect) # Column selection helpers

# Visualization and plotting
library(ggplot2)
library(ggpubr) # Publication-ready plots
library(patchwork) # Plot arrangement
library(gridExtra)
library(grid)

# Mapping
library(leaflet)
library(leafsync) # Synced leaflet maps
library(htmlwidgets)
library(htmltools)

# Statistical analysis
library(gmodels)
library(MASS) # Negative binomial models
library(pscl) # Zero-inflated models

# Other utilities
library(here) # Project paths
library(gt) # Tables
library(REDCapR) # REDCap data access
library(devtools)
library(anytime)
library(shiny)
library(DiagrammeR)
library(rsconnect)

setwd("~/Google Drive/My Drive/Typhoid R01/Data")
sites <- read.csv("Site_54.csv") # Maps the catchments out
phage <- read.csv("REDCapPhageData.csv") # REDCap phage assay data until January 2026
pcr <- read.csv("typhi_sitePos_20251231.csv") # PCR by site until 12/31/2025
monthly_cases <- read_excel("Site_monthlytyphoidcases_FINAL.xlsx") # Clinical cases by site and month until January 2026

#### DATA CLEANING ============================================================

# Create case_site from monthly_cases: aggregate total cases by site
case_site <- monthly_cases %>%
  mutate(
    # Sum all monthly case columns
    total_cases = rowSums(across(matches("^(2024|2025|2026)[_-]")), na.rm = TRUE),
    case_positivity_percent = round(100 * total_cases / Pop, 2)
  ) %>%
  dplyr::select(Site, Pop, total_cases, case_positivity_percent) %>%
  arrange(desc(case_positivity_percent))

# Create monthly_cases_long for time series analysis
monthly_cases_long <- monthly_cases %>%
  pivot_longer(
    cols = matches("^(2024|2025|2026)[_-]"), # All month columns
    names_to = "month",
    values_to = "cases"
  ) %>%
  mutate(
    month = ym(month), # Convert to date for proper ordering
    positivity_percent = round((cases / Pop) * 100, 2)
  )

# Sites with cases
case_y_site <- case_site %>% filter(total_cases > 0)

names(pcr)

pcr <- pcr %>%
  rename(Site = site, Positive = MS, Unknown = NA.) %>%
  dplyr::select(c("Site", "Positive", "Unknown")) # PCR positivity

phage <- phage %>%
  mutate(result_status = case_when(
    amplification_tw == 1 | direct_tw == 1 ~ "Positive",
    amplification_tw == 2 & direct_tw == 2 ~ "Negative",
    (amplification_tw == 2 & is.na(direct_tw)) | (direct_tw == 2 & is.na(amplification_tw)) ~ "Negative",
    TRUE ~ "Unknown"
  ))
table(phage$result_status)

freq_table_df <- phage %>% # Phage positivity
  rename(Site = site_id_tw) %>%
  dplyr::select(Site, result_status) %>%
  group_by(Site) %>%
  summarise(
    Positive = sum(result_status == "Positive"),
    Total = n(),
    `Positivity %` = round(100 * (Positive / Total), 2),
    .groups = "drop"
  )

#### TABLES AND FIGURES ####
# Amplification = Enrichment assay | 1 = +, 2 = -
# Direct = Direct plating assay | 1 = +, 2 = -

# Phage positivity by month - clean date and status
phage_month <- phage %>%
  mutate(
    date_sample = as.Date(ifelse(
      substr(date_results_tw, 1, 2) %in% c("29", "22") | substr(date_results_tw, 1, 4) == "2004", # key to correct years
      paste0("2024", substr(date_results_tw, 5, 10)),
      date_results_tw
    ), format = "%Y-%m-%d"),
    month = floor_date(date_sample, "month"),
    status = case_when(
      amplification_tw == 1 | direct_tw == 1 ~ "Positive",
      amplification_tw == 2 & direct_tw == 2 ~ "Negative",
      (amplification_tw == 2 & is.na(direct_tw)) |
        (is.na(amplification_tw) & direct_tw == 2) ~ "Negative",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(status), !is.na(site_id_tw))

# Table 3: Monthly summary of phage positives
phage_summary <- phage_month %>%
  group_by(month_year = month) %>%
  summarise(
    `Direct assay` = sum(direct_tw == 1, na.rm = TRUE),
    `Enrichment assay` = sum(amplification_tw == 1, na.rm = TRUE),
    Total = n(),
    .groups = "drop"
  ) %>%
  complete(
    month_year = seq.Date(
      from = min(phage_month$month, na.rm = TRUE),
      to = max(phage_month$month, na.rm = TRUE),
      by = "month"
    ),
    fill = list(`Direct assay` = 0, `Enrichment assay` = 0, Total = 0)
  )
sum(phage_summary$Total)

# Figure 1: Total Typhi Phage Positives by Month
phage_plot_data <- phage_summary %>%
  mutate(
    month_year = as.Date(month_year),
    total_positive = `Direct assay` + `Enrichment assay`,
    positivity_pct = ifelse(Total > 0, round(100 * total_positive / Total, 1), 0)
  ) %>%
  filter(month_year >= as.Date("2024-03-01") & month_year <= as.Date("2026-01-31"))


# 1b. Proportion Plot (%)
phage_proportion_plot <- ggplot(phage_plot_data, aes(x = month_year, y = positivity_pct)) +
  geom_col(fill = "#7260a7") +
  geom_text(aes(label = paste0(positivity_pct)), vjust = -0.4, size = 3.5) +
  scale_x_date(breaks = sort(unique(phage_plot_data$month_year)), date_labels = "%b %Y") +
  scale_y_continuous(limits = c(0, max(phage_plot_data$positivity_pct) * 1.2)) +
  theme_pubr(base_size = 14) +
  labs(
    x = "Month",
    y = expression(italic("S.") ~ "Typhi Phage Positivity (%)")
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

phage_proportion_plot
ggsave("Monthly_Phage_Proportions.png", phage_proportion_plot, width = 12, height = 4, dpi = 300)

# Phage data already cleaned above in phage_month (includes month variable)

# Aggregate to site-month level
phage_by_site_month <- phage_month %>%
  group_by(Site = site_id_tw, month) %>%
  summarise(
    total_samples = n(),
    positive_samples = sum(status == "Positive", na.rm = TRUE),
    total_plaques = sum(direct_counts_tw, na.rm = TRUE),
    phage_positivity = round(100 * positive_samples / total_samples, 2),
    .groups = "drop"
  )

# Aggregate to site level
phage_site <- phage_month %>%
  group_by(Site = site_id_tw) %>%
  summarise(
    total_samples = n(),
    positive_samples = sum(status == "Positive", na.rm = TRUE),
    total_plaques = sum(direct_counts_tw, na.rm = TRUE),
    phage_positivity = round(100 * positive_samples / total_samples, 2),
    .groups = "drop"
  )

# Final dataset
merged_df <- phage_by_site_month %>%
  left_join(monthly_cases_long, by = c("Site", "month")) %>%
  filter(!is.na(cases), !is.na(total_samples))
length(unique(merged_df$Site)) # should be 54 sites

#### Median positivity comparison ==============================================
data <- left_join(case_site, phage_site, by = c("Site"))

# Create a grouping variable based on whether there are any cases
data <- data %>%
  mutate(case_group = ifelse(total_cases > 0, "Typhoid cases", "No typhoid cases"))

data <- data %>%
  mutate(
    case_category = case_when(
      total_cases == 0 ~ "0 cases",
      total_cases >= 1 & total_cases <= 9 ~ "1-9 cases",
      total_cases >= 10 ~ "≥10 cases"
    ),
    # Make it a factor with ordered levels
    case_category = factor(case_category,
      levels = c("0 cases", "1-9 cases", "≥10 cases")
    )
  )

# Calculate positivity percentages for Phage
data <- data %>%
  mutate(phage_positivity_pct = (positive_samples / total_samples) * 100)

# Plot 1: Phage Positivity (Binary: Cases vs No Cases)
p_binary <- ggplot(data, aes(x = case_group, y = phage_positivity_pct, fill = case_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  theme_pubr(base_size = 14) +
  scale_fill_manual(values = c("Typhoid cases" = "#4361EE", "No typhoid cases" = "#BDBDBD")) +
  labs(
    x = "Catchment case status",
    y = expression(italic("S.") ~ "Typhi Phage Positivity (%)"),
    fill = "Site Group"
  ) +
  theme(legend.position = "none")

# Plot 2: Phage Positivity (Categorical breakdown)
p_categorical <- ggplot(data, aes(x = case_category, y = phage_positivity_pct, fill = case_category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  theme_pubr(base_size = 14) +
  scale_fill_manual(values = c("0 cases" = "#BDBDBD", "1-9 cases" = "#21908C", "≥10 cases" = "#440154")) +
  labs(
    x = "Catchment case count",
    y = expression(italic("S.") ~ "Typhi Phage Positivity (%)"),
    fill = "Burden"
  ) +
  theme(legend.position = "none")

# Arrange plots for poster
combined_boxplot_fig <- p_binary + p_categorical +
  plot_layout(widths = c(1, 1.5))

combined_boxplot_fig

ggsave("Phage_Positivity_Comparison.png", combined_boxplot_fig, width = 12, height = 4)

# Summary statistics for reference

    # Phage Positivity Summary (Binary)
    print(data %>%
      group_by(case_group) %>%
      summarise(
        n_sites = n(),
        median_pct = round(median(phage_positivity_pct, na.rm = TRUE), 1),
        mean_pct = round(mean(phage_positivity_pct, na.rm = TRUE), 1)
      ))

    #Phage Positivity Summary (Categorical Bins) 
    phage_summary_stats <- data %>%
      group_by(case_category) %>%
      summarise(
        n_sites = n(),
        median_pct = round(median(phage_positivity_pct, na.rm = TRUE), 1),
        mean_pct = round(mean(phage_positivity_pct, na.rm = TRUE), 1)
      )
    print(phage_summary_stats)

    # Exact p-values
      p_binary_val <- wilcox.test(phage_positivity_pct ~ case_group, data = data, exact = FALSE)$p.value
      cat("Binary (Cases vs No Cases): p =", format.pval(p_binary_val, digits = 3), "\n")

    # Kruskal-Wallis (Across all categories)
      k_val <- kruskal.test(phage_positivity_pct ~ case_category, data = data)$p.value
      cat("Kruskal-Wallis (Categorical): p =", format.pval(k_val, digits = 3), "\n")

#### MAIN MODEL ================================================================

# FORWARD DIRECTION: Do previous cases predict current phage positivity?
# Add lagged case variables for time series analysis
# Note: Cases are lagged (past), phage positivity is current (present)
merged_df <- merged_df %>%
  arrange(Site, month) %>%
  group_by(Site) %>%
  mutate(
    cases_lag1 = lag(cases, 1), # Cases from 1 month ago
    cases_lag2 = lag(cases, 2), # Cases from 2 months ago
    cases_lag3 = lag(cases, 3) # Cases from 3 months ago
  ) %>%
  ungroup()

# Model 1: Does presence of cases 1 month ago predict current phage positivity?
model_lag1 <- glm(
  cbind(positive_samples, total_samples - positive_samples) ~
    I(cases_lag1 > 0) + as.factor(Site),
  data = merged_df,
  family = binomial
)
exp(coef(model_lag1)) # Odds ratios
exp(confint(model_lag1)) # Confidence intervals

# Model 2: Does presence of cases 2 months ago predict current phage positivity?
model_lag2 <- glm(
  cbind(positive_samples, total_samples - positive_samples) ~
    I(cases_lag2 > 0) + as.factor(Site),
  data = merged_df,
  family = binomial
)
exp(coef(model_lag2))
exp(confint(model_lag2))

# Model 3: Does presence of cases 3 months ago predict current phage positivity?
model_lag3 <- glm(
  cbind(positive_samples, total_samples - positive_samples) ~
    I(cases_lag3 > 0) + as.factor(Site),
  data = merged_df,
  family = binomial
)
exp(coef(model_lag3)) # Odds ratios
exp(confint(model_lag3)) # Confidence intervals

# Plaques: Check distribution (majority are zeros, suggesting zero-inflated model)
mean(merged_df$total_plaques, na.rm = TRUE)
var(merged_df$total_plaques, na.rm = TRUE)
mean(merged_df$total_plaques == 0, na.rm = TRUE) # majority are 0's

# Standard neg binom model: Do case counts 1 month ago predict plaque counts?
model_nb <- glm.nb(
  total_plaques ~ cases_lag1 + as.factor(Site),
  data = merged_df
)
summary(model_nb)
exp(coef(model_nb))
# Use Wald confidence intervals (faster and more reliable with many site effects)
exp(confint.default(model_nb))

# Standard neg binom model: Do case counts 2 months ago predict plaque counts?
model_nb2 <- glm.nb(
  total_plaques ~ cases_lag2 + as.factor(Site),
  data = merged_df
)
summary(model_nb2)
exp(coef(model_nb2))
# Use Wald confidence intervals (faster and more reliable with many site effects)
5
# Zero-Inflated neg binom model: Do case counts 2 months ago predict plaque counts?
# Note: Using intercept-only zero-inflation to avoid overparameterization with Site effects
# Count part: cases_lag2 and Site effects; Zero-inflation: constant probability
model_zinb <- zeroinfl(
  total_plaques ~ cases_lag2 + as.factor(Site) | 1,
  dist = "negbin",
  data = merged_df
)

summary(model_zinb)
exp(coef(model_zinb))
exp(confint(model_zinb))

AIC(model_nb, model_zinb)

# REVERSE DIRECTION: Does previous phage presence predict future cases?
# Note: Phage presence is lagged (past), cases are current (present)
merged_df <- merged_df %>%
  arrange(Site, month) %>%
  group_by(Site) %>%
  mutate(
    phage_present_lag1 = lag(I(positive_samples > 0), 1), # Any positive samples 1 month ago (binary)
    phage_present_lag2 = lag(I(positive_samples > 0), 2), # Any positive samples 2 months ago (binary)
    plaques_lag1 = lag(total_plaques, 1), # Plaques from 1 month ago
    plaques_lag2 = lag(total_plaques, 2) # Plaques from 2 months ago
  ) %>%
  ungroup()

# Reverse Model 1: Does 1-month-ago phage positivity predict current cases?
model_cases_lag1 <- glm(
  I(cases > 0) ~ phage_present_lag1 + as.factor(Site),
  data = merged_df,
  family = binomial
)

exp(coef(model_cases_lag1)) # odds ratios
exp(confint(model_cases_lag1)) # CIs

# Reverse Model 2: Does 2-months-ago phage positivity predict current cases?
model_cases_lag2 <- glm(
  I(cases > 0) ~ phage_present_lag2 + as.factor(Site),
  data = merged_df,
  family = binomial
)

exp(coef(model_cases_lag2)) # odds ratios
exp(confint(model_cases_lag2)) # CIs

#### RESULTS SUMMARY FOR PRESENTATION ==========================================

# Helper function to extract key results (OR/IRR and CI) and hide Site effects
get_model_summary <- function(model, label, type = "OR") {
  # Extract coefficients and CIs
  if (inherits(model, "zeroinfl")) {
    items <- summary(model)$coefficients$count
    cis <- exp(confint(model))
    is_count <- grepl("count_", rownames(cis))
    cis <- cis[is_count, , drop = FALSE]
    rownames(cis) <- gsub("count_", "", rownames(cis))
  } else {
    items <- coef(summary(model))
    cis <- exp(confint.default(model))
  }

  # Filter out Site effects and Intercept
  # Using grepl ensures a logical vector of the same length as rownames
  main_effect_idx <- !grepl("Site", rownames(items), ignore.case = TRUE) &
    rownames(items) != "(Intercept)"

  # Check if we found any predictors to avoid row-mismatch errors
  if (sum(main_effect_idx) == 0) {
    message(paste("Warning: No main predictor found for model:", label))
    return(NULL)
  }

  res <- data.frame(
    Model = label,
    Predictor = rownames(items)[main_effect_idx],
    Type = type,
    Result = round(exp(items[main_effect_idx, 1]), 2),
    CI_Lower = round(cis[main_effect_idx, 1], 2),
    CI_Upper = round(cis[main_effect_idx, 2], 2),
    P_Value = format.pval(items[main_effect_idx, 4], digits = 3),
    stringsAsFactors = FALSE
  )
  return(res)
}

# Compile Table
summary_list <- list(
  get_model_summary(model_lag1, "Logistic (1-mo lag)", "OR"),
  get_model_summary(model_lag2, "Logistic (2-mo lag)", "OR"),
  get_model_summary(model_lag3, "Logistic (3-mo lag)", "OR"),
  get_model_summary(model_nb, "NegBinom (1-mo lag)", "IRR"),
  get_model_summary(model_nb2, "NegBinom (2-mo lag)", "IRR"),
  get_model_summary(model_cases_lag1, "Reverse Logistic (1-mo lag)", "OR"),
  get_model_summary(model_cases_lag2, "Reverse Logistic (2-mo lag)", "OR")
)

# Remove NULLs and combine
summary_table <- do.call(rbind, summary_list)

cat("\n--- KEY MODEL RESULTS FOR PRESENTATION ---\n")
print(summary_table)
cat("\nNote: All models adjusted for Catchment Site (fixed effects).\n")
cat("OR = Odds Ratio (Binary Positivity), IRR = Incidence Rate Ratio (Plaque Counts)\n")

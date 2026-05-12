#### TYPHOID R01 COVARIATE IMPORTANCE ANALYSIS ####
# Purpose: Quantify relative predictive strength of phage vs. environmental factors
# Outcome A: Clinical Case Presence (Binary: I(Cases > 0))
# Outcome B: Phage Positivity (Beta-Binomial)

rm(list = ls())


library(tidyverse)
library(readxl)
library(lubridate)
library(glmmTMB)
library(broom.mixed)
library(patchwork)

# 1. CONFIGURE PATHS
# Assuming running from repo root
data_dir <- "data"
manuscript_dir <- "data" 

# 2. LOAD DATA
# 2.1 Phage Data
phage_raw <- read.csv(file.path(data_dir, "REDCapPhageData.csv"))
phage_clean <- phage_raw %>%
  mutate(
    # Handle multiple date formats (M/D/YY from REDCap vs ISO YYYY-MM-DD)
    date_sample = as.Date(parse_date_time(date_results_tw, c("mdy", "ymd"))),
    # Additional month-year correction logic from original script
    date_sample = as.Date(ifelse(
      !is.na(date_sample) & (year(date_sample) == 2004 | year(date_sample) < 2024),
      paste0("2024", substr(as.character(date_sample), 5, 10)),
      as.character(date_sample)
    )),
    month = floor_date(date_sample, "month"),
    status = case_when(
      amplification_tw == 1 | direct_tw == 1 ~ "Positive",
      amplification_tw == 2 & direct_tw == 2 ~ "Negative",
      TRUE ~ "Negative"
    )
  ) %>%
  group_by(Site = site_id_tw, month) %>%
  summarise(
    total_samples = n(),
    positive_samples = sum(status == "Positive", na.rm = TRUE),
    total_plaques = sum(direct_counts_tw, na.rm = TRUE),
    phage_present = ifelse(positive_samples > 0, 1, 0),
    .groups = "drop"
  )

# 2.2 Clinical & Weather Data
cases_file <- file.path(manuscript_dir, "Site_monthlytyphoidcases_Jan 2026.xlsx")
case_sheets <- excel_sheets(cases_file)
cases_all <- map_df(case_sheets, function(s) {
  read_excel(cases_file, sheet = s) %>%
    pivot_longer(cols = matches("^(2024|2025|2026)[_-]"), names_to = "month_str", values_to = "clinical_cases") %>%
    mutate(month = ym(month_str))
}) %>%
  group_by(Site, month) %>%
  summarise(clinical_cases = max(clinical_cases, na.rm = TRUE), overlap = first(overlap), .groups = "drop") %>%
  filter(overlap == 1)

weather_raw <- read_excel(file.path(manuscript_dir, "weather_flow_hf183.xlsx"))
weather_clean <- weather_raw %>%
  mutate(
    Site = as.numeric(substr(as.character(Sample_ID), 1, 3)),
    month = floor_date(as.Date(date_collection_new), "month"),
    rainfall = ifelse(tolower(`x24_hours_rainfall_mm_recorded_from_0830_hrs_ist_of_yesterday_to_0830_hrs_ist_of_today`) %in% c("nil", "trace"), 0, 
                      as.numeric(`x24_hours_rainfall_mm_recorded_from_0830_hrs_ist_of_yesterday_to_0830_hrs_ist_of_today`)),
    avg_temp = (as.numeric(maximum_temp_o_c_recorded) + as.numeric(minimum_temp_o_c_recorded)) / 2,
    flow = as.numeric(g_per_hour_in_m3),
    humidity = as.numeric(relative_humidity_at_1730_hrs_percent_recorded),
    # HF183 Absolute Copy Number Logic: Using hf_183_tr (raw/continuous)
    hf183_abs = as.numeric(hf_183_tr)
  ) %>%
  group_by(Site, month) %>%
  summarise(avg_temp = mean(avg_temp, na.rm=T), total_rain = sum(rainfall, na.rm=T), 
            avg_flow = mean(flow, na.rm=T), 
            hf183_abs = mean(hf183_abs, na.rm=T),
            .groups = "drop")

# 2.3 Final Combined Analysis Set
analysis_df <- cases_all %>%
  inner_join(phage_clean, by = c("Site", "month")) %>%
  left_join(weather_clean, by = c("Site", "month"))

# 3. STATISTICAL MODELING
# We scale continuous predictors to improve model convergence
analysis_df <- analysis_df %>%
  mutate(across(c(avg_flow, total_rain, avg_temp, hf183_abs), scale, .names = "{.col}_scaled"))

run_importance <- function(outcome_f, family_obj, predictors, label) {
  res <- data.frame()
  
  # Multivariable Adjusted Model (Hierarchical)
  f_adj <- sprintf("%s ~ %s + (1|Site)", outcome_f, paste(predictors, collapse = " + "))
  m_adj <- glmmTMB(as.formula(f_adj), family = family_obj, data = analysis_df)
  tidy_adj <- tidy(m_adj, conf.int = TRUE, conf.method = "wald", exponentiate = TRUE, effects = "fixed") %>%
    filter(term != "(Intercept)") %>% 
    mutate(Model = "Adjusted", Outcome = label)
  
  # Individual Crude Models (Hierarchical)
  for (p in predictors) {
    f_crude <- sprintf("%s ~ %s + (1|Site)", outcome_f, p)
    m_crude <- glmmTMB(as.formula(f_crude), family = family_obj, data = analysis_df)
    tidy_p <- tidy(m_crude, conf.int = TRUE, conf.method = "wald", exponentiate = TRUE, effects = "fixed") %>%
      filter(grepl(p, term, fixed = TRUE)) %>%
      mutate(Model = "Crude", Outcome = label)
    res <- bind_rows(res, tidy_p)
  }
  bind_rows(res, tidy_adj)
}

env_preds_scaled_final <- c("avg_flow_scaled", "total_rain_scaled", "avg_temp_scaled", "hf183_abs_scaled")

# # A) Clinical Presence (Now using raw clinical cases as a predictor)
# importance_clinical <- run_importance(
#   "I(clinical_cases > 0)", binomial, c("phage_present", env_preds_scaled_final), "Clinical Case Presence (Y/N)"
# )

# B1) Phage Positivity (SD-Scaled Predictors for relative importance)
importance_phage_sd <- run_importance(
  "phage_present", binomial, c("clinical_cases", env_preds_scaled_final), "Phage Positivity (aOR - Standardized)"
)

# B2) Phage Positivity (Raw-Scale Predictors for absolute effect)
env_preds_raw <- c("avg_flow", "total_rain", "avg_temp", "hf183_abs")
importance_phage_raw <- run_importance(
  "phage_present", binomial, c("clinical_cases", env_preds_raw), "Phage Positivity (aOR - Raw Units)"
)

# 4. CONSOLIDATE AND OUTPUT
importance_results <- bind_rows(importance_clinical, importance_phage_sd, importance_phage_raw)

# Clean up term names for presentation
importance_results <- importance_results %>%
  mutate(Predictor = case_when(
    term == "phage_present" ~ "Phage Detection (Y/N)",
    term == "clinical_cases" ~ "Number of Clinical Cases",
    grepl("clinical_cases", term) ~ "Number of Clinical Cases",
    term == "avg_flow_scaled" ~ "Avg Flow Rate(1 SD)",
    term == "total_rain_scaled" ~ "Cumulative Rainfall (1 SD)",
    term == "avg_temp_scaled" ~ "Avg Temperature (1 SD)",
    term == "hf183_abs_scaled" ~ "Absolute copy numbers (1 SD)",
    term == "avg_flow" ~ "Avg Flow Rate (m3/h)",
    term == "total_rain" ~ "Cumulative Rainfall (mm)",
    term == "avg_temp" ~ "Avg Temperature (C)",
    term == "hf183_abs" ~ "Absolute copy numbers",
    TRUE ~ term
  ))

# Save Table
write.csv(importance_results %>% dplyr::select(Outcome, Predictor, Model, estimate, conf.low, conf.high, p.value), 
          "Typhoid_Covariate_Importance_Table.csv", row.names = FALSE)

# Generate Visualization (SD-Scaled Phage)
p_phage_sd <- importance_results %>%
  filter(Outcome == "Phage Positivity (aOR - Standardized)") %>%
  ggplot(aes(x = estimate, y = Predictor, color = Model)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(position = position_dodge(0.5), size = 4) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, position = position_dodge(0.5)) +
  geom_text(aes(label = round(estimate, 2)), vjust = -1.2, position = position_dodge(0.5), show.legend = FALSE, size = 4) +
  theme_minimal(base_size = 14) +
  labs(title = "Predictors of Phage Detection (Standardized)",
       subtitle = "Effect of 1 SD increase in each covariate",
       x = "Effect Size (Odds Ratio) and 95% CI", y = "") +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

ggsave("Phage_Positivity_Predictors_SD.png", p_phage_sd, width = 10, height = 7)

# Generate Visualization (Raw-Scale Phage)
p_phage_raw <- importance_results %>%
  filter(Outcome == "Phage Positivity (aOR - Raw Units)") %>%
  ggplot(aes(x = estimate, y = Predictor, color = Model)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(position = position_dodge(0.5), size = 4) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, position = position_dodge(0.5)) +
  geom_text(aes(label = round(estimate, 2)), vjust = -1.2, position = position_dodge(0.5), show.legend = FALSE, size = 4) +
  theme_minimal(base_size = 14) +
  labs(title = "Predictors of Phage Detection (Raw Units)",
       subtitle = "Effect of 1 unit increase (1mm rain, 1C temp, etc.)",
       x = "Effect Size (Odds Ratio) and 95% CI", y = "") +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

ggsave("Phage_Positivity_Predictors_Raw.png", p_phage_raw, width = 10, height = 7)

# # Also update the combined one for consistency
# p_combined <- ggplot(importance_results %>% filter(!grepl("Raw Units", Outcome)), 
#                     aes(x = estimate, y = Predictor, color = Model)) +
#   geom_vline(xintercept = 1, linetype = "dashed") +
#   geom_point(position = position_dodge(0.5), size = 3) +
#   geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, position = position_dodge(0.5)) +
#   geom_text(aes(label = round(estimate, 2)), vjust = -1.5, position = position_dodge(0.5), show.legend = FALSE, size = 3.5) +
#   facet_wrap(~Outcome, scales = "free_x") +
#   theme_minimal(base_size = 14) +
#   labs(x = "Estimate (OR/aOR) and 95% CI", y = "") +
#   theme(legend.position = "top", strip.text = element_text(face = "bold"))
# ggsave("Predictors.png", p_combined, width = 15, height = 8)

# 5. CORRELATION HEATMAP (TO SHOW COLLINEARITY)
# -------------------------------------------------------------------------
cor_data <- analysis_df %>%
  dplyr::select(clinical_cases, positive_samples, avg_flow, total_rain, avg_temp) %>%
  drop_na()

# Rename for heatmap
colnames(cor_data) <- c("Clinical Cases", "Phage Pos.", "Avg Flow", "Rainfall", "Temperature")
cor_matrix <- cor(cor_data)
cor_df <- as.data.frame(as.table(cor_matrix))

p_cor <- ggplot(cor_df, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = round(Freq, 2)), size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1)) +
  theme_minimal() +
  labs(title = "Covariate Correlation Heatmap", x = "", y = "", fill = "Correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Typhoid_Covariate_Correlations.png", p_cor, width = 10, height = 8)
p_cor

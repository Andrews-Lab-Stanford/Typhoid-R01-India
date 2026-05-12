#### TYPHOID R01 — EXPLORATORY BIWEEKLY MODELS ####
# Purpose: Test if collapsing weekly case counts into 2-week windows improves predictive signal.

library(tidyverse)
library(glmmTMB)
library(broom.mixed)
library(patchwork)

# ── 1. Load & Prepare Data ────────────────────────────────────────────────────
base_dir <- "/Users/esthjung/Library/CloudStorage/GoogleDrive-esthjung@stanford.edu/My Drive/Typhoid R01"
load(file.path(base_dir, "Data/analysis_ready_weekly.RData")) # merged_weekly, analysis_df_weekly

cat("=== EXPLORATORY BIWEEKLY ANALYSIS ===\n")

# Re-create analytic set with biweekly sums
# We use merged_weekly to get the unscaled raw lags
biweekly_df <- merged_weekly %>%
  filter(!is.na(temp), !is.na(hf183), !is.na(cases_lag8)) %>%
  mutate(
    cases_0_1 = cases_lag0 + cases_lag1,
    cases_2_3 = cases_lag2 + cases_lag3,
    cases_4_5 = cases_lag4 + cases_lag5,
    cases_6_7 = cases_lag6 + cases_lag7
  )

# ── 2. Model Function ─────────────────────────────────────────────────────────
fit_biweekly <- function(lag_label, lag_var) {
  cat("  Fitting Biweekly Lag", lag_label, "...\n")
  
  f_crude <- paste0(" ~ ", lag_var, " + (1 | Site)")
  f_adj   <- paste0(" ~ ", lag_var, " + flow + rainfall + temp + hf183 + pct_under_15 + pct_water_improved + pct_toilet_improved + pct_high_ses + (1 | Site)")
  
  tryCatch({
    # --- Q1 Positivity ---
    m_q1_crude <- glmmTMB(as.formula(paste0("cbind(positive_samples, total_samples - positive_samples)", f_crude)), family = betabinomial, data = biweekly_df)
    m_q1_adj   <- glmmTMB(as.formula(paste0("cbind(positive_samples, total_samples - positive_samples)", f_adj)), family = betabinomial, data = biweekly_df)
    
    # --- Q2 Abundance ---
    m_q2_crude <- glmmTMB(as.formula(paste0("total_plaques", f_crude)), ziformula = ~1, family = poisson, data = biweekly_df)
    m_q2_adj   <- glmmTMB(as.formula(paste0("total_plaques", f_adj)), ziformula = ~1, family = poisson, data = biweekly_df)
    
    # Extract
    res_q1_crude <- tidy(m_q1_crude, effects = "fixed", conf.int = TRUE) %>% filter(term == lag_var) %>% mutate(Outcome = "Q1 (Positivity)", Window = lag_label, Model = "Crude")
    res_q1_adj   <- tidy(m_q1_adj,   effects = "fixed", conf.int = TRUE) %>% filter(term == lag_var) %>% mutate(Outcome = "Q1 (Positivity)", Window = lag_label, Model = "Adjusted")
    res_q2_crude <- tidy(m_q2_crude, effects = "fixed", conf.int = TRUE) %>% filter(term == lag_var) %>% mutate(Outcome = "Q2 (Abundance)",  Window = lag_label, Model = "Crude")
    res_q2_adj   <- tidy(m_q2_adj,   effects = "fixed", conf.int = TRUE) %>% filter(term == lag_var) %>% mutate(Outcome = "Q2 (Abundance)",  Window = lag_label, Model = "Adjusted")
    
    bind_rows(res_q1_crude, res_q1_adj, res_q2_crude, res_q2_adj)
  }, error = function(e) { cat("  Error:", e$message, "\n"); NULL })
}

# ── 3. Run Models ─────────────────────────────────────────────────────────────
results <- bind_rows(
  fit_biweekly("Weeks 0-1", "cases_0_1"),
  fit_biweekly("Weeks 2-3", "cases_2_3"),
  fit_biweekly("Weeks 4-5", "cases_4_5"),
  fit_biweekly("Weeks 6-7", "cases_6_7")
)

# ── 4. Forest Plot ────────────────────────────────────────────────────────────
cat("\nGenerating Forest Plot...\n")

plot_data <- results %>%
  mutate(
    Window = factor(Window, levels = rev(c("Weeks 0-1", "Weeks 2-3", "Weeks 4-5", "Weeks 6-7"))),
    Model  = factor(Model, levels = c("Crude", "Adjusted")),
    Est    = exp(estimate),
    CI_lo  = exp(conf.low),
    CI_hi  = exp(conf.high)
  )

build_forest <- function(df, xlab, title, line_at = 1) {
  ggplot(df, aes(y = Window, x = Est, color = Model)) +
    geom_vline(xintercept = line_at, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.3, position = position_dodge(width = 0.5)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_text(aes(label = round(Est, 2)), 
              position = position_dodge(width = 0.5), 
              vjust = -1.2, size = 3, show.legend = FALSE) +
    scale_color_manual(values = c("Crude" = "#95a5a6", "Adjusted" = "#e67e22")) +
    labs(title = title, x = xlab, y = "Lag (biweekly)") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
}

p1 <- build_forest(plot_data %>% filter(Outcome == "Q1 (Positivity)"), 
                  "Odds Ratio (95% CI)", "A) Phage Positivity (OR)")

p2 <- build_forest(plot_data %>% filter(Outcome == "Q2 (Abundance)"), 
                  "Incidence Rate Ratio (95% CI)", "B) Phage Abundance (IRR)")

combined_fig <- (p1 | p2) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Biweekly Exploratory Analysis: Crude vs. Adjusted Estimates",
                  theme = theme(legend.position = "bottom"))

# Save
figures_dir <- file.path(base_dir, "Typhoid-R01-India/figures")
ggsave(file.path(figures_dir, "Exploratory_Biweekly_Forest_Plot.png"), combined_fig, width = 12, height = 7, dpi = 300)

# ── 5. Display Table ──────────────────────────────────────────────────────────
final_tab <- results %>%
  transmute(
    Outcome,
    Window,
    Model,
    Estimate = round(exp(estimate), 2),
    CI_Lower = round(exp(conf.low), 2),
    CI_Upper = round(exp(conf.high), 2),
    P_Value  = format.pval(p.value, digits = 3)
  )

print(as.data.frame(final_tab), right = FALSE)

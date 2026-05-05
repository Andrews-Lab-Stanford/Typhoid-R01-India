#### TYPHOID R01 — STEP 2: MAIN MODELS & MANUSCRIPT FIGURES ####
# Last Updated: May 2026
# Purpose: Run crude AND adjusted glmmTMB models across lags 0/1/2, generate figures.
#
# Crude model:    cases_lag + (1 | Site)  [site random effect only, no covariates]
# Adjusted model: cases_lag + weather + SES/WASH covariates + (1 | Site)
#
# INPUT:  Data/analysis_ready.RData  (run 1-data-cleaning/1_Data_Cleaning.R first)
# OUTPUT: figures/SAP_Slide_Q1.png       — OR  forest plot (betabinomial), crude + adjusted
#         figures/SAP_Slide_Q2.png       — IRR forest plot (ZIP Poisson), crude + adjusted
#         figures/SAP_Slide_Combined.png — Both panels side by side

library(tidyverse)
library(lubridate)
library(glmmTMB)
library(broom.mixed)
library(patchwork)

# ── Directories ────────────────────────────────────────────────────────────────
base_dir    <- "/Users/esthjung/Library/CloudStorage/GoogleDrive-esthjung@stanford.edu/My Drive/Typhoid R01"
data_dir    <- file.path(base_dir, "Data")
figures_dir <- file.path(base_dir, "Typhoid-R01-India/figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load Data ──────────────────────────────────────────────────────────────────
rdata_path <- file.path(data_dir, "analysis_ready.RData")
if (!file.exists(rdata_path)) {
  stop("Run 1-data-cleaning/1_Data_Cleaning.R first to generate: ", rdata_path)
}
load(rdata_path)   # loads: merged_data, res_table1, analysis_df
cat("=== TYPHOID R01: MAIN MODELS ===\n")
cat("Analytic dataset:", nrow(analysis_df), "site-months,",
    n_distinct(analysis_df$Site), "sites\n\n")

# ── Step 1 & 2: Fit crude + adjusted models using exact same N=787 dataset ────
# Note: We use analysis_df for ALL lags to ensure an apples-to-apples comparison.
# We use the predefined variables: clinical_cases (Lag 0), cases_lag1, cases_lag2
fit_models <- function(lag_n, lag_var) {
  cat("  Fitting Lag", lag_n, "...\n")
  
  # Crude: Predictor + random effect for site.
  f_crd <- paste0(" ~ ", lag_var, " + (1 | Site)")
  
  # Adjusted: Covariates + census variables + random effect for site.
  f_adj <- paste0(" ~ ", lag_var, " + flow + rainfall + temp + hf183 + pct_under_15 + pct_water_improved + pct_toilet_improved + pct_high_ses + (1 | Site)")

  
  tryCatch({
    list(
      lag = lag_n,
      lag_var = lag_var,
      bb_crd = glmmTMB(
        as.formula(paste0("cbind(positive_samples, total_samples - positive_samples)", f_crd)),
        family = betabinomial, data = analysis_df),
      bb_adj = glmmTMB(
        as.formula(paste0("cbind(positive_samples, total_samples - positive_samples)", f_adj)),
        family = betabinomial, data = analysis_df),
      zip_crd = glmmTMB(
        as.formula(paste0("total_plaques", f_crd)),
        ziformula = ~1, family = poisson, data = analysis_df),
      zip_adj = glmmTMB(
        as.formula(paste0("total_plaques", f_adj)),
        ziformula = ~1, family = poisson, data = analysis_df)
    )
  }, error = function(e) { cat("  ERROR at Lag", lag_n, ":", e$message, "\n"); NULL })
}

# ── Step 3: Run across lags 0, 1, 2 ──────────────────────────────────────────
all_fits <- list(
  fit_models(0, "clinical_cases"),
  fit_models(1, "cases_lag1"),
  fit_models(2, "cases_lag2")
)

# ── Step 4: Extract coefficients ──────────────────────────────────────────────
extract_est <- function(model, model_label, outcome_label, lag_n, lag_var) {
  tryCatch({
    tidy(model, effects = "fixed", conf.int = TRUE) %>%
      filter(term == lag_var) %>%
      transmute(
        Lag     = paste0("Lag ", lag_n),
        Model   = model_label,
        Outcome = outcome_label,
        Est     = exp(estimate),
        CI_lo   = exp(conf.low),
        CI_hi   = exp(conf.high)
      )
  }, error = function(e) NULL)
}

plot_data <- bind_rows(lapply(all_fits, function(r) {
  if (is.null(r)) return(NULL)
  bind_rows(
    extract_est(r$bb_crd,  "Crude",    "Q1", r$lag, r$lag_var),
    extract_est(r$bb_adj,  "Adjusted", "Q1", r$lag, r$lag_var),
    extract_est(r$zip_crd, "Crude",    "Q2", r$lag, r$lag_var),
    extract_est(r$zip_adj, "Adjusted", "Q2", r$lag, r$lag_var)
  )
})) %>%
  mutate(
    Lag   = factor(Lag,   levels = c("Lag 2", "Lag 1", "Lag 0")),
    # Reverse factor levels so Adjusted comes first (plots lower) and Crude plots higher
    Model = factor(Model, levels = c("Adjusted", "Crude"))
  )

cat("\nAll estimates (cases_lag coefficient):\n")
print(plot_data)

# ── Color palette: Crude = blue, Adjusted = red ────────────────────────────────
pal <- c("Crude" = "#2c7bb6", "Adjusted" = "#d7191c")

# ── Shared theme ──────────────────────────────────────────────────────────────
forest_theme <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.title       = element_text(face = "bold"),
    legend.position    = "bottom"
  )

# Helper: build one forest panel
build_forest <- function(data, xlab, show_y = TRUE) {
  p <- ggplot(data, aes(x = Est, y = Lag, color = Model)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi),
                   height = 0.2, linewidth = 0.8,
                   position = position_dodge(0.5)) +
    geom_point(size = 3.5,
               position = position_dodge(0.5)) +
    geom_text(aes(label = round(Est, 2)),
              vjust = -1.2, size = 3.8,
              position = position_dodge(0.5),
              show.legend = FALSE) +
    # Ensure legend still shows Crude first
    scale_color_manual(values = pal, name = "Model", breaks = c("Crude", "Adjusted")) +
    forest_theme +
    labs(x = xlab, y = if (show_y) "Lag (months)" else NULL)
  if (!show_y) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p
}



# ── Q1 Forest Plot: OR (Betabinomial) ─────────────────────────────────────────
fig_q1 <- build_forest(plot_data %>% filter(Outcome == "Q1"),
                        xlab = "Odds Ratio (95% CI)", show_y = TRUE)

# ── Q2 Forest Plot: IRR (ZIP Poisson) ─────────────────────────────────────────
fig_q2 <- build_forest(plot_data %>% filter(Outcome == "Q2"),
                        xlab = "Incidence Rate Ratio (95% CI)", show_y = FALSE)

# ── Combined — shared legend ────────────────────────────────────────────────────
fig_combined <- (fig_q1 | fig_q2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# ── Save ───────────────────────────────────────────────────────────────────────
ggsave(file.path(figures_dir, "SAP_Slide_Q1.png"),       fig_q1,       width = 8,  height = 5, dpi = 300)
ggsave(file.path(figures_dir, "SAP_Slide_Q2.png"),       fig_q2,       width = 8,  height = 5, dpi = 300)
ggsave(file.path(figures_dir, "SAP_Slide_Combined.png"), fig_combined, width = 14, height = 6, dpi = 300)

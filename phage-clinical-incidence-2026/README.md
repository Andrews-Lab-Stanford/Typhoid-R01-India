# Phage detection and clinical typhoid incidence, Vellore 2024–2026

Analysis code for the manuscript on *Salmonella* Typhi bacteriophage detection in wastewater and
blood culture-confirmed typhoid incidence. This is a separate analysis from the earlier work in
`1-data-cleaning/`, `2-analysis/` and `3-tables/`, which uses different time units, lag structure
and abundance models and does not reproduce the estimates reported here.

## Input data

Both scripts read `dataset_WWphage_R01.xlsx`, a single consolidated analysis file with a `Sheet1`
of one row per wastewater sample and a `data dictionary` sheet. **The dataset is not included in
this repository** and is excluded by `.gitignore`; place it at `data/dataset_WWphage_R01.xlsx` or
pass a path as the first argument.

The version used for the reported estimates contains 1,862 samples from the 38 non-overlapping
(catchment overlap level 1) catchments. Between-catchment analyses across all 54 catchments require
an export that also contains overlap levels 2 and 3.

## Scripts

Run from the repository root.

```bash
Rscript phage-clinical-incidence-2026/final_estimates.R [path/to/dataset.xlsx]
Rscript phage-clinical-incidence-2026/temporal_inference.R [path/to/dataset.xlsx] [n_reps]
```

| Script | Purpose | Runtime |
|:---|:---|:---|
| `final_estimates.R` | Every estimate reported in the Results. Writes `final_estimates.csv`. | a few seconds |
| `temporal_inference.R` | Supporting robustness work: lag profile and permutation null distributions. Writes `temporal_inference_results.csv`. | ~45 s at 1,000 replicates |

Outputs are written to the directory the script is run from and are not committed;
re-running regenerates them.

Requires R >= 4.1 with `readxl`, `survival`, `splines` and `lme4`. Tested on R 4.5.1 with
survival 3.8.3, lme4 1.1.37 and readxl 1.5.0.

## Specification

All within-catchment estimates come from one model: conditional logistic regression of sample-level
phage positivity, stratified by catchment, with calendar time modelled as a natural cubic spline in
fortnight (4 degrees of freedom). Stratification removes all between-catchment variation, so each
estimate is a within-catchment contrast, and the spline accounts for study-wide seasonal and secular
change in positivity without the 50 fortnight indicators that cause separation.

Exposure is the occurrence of at least one blood culture-confirmed case in the catchment within a
window of fortnights defined relative to the sample. `final_estimates.R` reports windows from the
sampled fortnight alone to the preceding 12 weeks, a window 8 to 16 weeks before the sample, a
window 8 weeks after it, the number of cases within the 8-week window, and models adjusted for log
HF183 concentration and measured flow. Sensitivity to the time adjustment is reported across splines
of 3 to 8 degrees of freedom and with no adjustment for calendar time. Variance between catchments
and over time is decomposed with a mixed-effects model; `lme4` does not converge reliably on the
adjusted models in these data and is not used for anything else.

`temporal_inference.R` is retained as supporting work rather than as the basis for any reported
estimate. It tests the primary association against three permutation nulls that preserve catchment
burden and the clustering of cases in time while destroying their alignment with phage timing.

## Guards on the input file

Both scripts stop rather than produce wrong numbers if the source file changes:

1. **`amplification_tw` and `direct_tw` must be coded `{1 = yes, 2 = no}`.** The data dictionary
   states 0/1 but the delivered file uses 1/2. If this is corrected upstream the assertion fires
   instead of silently inverting the outcome.
2. **`positive` must equal `amplification OR direct`.** Note that `direct_tw` is missing for the
   first nine fortnights, so `positive` reduces to amplification alone in that period.
3. **Catchment overlap scope is checked, not assumed**, so an export covering more than one overlap
   stratum is not silently analysed as if it were level 1.

Stata `.` missing values in `ln_hf183` and `hf_183_tr` are coerced to `NA` rather than being allowed
to turn the column into character.

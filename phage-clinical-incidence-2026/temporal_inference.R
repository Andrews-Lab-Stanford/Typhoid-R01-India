# ==============================================================================
# Typhoid R01 phage manuscript
# Temporal association, permutation-based inference, and negative controls
#
# Reads dataset_WWphage_R01.xlsx directly and reproduces every temporal number
# in the Results section dated 2026-10-12.
#
# Usage (from the repository root):
#   Rscript phage-clinical-incidence-2026/temporal_inference.R [path/to/dataset.xlsx] [n_reps]
# Default n_reps = 1000. Runtime ~3-5 min at 1000 reps.
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl); library(survival); library(lme4)
})

args    <- commandArgs(trailingOnly = TRUE)
XLSX    <- if (length(args) >= 1) args[1] else "data/dataset_WWphage_R01.xlsx"
NREP    <- if (length(args) >= 2) as.integer(args[2]) else 1000L
OUTDIR  <- "."
SEED    <- 20261012
set.seed(SEED)

stopifnot(file.exists(XLSX))
cat("Reading:", normalizePath(XLSX), "\n")
raw <- as.data.frame(read_excel(XLSX, sheet = "Sheet1", guess_max = 5000))
cat("Rows:", nrow(raw), " Cols:", ncol(raw), "\n\n")

# ------------------------------------------------------------------ recoding
# GUARD 1: amplification_tw / direct_tw are coded 2 = no, 1 = yes in the
# delivered file, NOT 0/1 as the data dictionary states. If the upstream file is
# ever corrected to 0/1 this assertion fires rather than silently inverting the
# outcome.
num <- function(x) suppressWarnings(as.numeric(as.character(x)))   # handles Stata "." 
chk <- function(v, nm) {
  u <- sort(unique(na.omit(v)))
  if (!identical(u, c(1, 2)))
    stop(sprintf("%s: expected values {1,2} (1=yes, 2=no) but found {%s}. ",
                 nm, paste(u, collapse = ",")),
         "The source file coding has changed - revisit the recode below.")
  invisible(TRUE)
}
raw$amplification_tw <- num(raw$amplification_tw)
raw$direct_tw        <- num(raw$direct_tw)
chk(raw$amplification_tw, "amplification_tw")
chk(raw$direct_tw,        "direct_tw")

d <- data.frame(
  site   = raw$ESsite,
  fn     = num(raw$fn_num),
  amp    = as.integer(raw$amplification_tw == 1),
  direct = ifelse(is.na(raw$direct_tw), NA_integer_, as.integer(raw$direct_tw == 1)),
  either = as.integer(num(raw$positive) == 1),
  cases  = num(raw$cases_fn),
  pop    = num(raw$catchment_pop),
  ovl    = num(raw$ca_overlap),
  ln183  = num(raw$ln_hf183),         # Stata "." -> NA
  flow   = num(raw$Flow_PerHourM3)
)

# GUARD 2: `positive` must equal (amp OR direct). Verified 1862/1862 on the
# 2026-10-12 file. Note this makes `positive` amplification-only wherever direct
# is missing (fortnights 2-11), i.e. differential ascertainment by period.
recon <- as.integer((d$amp == 1) | (!is.na(d$direct) & d$direct == 1))
if (any(recon != d$either))
  stop(sum(recon != d$either), " rows where `positive` != (amplification OR direct).")

# GUARD 3: scope. The 2026-10-12 file is non-overlapping catchments only.
cat("Catchment overlap levels present:", paste(sort(unique(d$ovl)), collapse = ", "), "\n")
if (length(unique(d$ovl)) == 1L)
  cat("NOTE: single overlap stratum (level ", unique(d$ovl),
      ") - between-catchment analyses across all 54 catchments are NOT computable ",
      "from this file.\n", sep = "")
cat(sprintf("Catchments: %d | samples: %d | population: %s | cases: %d\n",
            length(unique(d$site)), nrow(d),
            format(sum(tapply(d$pop, d$site, function(z) z[1])), big.mark = ","),
            sum(d$cases)))
cat(sprintf("Positivity -- amplification %d/%d (%.1f%%) | direct %d/%d (%.1f%%) | either %d/%d (%.1f%%)\n\n",
            sum(d$amp), nrow(d), 100 * mean(d$amp),
            sum(d$direct, na.rm = TRUE), sum(!is.na(d$direct)),
            100 * mean(d$direct, na.rm = TRUE),
            sum(d$either), nrow(d), 100 * mean(d$either)))

# ------------------------------------------- site x fortnight case count matrix
sites <- sort(unique(d$site)); fns <- min(d$fn):max(d$fn); nf <- length(fns)
C <- matrix(0, length(sites), nf, dimnames = list(sites, fns))
agg <- aggregate(cases ~ site + fn, d, max)
C[cbind(as.character(agg$site), as.character(agg$fn))] <- agg$cases
d$si <- match(d$site, sites); d$fi <- match(d$fn, fns)

# any case in the window [lo, hi] fortnights relative to the sample's fortnight
# (negative = before the sample, positive = after it)
win <- function(Cm, lo, hi) {
  out <- matrix(FALSE, nrow(Cm), ncol(Cm))
  for (j in seq_len(ncol(Cm))) {
    k <- (j + lo):(j + hi); k <- k[k >= 1 & k <= ncol(Cm)]
    if (length(k)) out[, j] <- rowSums(Cm[, k, drop = FALSE]) > 0
  }
  out
}
expose <- function(w) w[cbind(d$si, d$fi)]

OUTCOMES <- c(amplification = "amp", direct = "direct", either = "either")

# site-stratified conditional logistic: the within-catchment contrast, with no
# distributional assumption about the catchment effect
cl_fit <- function(oc, x) {
  dd <- data.frame(y = d[[oc]], x = x, st = d$site)
  dd <- dd[!is.na(dd$y), ]
  if (length(unique(dd$x)) < 2) return(c(NA, NA, NA, NA, NA))
  m <- try(clogit(y ~ x + strata(st), data = dd), silent = TRUE)
  if (inherits(m, "try-error")) return(c(NA, NA, NA, NA, NA))
  s <- summary(m)$coefficients
  c(beta = s[1, 1], or = exp(s[1, 1]),
    lo = exp(s[1, 1] - 1.96 * s[1, 3]), hi = exp(s[1, 1] + 1.96 * s[1, 3]), p = s[1, 5])
}
cl_beta <- function(oc, x) cl_fit(oc, x)[["beta"]]

results <- list()
add <- function(...) results[[length(results) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)

# ============================================================== 1. main estimates
cat("=== 1. PRIMARY ESTIMATES: any case in the past 8 weeks (fortnights -3..0) ===\n")
w_past <- win(C, -3, 0); x_past <- expose(w_past)
w_fut  <- win(C,  1, 4); x_fut  <- expose(w_fut)
cat(sprintf("exposed site-fortnights: past-8wk %d | future-8wk %d | overlap %d\n",
            sum(x_past), sum(x_fut), sum(x_past & x_fut)))
for (nm in names(OUTCOMES)) {
  r <- cl_fit(OUTCOMES[[nm]], x_past)
  cat(sprintf("  %-14s clogit OR %5.2f (%4.2f-%5.2f)  asymptotic p=%.5f\n",
              nm, r[["or"]], r[["lo"]], r[["hi"]], r[["p"]]))
  add(analysis = "past 8wk", outcome = nm, or = r[["or"]], lo = r[["lo"]], hi = r[["hi"]],
      p_asymptotic = r[["p"]], p_perm_indep = NA, p_perm_common = NA)
}
# mixed-effects version for comparability with the earlier draft
cat("\n  mixed-effects (catchment + fortnight random intercepts), for comparison:\n")
for (nm in names(OUTCOMES)) {
  dd <- data.frame(y = d[[OUTCOMES[[nm]]]], x = x_past,
                   st = factor(d$site), ft = factor(d$fn))
  dd <- dd[!is.na(dd$y), ]
  m <- try(glmer(y ~ x + (1 | st) + (1 | ft), data = dd, family = binomial,
                 control = glmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 3e5))), silent = TRUE)
  if (!inherits(m, "try-error")) {
    co <- summary(m)$coefficients
    conv <- if (is.null(m@optinfo$conv$lme4$code)) "" else "  [CONVERGENCE WARNING]"
    cat(sprintf("  %-14s glmer  OR %5.2f (%4.2f-%5.2f)%s\n", nm, exp(co[2, 1]),
                exp(co[2, 1] - 1.96 * co[2, 2]), exp(co[2, 1] + 1.96 * co[2, 2]), conv))
  }
}

# ================================================================ 2. lag profile
cat("\n=== 2. LAG PROFILE: case indicator at single-fortnight offset L ===\n")
cat("L<0 = cases before the sample; L>0 = cases after the sample\n")
cat(sprintf("%4s %-24s %-24s %-24s %8s\n", "L", "amplification", "direct", "either", "n_exposed"))
for (L in -5:5) {
  w <- win(C, L, L); x <- expose(w)
  rr <- lapply(OUTCOMES, cl_fit, x = x)
  cat(sprintf("%4d %s %8d\n", L,
      paste(sapply(rr, function(r) sprintf("%5.2f (%4.2f-%5.2f)", r[["or"]], r[["lo"]], r[["hi"]])),
            collapse = " "), sum(x)))
  for (nm in names(OUTCOMES)) {
    r <- rr[[nm]]
    add(analysis = paste0("lag L=", L), outcome = nm, or = r[["or"]], lo = r[["lo"]],
        hi = r[["hi"]], p_asymptotic = r[["p"]], p_perm_indep = NA, p_perm_common = NA)
  }
}

cat("\n=== 2b. 8-WEEK BLOCKS at increasing offset ===\n")
blocks <- list("past 8wk (L -3..0)" = c(-3, 0), "future 8wk (L 1..4)" = c(1, 4),
               "past 16-8wk (L -7..-4)" = c(-7, -4), "future 8-16wk (L 5..8)" = c(5, 8))
for (bn in names(blocks)) {
  b <- blocks[[bn]]; x <- expose(win(C, b[1], b[2]))
  rr <- lapply(OUTCOMES, cl_fit, x = x)
  cat(sprintf("  %-24s %s\n", bn,
      paste(sapply(rr, function(r) sprintf("%5.2f (%4.2f-%5.2f)", r[["or"]], r[["lo"]], r[["hi"]])),
            collapse = "  ")))
}

# ========================================================= 3. permutation nulls
# The asymptotic CI treats 1,862 samples as the unit of information. The actual
# information is ~47 case-fortnights in 16 catchments, and the case series is
# strongly clustered in time. These nulls preserve that structure and destroy
# only the feature under test.
#
#   A  independent circular time shift per catchment -- preserves each
#      catchment's burden and its within-catchment clustering exactly; destroys
#      alignment with phage timing AND with study-wide seasonality (liberal).
#   B  common circular shift applied to every catchment -- additionally
#      preserves cross-catchment synchrony, i.e. seasonality (conservative).
#      Only nf-1 distinct shifts exist, so p is granular to ~1/nf.
#   C  whole case series reassigned to a different catchment -- preserves every
#      temporal shape but breaks the catchment attribution.
shift_indep  <- function(Cm) t(apply(Cm, 1, function(r) { k <- sample.int(nf - 1, 1)
                                                          r[c((k + 1):nf, 1:k)] }))
shift_common <- function(Cm, k) Cm[, c((k + 1):nf, 1:k), drop = FALSE]
perm_site    <- function(Cm) Cm[sample.int(nrow(Cm)), , drop = FALSE]

run_null <- function(gen, nrep, lo, hi, label, obs) {
  M <- matrix(NA_real_, nrep, length(OUTCOMES), dimnames = list(NULL, names(OUTCOMES)))
  for (r in seq_len(nrep)) {
    x <- expose(win(gen(r), lo, hi))
    M[r, ] <- vapply(OUTCOMES, cl_beta, numeric(1), x = x)
  }
  cat("  ", label, "\n", sep = "")
  out <- setNames(numeric(length(OUTCOMES)), names(OUTCOMES))
  for (nm in names(OUTCOMES)) {
    v <- M[, nm]; v <- v[is.finite(v)]
    p <- (1 + sum(v >= obs[[nm]])) / (1 + length(v))
    out[[nm]] <- p
    cat(sprintf("    %-14s null median OR %.2f (95%% range %.2f-%.2f) | observed %.2f | p=%.4f\n",
                nm, exp(median(v)), exp(quantile(v, .025)), exp(quantile(v, .975)),
                exp(obs[[nm]]), p))
  }
  out
}

for (sp in list(list(nm = "PAST 8 WEEKS",   lo = -3, hi = 0),
                list(nm = "FUTURE 8 WEEKS", lo =  1, hi = 4))) {
  cat(sprintf("\n=== 3. PERMUTATION NULLS -- %s ===\n", sp$nm))
  obs <- lapply(OUTCOMES, cl_beta, x = expose(win(C, sp$lo, sp$hi)))
  pA <- run_null(function(r) shift_indep(C), NREP, sp$lo, sp$hi,
                 sprintf("A. independent circular shift per catchment (%d reps)", NREP), obs)
  pB <- run_null(function(r) shift_common(C, 1 + ((r - 1) %% (nf - 1))), nf - 1, sp$lo, sp$hi,
                 sprintf("B. common circular shift, all catchments (%d exact shifts)", nf - 1), obs)
  pC <- run_null(function(r) perm_site(C), NREP, sp$lo, sp$hi,
                 sprintf("C. case series reassigned to another catchment (%d reps)", NREP), obs)
  tag <- if (sp$lo < 0) "past 8wk" else "future 8wk"
  for (nm in names(OUTCOMES)) {
    r <- cl_fit(OUTCOMES[[nm]], expose(win(C, sp$lo, sp$hi)))
    add(analysis = paste(tag, "(permutation)"), outcome = nm, or = r[["or"]], lo = r[["lo"]],
        hi = r[["hi"]], p_asymptotic = r[["p"]], p_perm_indep = pA[[nm]], p_perm_common = pB[[nm]])
  }
}

# ================================================ 4. negative-control OUTCOME
# HF183 is a human faecal marker: it should track sewage strength and dilution
# but carries no information about typhoid specifically. If case periods were
# simply marking intervals of more concentrated sewage, they would predict it.
cat("\n=== 4. NEGATIVE-CONTROL OUTCOME: HF183 human faecal marker, and flow ===\n")
d$x_past <- x_past; d$x_fut <- x_fut
for (y in c("ln183", "flow")) for (ex in c("x_past", "x_fut")) {
  f <- as.formula(sprintf("%s ~ %s + (1|site) + (1|fn)", y, ex))
  m <- try(lmer(f, data = d), silent = TRUE)
  if (!inherits(m, "try-error")) {
    co <- summary(m)$coefficients
    cat(sprintf("  %-6s ~ %-7s beta %+8.3f (%+8.3f to %+8.3f) t=%+5.2f  %s\n",
                y, sub("x_", "", ex), co[2, 1], co[2, 1] - 1.96 * co[2, 2],
                co[2, 1] + 1.96 * co[2, 2], co[2, 3],
                if (abs(co[2, 3]) > 1.96) "<<< NOT NULL" else "null (as required)"))
  }
}
cat("\n  HF183 is nonetheless a live predictor of phage positivity:\n")
for (nm in names(OUTCOMES)) {
  dd <- d[!is.na(d[[OUTCOMES[[nm]]]]) & !is.na(d$ln183), ]
  co <- summary(glm(dd[[OUTCOMES[[nm]]]] ~ dd$ln183, family = binomial))$coefficients
  cat(sprintf("    %-14s OR %.3f per log unit (%.3f-%.3f) p=%.4f\n", nm,
              exp(co[2, 1]), exp(co[2, 1] - 1.96 * co[2, 2]),
              exp(co[2, 1] + 1.96 * co[2, 2]), co[2, 4]))
}
cat("\n  past-8wk association adjusted for HF183 and flow:\n")
for (nm in names(OUTCOMES)) {
  dd <- d[!is.na(d[[OUTCOMES[[nm]]]]) & !is.na(d$ln183) & !is.na(d$flow), ]
  m <- try(clogit(dd[[OUTCOMES[[nm]]]] ~ dd$x_past + dd$ln183 + dd$flow + strata(dd$site)),
           silent = TRUE)
  if (!inherits(m, "try-error")) {
    s <- summary(m)$coefficients
    cat(sprintf("    %-14s OR %5.2f (%4.2f-%5.2f)  [n=%d]\n", nm, exp(s[1, 1]),
                exp(s[1, 1] - 1.96 * s[1, 3]), exp(s[1, 1] + 1.96 * s[1, 3]), nrow(dd)))
  }
}

# ===================================================================== output
res <- do.call(rbind, results)
outfile <- file.path(OUTDIR, "temporal_inference_results.csv")
write.csv(res, outfile, row.names = FALSE)
cat("\nWrote", outfile, "(", nrow(res), "rows )\n")
cat("\nSeed:", SEED, "| permutation reps:", NREP, "\n")
cat("R:", R.version.string, "| survival", as.character(packageVersion("survival")),
    "| lme4", as.character(packageVersion("lme4")),
    "| readxl", as.character(packageVersion("readxl")), "\n")

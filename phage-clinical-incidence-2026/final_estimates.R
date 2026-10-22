# ==============================================================================
# Typhoid R01 phage manuscript — estimates reported in the Results
#
# One specification throughout: conditional logistic regression of sample-level
# phage positivity, stratified by catchment, with calendar time adjusted by a
# natural cubic spline in fortnight (4 df). Stratification removes all
# between-catchment variation, so every estimate is a within-catchment contrast.
#
# Usage (from the repository root):
#   Rscript phage-clinical-incidence-2026/final_estimates.R [path/to/dataset.xlsx]
# ==============================================================================
suppressPackageStartupMessages({library(readxl); library(survival); library(splines); library(lme4)})
a <- commandArgs(trailingOnly = TRUE)
XLSX <- if (length(a) >= 1) a[1] else "data/dataset_WWphage_R01.xlsx"
DF <- 4L
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
raw <- as.data.frame(read_excel(XLSX, sheet = "Sheet1", guess_max = 5000))

# amplification_tw / direct_tw are coded 2 = no, 1 = yes, NOT 0/1 as the data
# dictionary states. Stop rather than silently invert if this ever changes.
stopifnot(identical(sort(unique(na.omit(num(raw$amplification_tw)))), c(1, 2)),
          identical(sort(unique(na.omit(num(raw$direct_tw)))), c(1, 2)))
stopifnot(length(unique(num(raw$ca_overlap))) == 1L)   # level-1 catchments only

d <- data.frame(site = raw$ESsite, fn = num(raw$fn_num),
  amp    = as.integer(num(raw$amplification_tw) == 1),
  direct = ifelse(is.na(num(raw$direct_tw)), NA, as.integer(num(raw$direct_tw) == 1)),
  either = as.integer(num(raw$positive) == 1),
  cases  = num(raw$cases_fn), pop = num(raw$catchment_pop),
  ln183  = num(raw$ln_hf183), flow = num(raw$Flow_PerHourM3))
sites <- sort(unique(d$site)); fns <- min(d$fn):max(d$fn)
C <- matrix(0, length(sites), length(fns), dimnames = list(sites, fns))
ag <- aggregate(cases ~ site + fn, d, max)
C[cbind(as.character(ag$site), as.character(ag$fn))] <- ag$cases
d$si <- match(d$site, sites); d$fi <- match(d$fn, fns)
win <- function(lo, hi) { o <- matrix(FALSE, nrow(C), ncol(C))
  for (j in seq_len(ncol(C))) { k <- (j + lo):(j + hi); k <- k[k >= 1 & k <= ncol(C)]
    if (length(k)) o[, j] <- rowSums(C[, k, drop = FALSE]) > 0 }; o[cbind(d$si, d$fi)] }
d$w0 <- win(0,0); d$w4 <- win(-1,0); d$w8 <- win(-3,0); d$w12 <- win(-5,0)
d$pre816 <- win(-7,-4); d$next8 <- win(1,4)

rows <- list()
est <- function(oc, rhs, analysis, outcome, idx = 1) {
  vars <- c(oc, all.vars(as.formula(paste("~", rhs))), "site")
  dd <- d[stats::complete.cases(d[, vars]), ]
  m <- clogit(as.formula(sprintf("%s ~ %s + ns(fn, %d) + strata(site)", oc, rhs, DF)), data = dd)
  s <- summary(m)$coefficients
  r <- data.frame(analysis = analysis, outcome = outcome, n = nrow(dd),
                  or = exp(s[idx,1]), lo = exp(s[idx,1] - 1.96*s[idx,3]),
                  hi = exp(s[idx,1] + 1.96*s[idx,3]), p = s[idx,5])
  rows[[length(rows)+1]] <<- r
  cat(sprintf("  %-34s %-8s %5.2f (%4.2f-%5.2f)  p=%s\n", analysis, outcome,
      r$or, r$lo, r$hi, format.pval(r$p, digits = 2, eps = 1e-4)))
  invisible(r)
}
LAB <- c(amp = "amplif.", direct = "direct", either = "either")
cat(sprintf("38 non-overlapping catchments | %d samples | %d cases | %d catchments with cases | %d case-fortnights\n",
            nrow(d), sum(d$cases), sum(tapply(d$cases, d$site, sum) > 0), sum(C > 0)))
cat(sprintf("positivity: amplification %d/%d (%.1f%%), direct %d/%d (%.1f%%), either %d/%d (%.1f%%)\n\n",
  sum(d$amp), nrow(d), 100*mean(d$amp), sum(d$direct, na.rm=TRUE), sum(!is.na(d$direct)),
  100*mean(d$direct, na.rm=TRUE), sum(d$either), nrow(d), 100*mean(d$either)))

cat("=== Exposure window ===\n")
for (v in list(c("w0","same fortnight"), c("w4","case in past 4 weeks"),
               c("w8","case in past 8 weeks"), c("w12","case in past 12 weeks"),
               c("pre816","case 8-16 weeks before"), c("next8","case in next 8 weeks")))
  for (oc in names(LAB)) est(oc, v[1], v[2], LAB[[oc]])

cat("\n=== Case number WITHIN the 8-week window (gradient above one case) ===\n")
# number of cases in the same window that defines the primary exposure, so the
# reference group is fortnights with no case in the preceding 8 weeks
nwin <- matrix(0, nrow(C), ncol(C))
for (j in seq_len(ncol(C))) { k <- (j-3):j; k <- k[k >= 1 & k <= ncol(C)]
  nwin[, j] <- rowSums(C[, k, drop = FALSE]) }
d$ncase8 <- nwin[cbind(d$si, d$fi)]
cat("  samples by cases in window:", paste(sprintf("%s:%d", names(table(d$ncase8)),
    as.integer(table(d$ncase8))), collapse = " "), "\n")
d$n1 <- d$ncase8 == 1; d$n2 <- d$ncase8 == 2; d$n3 <- d$ncase8 >= 3
# each indicator contrasted against zero cases in the window
for (oc in c("amp","either")) for (v in c("n1","n2","n3")) {
  dd <- d[d$ncase8 == 0 | d[[v]], ]
  vars <- c(oc, v, "fn", "site"); dd <- dd[stats::complete.cases(dd[, vars]), ]
  m <- clogit(as.formula(sprintf("%s ~ %s + ns(fn,%d) + strata(site)", oc, v, DF)), data = dd)
  sm <- summary(m)$coefficients
  lbl <- c(n1 = "1 case in past 8 wk", n2 = "2 cases in past 8 wk", n3 = "3+ cases in past 8 wk")[[v]]
  rows[[length(rows)+1]] <- data.frame(analysis = lbl, outcome = LAB[[oc]], n = nrow(dd),
    or = exp(sm[1,1]), lo = exp(sm[1,1]-1.96*sm[1,3]), hi = exp(sm[1,1]+1.96*sm[1,3]), p = sm[1,5])
  cat(sprintf("  %-34s %-8s %5.2f (%4.2f-%5.2f)  p=%s\n", lbl, LAB[[oc]], exp(sm[1,1]),
      exp(sm[1,1]-1.96*sm[1,3]), exp(sm[1,1]+1.96*sm[1,3]), format.pval(sm[1,5], digits=2, eps=1e-4)))
}

cat("\n=== Past 8 weeks, adjusted ===\n")
for (oc in names(LAB)) est(oc, "w8 + ln183 + flow", "past 8 weeks, adj HF183 + flow", LAB[[oc]])
for (oc in names(LAB)) est(oc, "ln183", "log HF183, per unit", LAB[[oc]])

cat("\n=== Sensitivity of the past-8-week estimate to the time adjustment ===\n")
for (k in c(0, 3, 4, 5, 8)) for (oc in c("amp","either")) {
  rhs <- if (k == 0) "w8" else sprintf("w8 + ns(fn,%d)", k)
  m <- clogit(as.formula(sprintf("%s ~ %s + strata(site)", oc, rhs)), data = d[!is.na(d[[oc]]),])
  s <- summary(m)$coefficients
  cat(sprintf("  spline df=%-2s %-8s %5.2f (%4.2f-%5.2f)\n", if (k==0) "none" else k,
      LAB[[oc]], exp(s[1,1]), exp(s[1,1]-1.96*s[1,3]), exp(s[1,1]+1.96*s[1,3])))
}

cat("\n=== Variance decomposition (mixed-effects model, either assay) ===\n")
d$fnf <- factor(d$fn)
m <- glmer(either ~ w8 + (1|site) + (1|fnf), data = d, family = binomial,
           control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))
v <- as.data.frame(VarCorr(m)); tot <- sum(v$vcov)
for (i in seq_len(nrow(v)))
  cat(sprintf("  %-6s variance %.2f (%.0f%% of random variance)\n",
      ifelse(v$grp[i]=="site","catchment","fortnight"), v$vcov[i], 100*v$vcov[i]/tot))

res <- do.call(rbind, rows)
write.csv(res, "final_estimates.csv", row.names = FALSE)
cat("\nWrote final_estimates.csv (", nrow(res), "rows )\n")
cat("R:", R.version.string, "| survival", as.character(packageVersion("survival")),
    "| lme4", as.character(packageVersion("lme4")), "\n")

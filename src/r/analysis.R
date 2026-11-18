# analysis.R -- Option B: Bootstrapping
# Run from the project folder:  Rscript src/analysis.R

source("src/r/bootstrap.R")
set.seed(2026)

# ---- Data: complete cases --------------------------------------------------
nhefs <- read.csv("data/nhefs_2026.csv")
d <- nhefs[!is.na(nhefs$wt82_71), ]
for (v in c("sex", "race", "education", "exercise", "active")) d[[v]] <- factor(d[[v]])

fmt_ci <- function(ci) paste(sprintf("%.2f", ci), collapse = " to ")

# ---- Q1: difference in median weight change ---------------------------------
median_diff <- function(df) {
  median(df$wt82_71[df$qsmk == 1]) - median(df$wt82_71[df$qsmk == 0])
}

d2 <- d[, c("wt82_71", "qsmk")]   # only the columns we need (faster)
b_med     <- bootstrap(d2, median_diff, B = 10000)
b_med_str <- bootstrap(d2, median_diff, B = 10000, strata = "qsmk")

types <- c("normal", "basic", "percentile", "bca")
q1 <- data.frame(
  type       = types,
  ordinary   = sapply(types, function(t) fmt_ci(boot_ci(b_med, t))),
  stratified = sapply(types, function(t) fmt_ci(boot_ci(b_med_str, t))),
  row.names  = NULL
)

# ---- Q2: regression ----------------------------------------------------------
model <- wt82_71 ~ qsmk + sex + race + age + education + smokeintensity +
  smokeyrs + exercise + active + wt71
fit <- lm(model, data = d)

b_case  <- bootstrap(d, function(df) coef(lm(model, data = df))[["qsmk"]], B = 5000)
b_resid <- residual_bootstrap(fit, "qsmk", B = 5000)

q2 <- data.frame(
  method = c("Model-based", "Case bootstrap (random X)", "Residual bootstrap (fixed X)"),
  se     = c(summary(fit)$coefficients["qsmk", "Std. Error"], boot_se(b_case), boot_se(b_resid)),
  ci     = c(fmt_ci(confint(fit)["qsmk", ]),
             fmt_ci(boot_ci(b_case)),
             fmt_ci(boot_ci(b_resid)))
)

# residual spread in each exposure group
resid_sd <- tapply(resid(fit), d$qsmk, sd)

# ---- Q3: Monte Carlo error ---------------------------------------------------
# Repeat the Q1 bootstrap 20 times for each B and see how much results move.
Bs <- c(100, 500, 1000, 2000, 5000)
mc <- do.call(rbind, lapply(Bs, function(B) {
  t(replicate(20, {
    b <- bootstrap(d2, median_diff, B = B)
    c(B = B, se = boot_se(b), lower = boot_ci(b)[1], upper = boot_ci(b)[2])
  }))
}))
mc <- as.data.frame(mc)

q3 <- aggregate(cbind(se, lower, upper) ~ B, data = mc, FUN = sd)
names(q3) <- c("B", "sd_se", "sd_lower", "sd_upper")

# ---- Print results -----------------------------------------------------------
if (sys.nframe() == 0) {
  cat("Complete cases:", nrow(d), "of", nrow(nhefs), "\n\n")
  cat("Q1: median difference =", b_med$t0, " bootstrap SE =", round(boot_se(b_med), 3), "\n")
  print(q1)
  cat("\nQ2: qsmk coefficient =", round(coef(fit)[["qsmk"]], 3), "\n")
  print(q2)
  cat("\nResidual SD by qsmk:\n"); print(resid_sd)
  cat("\nQ3: Monte Carlo SD across 20 repeats\n"); print(q3)
}

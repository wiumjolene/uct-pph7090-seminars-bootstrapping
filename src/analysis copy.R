###############################################################################
# PPH7090S - Final assignment: Bootstrapping
# Jolene Wium | 18 November 2026
###############################################################################

library(boot)

# ----- Colour scheme ----------------------------------------------------------
uct_cols <- c(
  navy    = "#053259",  # primary: titles, axes
  blue    = "#1668A3",  # first series
  teal    = "#12999E",  # second series
  magenta = "#C32A66",  # callouts / estimates
  steel   = "#7196B7",  # neutral: reference lines, borders
  amber   = "#FCB028",  # highlights
  coral   = "#FA615C"   # alerts
)
par(col.main = uct_cols["navy"], col.lab = uct_cols["navy"],
    col.axis = uct_cols["navy"], fg = uct_cols["navy"])

# ----- Data -------------------------------------------------------------------
nhefs <- read.csv("data/nhefs_2026.csv")
dat <- subset(nhefs, !is.na(wt82_71))

n_total <- nrow(dat)
n_quit <- sum(dat$qsmk == 1)
n_cont <- sum(dat$qsmk == 0)

med_quit <- median(dat$wt82_71[dat$qsmk == 1])
med_cont <- median(dat$wt82_71[dat$qsmk == 0])
med_diff <- med_quit - med_cont

cat("Complete cases:", n_total, "\n")
cat("Quitters:", n_quit, "| Smokers:", n_cont, "\n")
cat("Difference in medians (kg):", round(med_diff, 3), "\n")

# ----- Missing data in the outcome (wt82_71) ----------------------------------
# Compare participants with and without a 1982 weight on baseline
# characteristics, to judge whether complete-case analysis is reasonable.

nhefs$wt_missing <- as.integer(is.na(nhefs$wt82_71))
n_missing <- sum(nhefs$wt_missing)
pct_missing <- 100 * mean(nhefs$wt_missing)
cat("Missing wt82_71:", n_missing, "of", nrow(nhefs),
    sprintf("(%.1f%%)", pct_missing), "\n")

# Continuous baseline variables: means, standardised mean difference, Welch t-test
miss_cont_vars <- c("age", "wt71", "ht", "smokeintensity", "smokeyrs",
                    "income", "cholesterol", "sbp", "dbp")

miss_cont <- do.call(rbind, lapply(miss_cont_vars, function(v) {
  x_miss <- nhefs[[v]][nhefs$wt_missing == 1]
  x_obs  <- nhefs[[v]][nhefs$wt_missing == 0]
  smd <- (mean(x_miss, na.rm = TRUE) - mean(x_obs, na.rm = TRUE)) /
    sqrt((var(x_miss, na.rm = TRUE) + var(x_obs, na.rm = TRUE)) / 2)
  data.frame(Variable = v,
             Mean_Missing = mean(x_miss, na.rm = TRUE),
             Mean_Observed = mean(x_obs, na.rm = TRUE),
             SMD = smd,
             p = t.test(x_miss, x_obs)$p.value)
}))
print(miss_cont, digits = 3, row.names = FALSE)

# Categorical baseline variables: % missing in each level, chi-squared test
miss_cat_vars <- c("qsmk", "sex", "race", "education", "exercise",
                   "active", "death")

miss_cat <- do.call(rbind, lapply(miss_cat_vars, function(v) {
  tab <- table(nhefs[[v]], nhefs$wt_missing)
  data.frame(Variable = v,
             Level = rownames(tab),
             N = as.vector(rowSums(tab)),
             N_Missing = as.vector(tab[, "1"]),
             Pct_Missing = as.vector(100 * tab[, "1"] / rowSums(tab)),
             p = suppressWarnings(chisq.test(tab)$p.value))
}))
print(miss_cat, digits = 3, row.names = FALSE)

# Fisher's exact test for the two most important binary comparisons
fisher_qsmk  <- fisher.test(table(nhefs$qsmk,  nhefs$wt_missing))$p.value
fisher_death <- fisher.test(table(nhefs$death, nhefs$wt_missing))$p.value
cat("Fisher p (quitting):", signif(fisher_qsmk, 3),
    "| Fisher p (death):", signif(fisher_death, 3), "\n")

# Share of the missing participants who died during follow-up
n_missing_died <- sum(nhefs$wt_missing == 1 & nhefs$death == 1)
cat("Missing participants who died:", n_missing_died, "of", n_missing, "\n")

# Joint model: which outcome-model covariates predict missingness?
miss_dat <- nhefs
miss_dat[c("sex", "race", "education", "exercise", "active")] <-
  lapply(miss_dat[c("sex", "race", "education", "exercise", "active")], factor)

miss_model <- glm(wt_missing ~ qsmk + sex + race + age + education +
                    smokeintensity + smokeyrs + exercise + active + wt71,
                  family = binomial, data = miss_dat)
print(summary(miss_model))

miss_or <- data.frame(
  Term = names(coef(miss_model)),
  OR = exp(coef(miss_model)),
  CI_Lower = exp(confint.default(miss_model)[, 1]),
  CI_Upper = exp(confint.default(miss_model)[, 2]),
  p = coef(summary(miss_model))[, "Pr(>|z|)"],
  row.names = NULL
)[-1, ]
print(miss_or, digits = 3, row.names = FALSE)

# Likelihood-ratio test against the intercept-only model
miss_lr_p <- with(miss_model,
                  pchisq(null.deviance - deviance, df.null - df.residual,
                         lower.tail = FALSE))
cat("LR test, covariates jointly predict missingness: p =",
    signif(miss_lr_p, 3), "\n")

# Plot: % missing by key groups
miss_plot <- c(
  "Continued smoking" = mean(nhefs$wt_missing[nhefs$qsmk == 0]),
  "Quit smoking"      = mean(nhefs$wt_missing[nhefs$qsmk == 1]),
  "Survived"          = mean(nhefs$wt_missing[nhefs$death == 0]),
  "Died"              = mean(nhefs$wt_missing[nhefs$death == 1]),
  "Age < 50"          = mean(nhefs$wt_missing[nhefs$age < 50]),
  "Age 50+"           = mean(nhefs$wt_missing[nhefs$age >= 50])
) * 100

barplot(miss_plot, horiz = TRUE, las = 1, xlim = c(0, max(miss_plot) * 1.15),
        col = rep(uct_cols[c("blue", "teal", "magenta")], each = 2),
        border = NA, xlab = "% missing wt82_71",
        main = "Missing 1982 weight by participant group")
abline(v = pct_missing, lty = 2, col = uct_cols["steel"], lwd = 2)

# ----- 1. Difference in medians -----------------------------------------------
median_diff_fn <- function(data, indices) {
  d <- data[indices, ]
  median(d$wt82_71[d$qsmk == 1]) -
    median(d$wt82_71[d$qsmk == 0])
}

set.seed(123)
median_boot <- boot(dat, median_diff_fn, R = 5000, strata = dat$qsmk)
median_se <- sd(median_boot$t[, 1])
median_ci <- boot.ci(median_boot, conf = 0.95, type = c("perc", "bca"))

print(median_boot)
cat("Bootstrap SE (median difference):", round(median_se, 3), "\n")
print(median_ci)

hist(median_boot$t[, 1], breaks = 40, col = uct_cols["blue"],
     border = "white", main = "Bootstrap distribution: median difference",
     xlab = "Difference in median weight change (kg)")
abline(v = median_boot$t0, col = uct_cols["magenta"], lwd = 2)
abline(v = 0, lty = 2, col = uct_cols["steel"])

# ----- 2. Regression ----------------------------------------------------------
cat_vars <- c("sex", "race", "education", "exercise", "active")
dat[cat_vars] <- lapply(dat[cat_vars], factor)

model_formula <- wt82_71 ~ qsmk + sex + race + age + education +
  smokeintensity + smokeyrs + exercise + active + wt71

model <- lm(model_formula, data = dat)
print(summary(model))

# Case bootstrap: resample participants and refit the model.
regression_coef_fn <- function(data, indices) {
  fit <- lm(model_formula, data = data[indices, ])
  unname(coef(fit)["qsmk"])
}

set.seed(123)
case_boot <- boot(dat, regression_coef_fn, R = 5000)
case_se <- sd(case_boot$t[, 1])
case_ci <- boot.ci(case_boot, type = "perc")$percent[4:5]

# Residual bootstrap: retain the predictors and resample centred residuals.
fitted_y <- fitted(model)
model_resid <- resid(model)
X <- model.matrix(model)

set.seed(456)
B <- 5000
resid_coef <- numeric(B)

for (b in seq_len(B)) {
  sampled_resid <- sample(model_resid - mean(model_resid), replace = TRUE)
  y_star <- fitted_y + sampled_resid
  fit_star <- lm.fit(x = X, y = y_star)
  resid_coef[b] <- fit_star$coefficients["qsmk"]
}

resid_se <- sd(resid_coef)
resid_ci <- quantile(resid_coef, c(0.025, 0.975))
ols_ci <- confint(model)["qsmk", ]

comparison <- data.frame(
  Method = c("Conventional OLS", "Case bootstrap", "Residual bootstrap"),
  Estimate = rep(unname(coef(model)["qsmk"]), 3),
  SE = c(coef(summary(model))["qsmk", "Std. Error"], case_se, resid_se),
  CI_Lower = c(ols_ci[1], case_ci[1], resid_ci[1]),
  CI_Upper = c(ols_ci[2], case_ci[2], resid_ci[2])
)
print(comparison, digits = 3, row.names = FALSE)

# Residual diagnostics
plot(fitted(model), resid(model),
     pch = 19, cex = 0.6, col = adjustcolor(uct_cols["blue"], alpha.f = 0.5),
     xlab = "Fitted weight change (kg)", ylab = "Residual (kg)",
     main = "Residuals versus fitted values")
abline(h = 0, col = uct_cols["magenta"], lty = 2, lwd = 2)

boxplot(resid(model) ~ dat$qsmk,
        names = c("Non-quitters", "Quitters"),
        col = uct_cols[c("blue", "teal")], border = uct_cols["navy"],
        ylab = "Residual (kg)", main = "Residuals by smoking status")
print(tapply(resid(model), dat$qsmk, sd))

# ----- 3. Monte Carlo error ---------------------------------------------------
# B_values <- c(200, 1000, 5000, 10000)
# n_runs <- 20
# mc_results <- data.frame()

# for (B in B_values) {
#   for (run in seq_len(n_runs)) {
#     set.seed(100000 + B + run)
#     result <- boot(dat, regression_coef_fn, R = B)
#     estimates <- result$t[, 1]
#     ci <- quantile(estimates, c(0.025, 0.975))

#     mc_results <- rbind(mc_results, data.frame(
#       B = B, Run = run, SE = sd(estimates),
#       CI_Lower = unname(ci[1]), CI_Upper = unname(ci[2])
#     ))
#   }
# }
# print(mc_results, digits = 4, row.names = FALSE)

# # Spread across the three runs for each value of B
# mc_summary <- aggregate(cbind(SE, CI_Lower, CI_Upper) ~ B,
#                         data = mc_results,
#                         FUN = function(x) diff(range(x)))
# names(mc_summary) <- c("B", "SE_Range", "Lower_CI_Range", "Upper_CI_Range")
# print(mc_summary, digits = 4, row.names = FALSE)

# plot(NA, xlim = range(B_values), ylim = range(mc_results$SE), log = "x",
#      xlab = "Bootstrap replicates (B)", ylab = "Bootstrap SE (kg)",
#      main = "Monte Carlo stability of the case-bootstrap SE")

# cols <- rep_len(uct_cols, n_runs)
# for (run in seq_len(n_runs)) {
#   d <- subset(mc_results, Run == run)
#   lines(d$B, d$SE, type = "b", pch = 19, col = cols[run])
# }
# legend("topright", legend = paste("Run", seq_len(n_runs)),
#        col = cols, lty = 1, pch = 19, bty = "n")
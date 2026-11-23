###############################################################################
# PPH7090S - Final assignment: Bootstrapping
# Jolene Wium | 18 November 2026
###############################################################################

library(boot)

# ----- Colour scheme ----------------------------------------------------------
uct_cols <- c(
  navy    = "#053259",
  blue    = "#1668A3",
  teal    = "#12999E",
  magenta = "#C32A66",
  steel   = "#7196B7",
  amber   = "#FCB028",
  coral   = "#FA615C"
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
B_values <- c(200, 1000, 5000, 10000)
n_runs <- 20
mc_results <- data.frame()

for (B in B_values) {
  for (run in seq_len(n_runs)) {
    set.seed(100000 + B + run)
    result <- boot(dat, regression_coef_fn, R = B)
    estimates <- result$t[, 1]
    ci <- quantile(estimates, c(0.025, 0.975))

    mc_results <- rbind(mc_results, data.frame(
      B = B, Run = run, SE = sd(estimates),
      CI_Lower = unname(ci[1]), CI_Upper = unname(ci[2])
    ))
  }
}
print(mc_results, digits = 4, row.names = FALSE)

# Spread across the 20 runs for each value of B
mc_summary <- aggregate(cbind(SE, CI_Lower, CI_Upper) ~ B,
                        data = mc_results,
                        FUN = function(x) diff(range(x)))
names(mc_summary) <- c("B", "SE_Range", "Lower_CI_Range", "Upper_CI_Range")
print(mc_summary, digits = 4, row.names = FALSE)

plot(NA, xlim = range(B_values), ylim = range(mc_results$SE), log = "x",
     xlab = "Bootstrap replicates (B)", ylab = "Bootstrap SE (kg)",
     main = "Monte Carlo stability of the case-bootstrap SE")

cols <- rep_len(uct_cols, n_runs)
for (run in seq_len(n_runs)) {
  d <- subset(mc_results, Run == run)
  lines(d$B, d$SE, type = "b", pch = 19, col = cols[run])
}
legend("topright", legend = paste("Run", seq_len(n_runs)),
       col = cols, lty = 1, pch = 19, bty = "n")
###############################################################################
# PPH7090S - Final assignment: Bootstrapping
# Missing data in the outcome (wt82_71)
# Jolene Wium | 22 November 2026
###############################################################################

# Compare participants with and without a 1982 weight on baseline
# characteristics, to judge whether complete-case analysis is reasonable.

nhefs <- read.csv("data/nhefs_2026.csv")
nhefs$wt_missing <- as.integer(is.na(nhefs$wt82_71))

n_missing <- sum(nhefs$wt_missing)
pct_missing <- 100 * mean(nhefs$wt_missing)
cat("Missing wt82_71:", n_missing, "of", nrow(nhefs),
    sprintf("(%.1f%%)", pct_missing), "\n")

# Missing participants who died during follow-up
n_missing_died <- sum(nhefs$wt_missing == 1 & nhefs$death == 1)
cat("Missing participants who died:", n_missing_died, "of", n_missing, "\n")

# Continuous baseline variables: means, standardised mean difference, Welch t-test
miss_cont_vars <- c("age", "wt71", "ht", "smokeintensity", "smokeyrs",
                    "cholesterol", "sbp", "dbp")

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

# Death during follow-up (binary): n (%), SMD for proportions, Fisher's exact test
p_died_miss <- mean(nhefs$death[nhefs$wt_missing == 1])
p_died_obs  <- mean(nhefs$death[nhefs$wt_missing == 0])
miss_death <- data.frame(
  n_Missing  = sum(nhefs$death[nhefs$wt_missing == 1]),
  n_Observed = sum(nhefs$death[nhefs$wt_missing == 0]),
  pct_Missing  = 100 * p_died_miss,
  pct_Observed = 100 * p_died_obs,
  SMD = (p_died_miss - p_died_obs) /
    sqrt((p_died_miss * (1 - p_died_miss) + p_died_obs * (1 - p_died_obs)) / 2),
  p = fisher.test(table(nhefs$wt_missing, nhefs$death))$p.value
)
print(miss_death, digits = 3, row.names = FALSE)

# Income (ordinal category code, not a continuous amount): median [IQR],
# Wilcoxon rank-sum test. No SMD, since a mean of category codes is not meaningful.
inc_miss <- nhefs$income[nhefs$wt_missing == 1]
inc_obs  <- nhefs$income[nhefs$wt_missing == 0]
miss_income <- data.frame(
  Median_Missing  = median(inc_miss, na.rm = TRUE),
  Q1_Missing      = quantile(inc_miss, 0.25, na.rm = TRUE, names = FALSE),
  Q3_Missing      = quantile(inc_miss, 0.75, na.rm = TRUE, names = FALSE),
  Median_Observed = median(inc_obs, na.rm = TRUE),
  Q1_Observed     = quantile(inc_obs, 0.25, na.rm = TRUE, names = FALSE),
  Q3_Observed     = quantile(inc_obs, 0.75, na.rm = TRUE, names = FALSE),
  p = wilcox.test(inc_miss, inc_obs, exact = FALSE)$p.value
)
print(miss_income, digits = 3, row.names = FALSE)

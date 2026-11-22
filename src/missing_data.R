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

# bootstrap.R -- a minimal bootstrap module (base R only)

# Draw B bootstrap samples of `data` and apply `statistic` to each.
# If `strata` is a column name, rows are resampled within each group.
bootstrap <- function(data, statistic, B = 2000, strata = NULL) {
  n <- nrow(data)
  draw <- function() {
    if (is.null(strata)) return(sample(n, n, replace = TRUE))
    unlist(lapply(split(1:n, data[[strata]]),
                  function(i) i[sample(length(i), replace = TRUE)]))
  }
  t_star <- replicate(B, statistic(data[draw(), ]))
  list(t0 = statistic(data), t_star = t_star, data = data, statistic = statistic)
}

boot_se <- function(b) sd(b$t_star)

# 95% (or other level) confidence intervals from a bootstrap() result.
boot_ci <- function(b, type = "percentile", level = 0.95) {
  a  <- (1 - level) / 2
  t0 <- b$t0
  ts <- b$t_star
  q  <- function(p) unname(quantile(ts, p))

  if (type == "normal") {
    bias <- mean(ts) - t0
    return(t0 - bias + c(-1, 1) * qnorm(1 - a) * sd(ts))
  }
  if (type == "percentile") return(q(c(a, 1 - a)))
  if (type == "basic")      return(2 * t0 - q(c(1 - a, a)))
  if (type == "bca") {
    # bias correction: share of replicates below the estimate
    z0 <- qnorm(mean(ts < t0) + 0.5 * mean(ts == t0))
    # acceleration: from the jackknife (leave one row out at a time)
    jack <- sapply(1:nrow(b$data), function(i) b$statistic(b$data[-i, ]))
    d <- mean(jack) - jack
    acc <- sum(d^3) / (6 * sum(d^2)^1.5)
    z <- qnorm(c(a, 1 - a))
    return(q(pnorm(z0 + (z0 + z) / (1 - acc * (z0 + z)))))
  }
  stop("type must be normal, percentile, basic or bca")
}

# Residual ("fixed-X") bootstrap for one coefficient of a linear model:
# keep X fixed, add resampled residuals to the fitted values and refit.
residual_bootstrap <- function(fit, term, B = 2000) {
  X   <- model.matrix(fit)
  yhat <- fitted(fit)
  e   <- resid(fit)
  t_star <- replicate(B, {
    y_star <- yhat + sample(e, replace = TRUE)
    lm.fit(X, y_star)$coefficients[term]
  })
  list(t0 = coef(fit)[[term]], t_star = unname(t_star))
}

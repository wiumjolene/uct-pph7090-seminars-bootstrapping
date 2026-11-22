#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: setup
#| include: false
#| cache: true
#| cache.extra: !expr tools::md5sum("src/analysis.R")
source("src/analysis.R")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: tbl-median-boot
#| tbl-cap: 'Bootstrap estimates for the difference in median weight change`#linebreak()`{=typst}`#text(size: 0.85em, style: "italic")[R = 5000, stratified by smoking status]`{=typst}'
#| echo: false

library(knitr)
library(kableExtra)

median_table <- data.frame(
  Method   = c("Percentile", "BCa"),
  Estimate = rep(median_boot$t0, 2),
  Bias     = rep(mean(median_boot$t[, 1]) - median_boot$t0, 2),
  SE       = rep(median_se, 2),
  CI_Lower = c(median_ci$percent[4], median_ci$bca[4]),
  CI_Upper = c(median_ci$percent[5], median_ci$bca[5])
)

median_table$CI <- sprintf("(%.2f, %.2f)", median_table$CI_Lower, median_table$CI_Upper)

tab <- kable(median_table[, c("Method", "Estimate", "Bias", "SE", "CI")],
             digits = 3, align = "lrrrr", booktabs = TRUE,
             col.names = c("CI method", "Estimate (kg)", "Bias", "SE", "95% CI"))

if (knitr::is_html_output() || knitr::is_latex_output()) {
  tab <- tab |>
    kable_styling(full_width = FALSE, latex_options = "HOLD_position") |>
    row_spec(0, bold = TRUE, color = "white", background = "#053259")
}

tab
#
#
#
#
#
#| label: fig-median-boot
#| fig-cap: 'Bootstrap distribution of the difference in median weight change`#linebreak()`{=typst}`#text(size: 0.85em, style: "italic")[R = 5000. Solid line: observed estimate; dashed line: zero]`{=typst}'
#| fig-width: 7
#| fig-height: 4.5
#| echo: false

par(col.main = unname(uct_cols["navy"]), col.lab = unname(uct_cols["navy"]),
    col.axis = unname(uct_cols["navy"]), fg = unname(uct_cols["navy"]))

hist(median_boot$t[, 1], breaks = 40, col = uct_cols["blue"],
     border = "white", main = "", xlim = range(c(0, median_boot$t[, 1])),
     xlab = "Difference in median weight change (kg)")
abline(v = median_boot$t0, col = uct_cols["magenta"], lwd = 2)
abline(v = 0, lty = 2, col = uct_cols["steel"], lwd = 2)
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: tbl-regression-boot
#| tbl-cap: 'Adjusted effect of quitting smoking on weight change`#linebreak()`{=typst}`#text(size: 0.85em, style: "italic")[Conventional vs bootstrap standard errors; R = 5000, percentile 95% CIs for bootstrap methods]`{=typst}'
#| echo: false

comparison_tab <- data.frame(
  Method   = comparison$Method,
  Estimate = sprintf("%.2f", comparison$Estimate),
  SE       = sprintf("%.3f", comparison$SE),
  CI       = sprintf("(%.2f, %.2f)", comparison$CI_Lower, comparison$CI_Upper)
)

knitr::kable(comparison_tab, align = "lrrr", row.names = FALSE,
             col.names = c("Method", "Estimate (kg)", "SE", "95% CI"))
```
#
#
#
#| label: fig-resid
#| fig-cap: "Residual diagnostics for the adjusted regression model."
#| fig-subcap:
#|   - "Residuals versus fitted values"
#|   - "Residuals by smoking status"
#| layout-ncol: 2
#| fig-width: 5
#| fig-height: 4.5
#| echo: false

par(col.main = uct_cols["navy"], col.lab = uct_cols["navy"],
    col.axis = uct_cols["navy"], fg = uct_cols["navy"])

plot(fitted(model), resid(model),
     pch = 19, cex = 0.6, col = adjustcolor(uct_cols["blue"], alpha.f = 0.5),
     xlab = "Fitted weight change (kg)", ylab = "Residual (kg)", main = "")
abline(h = 0, col = uct_cols["magenta"], lty = 2, lwd = 2)

par(col.main = uct_cols["navy"], col.lab = uct_cols["navy"],
    col.axis = uct_cols["navy"], fg = uct_cols["navy"])

boxplot(resid(model) ~ dat$qsmk,
        names = c("Non-quitters", "Quitters"),
        col = uct_cols[c("blue", "teal")], border = uct_cols["navy"],
        xlab = "", ylab = "Residual (kg)", main = "")

#
#
#
#
#
#
#| label: tbl-mc-summary
#| tbl-cap: 'Monte Carlo variability of the case bootstrap by number of replicates`#linebreak()`{=typst}`#text(size: 0.85em, style: "italic")[Range (max − min, kg) across 20 independent runs at each value of $B$]`{=typst}'
#| echo: false

mc_summary_tab <- data.frame(
  B        = format(mc_summary$B, big.mark = "\u00a0"),
  SE       = sprintf("%.4f", mc_summary$SE_Range),
  CI       = sprintf("%.3f, %.3f", mc_summary$Lower_CI_Range,
                     mc_summary$Upper_CI_Range)
)

knitr::kable(mc_summary_tab, align = "rrr", row.names = FALSE,
             col.names = c("$B$", "SE", "95% CI"))
```
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: tbl-baseline
#| tbl-cap: 'Baseline (1971) characteristics by smoking status`#linebreak()`{=typst}`#text(size: 0.85em, style: "italic")[Participants with observed weight change; mean (SD) or n (%)]`{=typst}'
#| echo: false

bl_quit <- dat$qsmk == 1
ind <- "   "

bl_mean <- function(x) {
  c(sprintf("%.1f (%.1f)", mean(x[!bl_quit]), sd(x[!bl_quit])),
    sprintf("%.1f (%.1f)", mean(x[bl_quit]), sd(x[bl_quit])))
}
bl_pct <- function(cond) {
  c(sprintf("%d (%.1f)", sum(cond[!bl_quit]), 100 * mean(cond[!bl_quit])),
    sprintf("%d (%.1f)", sum(cond[bl_quit]), 100 * mean(cond[bl_quit])))
}
bl_cat <- function(x, labels) {
  x <- as.character(x)
  t(sapply(names(labels), function(lv) bl_pct(x == lv)))
}

edu_lab <- c("1" = "8th grade or less", "2" = "Some high school",
             "3" = "High school", "4" = "Some college", "5" = "College or more")
exe_lab <- c("0" = "Much", "1" = "Moderate", "2" = "Little or none")
act_lab <- c("0" = "Very active", "1" = "Moderately active", "2" = "Inactive")

bl_rows <- rbind(
  c("Age (years)",              bl_mean(dat$age)),
  c("Female",                   bl_pct(as.character(dat$sex) == "1")),
  c("Black or other race",      bl_pct(as.character(dat$race) == "1")),
  c("Weight (kg)",              bl_mean(dat$wt71)),
  c("Cigarettes per day",       bl_mean(dat$smokeintensity)),
  c("Years smoked",             bl_mean(dat$smokeyrs)),
  c("Education",                "", ""),
  cbind(paste0(ind, edu_lab),   bl_cat(dat$education, edu_lab)),
  c("Recreational exercise",    "", ""),
  cbind(paste0(ind, exe_lab),   bl_cat(dat$exercise, exe_lab)),
  c("Daily physical activity",  "", ""),
  cbind(paste0(ind, act_lab),   bl_cat(dat$active, act_lab))
)

knitr::kable(as.data.frame(bl_rows), align = "lrr", row.names = FALSE,
             col.names = c("Characteristic",
                           sprintf("Non-quitters (n = %d)", n_cont),
                           sprintf("Quitters (n = %d)", n_quit)))
#
#
#
#
#
#
#
#
#
#| label: tbl-missingness
#| tbl-cap: 'Baseline characteristics of participants with missing and observed weight change`#linebreak()`{=typst}`#text(size: 0.85em, style: "italic")[SMD = standardised mean difference; p-values from Welch t-tests]`{=typst}'
#| echo: false

miss_labels <- c(
  age            = "Age (years)",
  wt71           = "Weight in 1971 (kg)",
  ht             = "Height (cm)",
  smokeintensity = "Cigarettes per day",
  smokeyrs       = "Years smoked",
  income         = "Income category",
  cholesterol    = "Cholesterol (mg/dL)",
  sbp            = "Systolic BP (mmHg)",
  dbp            = "Diastolic BP (mmHg)"
)

miss_cont_tab <- data.frame(
  Variable = miss_labels[miss_cont$Variable],
  Missing  = sprintf("%.1f", miss_cont$Mean_Missing),
  Observed = sprintf("%.1f", miss_cont$Mean_Observed),
  SMD      = sprintf("%.2f", miss_cont$SMD),
  p        = ifelse(miss_cont$p < 0.001, "<0.001", sprintf("%.3f", miss_cont$p))
)

knitr::kable(
  miss_cont_tab, align = "lrrrr", row.names = FALSE,
  col.names = c("Variable",
                sprintf("Missing (n = %d)", n_missing),
                sprintf("Observed (n = %d)", nrow(nhefs) - n_missing),
                "SMD", "p")
)
#
#
#
#
#
#
#
#
#
#
#| label: fig-mc-se
#| fig-cap: 'Monte Carlo stability of the case-bootstrap standard error`#linebreak()`{=typst}`#text(size: 0.85em, style: "italic")[Each line is one of 20 independent runs of the case bootstrap at each value of B]`{=typst}'
#| fig-width: 7
#| fig-height: 4.5
#| echo: false

par(col.main = unname(uct_cols["navy"]), col.lab = unname(uct_cols["navy"]),
    col.axis = unname(uct_cols["navy"]), fg = unname(uct_cols["navy"]))

plot(NA, xlim = range(B_values), ylim = range(mc_results$SE), log = "x",
     xaxt = "n", xlab = "Bootstrap replicates (B)", ylab = "Bootstrap SE (kg)",
     main = "")
axis(1, at = B_values, labels = format(B_values, big.mark = " "))

mc_line_cols <- rep_len(uct_cols, n_runs)
for (run in seq_len(n_runs)) {
  d <- subset(mc_results, Run == run)
  lines(d$B, d$SE, type = "b", pch = 19, cex = 0.7,
        col = adjustcolor(mc_line_cols[run], alpha.f = 0.7))
}
#
#
#
#
#
#
#
#
#| results: asis
for (f in c("src/analysis.R")) {
  cat("**", f, "**\n\n```r\n", sep = "")
  cat(readLines(f), sep = "\n")
  cat("\n```\n\n")
}
#
#
#
#

# Bootstrapping: weight change after quitting smoking

Final assignment for **PPH7090S – Seminars in Epidemiology & Biostatistics** (University of Cape Town).

**Author:** Jolene Wium (WMXJOL001)

This project uses bootstrap methods to estimate the uncertainty in the difference in weight change between people who quit smoking and those who continued, using data from the NHANES I Epidemiologic Follow-up Study (NHEFS), 1971–1982.

## Research questions

1. What is the difference in **median** weight change between quitters and non-quitters, and how uncertain is it? There is no simple formula for the standard error of a difference in medians, so the bootstrap is used.
2. After adjusting for baseline characteristics, does the conventional linear-regression standard error for the effect of quitting agree with **case** and **residual** bootstrap standard errors?
3. How many bootstrap replicates (*B*) are needed for stable results? (Monte Carlo error.)

## Methods at a glance

| Analysis | Approach |
|---|---|
| Difference in medians | Stratified bootstrap (by quit status), *B* = 5 000; percentile and BCa 95% CIs |
| Adjusted effect of quitting | Linear regression adjusted for sex, race, age, education, smoking intensity, years smoked, exercise, physical activity and baseline weight; standard errors compared across OLS, case bootstrap and residual bootstrap |
| Monte Carlo error | Case bootstrap repeated 20 times at *B* = 200, 1 000, 5 000 and 10 000 |
| Missing outcome | Participants with and without 1982 weight compared on baseline characteristics |

## Key results

- Quitters' median weight gain was **1.8 kg** higher than non-quitters' (95% BCa CI 0.9 to 2.7 kg).
- After adjustment, quitters gained on average **3.4 kg** more (95% CI 2.5 to 4.3 kg, case bootstrap).
- The case bootstrap standard error was slightly larger than the OLS standard error (0.477 vs 0.441 kg), consistent with mild heteroscedasticity.
- *B* = 5 000 was sufficient for stable standard errors and confidence limits.

See [`report/WMXJOL001.pdf`](report/WMXJOL001.pdf) for the full report.

## Repository structure

```
.
├── data/
│   └── nhefs_2026.csv        # NHEFS analysis dataset (1 629 participants)
├── src/
│   ├── analysis.R            # Median bootstrap, regression bootstraps, Monte Carlo error
│   └── missing_data.R        # Comparison of participants with missing vs observed outcome
├── report/
│   ├── WMXJOL001.qmd         # Quarto source for the report
│   └── WMXJOL001.pdf         # Rendered report
└── README.md
```

## Reproducing the report

### Requirements

- [R](https://www.r-project.org/) (4.x)
- [Quarto](https://quarto.org/) 1.4 or later (renders to PDF with Typst, so no LaTeX installation is needed)
- R packages:

```r
install.packages(c("boot", "knitr", "kableExtra", "rmarkdown"))
```

### Render

From the repository root:

```bash
quarto render report/WMXJOL001.qmd
```

The report sets the knitr working directory to the repository root, so the scripts in `src/` and the data in `data/` are found automatically. Random seeds are set in `src/analysis.R`, so bootstrap results are reproducible.

The full Monte Carlo analysis (4 values of *B* × 20 runs) takes several minutes. Results are cached after the first render and are recalculated only when a script in `src/` changes.

## Data

`wt82_71` (weight change in kg between 1971 and 1982) is the outcome and `qsmk` (quit smoking between baseline and 1982) is the exposure. Weight change is missing for 63 participants (3.9%), who are excluded from the main analyses; see Appendix B of the report.

The NHEFS data are publicly available and are widely used in Hernán & Robins, *Causal Inference: What If* (Chapman & Hall/CRC, 2020).

## Colour palette

All figures use the project colour scheme: `#053259`, `#1668A3`, `#12999E`, `#C32A66`, `#7196B7`, `#FCB028`, `#FA615C`.

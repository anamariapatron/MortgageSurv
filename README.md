# MortgageSurv <img src="man/figures/logo.png" align="right" height="139" alt="" />

> Survival Analysis Models for Mortgage Default Prediction

[![R package](https://img.shields.io/badge/R-package-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

**MortgageSurv** adapts survival models originally developed in biostatistics
to the context of the UK mortgage market. The event of interest is defined as
mortgage **origination** (birth), while **default** is reinterpreted as death.
This conceptual shift provides a novel framework for studying loan dynamics in
the period preceding default.

**Author:** Ana M Patrón Piñerez

The package implements:

| Module | Functions |
|--------|-----------|
| **Simulation** | `simulate_mortgage_data()`, `simLMPH()`, `find_lambda()`, `generate_landmarks()` |
| **Real Data** | `from_real_data()`, `validate_mortgage_data()`, `example_real_data()` |
| **Static Modelling** | `prepare_aft_data()`, `fit_static_models()`, `plot_cumhaz_static()` |
| **Dynamic Modelling** | `prepare_landmark_data()`, `fit_landmark_model()` |
| **Prediction** | `predict_risk_aft()`, `predict_risk_ph()` |
| **Evaluation** | `evaluate_landmark_risk()`, `get_percentiles()` |
| **Visualisation** | `plot_cumulative_deaths()`, `plot_death_times()`, `plot_default_dates()`, `plot_cumulative_hazard()`, `plot_landmark_densities()`, `plot_aft_densities()`, `plot_risk_comparison()` |

---

## Installation

```r
# From a local .tar.gz build:
install.packages("MortgageSurv_0.1.0.tar.gz", repos = NULL, type = "source")

# Or install directly from the source folder:
devtools::install("path/to/MortgageSurv")

# Or from GitHub:
devtools::install_github("anamariapatron/MortgageSurv")
```

### Dependencies

```r
install.packages(c("lubridate", "truncnorm", "ggplot2", "dplyr",
                   "survival", "Landmarking", "eha", "flexsurv", "tidyr"))
```

---

## Workflow

### Option A · Simulated data

```r
library(MortgageSurv)

sim <- simulate_mortgage_data(n = 1000, n_iterations = 7, seed = 1234)
```

### Option B · Real data

```r
library(MortgageSurv)

# See the expected format:
template <- example_real_data(n = 10, m = 3)

# Load your own data:
my_data <- read.csv("my_mortgages.csv")

sim <- from_real_data(
  data           = my_data,
  id_col         = "loan_id",
  time_col       = "months_to_event",
  status_col     = "defaulted",
  covariate_cols = c("interest_rate", "inflation", "lti", "age"),
  date_col       = "obs_date",
  start_date_col = "origination_date",
  iteration_col  = "period"
)
```

---

### 1 · Visualise default patterns

```r
plot_cumulative_deaths(sim$status_matrix)
plot_death_times(sim$time_matrix, sim$status_matrix)
plot_default_dates(sim$obsdate_matrix, sim$status_matrix)
```

### 2 · Fit static models (Cox, Weibull PH, Weibull AFT)

```r
df_aft <- prepare_aft_data(
  time_matrix        = sim$time_matrix,
  status_matrix      = sim$status_matrix,
  x1_matrix          = sim$x1_matrix,
  x2_matrix          = sim$x2_matrix,
  x3_matrix          = sim$x3_matrix,
  x4_matrix          = sim$x4_matrix,
  obsdate_matrix     = sim$obsdate_matrix,
  credit_start_dates = sim$credit_start_dates
)

fits <- fit_static_models(df_aft)
summary(fits$weibull_aft)
plot_cumhaz_static(fits$weibull_ph, fits$weibull_aft)
```

### 3 · Fit the dynamic Landmarking model

```r
landmarks <- c(54, 56, 58, 60, 62, 64, 66)   # must match n_iterations

df_long <- prepare_landmark_data(
  time_matrix        = sim$time_matrix,
  status_matrix      = sim$status_matrix,
  x1_matrix          = sim$x1_matrix,
  x2_matrix          = sim$x2_matrix,
  x3_matrix          = sim$x3_matrix,
  x4_matrix          = sim$x4_matrix,
  obsdate_matrix     = sim$obsdate_matrix,
  credit_start_dates = sim$credit_start_dates,
  landmarks          = landmarks
)

lm_fit <- fit_landmark_model(df_long, x_L = landmarks, x_hor = landmarks + 12)
```

### 4 · Visualise predictions

```r
plot_cumulative_hazard(lm_fit, person_id = 3)
plot_landmark_densities(lm_fit, landmarks)
```

### 5 · Compare AFT vs. Landmarking

```r
risk_aft <- predict_risk_aft(fits$weibull_aft, newdata = df_aft,
                              landmarks = landmarks)

plot_risk_comparison(lm_fit, risk_aft, landmarks,
                     output_pdf = "comparison.pdf")
```

### 6 · Pseudo-evaluation metric

```r
ev <- evaluate_landmark_risk(
  matrices = list(sim$x1_matrix, sim$x2_matrix,
                  sim$x3_matrix, sim$x4_matrix),
  betas    = matrix(c(
     50,  0.9,  0.8, -0.2,
     51, 0.92, 0.81, -0.21,
     49, 0.88, 0.79, -0.19,
     52, 0.95, 0.82, -0.22,
     48, 0.89, 0.78, -0.18,
      5,  0.9,  0.8, -0.2,
     51, 0.91, 0.81, -0.21
  ), nrow = 7, byrow = TRUE),
  thetas   = matrix(c(
    0.5, 22.25, 0.5, 22.49, 0.5, 22.50,
    0.5, 22.73, 0.5, 22.96, 0.5, 23.19, 0.5, 23.41
  ), nrow = 7, byrow = TRUE)
)
ev$result
```

---

## Model description

### Covariates

| Variable | Description | Distribution |
|----------|-------------|--------------|
| `x1` | UK Bank Rate (interest rate) | Truncated Normal (0.001, 0.06) |
| `x2` | CPI inflation | Truncated Normal (0, 0.10) |
| `x3` | Loan-to-income ratio (LTI) | Truncated Normal (1, 6) |
| `x4` | Borrower age | Truncated Normal (20, 70) |

### Landmark PH Weibull simulation

Event times are generated via a piecewise Weibull model where the baseline
hazard is recalibrated at each landmark:

$$h(t | x, k) = h_0^{(k)}(t) \cdot \exp(\beta_k^\top x), \quad t \in [L_k, L_{k+1})$$

### Dynamic Landmarking

The Landmarking (LOCF) approach fits a separate Cox model at each landmark
$L_k$, using only individuals still at risk and their last-observation-carried-
forward covariate values.

---

## Citation

If you use this package in your research, please cite:

> Patrón Piñerez, A.M. (2025). *MortgageSurv: Survival Analysis Models for Mortgage
> Default Prediction*. R package version 0.1.0.

---

## License

MIT © 2025 Ana M Patrón Piñerez

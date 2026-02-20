# Cohort Survival Simulation with simulate_survivorship()

## Overview

Understanding population survivorship is fundamental to conservation
assessments, harvest management, and life history comparisons. The
[`simulate_survivorship()`](https://brian-j-moe.github.io/vitalBayes/reference/simulate_survivorship.md)
function performs Monte Carlo simulation of cohort survival using
age-specific mortality schedules from
[`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md).
By integrating over mortality schedules that reflect uncertainty in
growth and maturity parameters, the function provides fully propagated
uncertainty bounds on survival metrics.

This vignette covers:

- Mathematical framework for survival simulation
- Integration with vitalBayes mortality estimation
- Key survival metrics and interpretation
- Visualization and customization

## Mathematical Framework

### From Mortality to Survival

Cumulative survival to age \\t\\ given an age-specific mortality
schedule \\M(a)\\ is:

\\S(t) = \exp\left(-\int_0^t M(a) \\ da\right)\\

The function evaluates this integral numerically via trapezoidal
quadrature.

### Discrete Cohort Simulation

Demographic stochasticity is incorporated via binomial transitions:

\\N\_{t+1} \sim \text{Binomial}\left(N_t, \frac{S(t+1)}{S(t)}\right)\\

### Derived Metrics

**Mean age at death:** \\\bar{t}\_{death} = \frac{\sum_t t \cdot
d_t}{\sum_t d_t}\\

**Survival to maturity:** \\S(t\_{mat})\\

**Survival to maximum age:** \\S(t\_{max})\\

## Basic Usage

### With Stanfit Objects

``` r
library(vitalBayes)
library(data.table)

# Fit models in the vitalBayes workflow
# growth_fit <- fit_bayesian_growth(..., k_based = FALSE)

# Generate mortality schedules with full correlation preservation
mort <- get_stochastic_mortality(
  method     = "CW",
  growth_fit = growth_fit,
  sex        = 1,
  iter       = 2000,
  scaled     = TRUE,
  p          = 0.001,
  print_plot = FALSE
)

# Simulate survivorship
surv <- simulate_survivorship(
  mc_object = mort,
  n         = 50000,
  n_iter    = 2000,
  mode      = "random"
)

# Examine results
surv$Aggregate$Age_of_Death
surv$Aggregate$Survival_to_tmat
surv$Aggregate$Survival_to_tmax
```

### With Manual Parameters

``` r
# Generate mortality with specified correlation
mort <- get_stochastic_mortality(
  method       = "CW",
  Linf         = c(126, 10),
  L0           = c(35, 3),
  Lmat         = c(83, 5),
  tmat         = c(47, 4),
  maturity_cor = 0.5,
  growth_model = "vb",
  iter         = 2000,
  print_plot   = FALSE
)

surv <- simulate_survivorship(
  mc_object = mort,
  n         = 50000,
  n_iter    = 2000
)
```

## Output Structure

``` r
names(surv)
#> [1] "Aggregate" "Per_Set"   "Raw"       "Plot"

# Aggregate summaries
names(surv$Aggregate)
#> [1] "Survival" "Age_of_Death" "Survival_to_tmat" "Survival_to_tmax"

# Per-set summaries
names(surv$Per_Set)

# Raw trajectory data
head(surv$Raw)
#>      age survivors prop_surv  iter set_id
#>    <int>     <int>     <num> <int>  <int>
```

## Simulation Modes

### Random Sampling

``` r
surv <- simulate_survivorship(
  mc_object = mort,
  n         = 50000,
  n_iter    = 5000,
  mode      = "random"
)
```

### Per-Set Mode

``` r
surv_perset <- simulate_survivorship(
  mc_object = mort,
  n         = 50000,
  n_iter    = 100,
  mode      = "per_set"
)

# Examine variation across parameter sets
surv_perset$Per_Set$Age_of_Death
```

## Interpreting Metrics

### Mean Age at Death

``` r
surv$Aggregate$Age_of_Death
#>     mean    sd  lower  upper
#>    <num> <num>  <num>  <num>
#> 1:  42.1  8.92  26.4   58.3
```

### Survival to Maturity

``` r
surv$Aggregate$Survival_to_tmat
#>     mean      sd   lower   upper
#>    <num>   <num>   <num>   <num>
#> 1: 0.612  0.0891   0.448   0.782
```

### Survival to Maximum Age

``` r
surv$Aggregate$Survival_to_tmax
#>       mean      sd    lower    upper
#>      <num>   <num>    <num>    <num>
#> 1: 0.00092 0.00041 0.000341 0.00172
```

## Visualization

### Default Plot

``` r
surv <- simulate_survivorship(
  mc_object = mort,
  n         = 50000,
  n_iter    = 2000,
  palette   = "synthwave"
)
```

### Custom Aesthetics

``` r
surv <- simulate_survivorship(
  mc_object    = mort,
  n            = 50000,
  n_iter       = 2000,
  palette      = "viridis",
  bar_alpha    = 0.7,
  vline_color  = "#C0392B",
  title        = "Cohort Survivorship",
  base_size    = 12
)
```

### Custom X-Axis Breaks

``` r
surv <- simulate_survivorship(
  mc_object = mort,
  n_iter    = 2000,
  x_breaks  = function(tmax) seq(0, ceiling(tmax), by = 10)
)
```

## Sex-Specific Comparisons

``` r
# Females
mort_f <- get_stochastic_mortality(
  method = "CW", growth_fit = growth_fit, sex = 1,
  iter = 2000, print_plot = FALSE
)
surv_f <- simulate_survivorship(mort_f, n_iter = 2000, print_plot = FALSE)
surv_f$Aggregate$Survival[, sex := "Female"]

# Males
mort_m <- get_stochastic_mortality(
  method = "CW", growth_fit = growth_fit, sex = 2,
  iter = 2000, print_plot = FALSE
)
surv_m <- simulate_survivorship(mort_m, n_iter = 2000, print_plot = FALSE)
surv_m$Aggregate$Survival[, sex := "Male"]

# Combined plot
combined <- rbind(surv_f$Aggregate$Survival, surv_m$Aggregate$Survival)

library(ggplot2)
ggplot(combined, aes(x = age, y = mean_surv, color = sex, fill = sex)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  labs(x = "Age (years)", y = "Cumulative Survival",
       title = "Sex-Specific Survivorship Curves") +
  theme_bw()
```

## Practical Considerations

### Cohort Size

``` r
# Small cohort: More demographic stochasticity
surv_small <- simulate_survivorship(mc_object = mort, n = 1000, n_iter = 2000,
                                    print_plot = FALSE)

# Large cohort: Nearly deterministic
surv_large <- simulate_survivorship(mc_object = mort, n = 100000, n_iter = 2000,
                                    print_plot = FALSE)
```

For most applications, `n = 50000` provides a good balance.

### Iteration Count

``` r
# Quick exploratory
surv_quick <- simulate_survivorship(mc_object = mort, n_iter = 500,
                                    print_plot = FALSE)

# Final analysis
surv_final <- simulate_survivorship(mc_object = mort, n_iter = 5000,
                                    print_plot = FALSE)
```

## Sensitivity Analysis

``` r
per_set_data <- surv$Per_Set$Age_of_Death
per_set_merged <- merge(per_set_data, mort$Parameters, by = "set_id")

cor_matrix <- cor(per_set_merged[, .(mean_age, Linf, L0, Lmat, tmat, tmax)],
                  use = "complete.obs")
cor_matrix["mean_age", ]
```

## Troubleshooting

| Issue                    | Cause                 | Solution                      |
|--------------------------|-----------------------|-------------------------------|
| Slow computation         | Too many iterations   | Reduce n_iter for exploration |
| Survival = 0 or 1        | Age beyond schedule   | Extend age_seq in mortality   |
| Missing Survival_to_tmat | No tmat in Parameters | Ensure tmat provided          |
| NaN values               | Negative mortality    | Check mortality schedules     |

## See Also

- [`vignette("mortality_estimation")`](https://brian-j-moe.github.io/vitalBayes/articles/mortality_estimation.md) -
  Generating mortality schedules
- [`vignette("fit_bayesian_growth")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_growth.md) -
  Growth model fitting

## References

Caswell, H. (2001). *Matrix Population Models* (2nd ed.). Sinauer
Associates.

Cortes, E. (2002). Incorporating uncertainty into demographic modeling.
*Conservation Biology*, 16(4), 1048-1062.

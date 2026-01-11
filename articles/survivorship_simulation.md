# Cohort Survival Simulation with simulate_survivorship()

## Overview

Understanding population survivorship is fundamental to conservation
assessments, harvest management, and life history comparisons. While
natural mortality rates describe instantaneous hazards, survivorship
curves translate these into intuitive metrics: what proportion of a
cohort survives to reproductive maturity? What is the expected lifespan?
How does survival uncertainty propagate through to population-level
predictions?

The
[`simulate_survivorship()`](https://brian-j-moe.github.io/vitalBayes/reference/simulate_survivorship.md)
function performs Monte Carlo simulation of cohort survival using
age-specific mortality schedules from
[`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md).
By integrating over mortality schedules that themselves reflect
uncertainty in growth parameters, the function provides fully propagated
uncertainty bounds on survival metrics—bounds that honestly reflect our
knowledge (or lack thereof) about the population.

This vignette covers:

- The mathematical framework for survival simulation
- Integration with vitalBayes mortality estimation
- Key survival metrics and their interpretation
- Visualization options for publication
- Practical considerations for simulation design

## Mathematical Framework

### From Mortality to Survival

The relationship between instantaneous mortality rate \\M(t)\\ and
cumulative survival \\S(t)\\ is fundamental to demographic analysis. For
a cohort exposed to age-specific mortality, cumulative survival to age
\\t\\ is:

\\S(t) = \exp\left(-\int_0^t M(a) \\ da\right)\\

The function implements this relationship via trapezoidal numerical
integration over the mortality schedule. For each simulation iteration,
the integral is computed to yield a cumulative survival curve, and then
discrete cohort transitions are simulated using binomial survival.

### Discrete Cohort Simulation

While the cumulative survival curve is deterministic given a mortality
schedule, actual populations experience demographic stochasticity—random
variation in who survives each year. The function simulates this by
treating each age transition as a binomial process:

\\N\_{t+1} \sim \text{Binomial}\left(N_t, \frac{S(t+1)}{S(t)}\right)\\

where \\N_t\\ is the number of individuals alive at age \\t\\, and
\\S(t+1)/S(t)\\ is the conditional probability of surviving from age
\\t\\ to age \\t+1\\.

This stochastic formulation means that even with identical mortality
schedules, different simulation runs will produce slightly different
survival trajectories. This variation is appropriate when projecting
small populations where demographic stochasticity matters, but averaging
over many iterations recovers the expected (deterministic) survival
curve.

### Derived Metrics

Several biologically meaningful metrics emerge from the survival
simulations:

**Mean age at death** is computed by weighting ages by the number of
deaths at each age: \\\bar{t}\_{death} = \frac{\sum_t t \cdot
d_t}{\sum_t d_t}\\ where \\d_t\\ is the number of individuals dying at
age \\t\\.

**Survival to maturity** is the proportion of the birth cohort surviving
to age \\t\_{mat}\\, interpolated from the survival curve: \\S(t\_{mat})
= \text{interpolate}(t\_{mat}, \\t_i, S_i\\)\\

**Survival to maximum age** is analogously computed at \\t\_{max}\\,
representing the probability that any individual in the birth cohort
reaches the estimated maximum lifespan.

## Basic Usage

### Standard Workflow

The typical workflow begins with mortality estimation, then feeds those
schedules into survival simulation:

``` r
library(vitalBayes)
library(data.table)

# Step 1: Fit growth model (see fit_bayesian_growth vignette)
# growth_fit <- fit_bayesian_growth(...)
# maturity_fit <- fit_bayesian_maturity(...)

# Step 2: Generate mortality schedules
mort <- get_stochastic_mortality(
  method       = "CW",
  growth_fit   = growth_fit,
  maturity_fit = maturity_fit,
  sex          = 2,  # Males
  iter         = 2000,
  scaled       = TRUE,
  p            = 0.001,
  print_plot   = FALSE  # Suppress plot for now
)

# Step 3: Simulate survivorship
surv <- simulate_survivorship(
  mc_object = mort,
  n         = 50000,  # Starting cohort size
  n_iter    = 2000,   # Number of simulation iterations
  mode      = "random"
)

# View the plot (automatically printed by default)
# surv$Plot

# Examine key metrics
surv$Aggregate$Age_of_Death
#>     mean    sd  lower  upper
#>    <num> <num>  <num>  <num>
#> 1:  8.42  2.31   4.89  13.21
```

### Understanding the Output Structure

The function returns a list with four components:

``` r
# Aggregate summaries across all simulations
names(surv$Aggregate)
#> [1] "Survival"         "Age_of_Death"     "Survival_to_tmat" "Survival_to_tmax"

# Per-set summaries (one row per parameter set)
names(surv$Per_Set)
#> [1] "Survival"         "Age_of_Death"     "Survival_to_tmat" "Survival_to_tmax"

# Raw trajectory data
head(surv$Raw)
#>      age survivors prop_surv  iter set_id
#>    <int>     <int>     <num> <int>  <int>
#> 1:     0     50000    1.0000     1    847
#> 2:     1     34218    0.6844     1    847
#> 3:     2     27893    0.5579     1    847
#> ...

# The plot object
class(surv$Plot)
#> [1] "gg" "ggplot"
```

The `Aggregate` component provides the metrics most commonly needed for
reporting: population-level means and credible intervals. The `Per_Set`
component is useful for examining how survival varies across different
parameter combinations. The `Raw` data enables custom analyses and
alternative visualizations.

## Simulation Modes

### Random Sampling Mode

In `mode = "random"`, each simulation iteration draws a random parameter
set from those generated by
[`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md).
This mode is appropriate when you want overall population estimates that
integrate over parameter uncertainty:

``` r
surv_random <- simulate_survivorship(
  mc_object = mort,
  n         = 50000,
  n_iter    = 5000,  # More iterations for stable estimates
  mode      = "random"
)

# The number of iterations controls how well we estimate
# the expectation over parameter uncertainty
```

The number of iterations (`n_iter`) controls how well you approximate
the expected survival metrics under parameter uncertainty. For final
analyses, 2000-5000 iterations typically provide stable estimates. For
exploratory work, 500-1000 may suffice.

### Per-Set Mode

In `mode = "per_set"`, the function runs `n_iter` simulations for *each*
parameter set, enabling examination of how survival varies across the
parameter space:

``` r
# Run multiple simulations per parameter set
# Warning: This produces n_iter × n_sets total simulations
surv_perset <- simulate_survivorship(
  mc_object = mort,
  n         = 50000,
  n_iter    = 100,   # 100 simulations per parameter set
  mode      = "per_set"
)

# Examine variation across parameter sets
surv_perset$Per_Set$Age_of_Death
#>    set_id mean_age   sd_age  Linf      k     t0   tmat   tmax
#>     <int>    <num>    <num> <num>  <num>  <num>  <num>  <num>
#> 1:      1     7.89    0.234  95.2 0.0891  -1.42   8.44  32.45
#> 2:      2     9.12    0.198  102.1 0.0754 -1.56   9.21  38.12
#> ...
```

This mode is computationally expensive (total simulations = `n_iter` ×
number of parameter sets), but it’s valuable for:

- Understanding which parameters most strongly influence survival
- Identifying parameter regions with unusual behavior
- Building intuition about the model’s sensitivity

## Visualization

### Default Plot

The default plot shows a bar chart of mean survival by age with error
bars representing 95% credible intervals. A vertical line and shaded
band mark the mean age at death and its uncertainty. The caption reports
key survival metrics:

``` r
surv <- simulate_survivorship(
  mc_object = mort,
  n         = 50000,
  n_iter    = 2000,
  title     = "Cohort Survivorship: Male Spiny Dogfish",
  subtitle  = "Chen-Watanabe mortality with Gompertz senescence"
)
```

### Color Palette Options

Like other vitalBayes functions,
[`simulate_survivorship()`](https://brian-j-moe.github.io/vitalBayes/reference/simulate_survivorship.md)
supports multiple color palettes:

``` r
# Synthwave (default)
surv_synth <- simulate_survivorship(mc_object = mort, palette = "synthwave")

# Okabe-Ito (colorblind-friendly)
surv_okabe <- simulate_survivorship(mc_object = mort, palette = "okabe")

# Viridis family
surv_plasma <- simulate_survivorship(mc_object = mort, palette = "plasma")
```

### Custom Colors

For complete control, specify individual colors:

``` r
surv_custom <- simulate_survivorship(
  mc_object   = mort,
  bar_fill    = "#1B9E77",   # Teal bars
  bar_color   = "#117755",   # Darker teal outline
  vline_color = "#D95F02",   # Orange mean age line
  rect_fill   = "#D95F02",   # Orange uncertainty band
  rect_alpha  = 0.3
)
```

### Customizing the ggplot Object

The returned plot is a standard ggplot2 object, so you can modify it:

``` r
library(ggplot2)

surv$Plot +
  # Add reference lines
  geom_vline(xintercept = c(8, 12), linetype = "dashed", color = "gray50") +
  annotate("text", x = 8.2, y = 0.8, label = "Maturity", hjust = 0, size = 3) +
  annotate("text", x = 12.2, y = 0.8, label = "Peak reproduction", hjust = 0, size = 3) +
  # Customize labels
  labs(
    title = "Predicted Survivorship Curve",
    caption = NULL  # Remove default caption if desired
  ) +
  # Modify theme
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = "gray90")
  )
```

### Alternative Visualization: Survival Curves

For some applications, a continuous survival curve is more appropriate
than the discrete bar chart:

``` r
library(ggplot2)

# Use the aggregate survival data
ggplot(surv$Aggregate$Survival, aes(x = age)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              fill = vital_palette(1), alpha = 0.3) +
  geom_line(aes(y = mean_surv), color = vital_palette(1, "synthwave")[3],
            linewidth = 1.2) +
  # Add key ages
  geom_point(data = data.frame(
    age = c(0, 8, 12, 30),
    surv = surv$Aggregate$Survival[age %in% c(0, 8, 12, 30), mean_surv]
  ), aes(y = surv), size = 3, color = vital_palette(1, "synthwave")[1]) +
  labs(
    x = "Age (years)",
    y = "Cumulative Survival",
    title = "Type I Survivorship Curve"
  ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  theme_bw(base_size = 12)
```

## Interpreting Survival Metrics

### Mean Age at Death

The mean age at death represents the average lifespan of individuals in
the cohort. For populations with high juvenile mortality (typical of
elasmobranchs), this value will be substantially lower than the maximum
age:

``` r
surv$Aggregate$Age_of_Death
#>     mean    sd  lower  upper
#>    <num> <num>  <num>  <num>
#> 1:  8.42  2.31   4.89  13.21

# Interpretation: On average, individuals live about 8.4 years,
# with 95% of simulation runs producing mean lifespans between
# 4.9 and 13.2 years
```

The uncertainty here reflects both parameter uncertainty (from the
growth model posterior) and the scaling relationship for mortality. Wide
credible intervals suggest that survival metrics are sensitive to
uncertainty in growth parameters.

### Survival to Maturity

This metric directly addresses a key conservation question: what
proportion of individuals survive long enough to reproduce?

``` r
surv$Aggregate$Survival_to_tmat
#>     mean      sd   lower   upper
#>    <num>   <num>   <num>   <num>
#> 1: 0.312  0.0891   0.168   0.512

# Interpretation: About 31% of individuals survive to maturity,
# with substantial uncertainty (95% CI: 17% - 51%)
```

For populations where maturity occurs late relative to early mortality,
this proportion can be surprisingly low. Values below ~20% may indicate
populations vulnerable to recruitment failure under increased mortality
pressure.

### Survival to Maximum Age

This metric typically yields very small probabilities, as expected:

``` r
surv$Aggregate$Survival_to_tmax
#>       mean      sd    lower    upper
#>      <num>   <num>    <num>    <num>
#> 1: 0.00089 0.00041 0.000312 0.00178

# Interpretation: About 0.09% (roughly 1 in 1,100) individuals 
# survive to estimated maximum age
```

This metric is primarily useful as a sanity check. If `p` was set to
0.001 in
[`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md),
you’d expect survival to tmax to be near 0.001—substantial deviation
suggests the mortality schedule shape is inconsistent with the scaling
assumption.

## Practical Considerations

### Cohort Size Selection

The `n` parameter controls the starting cohort size for each simulation.
Larger cohorts reduce demographic stochasticity:

``` r
# Small cohort: More demographic stochasticity
surv_small <- simulate_survivorship(mc_object = mort, n = 1000, n_iter = 2000)

# Large cohort: Nearly deterministic given the mortality schedule
surv_large <- simulate_survivorship(mc_object = mort, n = 100000, n_iter = 2000)

# Compare SD of mean survival at age 10
surv_small$Aggregate$Survival[age == 10, se_surv]  # Larger SE
surv_large$Aggregate$Survival[age == 10, se_surv]  # Smaller SE
```

For most applications, `n = 50000` provides a good balance. Use larger
values when you want to minimize demographic noise and isolate parameter
uncertainty. Use smaller values when demographic stochasticity is itself
of interest (e.g., for small population projections).

### Iteration Count

The number of iterations affects how well you characterize the
uncertainty distribution:

``` r
# Quick exploratory analysis
surv_quick <- simulate_survivorship(mc_object = mort, n_iter = 500)

# Final analysis for publication
surv_final <- simulate_survivorship(mc_object = mort, n_iter = 5000)
```

As a rule of thumb, credible interval bounds stabilize faster than
means, and survival probabilities at extreme ages require more
iterations to estimate precisely.

### Computational Efficiency

Survival simulation is computationally intensive because it involves
numerical integration for each parameter set. Some strategies for
managing computation time:

``` r
# Reduce parameter sets in mortality estimation
mort_fast <- get_stochastic_mortality(
  method     = "CW",
  growth_fit = growth_fit,
  sex        = 2,
  iter       = 500,  # Fewer parameter draws
  print_plot = FALSE
)

# Fewer survival iterations
surv_fast <- simulate_survivorship(
  mc_object = mort_fast,
  n         = 20000,  # Smaller cohort
  n_iter    = 500     # Fewer iterations
)
```

For exploratory work, 500-1000 parameter sets and survival iterations
are usually sufficient to identify qualitative patterns.

## Sex-Specific Comparisons

To compare survival between sexes, run the workflow separately for each:

``` r
# Females
mort_f <- get_stochastic_mortality(
  method = "CW", growth_fit = growth_fit, maturity_fit = maturity_fit,
  sex = 1, iter = 2000, print_plot = FALSE
)
surv_f <- simulate_survivorship(mort_f, n = 50000, n_iter = 2000, print_plot = FALSE)
surv_f$Aggregate$Survival[, sex := "Female"]

# Males
mort_m <- get_stochastic_mortality(
  method = "CW", growth_fit = growth_fit, maturity_fit = maturity_fit,
  sex = 2, iter = 2000, print_plot = FALSE
)
surv_m <- simulate_survivorship(mort_m, n = 50000, n_iter = 2000, print_plot = FALSE)
surv_m$Aggregate$Survival[, sex := "Male"]

# Combined plot
combined <- rbind(surv_f$Aggregate$Survival, surv_m$Aggregate$Survival)

library(ggplot2)
ggplot(combined, aes(x = age, y = mean_surv, color = sex, fill = sex)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = vital_palette(2)) +
  scale_fill_manual(values = vital_palette(2)) +
  labs(x = "Age (years)", y = "Cumulative Survival",
       title = "Sex-Specific Survivorship Curves") +
  theme_bw() +
  theme(legend.position = "top")
```

## Sensitivity Analysis

The `Per_Set` output enables informal sensitivity analysis:

``` r
# Examine correlation between parameters and survival metrics
per_set_data <- surv$Per_Set$Age_of_Death

# Which parameters most strongly predict mean age at death?
cor_matrix <- cor(per_set_data[, .(mean_age, Linf, k, t0, tmat, tmax)])
cor_matrix["mean_age", ]
#>   mean_age      Linf         k        t0      tmat      tmax 
#>     1.0000    0.3241   -0.5892    0.2134    0.4521    0.7823

# Interpretation: tmax is the strongest predictor (positively correlated),
# followed by k (negatively correlated) and tmat (positively correlated)
```

For formal sensitivity analysis, consider varying individual parameters
while holding others fixed, or using variance-based methods like Sobol
indices.

## Troubleshooting

| Issue                                    | Possible Cause                      | Solution                                                                                                                                 |
|------------------------------------------|-------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------|
| Very slow computation                    | Too many iterations or large cohort | Reduce `n_iter` or `n` for exploratory work                                                                                              |
| Survival probabilities of exactly 0 or 1 | Age beyond mortality schedule range | Extend `age_seq` in [`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md)       |
| Missing Survival_to_tmat                 | No tmat in Parameters               | Provide `maturity_fit` to [`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md) |
| NaN values in survival                   | Negative mortality values           | Check mortality schedules for numerical issues                                                                                           |
| Unrealistic survival rates               | Mortality scaling mismatch          | Verify `p` or `M_target` in mortality estimation                                                                                         |

## See Also

- [`vignette("mortality_estimation")`](https://brian-j-moe.github.io/vitalBayes/articles/mortality_estimation.md)
  — Generating mortality schedules for survival simulation
- [`vignette("fit_bayesian_growth")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_growth.md)
  — Growth model fitting that feeds the workflow
- [`vignette("fit_bayesian_maturity")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_maturity.md)
  — Maturity estimation for tmat and Lmat
- [Statistical
  Background](https://brian-j-moe.github.io/vitalBayes/articles/Understanding_vitalBayes.html#survival)
  — Mathematical foundations

## References

Caswell, H. (2001). *Matrix Population Models: Construction, Analysis,
and Interpretation* (2nd ed.). Sinauer Associates.

Cortes, E. (2002). Incorporating uncertainty into demographic modeling:
Application to shark populations and their conservation. *Conservation
Biology*, 16(4), 1048-1062.

Simpfendorfer, C. A. (1999). Mortality estimates and demographic
analysis for the Australian sharpnose shark, *Rhizoprionodon taylori*,
from northern Australia. *Fishery Bulletin*, 97(4), 978-986.

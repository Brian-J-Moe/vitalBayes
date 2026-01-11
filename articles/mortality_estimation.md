# Age-Specific Natural Mortality with get_stochastic_mortality()

## Overview

Natural mortality (\\M\\) is among the most difficult parameters to
estimate in fisheries science, yet it profoundly influences population
dynamics, sustainable harvest levels, and conservation assessments. For
elasmobranchs, obtaining direct mortality estimates is particularly
challenging: long lifespans make mark-recapture studies impractical,
populations are often too sparse for reliable catch-curve analysis, and
the absence of otoliths necessitates alternative aging structures with
their own uncertainties.

The
[`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md)
function addresses this challenge by implementing three established
life-history-based mortality models within a Monte Carlo framework that
properly propagates uncertainty in von Bertalanffy growth parameters.
When integrated with vitalBayes growth model outputs, the function
preserves the joint posterior distribution of life history parameters,
yielding mortality schedules with realistic, honest uncertainty bounds.

This vignette covers:

- The biological and mathematical foundations of each mortality model
- Integration with vitalBayes growth model outputs
- Scaling mortality schedules to match empirical relationships
- Visualization with the vitalBayes color system
- Best practices for model selection and interpretation

## Mathematical Foundations

### Chen-Watanabe Model (1989)

The Chen-Watanabe model derives age-specific mortality from von
Bertalanffy growth parameters under the biological reasoning that
metabolic demands—and thus vulnerability to mortality—scale with growth
rate. The core equation expresses instantaneous mortality as:

\\M(t) = \frac{k}{1 - e^{-k(t - t_0)}}\\

where \\k\\ is the von Bertalanffy growth coefficient and \\t_0\\ is the
theoretical age at length zero. This formulation produces high juvenile
mortality that declines rapidly during the growth phase, stabilizing as
individuals approach asymptotic size.

For long-lived species like elasmobranchs, the single-phase model may
underestimate senescence mortality. The optional two-phase extension
adds a late-life component that increases mortality after a transition
age \\t_m\\. vitalBayes implements both Gompertz and logistic senescence
functions, with smooth blending to avoid discontinuities:

**Gompertz senescence:** \\M\_{late}(t) = A \cdot e^{B \cdot t}\\

**Logistic senescence:** \\M\_{late}(t) = \frac{K}{1 + e^{-r(t -
t\_{max})}}\\

The transition is smoothed using a sigmoid weight function: \\w(t) =
\frac{1}{1 + e^{-2(t - t_m)/\sigma_w}}\\

where \\\sigma_w\\ controls the width of the transition zone.

### Peterson-Wroblewski Model (1984)

The Peterson-Wroblewski model takes a purely allometric approach,
expressing mortality as a power function of body weight:

\\M(W) = 1.92 \cdot W^{-0.25}\\

where \\W\\ is body weight in grams. This model captures the
well-documented pattern that smaller individuals experience higher
mortality across taxa, though it does not explicitly incorporate
age-dependent processes beyond the weight-age relationship implied by
growth.

To apply this model, you must provide a length-weight function that
converts von Bertalanffy predicted lengths to body mass. The exponent
(-0.25) represents the metabolic scaling relationship and is remarkably
consistent across fish species.

### Lorenzen Model (1996, 2022)

Lorenzen’s framework offers two formulations that capture size-dependent
mortality:

**Weight-based formulation:** \\M(W) = \alpha \cdot W^{\beta}\\

where \\\alpha \sim \mathcal{N}(3.69, 0.502)\\ and \\\beta \sim
\mathcal{N}(-0.305, 0.029)\\. The parameter distributions incorporate
uncertainty from the meta-analysis underlying the model.

**Growth-based formulation (Lorenzen 2022):** \\\ln M = 0.28 - 1.30
\ln\left(\frac{L}{L\_\infty}\right) + 1.08 \ln(k)\\

This newer formulation directly incorporates von Bertalanffy parameters
without requiring a length-weight relationship, making it particularly
useful when allometric parameters are uncertain or unavailable. The
coefficients include uncertainty: \\0.28 \pm 0.105\\, \\-1.30 \pm
0.059\\, and \\1.08 \pm 0.082\\.

## Integration with vitalBayes

### Why Joint Posterior Sampling Matters

When growth parameters are estimated via Bayesian methods, the posterior
distribution typically exhibits substantial correlations. In von
Bertalanffy models, \\L\_\infty\\ and \\k\\ are almost always negatively
correlated—high asymptotic size tends to associate with slower growth,
and vice versa. The correlation between \\k\\ and \\t_0\\ is often
positive.

The traditional approach of specifying parameters as independent
`c(mean, sd)` vectors ignores these correlations, occasionally
generating biologically implausible combinations (e.g., high
\\L\_\infty\\ with high \\k\\) that never appeared in the original
posterior. Since Chen-Watanabe mortality is directly proportional to
\\k\\, and all models depend on the growth trajectory, these implausible
combinations propagate into the mortality schedules, potentially
inflating uncertainty bounds.

[`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md)
addresses this by accepting vitalBayes fit objects directly, drawing
correlated parameter samples from the joint posterior. This yields
mortality schedules with *narrower but more honest* uncertainty
intervals that reflect genuine biological uncertainty rather than
sampling artifacts.

### Extracting Parameters from vitalBayes Fits

The workhorse function for parameter extraction is
[`extract_lh_params()`](https://brian-j-moe.github.io/vitalBayes/reference/extract_lh_params.md),
which pulls posterior draws while preserving correlations:

``` r
library(vitalBayes)
library(data.table)

# After fitting growth and maturity models
# growth_fit <- fit_bayesian_growth(...)
# maturity_fit <- fit_bayesian_maturity(...)

# Extract full posterior draws for males
params_male <- extract_lh_params(
  growth_fit   = growth_fit,
  maturity_fit = maturity_fit,
  sex          = 2,  # 1 = female, 2 = male
  format       = "draws",
  n_draws      = 2000
)

# View structure
head(params_male)
#>     Linf      k      L0       t0    tmat .draw sex
#>    <num>  <num>   <num>    <num>   <num> <int> <char>
#> 1:  98.2 0.0823  33.1    -1.42     8.92     1     M
#> 2:  95.8 0.0891  34.0    -1.35     8.44     2     M
#> ...

# Check correlations are preserved
cor(params_male[, .(Linf, k, t0)])
#>        Linf      k     t0
#> Linf  1.000 -0.712  0.234
#> k    -0.712  1.000 -0.456
#> t0    0.234 -0.456  1.000
```

The correlation structure you see here should match what you’d compute
from the original `growth_fit` posterior—that’s the key benefit of joint
sampling.

## Basic Usage

### Using vitalBayes Fits (Recommended)

The most statistically rigorous approach passes fit objects directly:

``` r
# Chen-Watanabe with two-phase senescence
mort_cw <- get_stochastic_mortality(
  method       = "CW",
  growth_fit   = growth_fit,
  maturity_fit = maturity_fit,
  sex          = 2,  # Males
  iter         = 2000,
  scaled       = TRUE,
  M_target     = NULL,  # Derive from survival probability
  p            = 0.001, # 0.1% survive to tmax
  two_phase    = TRUE,
  late_model   = "gompertz"
)

# Examine output
mort_cw$Summary
#>    age_round M_median   M_mean  M_lower  M_upper
#>        <num>    <num>    <num>    <num>    <num>
#> 1:      0.00   0.4521   0.4892   0.2891   0.7124
#> 2:      0.02   0.3982   0.4245   0.2654   0.6234
#> ...

# View the plot
mort_cw$Plot
```

### Manual Parameter Specification

When vitalBayes fits aren’t available (e.g., working from published
literature values), you can specify parameters as `c(mean, sd)` vectors:

``` r
# Parameters from published growth study
# Note: This approach ignores correlations

mort_manual <- get_stochastic_mortality(
  method = "CW",
  Linf   = c(120, 8),     # cm, mean and SD
  k      = c(0.08, 0.015),
  t0     = c(-1.5, 0.3),
  tmat   = c(12, 1.5),    # Age at maturity for CW transition
  iter   = 2000,
  scaled = TRUE,
  p      = 0.001
)
```

## Scaling Mortality Schedules

### The Scaling Problem

Life-history-based mortality models produce *relative* schedules that
capture age-dependent patterns but may not match the overall mortality
level expected for a population. Scaling adjusts the mean mortality to
match an external target while preserving the age-specific shape.

The scaling transformation is: \\M\_{scaled}(t) =
\frac{M\_{raw}(t)}{\bar{M}\_{raw}} \times M\_{target}\\

### Scaling Options

**1. Survival probability to maximum age:**

When no external mortality information is available, you can derive
\\M\_{target}\\ from an assumed probability of surviving to maximum age:

\\M\_{target} = \frac{-\ln(p)}{t\_{max}}\\

``` r
mort_scaled <- get_stochastic_mortality(
  method   = "CW",
  growth_fit = growth_fit,
  sex      = 2,
  scaled   = TRUE,
  M_target = NULL,  # Triggers survival-based scaling
  p        = 0.001  # 0.1% survive to tmax
)
```

A common choice is \\p = 0.001\\ (0.1%), which implies that the
instantaneous mortality rate, averaged over the lifespan, equals
approximately \\-\ln(0.001)/t\_{max}\\. This is a conservative
assumption for many elasmobranchs.

**2. Hoenig-type empirical relationships:**

Hoenig (1983) and Then et al. (2015) provide regression relationships
between maximum age and natural mortality. You can pass these as
functions:

``` r
# Hoenig (1983) for fish
hoenig_fish <- function(tmax) exp(1.46 - 1.01 * log(tmax))

# Then et al. (2015) updated relationship
then_2015 <- function(tmax) 4.899 * tmax^(-0.916)

# With uncertainty in the relationship itself
then_stochastic <- function(tmax) {
  rnorm(1, 4.899, 0.327) * tmax^rnorm(1, -0.916, 0.043)
}

mort_hoenig <- get_stochastic_mortality(
  method   = "CW",
  growth_fit = growth_fit,
  sex      = 2,
  scaled   = TRUE,
  M_target = then_stochastic  # Function evaluated per iteration
)
```

**3. Fixed target value:**

If you have an independent mortality estimate (e.g., from catch curves
or tagging), you can scale directly to that value:

``` r
mort_fixed <- get_stochastic_mortality(
  method   = "L",
  growth_fit = growth_fit,
  sex      = 1,  # Females
  scaled   = TRUE,
  M_target = 0.15  # Fixed value from catch-curve analysis
)
```

## Model Selection Guidelines

### When to Use Each Model

**Chen-Watanabe** is typically preferred for elasmobranchs because:

- It explicitly incorporates age-dependent processes through the growth
  function
- The two-phase extension captures senescence patterns observed in
  long-lived species
- The transition age can be informed by maturity data
- It requires only von Bertalanffy parameters (no length-weight
  relationship needed)

**Peterson-Wroblewski** is useful when:

- You want pure size-based mortality without age assumptions
- Reliable length-weight parameters are available
- You’re comparing across species with different growth patterns

**Lorenzen** is advantageous when:

- You want to incorporate uncertainty in the allometric relationship
  itself
- The growth-based formulation (2022) is attractive when length-weight
  data are unavailable
- Meta-analytic uncertainty is important for your application

### Practical Decision Framework

| Scenario                                   | Recommended Model       | Rationale                                       |
|--------------------------------------------|-------------------------|-------------------------------------------------|
| Long-lived elasmobranch with maturity data | CW two-phase            | Captures both juvenile and senescence mortality |
| Unknown age structure, good size data      | Lorenzen (weight-based) | Size-dependent without age assumptions          |
| Limited data, literature growth params     | CW single-phase         | Minimal parameter requirements                  |
| Comparing across diverse taxa              | Peterson-Wroblewski     | Standardized allometric approach                |
| No length-weight relationship              | Lorenzen (growth-based) | Uses only VB parameters                         |

## Visualization

### Color Palette Options

vitalBayes provides several color palettes through
[`vital_palette()`](https://brian-j-moe.github.io/vitalBayes/reference/vital_palette.md):

``` r
# Synthwave (default) - vitalBayes hex sticker colors
mort_synth <- get_stochastic_mortality(
  method = "CW", growth_fit = growth_fit, sex = 2,
  palette = "synthwave"
)

# Colorblind-friendly: Okabe-Ito palette
mort_okabe <- get_stochastic_mortality(
  method = "CW", growth_fit = growth_fit, sex = 2,
  palette = "okabe"
)

# Viridis family
mort_viridis <- get_stochastic_mortality(
  method = "CW", growth_fit = growth_fit, sex = 2,
  palette = "viridis"
)

# Custom colors
mort_custom <- get_stochastic_mortality(
  method = "CW", growth_fit = growth_fit, sex = 2,
  fill_color = "#2E86AB",
  line_color = "#A23B72"
)
```

### Customizing Plots

The returned plot is a ggplot2 object, so you can modify it further:

``` r
library(ggplot2)

# Add annotations
mort_cw$Plot +
  geom_vline(xintercept = 12, linetype = "dashed", color = "gray40") +
  annotate("text", x = 12.5, y = 0.4, label = "Age at maturity",
           hjust = 0, size = 3) +
  labs(title = "Spiny Dogfish Natural Mortality",
       subtitle = "Chen-Watanabe model with Gompertz senescence")
```

### Comparing Models

``` r
# Fit all three models
mort_cw <- get_stochastic_mortality(
  method = "CW", growth_fit = growth_fit, sex = 2,
  print_plot = FALSE
)

mort_pw <- get_stochastic_mortality(
  method = "PW", growth_fit = growth_fit, sex = 2,
  lw_fun = function(L) 0.0001 * L^3.1,  # Example L-W function
  print_plot = FALSE
)

mort_lor <- get_stochastic_mortality(
  method = "L", growth_fit = growth_fit, sex = 2,
  weight_based = FALSE,
  print_plot = FALSE
)

# Combine for comparison plot
library(ggplot2)

combined <- rbind(
  mort_cw$Summary[, model := "Chen-Watanabe"],
  mort_pw$Summary[, model := "Peterson-Wroblewski"],
  mort_lor$Summary[, model := "Lorenzen"]
)

ggplot(combined, aes(x = age_round)) +
  geom_ribbon(aes(ymin = M_lower, ymax = M_upper, fill = model),
              alpha = 0.3) +
  geom_line(aes(y = M_median, color = model), linewidth = 1) +
  scale_color_manual(values = vital_palette(3)) +
  scale_fill_manual(values = vital_palette(3)) +
  labs(x = "Age (years)", y = "Instantaneous Mortality (M)",
       title = "Comparison of Mortality Models") +
  theme_bw() +
  theme(legend.position = "top")
```

## Modular Mortality Functions

For applications where you need mortality at specific ages without the
full Monte Carlo machinery, vitalBayes exports the individual model
functions:

``` r
# Chen-Watanabe at specific ages
ages <- c(0, 1, 5, 10, 20, 30)
M_cw <- M_chen_watanabe(
  age   = ages,
  Linf  = 100,
  k     = 0.08,
  t0    = -1.5,
  tmax  = 35,
  tmat  = 12,
  two_phase   = TRUE,
  late_model  = "gompertz"
)

# Peterson-Wroblewski
lw_func <- function(L) 0.0001 * L^3.1
M_pw <- M_peterson_wroblewski(
  age    = ages,
  Linf   = 100,
  k      = 0.08,
  t0     = -1.5,
  lw_fun = lw_func
)

# Lorenzen (growth-based)
M_lor <- M_lorenzen(
  age          = ages,
  Linf         = 100,
  k            = 0.08,
  t0           = -1.5,
  weight_based = FALSE,
  sample_params = FALSE  # Use mean values, not sampled
)

# View results
data.table(age = ages, CW = M_cw, PW = M_pw, Lorenzen = M_lor)
```

The
[`scale_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/scale_mortality.md)
helper applies the scaling transformation:

``` r
M_raw <- M_chen_watanabe(ages, Linf = 100, k = 0.08, t0 = -1.5, tmax = 35)

# Scale to Hoenig target
M_scaled <- scale_mortality(M_raw, M_target = then_2015, tmax = 35)
```

## Downstream Integration: Survival Simulation

The mortality schedules from
[`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md)
feed directly into
[`simulate_survivorship()`](https://brian-j-moe.github.io/vitalBayes/reference/simulate_survivorship.md)
for population projections:

``` r
# Generate mortality schedules
mort <- get_stochastic_mortality(
  method     = "CW",
  growth_fit = growth_fit,
  maturity_fit = maturity_fit,
  sex        = 2,
  iter       = 2000,
  scaled     = TRUE,
  p          = 0.001,
  print_plot = FALSE
)

# Simulate cohort survival
surv <- simulate_survivorship(
  mc_object = mort,
  n         = 50000,   # Cohort size
  n_iter    = 2000,    # Simulation iterations
  mode      = "random" # Sample parameter sets randomly
)

# View results
surv$Aggregate$Age_of_Death
#>     mean    sd  lower  upper
#>    <num> <num>  <num>  <num>
#> 1:  8.42  2.31   4.89  13.21

surv$Aggregate$Survival_to_tmax
#>       mean      sd    lower    upper
#>      <num>   <num>    <num>    <num>
#> 1: 0.00089 0.00041 0.000312 0.00178
```

See
[`vignette("survivorship_simulation")`](https://brian-j-moe.github.io/vitalBayes/articles/survivorship_simulation.md)
for comprehensive coverage of survival analysis.

## Troubleshooting

| Issue                       | Possible Cause                               | Solution                                                                   |
|-----------------------------|----------------------------------------------|----------------------------------------------------------------------------|
| Negative mortality values   | CW two-phase Taylor expansion breakdown      | Function automatically caps; consider single-phase or different late_model |
| Very wide uncertainty bands | Low correlation in posterior                 | Check growth model convergence; correlations should be substantial         |
| Unrealistic tmax estimates  | Linf_factor too extreme                      | Reduce from 0.999 to 0.99 or 0.95                                          |
| Missing tmat error          | CW needs age-at-maturity                     | Provide maturity_fit or manual tmat parameter                              |
| lw_fun error                | PW/Lorenzen weight-based needs length-weight | Provide function: `lw_fun = function(L) a * L^b`                           |

## See Also

- [`vignette("survivorship_simulation")`](https://brian-j-moe.github.io/vitalBayes/articles/survivorship_simulation.md)
  — Cohort survival analysis using mortality schedules
- [`vignette("fit_bayesian_growth")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_growth.md)
  — Growth model fitting that feeds mortality estimation
- [`vignette("partial_pooling")`](https://brian-j-moe.github.io/vitalBayes/articles/partial_pooling.md)
  — Hierarchical models for sex-specific parameters
- [Statistical
  Background](https://brian-j-moe.github.io/vitalBayes/articles/Understanding_vitalBayes.html#mortality)
  — Full mathematical derivations

## References

Chen, S., & Watanabe, S. (1989). Age dependence of natural mortality
coefficient in fish population dynamics. *Nippon Suisan Gakkaishi*,
55(2), 205-208.

Hoenig, J. M. (1983). Empirical use of longevity data to estimate
mortality rates. *Fishery Bulletin*, 82(1), 898-903.

Lorenzen, K. (1996). The relationship between body weight and natural
mortality in juvenile and adult fish: a comparison of natural ecosystems
and aquaculture. *Journal of Fish Biology*, 49(4), 627-642.

Lorenzen, K. (2022). Size- and age-dependent natural mortality in fish
populations: Biology, models, implications, and a generalized
length-inverse model. *Fisheries Research*, 255, 106454.

Peterson, I., & Wroblewski, J. S. (1984). Mortality rate of fishes in
the pelagic ecosystem. *Canadian Journal of Fisheries and Aquatic
Sciences*, 41(7), 1117-1120.

Then, A. Y., Hoenig, J. M., Hall, N. G., & Hewitt, D. A. (2015).
Evaluating the predictive performance of empirical estimators of natural
mortality rate using information on over 200 fish species. *ICES Journal
of Marine Science*, 72(1), 82-92.

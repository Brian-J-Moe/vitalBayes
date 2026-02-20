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
properly propagates uncertainty. The function accepts either vitalBayes
stanfit objects (preserving posterior correlations) or manual
`c(mean, sd)` specifications, with explicit handling of covariance
between maturity parameters.

This vignette covers:

- The unified mortality framework: \\M(t) = M\_\infty / G(t)\\
- Parameter input options: stanfit objects vs. manual specification
- Covariance handling for maturity parameters
- The three supported mortality models
- Integration with survivorship simulation

## The Unified Mortality Framework

### Core Formulation

The vitalBayes mortality framework expresses age-specific Chen-Watanabe
mortality as:

\\M(t) = \frac{M\_\infty}{G(t)}\\

where \\M\_\infty\\ is the asymptotic mortality rate and \\G(t) =
L(t)/L\_\infty\\ is the relative size at age \\t\\.

### The Unified \\M\_\infty\\ Anchor

The asymptotic mortality rate is computed using the von Bertalanffy
formula regardless of which growth model is used:

\\M\_\infty = \frac{1}{t\_{mat}} \ln\left(\frac{L\_\infty -
L_0}{L\_\infty - L\_{mat}}\right)\\

This provides a consistent mortality anchor across growth models.

### Model-Specific Native \\k\\ Derivations

Each growth model requires its own \\k\\ derivation:

**Von Bertalanffy:** \\k\_{VB} = \frac{1}{t\_{mat}}
\ln\left(\frac{L\_\infty - L_0}{L\_\infty - L\_{mat}}\right)\\

**Gompertz:** \\k\_{Gomp} = -\frac{1}{t\_{mat}}
\ln\left(\frac{\ln(L\_\infty / L\_{mat})}{\ln(L\_\infty / L_0)}\right)\\

**Logistic:** \\k\_{Log} = -\frac{1}{t\_{mat}}
\ln\left(\frac{L\_\infty/L\_{mat} - 1}{L\_\infty/L_0 - 1}\right)\\

## Parameter Input Options

### Option 1: vitalBayes Stanfit Objects (Recommended)

When vitalBayes fit objects are provided, posterior correlations between
parameters are preserved:

``` r
library(vitalBayes)
library(data.table)

# After fitting models in the vitalBayes workflow:
# birth_fit <- fit_bayesian_birth(...)
# L50_fit <- fit_bayesian_maturity(maturity = "mat", lt = "fl", ...)
# t50_fit <- fit_bayesian_maturity(maturity = "mat", age = "age", ...)
# growth_fit <- fit_bayesian_growth(..., k_based = FALSE)  # Maturity-based

# Full correlation preservation (all from growth_fit)
mort <- get_stochastic_mortality(
  method     = "CW",
  growth_fit = growth_fit,   # Provides Linf, L0, Lmat, tmat
  sex        = 1,            # Female
  growth_model = "vb",
  iter       = 2000
)
```

### Option 2: Multiple Stanfit Objects

When parameters come from different model fits:

``` r
# Linf and L0 from growth fit, maturity from separate fits
mort <- get_stochastic_mortality(
  method              = "CW",
  growth_fit          = growth_fit,           # Provides Linf, L0
  length_maturity_fit = L50_fit,              # Provides Lmat
  age_maturity_fit    = t50_fit,              # Provides tmat
  maturity_cor        = 0.5,                  # Assumed Lmat-tmat correlation
  sex                 = 1,
  growth_model        = "vb",
  iter                = 2000
)
```

### Option 3: Manual Parameter Specification

When stanfit objects aren’t available:

``` r
mort <- get_stochastic_mortality(
  method       = "CW",
  Linf         = c(126, 10),   # c(mean, sd)
  L0           = c(35, 3),
  Lmat         = c(83, 5),
  tmat         = c(47, 4),
  maturity_cor = 0.5,          # Lmat-tmat correlation
  growth_model = "vb",
  iter         = 2000
)
```

### Option 4: Mixed Sources

You can combine stanfit objects with manual specifications:

``` r
# Growth parameters from fit, maturity manual
mort <- get_stochastic_mortality(
  method     = "CW",
  growth_fit = growth_fit,     # Provides Linf, L0
  Lmat       = c(83, 5),       # Manual
  tmat       = c(47, 4),       # Manual
  maturity_cor = 0.5,
  sex        = 1
)

# Or use birth_fit for L0
mort <- get_stochastic_mortality(
  method    = "CW",
  Linf      = c(126, 10),
  birth_fit = birth_fit,       # Provides L0 from b50
  Lmat      = c(83, 5),
  tmat      = c(47, 4)
)
```

### Parameter Priority

When multiple sources are available:

| Parameter | Priority                                    |
|-----------|---------------------------------------------|
| Linf      | growth_fit \> manual                        |
| L0        | growth_fit \> birth_fit \> manual           |
| Lmat      | growth_fit \> length_maturity_fit \> manual |
| tmat      | growth_fit \> age_maturity_fit \> manual    |

## Covariance Handling

### Why Covariance Matters

Life history parameters are typically correlated. For maturity
parameters, \\L\_{mat}\\ and \\t\_{mat}\\ are positively
correlated—larger individuals tend to mature at older ages. Ignoring
these correlations can generate biologically implausible parameter
combinations and artificially inflate uncertainty bounds.

### Correlation Scenarios

**1. Full Correlation Preservation**

When all parameters come from the same `growth_fit`:

``` r
mort <- get_stochastic_mortality(
  method     = "CW",
  growth_fit = growth_fit,  # Has Linf, L0, Lmat, tmat
  sex        = 1
)
# Message: "All parameters from growth_fit - correlations fully preserved."
```

**2. Bivariate Normal Sampling**

When \\L\_{mat}\\ and \\t\_{mat}\\ come from separate sources:

``` r
mort <- get_stochastic_mortality(
  method              = "CW",
  growth_fit          = growth_fit,
  length_maturity_fit = L50_fit,
  age_maturity_fit    = t50_fit,
  maturity_cor        = 0.6,
  sex                 = 1
)
# Message: "Lmat and tmat from separate fits - using specified correlation (rho = 0.60)."
```

**3. Independent Sampling**

Set `maturity_cor = 0` for independent sampling:

``` r
mort <- get_stochastic_mortality(
  method       = "CW",
  Linf         = c(126, 10),
  L0           = c(35, 3),
  Lmat         = c(83, 5),
  tmat         = c(47, 4),
  maturity_cor = 0,
  growth_model = "vb"
)
```

### Choosing `maturity_cor`

The default `maturity_cor = 0.5` reflects the biological expectation
that larger individuals mature later. Alternatives:

- **Estimate from data**: If you have paired observations
- **Literature values**: Species-specific correlations from published
  studies
- **Set to NA**: Function attempts to estimate from aligned posterior
  draws

## Mortality Models

### Chen-Watanabe (CW)

``` r
mort_cw <- get_stochastic_mortality(
  method     = "CW",
  growth_fit = growth_fit,
  sex        = 1,
  scaled     = TRUE,
  p          = 0.001
)

# With two-phase senescence
mort_cw_2p <- get_stochastic_mortality(
  method     = "CW",
  growth_fit = growth_fit,
  sex        = 1,
  two_phase  = TRUE,
  late_model = "gompertz"
)
```

### Peterson-Wroblewski (PW)

``` r
lw_func <- function(L) 0.0001 * L^3.1

mort_pw <- get_stochastic_mortality(
  method     = "PW",
  growth_fit = growth_fit,
  sex        = 1,
  lw_fun     = lw_func
)
```

### Lorenzen (L)

``` r
# Growth-based
mort_lor <- get_stochastic_mortality(
  method       = "L",
  growth_fit   = growth_fit,
  sex          = 1,
  weight_based = FALSE
)

# Weight-based
mort_lor_wt <- get_stochastic_mortality(
  method       = "L",
  growth_fit   = growth_fit,
  sex          = 1,
  weight_based = TRUE,
  lw_fun       = lw_func
)
```

## Output Structure

``` r
names(mort)
#> [1] "Schedules"  "Parameters" "Summary"    "Plot"

# Schedules: age-specific M for each parameter set
head(mort$Schedules)

# Parameters: sampled values with derived quantities
head(mort$Parameters)

# Summary: age-wise statistics
head(mort$Summary)
```

The plot caption indicates the correlation status and key model choices.

## Scaling Mortality

``` r
# Survival probability (default)
mort <- get_stochastic_mortality(
  method = "CW", growth_fit = growth_fit, sex = 1,
  scaled = TRUE, p = 0.001
)

# Empirical relationship
then_2015 <- function(tmax) 4.899 * tmax^(-0.916)
mort <- get_stochastic_mortality(
  method = "CW", growth_fit = growth_fit, sex = 1,
  scaled = TRUE, M_target = then_2015
)

# Fixed target
mort <- get_stochastic_mortality(
  method = "CW", growth_fit = growth_fit, sex = 1,
  scaled = TRUE, M_target = 0.15
)
```

## Integration with Survivorship

``` r
mort <- get_stochastic_mortality(
  method     = "CW",
  growth_fit = growth_fit,
  sex        = 1,
  iter       = 2000,
  print_plot = FALSE
)

surv <- simulate_survivorship(
  mc_object = mort,
  n         = 50000,
  n_iter    = 2000
)

surv$Aggregate$Age_of_Death
surv$Aggregate$Survival_to_tmat
```

## Reporting Results

Document parameter sources, correlation handling, and model choices:

> “Natural mortality was estimated using the Chen-Watanabe model within
> the vitalBayes unified framework. Growth parameters (\\L\_\infty\\,
> \\L_0\\) were drawn from the joint posterior of a von Bertalanffy
> growth model, while maturity parameters (\\L\_{mat}\\, \\t\_{mat}\\)
> were obtained from separate maturity ogive fits with assumed
> correlation \\\rho = 0.5\\. Mortality schedules were scaled to achieve
> 0.1% survival to maximum age (n = 2,000 Monte Carlo iterations).”

## See Also

- [`vignette("survivorship_simulation")`](https://brian-j-moe.github.io/vitalBayes/articles/survivorship_simulation.md) -
  Cohort survival analysis
- [`vignette("fit_bayesian_growth")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_growth.md) -
  Growth model fitting
- [`vignette("chen_watanabe_reparameterization")`](https://brian-j-moe.github.io/vitalBayes/articles/chen_watanabe_reparameterization.md) -
  Mathematical framework

## References

Chen, S., & Watanabe, S. (1989). Age dependence of natural mortality
coefficient in fish population dynamics. *Nippon Suisan Gakkaishi*,
55(2), 205-208.

Lorenzen, K. (2022). Size- and age-dependent natural mortality in fish
populations. *Fisheries Research*, 255, 106454.

Peterson, I., & Wroblewski, J. S. (1984). Mortality rate of fishes in
the pelagic ecosystem. *Canadian Journal of Fisheries and Aquatic
Sciences*, 41(7), 1117-1120.

Then, A. Y., et al. (2015). Evaluating the predictive performance of
empirical estimators of natural mortality rate. *ICES Journal of Marine
Science*, 72(1), 82-92.

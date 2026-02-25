# Maturity Ogive Estimation with fit_bayesian_maturity()

## Overview

The
[`fit_bayesian_maturity()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md)
function estimates length-at-maturity (\\L\_{50}\\) and/or
age-at-maturity (\\t\_{50}\\) using Bayesian probit regression. It
supports single-sex and two-sex models with optional **partial pooling**
between sexes.

## The Model

\\P(\text{mature}\_i) = \Phi\[\beta \cdot (x_i - x\_{50})\]\\

where \\x\\ is either length (\\L\\) or age (\\t\\), and \\x\_{50}\\ is
the value at which 50% of individuals are mature.

**Derived quantities:**

- \\x\_{05} = x\_{50} - 1.645/\beta\\ (5% maturity)
- \\x\_{95} = x\_{50} + 1.645/\beta\\ (95% maturity)
- Transition width \\= x\_{95} - x\_{05} = 3.29/\beta\\

## Basic Usage: Single-Sex Model

``` r
library(vitalBayes)
library(data.table)

# Load simulated data
data(growth_data)

# Filter to free-swimming females only
female_data <- growth_data[embryo == FALSE & sex == "female"]

# Fit length-at-maturity
L50_fit <- fit_bayesian_maturity(
 maturity = "mat",
 lt       = "fl",
 data     = female_data,
 chains   = 4,
 parallel = TRUE
)

# View results
L50_fit$summary(c("L50", "slope", "L05", "L95", "transition_width"))
```

## Two-Sex Models

When sex data are available, fit both sexes simultaneously:

``` r
# Prepare data (both sexes)
mat_data <- growth_data[embryo == FALSE & !is.na(mat)]

# Two-sex model with partial pooling (recommended)
L50_fit_2sex <- fit_bayesian_maturity(
 maturity    = "mat",
 lt          = "fl",
 sex         = "sex",
 data        = mat_data,
 use_pooling = TRUE,   # Enable partial pooling
 prior_tau   = 0.5     # Half-normal scale for between-sex SD
)

# Sex-specific estimates
L50_fit_2sex$summary("L50")

# Difference between sexes
L50_fit_2sex$summary("L50_diff")
```

### Why Partial Pooling?

Elasmobranch datasets often have imbalanced sex ratios. With 150 females
but only 34 males, the male estimate would have very wide credible
intervals if fit independently.

**Partial pooling** addresses this by modeling sex-specific parameters
as draws from a common distribution:

\\\log(L\_{50,s}) = \mu + \tau \cdot \eta_s, \quad \eta_s \sim
\mathcal{N}(0,1)\\

The model *learns* the appropriate degree of pooling (\\\tau\\) from the
data:

- Large \\\tau\\ → sexes are different → estimates stay separated
- Small \\\tau\\ → sexes are similar → estimates shrink together
- Sparse sex → borrows strength from the other sex

For a comprehensive treatment, see `vignette("partial_pooling")` or the
[Partial Pooling
section](https://brian-j-moe.github.io/vitalBayes/articles/Understanding_vitalBayes.html#maturity)
in the Statistical Methods guide.

## Fitting Both Length and Age Maturity

``` r
# Filter to individuals with both maturity and age data
mat_aged <- growth_data[embryo == FALSE & !is.na(mat) & !is.na(age)]

# Fit both models simultaneously
mat_fits <- fit_bayesian_maturity(
 maturity = "mat",
 lt       = "fl",
 age      = "age",
 sex      = "sex",
 data     = mat_aged,
 use_pooling = TRUE
)

# Access individual fits
mat_fits$length$summary("L50")
mat_fits$age$summary("t50")
```

## Sex Coding

The function auto-detects common sex coding conventions:

``` r
# All of these work automatically:
# sex = c("F", "M")           # English
# sex = c("Female", "Male")   # English full
# sex = c("female", "male")   # lowercase (used in simulated data)
# sex = c(1, 2)               # Numeric (1=female, 2=male)
# sex = c("Hembra", "Macho")  # Spanish
# sex = c("Femelle", "Mâle")  # French
# sex = c("メス", "オス")      # Japanese

# For non-standard coding, specify explicitly:
fit <- fit_bayesian_maturity(
 maturity = "mat",
 lt       = "fl",
 sex      = "sex",
 female   = "A",   # Your female code
 male     = "B",   # Your male code
 data     = your_data
)
```

## Prior Specification

Priors use a CV-based approach for intuitive specification:

``` r
# Custom priors
L50_fit <- fit_bayesian_maturity(
 maturity   = "mat",
 lt         = "fl",
 data       = mat_data,
 mean_L50   = 75,     # Prior mean (cm)
 cv_L50     = 0.2,    # 20% CV → SD = 15 cm
 mean_slope = 0,      # Log-scale prior for slope
 sd_slope   = 1
)
```

## Visualization

``` r
# Basic maturity ogive
plot_maturity_ogive(
 fit  = L50_fit_2sex,
 type = "length",
 data = mat_data,
 x_col = "fl",
 maturity_col = "mat",
 sex_col = "sex"
)

# With French labels
plot_maturity_ogive(
 fit        = L50_fit_2sex,
 type       = "length",
 data       = mat_data,
 sex_labels = c("female" = "Femelle", "male" = "Mâle"),
 x_lab      = "Longueur (cm)",
 y_lab      = "Probabilité de maturité"
)
```

## Posterior Predictive Checks

``` r
# Classification accuracy
L50_fit_2sex$summary("prop_correct_rep")

# Mean predicted probability by actual status
L50_fit_2sex$summary(c("mean_p_mature_f", "mean_p_mature_m",
                      "mean_p_immature_f", "mean_p_immature_m"))
```

## Using Output in Growth Models

The \\L\_{50}\\ and \\t\_{50}\\ estimates become informative priors for
the maturity-based growth parameterization:

``` r
# Fit t50 for age-at-maturity
t50_fit <- fit_bayesian_maturity(
 maturity = "mat",
 age      = "age",
 sex      = "sex",
 data     = mat_aged
)

# Pass to growth model
growth_fit <- fit_bayesian_growth(
 lt      = "fl",
 age     = "age",
 sex     = "sex",
 data    = growth_data[embryo == FALSE & !is.na(age)],
 k_based = FALSE,           # Use maturity-based parameterization
 L50_fit = L50_fit_2sex,    # Length-at-maturity fit
 t50_fit = t50_fit          # Age-at-maturity fit
)
```

## Troubleshooting

| Issue                       | Solution                                                |
|-----------------------------|---------------------------------------------------------|
| Divergent transitions       | Increase `adapt_delta`; check for separation in data    |
| Very wide \\\tau\\ estimate | Data may not support pooling; try `use_pooling = FALSE` |
| Poor classification         | Check maturity coding (must be 0/1)                     |
| Sex detection fails         | Use explicit `female`/`male` arguments                  |

## See Also

- `vignette("partial_pooling")` — Detailed treatment of hierarchical
  modeling for imbalanced data
- [Statistical Methods: Maturity
  Estimation](https://brian-j-moe.github.io/vitalBayes/articles/Understanding_vitalBayes.html#maturity)
  — Full mathematical derivation
- [Statistical Methods: Probit
  Link](https://brian-j-moe.github.io/vitalBayes/articles/Understanding_vitalBayes.html#probit)
  — Why probit over logit
- [`plot_maturity_ogive()`](https://brian-j-moe.github.io/vitalBayes/reference/plot_maturity_ogive.md)
  — Visualization
- [`fit_bayesian_growth()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md)
  — Using maturity estimates in growth models

------------------------------------------------------------------------

*This document is part of the vitalBayes R package. For bug reports,
feature requests, or questions, please visit the GitHub repository.*

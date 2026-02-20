# Growth Modeling with fit_bayesian_growth()

## Overview

The
[`fit_bayesian_growth()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md)
function fits von Bertalanffy, Gompertz, or Logistic growth models using
Bayesian methods. A central methodological feature is the
**maturity-based parameterization**, which derives the growth
coefficient \\k\\ from observable maturity metrics rather than
estimating it directly.

## Growth Model Equations

**von Bertalanffy** (`model = "v"`): \\L(t) = L\_\infty - (L\_\infty -
L_0) e^{-kt}\\

**Gompertz** (`model = "g"`): \\L(t) = L\_\infty
\exp\left\[-\ln\left(\frac{L\_\infty}{L_0}\right) e^{-kt}\right\]\\

**Logistic** (`model = "l"`): \\L(t) = \frac{L\_\infty}{1 +
\left(\frac{L\_\infty}{L_0} - 1\right) e^{-kt}}\\

## Two Parameterization Approaches

### 1. K-Based (Traditional)

Directly estimates the growth coefficient \\k\\:

``` r
library(vitalBayes)
library(data.table)

# Load simulated data
data(growth_data)

# Filter to non-embryos with age data
gdata <- growth_data[embryo == FALSE & !is.na(age)]

# Traditional k-based von Bertalanffy
growth_k <- fit_bayesian_growth(
  lt      = "fl",
  age     = "age",
  sex     = "sex",
  data    = gdata,
  model   = "v",
  k_based = TRUE,
  CV_k    = 0.5,        # 50% CV on k prior (high uncertainty)
  parallel = TRUE
)

growth_k$summary(c("Linf", "L0", "k", "sigma"))
```

### 2. Maturity-Based (Recommended)

Derives \\k\\ from maturity parameters (\\L\_{mat}\\, \\t\_{mat}\\),
ensuring biological consistency. The derivation differs by growth model:

**von Bertalanffy:** \\k = \frac{1}{t\_{mat}} \ln\left(\frac{L\_\infty -
L_0}{L\_\infty - L\_{mat}}\right)\\

**Gompertz:** \\k = \frac{1}{t\_{mat}} \ln\left(\frac{\ln(L\_\infty /
L_0)}{\ln(L\_\infty / L\_{mat})}\right)\\

**Logistic:** \\k = \frac{1}{t\_{mat}}
\ln\left(\frac{L\_{mat}(L\_\infty - L_0)}{L_0(L\_\infty -
L\_{mat})}\right)\\

Each derivation ensures the growth curve passes exactly through the
point \\(t\_{mat}, L\_{mat})\\, anchoring the model to an observable
biological milestone.

``` r
# First, fit maturity models
mat_data <- growth_data[embryo == FALSE & !is.na(mat)]

L50_fit <- fit_bayesian_maturity(
  maturity = "mat", lt = "fl", sex = "sex",
  data = mat_data,
  use_pooling = TRUE
)

t50_fit <- fit_bayesian_maturity(
  maturity = "mat", age = "age", sex = "sex",
  data = mat_data[!is.na(age)],
  use_pooling = TRUE
)

# Optional: fit birth model for L0 prior
birth_fit <- fit_bayesian_birth(
  embryo_lts = growth_data[embryo == TRUE, fl],
  free_swimming_lts = growth_data[embryo == FALSE, fl]
)

# Maturity-based growth model
growth_mat <- fit_bayesian_growth(
  lt        = "fl",
  age       = "age",
  sex       = "sex",
  data      = gdata,
  model     = "v",
  k_based   = FALSE,                    # Use maturity-based parameterization
  length.mature_stanfit = L50_fit,      # Provides Lmat prior
  age.mature_stanfit    = t50_fit,      # Provides tmat prior
  birth_stanfit         = birth_fit,    # Provides L0 prior
  parallel  = TRUE
)

# k is now a derived quantity
growth_mat$summary(c("Linf", "L0", "Lmat", "tmat", "k", "sigma"))
```

**Why maturity-based?**

1.  **Reduced correlation**: \\(L\_\infty, k)\\ are notoriously
    correlated; maturity parameters break this
2.  **Observable anchoring**: Maturity milestones fall within the data
    range, unlike \\L\_\infty\\
3.  **Prior propagation**: Uncertainty from maturity models flows into
    growth estimates
4.  **Biological coherence**: The growth curve *must* pass through
    \\(t\_{mat}, L\_{mat})\\

For the full derivation, see the [Maturity-Based Parameterization
section](https://brian-j-moe.github.io/vitalBayes/articles/Understanding_vitalBayes.html#maturity-param)
in the Statistical Methods guide.

## The \\L\_\infty\\ Constraint

A critical issue: unconstrained \\L\_\infty\\ often converges to values
*below* the largest observed individuals—biologically impossible.

**vitalBayes enforces** \\L\_\infty \> L\_{max}\\ (the maximum observed
length). For the biological and statistical rationale, see [Why
\\L\_\infty\\ Must Exceed Maximum Observed
Length](https://brian-j-moe.github.io/vitalBayes/articles/Understanding_vitalBayes.html#linf)
in the Statistical Methods guide.

``` r
# Lmax is auto-detected from data
growth_fit <- fit_bayesian_growth(
  lt   = "fl",
  age  = "age",
  data = gdata
)
# Message: "Lmax from data: 98.5 cm"

# Or specify manually (e.g., if you have length data without age)
growth_fit <- fit_bayesian_growth(
  lt   = "fl",
  age  = "age",
  data = gdata,
  Lmax = c(100, 95)  # Female, Male
)
```

## Two-Sex Models: Pooling Strategies

When sample sizes are imbalanced between sexes (common in elasmobranch
research), partial pooling borrows strength across sexes to reduce
uncertainty for the sparse group.

### The Double-Pooling Problem

A subtle issue arises when using **maturity-based parameterization**
with **partial pooling**: if the upstream maturity models
([`fit_bayesian_maturity()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md))
were themselves fit with `use_pooling = TRUE`, the maturity parameters
(\\L\_{mat}\\, \\t\_{mat}\\) already contain pooled estimates. Pooling
them *again* in the growth model can:

1.  Over-shrink sex differences toward the population mean
2.  In extreme cases, reverse genuine biological dimorphism
3.  Artificially tighten credible intervals

### Selective Pooling (Default)

The `pool_maturity` argument controls whether maturity parameters enter
the hierarchical structure:

``` r
# When both maturity fits are CmdStanMCMC objects from vitalBayes,
# pool_maturity auto-detects to FALSE (selective pooling)
growth_2sex <- fit_bayesian_growth(
  lt          = "fl",
  age         = "age",
  sex         = "sex",
  data        = gdata,
  model       = "v",
  k_based     = FALSE,
  length.mature_stanfit = L50_fit,
  age.mature_stanfit    = t50_fit,
  use_pooling = TRUE,     # Partial pooling enabled
  # pool_maturity = NULL  # Auto-detects to FALSE
  parallel    = TRUE
)
# Message: "Auto-detected vitalBayes maturity fits: using selective pooling"
```

**Under selective pooling (`pool_maturity = FALSE`):**

| Parameter     | Pooled? | Prior Source                        |
|---------------|---------|-------------------------------------|
| \\L\_\infty\\ | Yes     | Data-derived (needs regularization) |
| \\L_0\\       | Yes     | Birth model or default              |
| \\L\_{mat}\\  | No      | Direct from maturity fit            |
| \\t\_{mat}\\  | No      | Direct from maturity fit            |

This ensures \\L\_{mat}\\ and \\t\_{mat}\\ preserve their sex-specific
biological signal while \\L\_\infty\\ and \\L_0\\ benefit from
hierarchical shrinkage.

### Full Pooling with Anchoring

If you prefer full pooling (or must use it with manual priors), the
function uses **widened anchoring priors** (3x original SD) to prevent
over-constraint:

``` r
# Force full pooling explicitly
growth_full <- fit_bayesian_growth(
  lt          = "fl",
  age         = "age",
  sex         = "sex",
  data        = gdata,
  model       = "v",
  k_based     = FALSE,
  length.mature_stanfit = L50_fit,
  age.mature_stanfit    = t50_fit,
  use_pooling   = TRUE,
  pool_maturity = TRUE,   # Override auto-detection
  parallel      = TRUE
)
# Note: "pool_maturity = TRUE with vitalBayes maturity fits may cause 
#        double-pooling. Using widened priors (3x SD) to mitigate."
```

### Manual Priors (Auto-Detects Full Pooling)

When providing manual priors instead of vitalBayes fits, full pooling is
the default since there’s no prior pooling to double:

``` r
# Manual priors: pool_maturity defaults to TRUE
growth_manual <- fit_bayesian_growth(
  lt          = "fl",
  age         = "age", 
  sex         = "sex",
  data        = gdata,
  model       = "v",
  k_based     = FALSE,
  prior_Lmat  = rbind(c(72, 8), c(68, 8)),   # Female, Male: mean, SD
  prior_tmat  = rbind(c(13, 2), c(11, 2)),
  use_pooling = TRUE
  # pool_maturity auto-detects to TRUE
)
```

### Decision Guide

| Scenario                                                     | `pool_maturity`   | Rationale                                         |
|--------------------------------------------------------------|-------------------|---------------------------------------------------|
| vitalBayes maturity fits + `use_pooling = TRUE` in maturity  | `FALSE` (auto)    | Avoid double-pooling                              |
| vitalBayes maturity fits + `use_pooling = FALSE` in maturity | Could use `TRUE`  | Single pooling stage is safe                      |
| Manual priors                                                | `TRUE` (auto)     | No prior pooling to compound                      |
| Want maximum shrinkage                                       | `TRUE` (explicit) | Accept tighter CIs, check for dimorphism reversal |

## Comparing Growth Models

``` r
# Fit all three models
vb_fit <- fit_bayesian_growth(
  lt = "fl", age = "age", sex = "sex", data = gdata,
  model = "v", k_based = FALSE,
  length.mature_stanfit = L50_fit, age.mature_stanfit = t50_fit
)

gomp_fit <- fit_bayesian_growth(
  lt = "fl", age = "age", sex = "sex", data = gdata,
  model = "g", k_based = FALSE,
  length.mature_stanfit = L50_fit, age.mature_stanfit = t50_fit
)

logis_fit <- fit_bayesian_growth(
  lt = "fl", age = "age", sex = "sex", data = gdata,
  model = "l", k_based = FALSE,
  length.mature_stanfit = L50_fit, age.mature_stanfit = t50_fit
)

# Compare via LOO-CV
loo_vb <- compute_loo(vb_fit)
loo_gomp <- compute_loo(gomp_fit)
loo_logis <- compute_loo(logis_fit)

compare_loo(
  "von Bertalanffy" = loo_vb,
  "Gompertz" = loo_gomp,
  "Logistic" = loo_logis
)
```

## CV-Based Prior Specification

Priors are specified via coefficient of variation for intuitive scaling:

``` r
growth_fit <- fit_bayesian_growth(
  lt   = "fl",
  age  = "age",
  data = gdata,
  
  # Prior CVs (proportion of mean)
  CV_Linf = 0.20,    # 20% uncertainty on Linf
  CV_L0   = 0.30,    # 30% uncertainty on L0
  CV_k    = 0.50,    # 50% uncertainty on k (if k_based = TRUE)
  CV_Lmat = 0.20,    # 20% uncertainty on Lmat
  CV_tmat = 0.30,    # 30% uncertainty on tmat
  
  # Linf prior mean = 1.05 * Lmax by default
  Linf_multiplier = 1.05
)
```

## Visualization

``` r
# Basic growth curve
plot_growth_curve(
  fit        = growth_2sex,
  data       = gdata,
  age_col    = "age",
  length_col = "fl",
  sex_col    = "sex"
)

# Multilingual support
plot_growth_curve(
  fit        = growth_2sex,
  data       = gdata,
  sex_labels = c("female" = "Hembra", "male" = "Macho"),
  x_lab      = "Edad (años)",
  y_lab      = "Longitud (cm)"
)

# Compare models visually
compare_growth_models(
  "von Bertalanffy" = vb_fit,
  "Gompertz" = gomp_fit,
  "Logistic" = logis_fit,
  data = gdata,
  age_col = "age",
  length_col = "fl"
)
```

## Posterior Predictive Checks

``` r
# Built-in PPC metrics
growth_2sex$summary(c("rmse_f", "rmse_m", "mean_residual_f", "mean_residual_m"))

# Residual diagnostics
plot_residuals(
  fit        = growth_2sex,
  data       = gdata,
  age_col    = "age",
  length_col = "fl",
  type       = "all"
)
```

## Complete Workflow Example

``` r
# Load data
data(growth_data)

# ---- Stage 1: Birth ----
birth_fit <- fit_bayesian_birth(
  embryo_lts = growth_data[embryo == TRUE, fl],
  free_swimming_lts = growth_data[embryo == FALSE, fl]
)

# ---- Stage 2: Maturity ----
mat_data <- growth_data[embryo == FALSE & !is.na(mat)]

L50_fit <- fit_bayesian_maturity(
  maturity = "mat", lt = "fl", sex = "sex",
  data = mat_data, 
  use_pooling = TRUE
)

t50_fit <- fit_bayesian_maturity(
  maturity = "mat", age = "age", sex = "sex",
  data = mat_data[!is.na(age)], 
  use_pooling = TRUE
)

# ---- Stage 3: Growth ----
# Note: pool_maturity auto-detects to FALSE (selective pooling)
# since L50_fit and t50_fit are CmdStanMCMC objects
growth_fit <- fit_bayesian_growth(
  lt        = "fl",
  age       = "age",
  sex       = "sex",
  data      = growth_data[embryo == FALSE & !is.na(age)],
  model     = "v",
  k_based   = FALSE,
  birth_stanfit         = birth_fit,
  length.mature_stanfit = L50_fit,
  age.mature_stanfit    = t50_fit,
  use_pooling = TRUE
)

# ---- Summary ----
create_parameter_table(
  birth = birth_fit,
  L50 = L50_fit,
  t50 = t50_fit,
  growth = growth_fit
)
```

## Troubleshooting

| Issue                         | Solution                                                              |
|-------------------------------|-----------------------------------------------------------------------|
| Divergent transitions         | Increase `adapt_delta` (0.95 to 0.99)                                 |
| \\k\\ hitting boundaries      | Check that \\L\_{mat} \< L\_\infty\\ and \\L_0 \< L\_{mat}\\          |
| \\L\_\infty\\ too low         | Increase `Lmax` or `Linf_multiplier`                                  |
| Poor fit at young ages        | Consider different growth model (Gompertz often better for juveniles) |
| Sex differences reversed      | Check `pool_maturity`; try `pool_maturity = FALSE`                    |
| Over-tight credible intervals | May indicate double-pooling; use selective pooling                    |

## See Also

- [`vignette("partial_pooling")`](https://brian-j-moe.github.io/vitalBayes/articles/partial_pooling.md)
  — When and how to use hierarchical structure for imbalanced sex ratios
- [Statistical Methods: Growth
  Models](https://brian-j-moe.github.io/vitalBayes/articles/Understanding_vitalBayes.html#growth)
  — Full mathematical derivation
- [Statistical Methods: Maturity-Based
  Parameterization](https://brian-j-moe.github.io/vitalBayes/articles/Understanding_vitalBayes.html#maturity-param)
  — Methodological foundation
- [Statistical Methods: CV-Based
  Priors](https://brian-j-moe.github.io/vitalBayes/articles/Understanding_vitalBayes.html#cv-priors)
  — How priors are specified
- [`plot_growth_curve()`](https://brian-j-moe.github.io/vitalBayes/reference/plot_growth_curve.md),
  [`compare_growth_models()`](https://brian-j-moe.github.io/vitalBayes/reference/compare_growth_models.md)
  — Visualization
- [`compute_loo()`](https://brian-j-moe.github.io/vitalBayes/reference/compute_loo.md),
  [`compare_loo()`](https://brian-j-moe.github.io/vitalBayes/reference/compare_loo.md)
  — Model comparison

------------------------------------------------------------------------

*This document is part of the vitalBayes R package. For bug reports,
feature requests, or questions, please visit the GitHub repository.*

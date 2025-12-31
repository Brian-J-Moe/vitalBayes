# Partial Pooling for Imbalanced Sex Ratios

## The Problem: Imbalanced Sex Ratios

Elasmobranch datasets frequently have unequal sample sizes between
sexes. The `imbalanced_data` example dataset illustrates this common
scenario with 150 females but only 34 males. Fitting separate models for
each sex leads to:

- **Wide credible intervals** for the sparse sex (males)
- **Unstable estimates** driven by a few influential observations
- **Inefficient use of information** — we ignore the fact that both
  sexes are *the same species*

``` r
library(vitalBayes)
library(data.table)

# Load the imbalanced dataset (150F, 34M, 13 embryos)
data(imbalanced_data)

# Check the sex ratio
imbalanced_data[embryo == FALSE, .N, by = sex]
#    sex   N
# 1: female  150
# 2: male     34
```

## Three Estimation Strategies

Consider estimating a parameter \\\theta\\ (e.g., \\L\_{50}\\) for each
sex:

### 1. Complete Pooling (No Sex Effect)

Assume both sexes share identical parameters:

\\\theta\_{\text{female}} = \theta\_{\text{male}} = \theta\\

**Problem**: Ignores real biological differences between sexes.

### 2. No Pooling (Fully Separate)

Estimate each sex independently:

\\\theta\_{\text{female}} \perp \theta\_{\text{male}}\\

**Problem**: Ignores that both sexes are the same species. The sparse
sex gets unreliable estimates.

### 3. Partial Pooling (Hierarchical)

Model sex-specific parameters as draws from a common distribution:

\\\theta_s \sim \mathcal{N}(\mu, \tau^2)\\

where \\\mu\\ is the species-level mean and \\\tau\\ controls
between-sex variation.

**Key insight**: The model *learns* \\\tau\\ from the data:

- If sexes are similar → small \\\tau\\ → estimates shrink together
- If sexes are different → large \\\tau\\ → estimates stay separated
- Sparse sex → borrows strength from data-rich sex

## The vitalBayes Implementation

### Non-Centered Parameterization

For numerical stability, vitalBayes uses a non-centered
parameterization:

\\\log(\theta_s) = \mu + \tau \cdot \eta_s, \quad \eta_s \sim
\mathcal{N}(0, 1)\\

This separates the global mean (\\\mu\\) from sex-specific deviations
(\\\eta_s\\), which dramatically improves MCMC sampling when \\\tau\\ is
small.

### Prior on \\\tau\\

The between-sex standard deviation uses a half-normal prior:

\\\tau \sim \text{Half-Normal}(0, \sigma\_\tau)\\

where \\\sigma\_\tau\\ is set via the `prior_tau` argument. This prior:

- Allows \\\tau = 0\\ (complete pooling) when data support it
- Permits large \\\tau\\ when sexes genuinely differ
- Avoids the heavy tails of half-Cauchy that can cause divergences

## Practical Usage

### Enabling Partial Pooling

``` r
# Prepare maturity data from the imbalanced dataset
mat_data <- imbalanced_data[embryo == FALSE & !is.na(mat)]

# With partial pooling (recommended for imbalanced data)
L50_pooled <- fit_bayesian_maturity(
 maturity    = "mat",
 lt          = "fl",
 sex         = "sex",
 data        = mat_data,
 use_pooling = TRUE,    # Enable hierarchical structure
 prior_tau   = 0.5      # Half-normal scale (on log scale)
)

# Without partial pooling (separate estimation)
L50_unpooled <- fit_bayesian_maturity(
 maturity    = "mat",
 lt          = "fl",
 sex         = "sex",
 data        = mat_data,
 use_pooling = FALSE
)
```

### Comparing Results

``` r
# Pooled estimates
L50_pooled$summary("L50")

# Unpooled estimates
L50_unpooled$summary("L50")

# The pooled male estimate will typically have:
# 1. Similar median to unpooled
# 2. Narrower credible interval
# 3. Slight shrinkage toward the female estimate
```

### Formal Comparison

Use
[`compare_pooling()`](https://brian-j-moe.github.io/vitalBayes/reference/compare_pooling.md)
to quantify the difference:

``` r
compare_pooling(
 pooled   = L50_pooled,
 unpooled = L50_unpooled,
 params   = "L50"
)

# Output shows:
# - Point estimates for each approach
# - CI widths (pooled typically narrower for sparse sex)
# - Shrinkage magnitude
```

## Understanding Shrinkage

Partial pooling produces **shrinkage**: estimates for the sparse group
are pulled toward the overall mean. The amount of shrinkage depends on:

1.  **Sample size ratio**: More imbalance → more shrinkage for sparse
    sex
2.  **Observed difference**: Large apparent differences → less shrinkage
3.  **Within-group variance**: High variance → more shrinkage

### Visualizing Shrinkage

``` r
library(ggplot2)

# Extract posteriors
draws_pooled <- L50_pooled$draws("L50", format = "df")
draws_unpooled <- L50_unpooled$draws("L50", format = "df")

# Compare male estimates (the sparse sex)
ggplot() +
 geom_density(data = draws_unpooled, aes(x = `L50[2]`, fill = "Unpooled"), 
              alpha = 0.5) +
 geom_density(data = draws_pooled, aes(x = `L50[2]`, fill = "Pooled"), 
              alpha = 0.5) +
 labs(x = "Male L50 (cm)", y = "Density", 
      title = "Effect of Partial Pooling on Male L50 Estimate",
      subtitle = "Note the narrower credible interval with pooling") +
 theme_vital()
```

## When to Use Partial Pooling

### Recommended (use_pooling = TRUE)

- **Imbalanced sample sizes** between sexes (like the `imbalanced_data`
  example)
- **Moderate expected differences** between sexes
- **Sparse data** overall
- **Downstream use** of estimates (growth models, population models)

### Optional (use_pooling = FALSE)

- **Balanced sample sizes** and adequate data for both sexes (like the
  `growth_data` example)
- **Strong prior belief** that sexes are very different
- **Exploratory analysis** to assess sex differences without shrinkage

## Partial Pooling in Growth Models

The same principles apply to growth parameters:

``` r
# Prepare growth data from imbalanced dataset
gdata <- imbalanced_data[embryo == FALSE & !is.na(age)]

# Partial pooling on Linf, L0, k, Lmat, tmat
growth_pooled <- fit_bayesian_growth(
 lt          = "fl",
 age         = "age",
 sex         = "sex",
 data        = gdata,
 use_pooling = TRUE,
 prior_tau   = 0.2     # Tighter for growth (less between-sex variation expected)
)

# Sex-specific estimates with uncertainty reduction
growth_pooled$summary(c("Linf", "k"))

# Sex differences (still estimable!)
growth_pooled$summary(c("Linf_diff", "k_diff"))
```

## Choosing `prior_tau`

The `prior_tau` argument controls the half-normal scale for between-sex
SD:

| Value | Interpretation                 | Use Case                       |
|-------|--------------------------------|--------------------------------|
| 0.1   | Expect very similar sexes      | Growth rate, measurement error |
| 0.2   | Expect modest differences      | Linf, L0                       |
| 0.5   | Expect moderate differences    | L50, t50 (default)             |
| 1.0   | Expect substantial differences | Rare; consider no pooling      |

**Rule of thumb**: Since parameters are on log scale, `prior_tau = 0.5`
corresponds to roughly 50% expected variation between sexes (before
seeing data).

## Diagnosing Pooling Behavior

### Check the Estimated \\\tau\\

``` r
# Large tau = sexes are different
# Small tau = sexes are similar (more shrinkage)
L50_pooled$summary("tau_L50")
```

### Compare LOO-CV

``` r
loo_pooled <- compute_loo(L50_pooled)
loo_unpooled <- compute_loo(L50_unpooled)

compare_loo(
 "Partial Pooling" = loo_pooled,
 "No Pooling" = loo_unpooled
)

# Pooled model often has better (higher) elpd when:
# - Sample sizes are imbalanced
# - True sex differences are modest
```

## Common Questions

### Does pooling bias my estimates?

No. Partial pooling produces **regularized** estimates that trade a
small amount of bias for substantial variance reduction. For the sparse
sex, this trade-off almost always improves mean squared error.

### Can I still detect sex differences?

Yes! The model estimates `L50_diff`, `Linf_diff`, etc. — the posterior
difference between sexes. If sexes truly differ, the credible interval
for these differences will exclude zero.

### What if sexes are truly very different?

The model will estimate a large \\\tau\\, which reduces shrinkage. In
the limit of very large \\\tau\\, partial pooling converges to no
pooling.

### Should I always use pooling?

For two-sex elasmobranch models with typical sample sizes, yes. The only
case to avoid pooling is when you have abundant, balanced data for both
sexes (like the `growth_data` example with 189F, 176M) AND strong prior
belief in very different parameters.

## Example: Full Workflow with Pooling

``` r
# ---- Data Prep ----
# Use the imbalanced dataset to demonstrate pooling benefits
data(imbalanced_data)
mat_data <- imbalanced_data[embryo == FALSE & !is.na(mat)]
gdata <- imbalanced_data[embryo == FALSE & !is.na(age)]

# ---- Maturity with Pooling ----
L50_fit <- fit_bayesian_maturity(
 maturity = "mat", lt = "fl", sex = "sex",
 data = mat_data,
 use_pooling = TRUE,
 prior_tau = 0.5
)

t50_fit <- fit_bayesian_maturity(
 maturity = "mat", age = "age", sex = "sex",
 data = mat_data[!is.na(age)],
 use_pooling = TRUE,
 prior_tau = 0.5
)

# ---- Growth with Pooling ----
growth_fit <- fit_bayesian_growth(
 lt = "fl", age = "age", sex = "sex",
 data = gdata,
 k_based = FALSE,
 L50_fit = L50_fit,
 t50_fit = t50_fit,
 use_pooling = TRUE,
 prior_tau = 0.2
)

# ---- Results ----
# Sex-specific estimates with appropriate uncertainty
growth_fit$summary(c("Linf", "k"))

# Credible sex differences
growth_fit$summary(c("Linf_diff", "k_diff"))
```

## Comparing Datasets: Balanced vs Imbalanced

The package includes datasets that illustrate when pooling matters most:

``` r
# Imbalanced data (150F, 34M) - pooling helps
data(imbalanced_data)
imbalanced_data[embryo == FALSE, .N, by = sex]

# Balanced data (189F, 176M) - pooling still works but less critical
data(growth_data)
growth_data[embryo == FALSE, .N, by = sex]

# Limited data (24F, 18M) - pooling essential for both sexes
data(limited_data)
limited_data[embryo == FALSE, .N, by = sex]
```

## Summary

| Aspect          | No Pooling    | Partial Pooling           |
|-----------------|---------------|---------------------------|
| Sparse sex CI   | Wide          | Narrower                  |
| Point estimates | Unstable      | Regularized               |
| Shrinkage       | None          | Data-adaptive             |
| Sex differences | Direct        | Via difference parameters |
| Recommended for | Balanced data | Imbalanced data           |

For most elasmobranch life history analyses, **partial pooling is the
recommended default** when fitting two-sex models.

## See Also

- [Statistical Methods: Partial
  Pooling](https://brian-j-moe.github.io/vitalBayes/articles/Understanding_vitalBayes.html#maturity)
  — Full mathematical derivation with non-centered parameterization
  details
- [Statistical Methods: CV-Based
  Priors](https://brian-j-moe.github.io/vitalBayes/articles/Understanding_vitalBayes.html#cv-priors)
  — How prior_tau works
- [`vignette("fit_bayesian_maturity")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_maturity.md)
  — Maturity model fitting
- [`vignette("fit_bayesian_growth")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_growth.md)
  — Growth model fitting
- [`vignette("model_diagnostics")`](https://brian-j-moe.github.io/vitalBayes/articles/model_diagnostics.md)
  — Comparing pooled vs unpooled via LOO-CV

------------------------------------------------------------------------

*This document is part of the vitalBayes R package. For bug reports,
feature requests, or questions, please visit the GitHub repository.*

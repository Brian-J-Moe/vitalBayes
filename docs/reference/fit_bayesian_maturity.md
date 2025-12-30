# Fit Bayesian Maturity Models

Fits binomial regression(s) with probit link function for estimating
length-at-50%-maturity (L50) and/or age-at-50%-maturity (t50). Supports
both single-sex and hierarchical two-sex models with optional partial
pooling between sexes.

The probit link is chosen for consistency with the birth model and for
its threshold-crossing interpretation: latent developmental readiness is
assumed normally distributed, and an individual matures when readiness
exceeds a threshold. This interpretation is more biologically intuitive
than log-odds, which are not standard reporting metrics in elasmobranch
research.

## Usage

``` r
fit_bayesian_maturity(
  maturity,
  lt = NULL,
  age = NULL,
  sex = NULL,
  female = NULL,
  male = NULL,
  data = NULL,
  mean_L50 = NULL,
  mean_t50 = NULL,
  cv_L50 = 0.3,
  cv_t50 = 0.3,
  mean_slope = 0,
  sd_slope = 1,
  use_pooling = TRUE,
  prior_tau = 0.5,
  parallel = TRUE,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 500,
  seed = 1234,
  ...
)
```

## Arguments

- maturity:

  Column name in `data` or numeric vector of binary maturity indicators
  (0 = immature, 1 = mature).

- lt:

  Optional column name or vector of observed lengths (\> 0). Required
  for L50.

- age:

  Optional column name or vector of observed ages (\>= 0). Required for
  t50.

- sex:

  Optional column name or vector for sex. If provided, fits two-sex
  model. Supports auto-detection of common coding conventions including:
  F/M, Female/Male, 1/2 (numeric), and equivalents in Spanish,
  Portuguese, French, German, Italian, Japanese, Chinese, and Russian.
  See Details.

- female:

  Optional. Explicit specification of how females are coded in the data.
  If `NULL` (default), auto-detection is attempted. Use with `male` for
  non-standard coding schemes.

- male:

  Optional. Explicit specification of how males are coded in the data.
  Must be provided together with `female`.

- data:

  A data.frame or data.table containing referenced columns.

- mean_L50, mean_t50:

  Prior means on natural scale. If `NULL` (default), computed from data
  as midpoint of transition zone. For two-sex models, can be length-2
  vectors for sex-specific priors.

- cv_L50, cv_t50:

  Coefficient of variation for priors. Default 0.3.

- mean_slope, sd_slope:

  Prior mean and SD for slope on log scale. Default `mean_slope = 0`,
  `sd_slope = 1`.

- use_pooling:

  Logical. For two-sex models, use partial pooling? Default `TRUE`.

- prior_tau:

  Scale for half-normal prior on between-sex SD. Default 0.5.

- parallel:

  Logical. Run chains in parallel? Default `TRUE`.

- chains:

  Integer. Number of MCMC chains. Default 4.

- iter_warmup:

  Integer. Warmup iterations per chain. Default 1000.

- iter_sampling:

  Integer. Sampling iterations per chain. Default 1000.

- refresh:

  Integer. Progress update frequency. Default 500.

- seed:

  Integer. Random seed. Default 1234.

- ...:

  Additional arguments passed to `$sample()`.

## Value

For single model (length or age only): A `CmdStanMCMC` object. For both
models: A named list with `$length` and `$age` fits.

## Details

The linear predictor is parameterized directly in terms of x50: \$\$p_i
= \Phi\[\beta \times (x_i - x\_{50})\]\$\$

where \\\Phi\\ is the standard normal CDF, and x represents either
length or age.

## Prior Specification

Priors are constructed using a coefficient of variation approach. The
prior mean for L50/t50 is automatically computed as the midpoint between
the smallest mature individual and largest immature individual. The
prior SD is `mean * cv`.

## Initialization

Initial values are set to the data-derived midpoint estimates, ensuring
chains start near the high-probability region of the posterior.

## See also

[`vignette("fit_bayesian_maturity")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_maturity.md)
for usage examples.

[`vignette("partial_pooling")`](https://brian-j-moe.github.io/vitalBayes/articles/partial_pooling.md)
for detailed explanation of hierarchical modeling for imbalanced sex
ratios.

[Statistical Methods: Maturity
Estimation](https://brian-j-moe.github.io/vitalBayes/doc/Understanding_vitalBayes.html#maturity)
for mathematical derivation including partial pooling.

[Statistical Methods: Probit
Link](https://brian-j-moe.github.io/vitalBayes/doc/Understanding_vitalBayes.html#probit)
for justification of probit over logit.

[`fit_bayesian_birth`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_birth.md),
[`fit_bayesian_growth`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md),
[`plot_maturity_ogive`](https://brian-j-moe.github.io/vitalBayes/reference/plot_maturity_ogive.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Load simulated data
data(growth_data)

# Single-sex length-at-maturity (females only)
L50_fit <- fit_bayesian_maturity(
  maturity = "mat",
  lt       = "fl",
  data     = growth_data[sex == "female" & embryo == FALSE]
)

# Two-sex model with partial pooling
L50_fit_2sex <- fit_bayesian_maturity(
  maturity    = "mat",
  lt          = "fl",
  sex         = "sex",
  data        = growth_data[embryo == FALSE],
  use_pooling = TRUE
)

# Example with imbalanced data (150F, 34M) - pooling is especially helpful
data(imbalanced_data)
L50_imbalanced <- fit_bayesian_maturity(
  maturity    = "mat",
  lt          = "fl",
  sex         = "sex",
  data        = imbalanced_data[embryo == FALSE],
  use_pooling = TRUE  # Borrows strength across sexes
)

# Fit both length and age maturity models
mat_fits <- fit_bayesian_maturity(
  maturity = "mat",
  lt       = "fl",
  age      = "age",
  sex      = "sex",
  data     = growth_data[embryo == FALSE]
)
} # }
```

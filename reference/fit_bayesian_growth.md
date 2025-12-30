# Fit Bayesian Individual Growth Models

Fits von Bertalanffy, Gompertz, or Logistic growth models using
precompiled Stan models via the instantiate package. Supports either
traditional k-based or maturity-based parameterization, and both
single-sex and hierarchical two-sex models with optional partial pooling
between sexes.

## Usage

``` r
fit_bayesian_growth(
  lt,
  age,
  sex = NULL,
  female = NULL,
  male = NULL,
  data = NULL,
  model = c("v", "g", "l"),
  k_based = FALSE,
  birth_stanfit = NULL,
  length.mature_stanfit = NULL,
  age.mature_stanfit = NULL,
  Lmax = NULL,
  Linf_multiplier = 1.05,
  CV_Linf = 0.2,
  CV_k = 0.5,
  CV_L0 = 0.3,
  CV_Lmat = 0.2,
  CV_tmat = 0.3,
  prior_L0 = NULL,
  prior_Lmat = NULL,
  prior_tmat = NULL,
  prior_k = NULL,
  use_pooling = TRUE,
  prior_tau = 0.2,
  robust = FALSE,
  loc_sig = 0,
  scale_sig = 1,
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

- lt:

  Column name in `data` or numeric vector of observed lengths (\> 0).

- age:

  Column name in `data` or numeric vector of observed ages (\>= 0).

- sex:

  Optional column name or vector for sex. If provided, fits two-sex
  model. Supports auto-detection of common coding conventions including:
  F/M, Female/Male, 1/2 (numeric), and equivalents in multiple
  languages. See
  [`fit_bayesian_maturity`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md)
  for full list.

- female:

  Optional. Explicit specification of how females are coded in the data.
  If `NULL` (default), auto-detection is attempted.

- male:

  Optional. Explicit specification of how males are coded in the data.
  Must be provided together with `female`.

- data:

  A data.frame or data.table containing referenced columns. May include
  incomplete cases (length without age) which are used only to determine
  Lmax.

- model:

  Character. Growth model type: `"v"` (von Bertalanffy), `"g"`
  (Gompertz), or `"l"` (Logistic). Default `"v"`.

- k_based:

  Logical. If `TRUE`, uses traditional k-based parameterization. If
  `FALSE` (default), derives k from maturity parameters (Lmat, tmat).

- birth_stanfit:

  Optional CmdStanMCMC from
  [`fit_bayesian_birth`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_birth.md).
  Used to set informative priors for L0.

- length.mature_stanfit:

  Optional CmdStanMCMC from
  [`fit_bayesian_maturity`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md)
  with L50 parameter. Used to set priors for Lmat.

- age.mature_stanfit:

  Optional CmdStanMCMC from
  [`fit_bayesian_maturity`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md)
  with t50 parameter. Used to set priors for tmat.

- Lmax:

  Numeric. Maximum observed length by sex for setting Linf prior and
  lower bound. If `NULL` (default), computed from all length data in
  `data` (including rows with missing age). Can be a single value
  (applied to both sexes) or length-2 vector for sex-specific values
  c(female, male).

- Linf_multiplier:

  Numeric. Multiplier for Lmax to get Linf prior mean. Default 1.05
  (i.e., 5% larger than observed max).

- CV_Linf:

  Coefficient of variation for Linf prior. Default 0.2.

- CV_k:

  Coefficient of variation for k prior. Default 0.5.

- CV_L0:

  Coefficient of variation for L0 prior (if not from birth fit). Default
  0.3.

- CV_Lmat:

  Coefficient of variation for Lmat prior (if not from maturity fit).
  Default 0.2.

- CV_tmat:

  Coefficient of variation for tmat prior (if not from maturity fit).
  Default 0.3.

- prior_L0:

  Manual prior for L0 on natural scale: `c(mean, sd)`. Ignored if
  `birth_stanfit` provided.

- prior_Lmat:

  Manual prior for Lmat on natural scale. Only for maturity-based.

- prior_tmat:

  Manual prior for tmat on natural scale. Only for maturity-based.

- prior_k:

  Manual prior for k on natural scale. Only for k-based models.

- use_pooling:

  Logical. For two-sex models, use partial pooling on location
  parameters? Default `TRUE`. Slope parameters are never pooled.

- prior_tau:

  Scale for half-normal priors on between-sex SD. Default 0.2.

- robust:

  Logical. Use Student-t observation model? Default `FALSE`.

- loc_sig, scale_sig:

  Location and scale for Cauchy prior on sigma. Defaults: `loc_sig = 0`,
  `scale_sig = 1`.

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

A `CmdStanMCMC` object.

## Details

### Growth Model Equations

**von Bertalanffy (model = "v"):** \$\$L(t) = L\_\infty - (L\_\infty -
L_0) e^{-kt}\$\$

**Gompertz (model = "g"):** \$\$L(t) = L\_\infty
\exp\left(-\ln\left(\frac{L\_\infty}{L_0}\right) e^{-kt}\right)\$\$

**Logistic (model = "l"):** \$\$L(t) = \frac{L\_\infty}{1 +
\left(\frac{L\_\infty}{L_0} - 1\right) e^{-kt}}\$\$

### Prior Specification

Priors are constructed using coefficient of variation (CV) arguments:

**Linf:** Prior mean defaults to `1.05 * Lmax`, where Lmax is the
maximum observed length (including incomplete cases without age data).
The prior SD is `mean * CV_Linf`. Users can override Lmax via the `Lmax`
argument.

**k:** Prior mean is estimated from the data by solving the growth
equation for k at each observation and averaging. Prior SD is
`mean * CV_k`.

**L0, Lmat, tmat:** Can be informed by prior birth/maturity model fits,
or specified manually.

### Initialization

Initial values are set to data-derived estimates. For k, the
Gulland-Holt linearization method is used to provide a reasonable
starting value.

## See also

[`vignette("fit_bayesian_growth")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_growth.md)
for usage examples.

[`vignette("partial_pooling")`](https://brian-j-moe.github.io/vitalBayes/articles/partial_pooling.md)
for hierarchical modeling of imbalanced sex ratios.

[Statistical Methods: Growth
Models](https://brian-j-moe.github.io/vitalBayes/doc/Understanding_vitalBayes.html#growth)
for equations and model comparison.

[Statistical Methods: Maturity-Based
Parameterization](https://brian-j-moe.github.io/vitalBayes/doc/Understanding_vitalBayes.html#maturity-param)
for the derivation of k from maturity parameters.

[Statistical Methods: L-infinity
Constraint](https://brian-j-moe.github.io/vitalBayes/doc/Understanding_vitalBayes.html#linf)
for why Linf must exceed maximum observed length.

[`fit_bayesian_birth`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_birth.md),
[`fit_bayesian_maturity`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md),
[`plot_growth_curve`](https://brian-j-moe.github.io/vitalBayes/reference/plot_growth_curve.md),
[`compare_growth_models`](https://brian-j-moe.github.io/vitalBayes/reference/compare_growth_models.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Load simulated data
data(growth_data)

# Workflow: birth -> maturity -> growth

# Step 1: Fit birth model
birth_fit <- fit_bayesian_birth(
  embryo_lts = growth_data[embryo == TRUE, fl],
  free_swimming_lts = growth_data[embryo == FALSE, fl]
)

# Step 2: Fit maturity models
mat_data <- growth_data[embryo == FALSE & !is.na(mat)]
L50_fit <- fit_bayesian_maturity(
  maturity = "mat", lt = "fl", sex = "sex",
  data = mat_data, use_pooling = TRUE
)
t50_fit <- fit_bayesian_maturity(
  maturity = "mat", age = "age", sex = "sex",
  data = mat_data[!is.na(age)], use_pooling = TRUE
)

# Step 3: Fit growth model (maturity-based, pooled)
growth_fit <- fit_bayesian_growth(
  lt = "fl", age = "age", sex = "sex",
  data = growth_data[embryo == FALSE & !is.na(age)],
  model = "vb",
  k_based = FALSE,
  birth_fit = birth_fit,
  L50_fit = L50_fit,
  t50_fit = t50_fit,
  use_pooling = TRUE
)

# View key parameters
growth_fit$summary(c("Linf", "L0", "k", "Lmat", "tmat"))

# Example with imbalanced data
data(imbalanced_data)
gdata <- imbalanced_data[embryo == FALSE & !is.na(age)]
growth_imbal <- fit_bayesian_growth(
  lt = "fl", age = "age", sex = "sex",
  data = gdata,
  use_pooling = TRUE  # Important for sparse male data
)
} # }
```

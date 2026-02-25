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
  CV_delta = 0.5,
  CV_k = 0.5,
  CV_L0 = 0.3,
  CV_Lmat = 0.2,
  CV_tmat = 0.3,
  prior_L0 = NULL,
  prior_Lmat = NULL,
  prior_tmat = NULL,
  prior_k = NULL,
  use_pooling = TRUE,
  pool_maturity = NULL,
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

- CV_delta:

  Coefficient of variation for the \\\delta_L\\ excess in the
  delta-gamma \\L\_\infty\\ prior. Controls the shape of the
  \\\text{Gamma}(\alpha, \beta)\\ distribution on \\\delta_L =
  L\_\infty - L\_{max}\\. A value of 0.50 means "50\\ uncertainty about
  how far above \\L\_{max}\\ the asymptote sits." Default 0.50.

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

- pool_maturity:

  Logical or `NULL`. For two-sex maturity-based models with
  `use_pooling = TRUE`, controls whether maturity parameters (Lmat,
  tmat) are included in the hierarchical pooling structure.

  - `NULL` (default): Auto-detect based on maturity fit source. Set to
    `FALSE` if both `length.mature_stanfit` and `age.mature_stanfit` are
    CmdStanMCMC objects (vitalBayes fits); otherwise `TRUE`.

  - `FALSE`: Selective pooling. Only Linf and L0 are pooled; Lmat and
    tmat use direct sex-specific priors. Recommended when maturity
    priors come from vitalBayes fits (avoids double-pooling).

  - `TRUE`: Full pooling with anchoring. All parameters are pooled, but
    maturity parameters use widened priors (3x SD) to prevent
    over-constraint.

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

A `CmdStanMCMC` object with the following key parameters:

- Linf:

  Asymptotic length(s) by sex

- L0:

  Length at age 0 (birth size) by sex

- k:

  Growth coefficient(s) by sex (derived if maturity-based)

- Lmat, tmat:

  Maturity parameters (if maturity-based)

- sigma:

  Observation error SD by sex

- \*\_diff:

  Sex differences (female - male) for each parameter

- log_lik:

  Pointwise log-likelihood for LOO-CV

- y_pred, y_rep:

  Posterior predictions and replications

## Details

### Growth Model Equations

**von Bertalanffy (model = "v"):** \$\$L(t) = L\_\infty - (L\_\infty -
L_0) e^{-kt}\$\$

**Gompertz (model = "g"):** \$\$L(t) = L\_\infty
\exp\left(-\ln\left(\frac{L\_\infty}{L_0}\right) e^{-kt}\right)\$\$

**Logistic (model = "l"):** \$\$L(t) = \frac{L\_\infty}{1 +
\left(\frac{L\_\infty}{L_0} - 1\right) e^{-kt}}\$\$

### Maturity-Based Parameterization

When `k_based = FALSE`, the growth coefficient \\k\\ is derived from
maturity parameters rather than estimated directly:

\$\$k = \frac{1}{t\_{mat}} \ln\left(\frac{L\_\infty - L_0}{L\_\infty -
L\_{mat}}\right)\$\$

This parameterization breaks the notorious \\(L\_\infty, k)\\
correlation and anchors the growth curve to observable maturity
milestones.

### Selective Pooling (Two-Sex Models)

When fitting two-sex models with maturity-based parameterization and
partial pooling enabled, the function implements **selective pooling**
to avoid double-pooling maturity parameters that already carry
information from upstream maturity model fits.

The `pool_maturity` argument controls this behavior:

- `pool_maturity = FALSE` (default when vitalBayes maturity fits
  provided): \\L\_\infty\\ uses soft pairwise shrinkage on the log scale
  and \\L_0\\ is pooled hierarchically. Maturity parameters
  (\\L\_{mat}\\, \\t\_{mat}\\) use direct sex-specific priors from the
  upstream fits, preserving biological signal without risk of
  double-pooling.

- `pool_maturity = TRUE` (default when manual priors provided): All
  parameters are pooled, but maturity parameters use widened anchoring
  priors (3x original SD) to prevent over-constraint while still
  allowing some regularization.

The auto-detection logic sets `pool_maturity = FALSE` when
`length.mature_stanfit` and `age.mature_stanfit` are both CmdStanMCMC
objects (i.e., vitalBayes maturity model fits), since these already
contain pooled estimates when `use_pooling = TRUE` was used in the
maturity models.

### Prior Specification

Priors are constructed using coefficient of variation (CV) arguments:

**Linf:** Uses a delta-gamma parameterization to avoid boundary pile-up
near Lmax. Internally, \\\delta = L\_\infty - L\_{max}\\ is estimated
with \\\delta \sim \text{Gamma}(\alpha, \beta)\\, where \\\alpha =
1/\text{CV\\delta}^2\\ and \\\beta = 1/(\mu\_\delta \cdot
\text{CV\\delta}^2)\\. The prior mean for \\\delta\\ defaults to
`Lmax * (Linf_multiplier - 1)`. `CV_delta` controls uncertainty about
the excess above \\L\_{max}\\, not about \\L\_\infty\\ itself (see
vignette for details). For two-sex models with pooling, between-sex
shrinkage is applied via a soft penalty on \\\log(L\_\infty)\\,
preserving the proportional interpretation of `prior_tau`.

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
for usage examples with gulper shark data.

`vignette("partial_pooling")` for hierarchical modeling of imbalanced
sex ratios.

[Statistical Methods: Growth
Models](https://brian-j-moe.github.io/vitalBayes/doc/vitalBayes_stats_explained.html#growth)
for equations and model comparison.

[Statistical Methods: Maturity-Based
Parameterization](https://brian-j-moe.github.io/vitalBayes/doc/vitalBayes_stats_explained.html#maturity-param)
for the derivation of k from maturity parameters.

[Statistical Methods: L-infinity
Constraint](https://brian-j-moe.github.io/vitalBayes/doc/vitalBayes_stats_explained.html#linf)
for why Linf must exceed maximum observed length.

[`fit_bayesian_birth`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_birth.md),
[`fit_bayesian_maturity`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md),
[`plot_growth_curve`](https://brian-j-moe.github.io/vitalBayes/reference/plot_growth_curve.md),
[`compare_growth_models`](https://brian-j-moe.github.io/vitalBayes/reference/compare_growth_models.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Workflow: birth -> maturity -> growth

# Step 1: Fit birth model
birth_fit <- fit_bayesian_birth(
  embryo_lts = sharks[embryo == TRUE, length],
  free_swimming_lts = sharks[embryo == FALSE, length]
)

# Step 2: Fit maturity models (with pooling)
mat_length_fit <- fit_bayesian_maturity(
  maturity = mat, lt = length, sex = sex,
  data = sharks[embryo == FALSE],
  use_pooling = TRUE
)

mat_age_fit <- fit_bayesian_maturity(
  maturity = mat, age = age, sex = sex,
  data = sharks[embryo == FALSE & !is.na(age)],
  use_pooling = TRUE
)

# Step 3: Fit growth model (maturity-based, selective pooling)
# pool_maturity auto-detects to FALSE since maturity fits are CmdStanMCMC
growth_fit <- fit_bayesian_growth(
  lt = length, age = age, sex = sex,
  data = sharks,
  model = "v",
  k_based = FALSE,
  birth_stanfit = birth_fit,
  length.mature_stanfit = mat_length_fit,
  age.mature_stanfit = mat_age_fit,
  use_pooling = TRUE,
  robust = TRUE
)

# Explicit control: force full pooling even with vitalBayes fits
growth_fit_full <- fit_bayesian_growth(
  ...,
  pool_maturity = TRUE  # Override auto-detection
)

# Alternative: Manual priors (auto-detects pool_maturity = TRUE)
growth_fit_manual <- fit_bayesian_growth(
  lt = length, age = age, sex = sex,
  data = sharks,
  k_based = FALSE,
  prior_Lmat = rbind(c(70, 10), c(65, 10)),  # Female, Male
  prior_tmat = rbind(c(12, 2), c(10, 2)),
  use_pooling = TRUE
)
} # }
```

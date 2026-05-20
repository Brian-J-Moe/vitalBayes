# Growth Modeling with fit_bayesian_growth()

## Overview

The
[`fit_bayesian_growth()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md)
function fits von Bertalanffy, Gompertz, or Logistic growth models using
Bayesian methods. A central design choice is the **maturity-based
parameterization**, which derives the growth coefficient \\k\\ from
observable maturity metrics rather than estimating it directly. This
anchors the growth curve to biologically interpretable quantities and
propagates uncertainty from upstream maturity estimates into the growth
posterior.

## Growth Model Equations

All three growth models are implemented in \\L_0\\-based form, where
\\L_0\\ is length at birth — a directly observable quantity, unlike the
VB parameter \\t_0\\ (the theoretical age at length zero). This ensures
that priors can be set from embryo or neonate measurements.

### Von Bertalanffy

\\L(t) = L\_\infty - (L\_\infty - L_0) \\ e^{-kt}\\

The von Bertalanffy (VB) model describes growth as a constant
exponential approach to \\L\_\infty\\. The absolute growth rate
\\dL/dt\\ is highest at birth and declines monotonically, with an
inflection point at \\0.632 \\ L\_\infty\\. This pattern suits species
where somatic growth is fastest in early life and steadily decelerates
as the organism approaches asymptotic size.

The VB model is the most widely used growth function in fisheries
science and provides the baseline parameterization for the Chen-Watanabe
natural mortality model (see
`vignette("chen_watanabe_reparameterization")`). However, it is
sensitive to data coverage at extreme ages: sparse adult observations
can produce biologically implausible \\L\_\infty\\ estimates, and the
strong \\L\_\infty\\-\\k\\ correlation under traditional
parameterization can lead to poorly identified posteriors.

### Gompertz

\\L(t) = L\_\infty \exp\\\left\[-\ln\\\left(\frac{L\_\infty}{L_0}\right)
e^{-kt}\right\]\\

The Gompertz model describes growth where the *rate of deceleration*
itself decelerates — an exponentially declining growth rate rather than
a linearly declining one. Growth is initially rapid but slows earlier
and more abruptly than under VB dynamics. The inflection point occurs at
\\L\_\infty / e \approx 0.368 \\ L\_\infty\\, much earlier in ontogeny
than the VB inflection.

The early inflection makes the Gompertz particularly suitable for
species exhibiting rapid juvenile growth followed by pronounced
deceleration. Many small coastal elasmobranchs (bonnetheads, Atlantic
sharpnose sharks, some skate species) show this pattern, with growth
effectively ceasing near reproductive maturity. The Gompertz model also
tends to produce more stable \\L\_\infty\\ estimates than the VB when
adult data are sparse, since it doesn’t require the data to resolve the
late-asymptotic curvature as precisely.

### Logistic

\\L(t) = \frac{L\_\infty}{1 + \left(\frac{L\_\infty}{L_0} - 1\right)
e^{-kt}}\\

The Logistic model is symmetric around its inflection point at
\\L\_\infty / 2\\, with growth accelerating before the midpoint and
decelerating after. This sigmoidal trajectory is less commonly applied
in fisheries, but it can be appropriate for species where early juvenile
growth is initially slow (due to nutritional constraints, habitat
transitions, or ontogenetic diet shifts) before accelerating during a
rapid growth phase and then decelerating toward asymptotic size.

In practice, the Logistic model provides the tightest approach to
\\L\_\infty\\ — it reaches asymptotic size faster than either VB or
Gompertz for equivalent parameter values. This can sometimes produce
\\L\_\infty\\ estimates uncomfortably close to the maximum observed
length.

## Two Parameterization Approaches

### K-Based (Traditional)

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

### Maturity-Based (Recommended)

The maturity-based parameterization derives \\k\\ from the maturity
milestones \\(L\_{mat}, t\_{mat})\\, ensuring that the growth curve
passes through the maturity point and that prior information from
maturity models propagates into the growth posterior. The derivation is
model-specific: substituting \\L(t\_{mat}) = L\_{mat}\\ into each growth
equation and solving for \\k\\ yields:

**Von Bertalanffy:**

\\k\_{VB} = \frac{1}{t\_{mat}} \ln\\\left(\frac{L\_\infty -
L_0}{L\_\infty - L\_{mat}}\right)\\

**Gompertz:**

\\k_g = \frac{1}{t\_{mat}} \ln\\\left(\frac{\ln(L\_\infty /
L_0)}{\ln(L\_\infty / L\_{mat})}\right)\\

**Logistic:**

\\k_l = \frac{1}{t\_{mat}} \ln\\\left(\frac{L\_{mat}(L\_\infty -
L_0)}{L_0(L\_\infty - L\_{mat})}\right)\\

In each case, \\k\\ is a deterministic function of \\(L\_\infty, L_0,
L\_{mat}, t\_{mat})\\ — parameters that are either directly estimated in
the Stan model or informed by upstream birth and maturity fits. This
eliminates the need to specify a prior on \\k\\ directly and breaks the
\\L\_\infty\\-\\k\\ posterior correlation that plagues traditional
parameterizations. Note that \\k\\ has different meanings across the
three models and should not be directly compared (e.g., Gompertz \\k_g\\
is generally larger than VB \\k\_{VB}\\ for the same species). What
matters for biological inference is the resulting growth trajectory, not
the numerical value of \\k\\.

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
  length.mature_stanfit = L50_fit,       # Provides Lmat prior
  age.mature_stanfit    = t50_fit,       # Provides tmat prior
  birth_stanfit         = birth_fit,     # Provides L0 prior
  parallel  = TRUE
)

# k is now a derived quantity
growth_mat$summary(c("Linf", "L0", "Lmat", "tmat", "k", "sigma"))
```

**Why maturity-based?** The maturity-based parameterization offers
several advantages. First, it reduces posterior correlation:
\\(L\_\infty, k)\\ are typically heavily correlated under traditional
parameterization, and replacing \\k\\ with \\(L\_{mat}, t\_{mat})\\
breaks this collinearity. Second, it provides observable anchoring:
maturity milestones fall within the observed data range, unlike
\\L\_\infty\\ which is an extrapolation beyond the largest individuals.
Third, prior information from upstream maturity models propagates
directly into growth estimates, providing informative regularization for
both sexes. Fourth, it enforces biological coherence: the fitted growth
curve necessarily passes through \\(t\_{mat}, L\_{mat})\\, which is a
real biological constraint rather than a mathematical convenience.

For the full derivation, see the model equations reference or the Stan
source code.

## The \\L\_\infty\\ Constraint

A critical issue in growth modeling: unconstrained \\L\_\infty\\ often
converges to values *below* the largest observed individuals —
biologically impossible. vitalBayes enforces \\L\_\infty \> L\_{max}\\
(the maximum observed length in the data), with the prior mean set at
`Linf_multiplier` \\\times\\ \\L\_{max}\\ (default 1.05, i.e., 5% above
maximum observed length).

``` r

# Lmax is auto-detected from data (including rows without age)
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

Note that the `data` argument can include incomplete cases (length
without age), which are used only to determine \\L\_{max}\\. This is
useful when your dataset contains measured individuals that were not
aged — their lengths still inform the plausible range of \\L\_\infty\\.

## Observation Model

By default, vitalBayes models observation error as lognormal — the log
of predicted length plus Gaussian noise. This ensures positive
predictions, accommodates the typically multiplicative nature of growth
measurement error (larger individuals have proportionally larger
errors), and produces well-behaved likelihoods.

For datasets with outliers or heavy-tailed residuals, the
`robust = TRUE` option switches to a Student-t observation model, which
downweights extreme observations:

``` r

growth_robust <- fit_bayesian_growth(
  lt = "fl", age = "age", sex = "sex",
  data = gdata,
  robust = TRUE  # Student-t instead of lognormal
)
```

## Two-Sex Models: Pooling Strategies

When sample sizes are imbalanced between sexes (common in elasmobranch
research), partial pooling borrows strength across sexes to reduce
uncertainty for the sparse group. For the general theory of partial
pooling, see `vignette("partial_pooling")`.

### The Double-Pooling Problem

A subtle issue arises when using **maturity-based parameterization**
with **partial pooling**: if the upstream maturity models
([`fit_bayesian_maturity()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md))
were themselves fit with `use_pooling = TRUE`, the maturity parameters
(\\L\_{mat}\\, \\t\_{mat}\\) already contain pooled estimates. Pooling
them *again* in the growth model can over-shrink sex differences toward
the population mean, in extreme cases reversing genuine biological
dimorphism, and can artificially tighten credible intervals.

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
function uses **widened anchoring priors** (3\\\times\\ original SD) to
prevent over-constraint:

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

When providing manual priors instead of vitalBayes fits, pooling across
all parameters is the default since there’s no prior pooling to double:

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

| Scenario | `pool_maturity` | Rationale |
|----|----|----|
| vitalBayes maturity fits + `use_pooling = TRUE` in maturity | `FALSE` (auto) | Avoid double-pooling |
| vitalBayes maturity fits + `use_pooling = FALSE` in maturity | Could use `TRUE` | Single pooling stage is safe |
| Manual priors | `TRUE` (auto) | No prior pooling to compound |
| Want maximum shrinkage | `TRUE` (explicit) | Accept tighter CIs, check for dimorphism reversal |

## Comparing Growth Models

Selecting among growth models is a standard part of the vitalBayes
workflow. All three models can be fit with the same data and prior
structure, then compared via LOO-CV:

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

The selected model can then be passed to
[`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md)
for mortality estimation, regardless of which growth model was chosen
(see `vignette("chen_watanabe_reparameterization")`).

## CV-Based Prior Specification

Priors are specified via coefficient of variation for intuitive,
scale-invariant control:

``` r

growth_fit <- fit_bayesian_growth(
  lt   = "fl",
  age  = "age",
  data = gdata,

  # Prior CVs (proportion of mean)
  CV_delta = 0.50,   # 50% uncertainty on delta (excess above Lmax)
  CV_L0    = 0.30,   # 30% uncertainty on L0
  CV_k     = 0.50,   # 50% uncertainty on k (if k_based = TRUE)
  CV_Lmat  = 0.20,   # 20% uncertainty on Lmat
  CV_tmat  = 0.30,   # 30% uncertainty on tmat

  # Linf prior mean = 1.05 * Lmax by default
  Linf_multiplier = 1.05
)
```

Most CVs operate on the parameter itself (e.g., `CV_L0 = 0.30` means 30%
relative uncertainty on \\L_0\\). The exception is `CV_delta`, which
controls uncertainty about the *excess* \\\delta_L = L\_\infty -
L\_{max}\\ rather than about \\L\_\infty\\ directly. At the default
`CV_delta = 0.50`, this produces a gamma prior on \\\delta_L\\ with
shape \\\alpha = 4\\, giving a proper mode and well-behaved HMC
geometry. The `Linf_multiplier` controls the prior *mean* of the excess
(default: 5% above \\L\_{max}\\), while `CV_delta` controls how tightly
concentrated the prior is around that mean.

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

| Issue | Solution |
|----|----|
| Divergent transitions | Increase `adapt_delta` (0.95 → 0.99) |
| \\k\\ hitting boundaries | Check that \\L\_{mat} \< L\_\infty\\ and \\L_0 \< L\_{mat}\\ |
| \\L\_\infty\\ too low | Increase `Lmax` or `Linf_multiplier` |
| \\L\_\infty\\ boundary pile-up (posterior median at \\L\_{max}\\) | Increase `CV_delta` or try a different growth model (logistic naturally prefers smaller excess) |
| Poor fit at young ages | Consider different growth model (Gompertz often better for juveniles) |
| Sex differences reversed | Check `pool_maturity`; try `pool_maturity = FALSE` |
| Over-tight credible intervals | May indicate double-pooling; use selective pooling |

## See Also

- `vignette("partial_pooling")` — When and how to use hierarchical
  structure for imbalanced sex ratios
- `vignette("chen_watanabe_reparameterization")` — How growth parameters
  feed into CW mortality estimation
- [`vignette("mortality_estimation")`](https://brian-j-moe.github.io/vitalBayes/articles/mortality_estimation.md)
  — Natural mortality from growth posteriors
- [`vignette("model_diagnostics")`](https://brian-j-moe.github.io/vitalBayes/articles/model_diagnostics.md)
  — LOO-CV model comparison and convergence checks
- [`plot_growth_curve()`](https://brian-j-moe.github.io/vitalBayes/reference/plot_growth_curve.md),
  [`compare_growth_models()`](https://brian-j-moe.github.io/vitalBayes/reference/compare_growth_models.md)
  — Visualization
- [`compute_loo()`](https://brian-j-moe.github.io/vitalBayes/reference/compute_loo.md),
  [`compare_loo()`](https://brian-j-moe.github.io/vitalBayes/reference/compare_loo.md)
  — Model comparison

------------------------------------------------------------------------

*This document is part of the vitalBayes R package. For bug reports,
feature requests, or questions, please visit the GitHub repository.*

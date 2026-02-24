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
function addresses these challenges by generating Monte Carlo
age-specific mortality schedules from growth model posteriors. By
drawing parameter sets from the joint posterior distribution of a fitted
growth model, the function propagates parameter uncertainty through to
mortality estimates, producing credible intervals that honestly reflect
our knowledge (or lack thereof) about the population. Three mortality
models are available, each with distinct theoretical foundations, and
all are compatible with any growth model fitted via
[`fit_bayesian_growth()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md).

This vignette covers:

- Equations and biological rationale for each mortality model
- Scaling mortality schedules to empirical targets
- The Monte Carlo workflow and parameter handling
- Model comparison and sensitivity analysis
- Integration with downstream survival simulation

## Mortality Models

### Chen-Watanabe (1989)

The Chen-Watanabe model derives age-specific mortality from the
relationship between body growth and mortality risk under von
Bertalanffy dynamics. In the \\L_0\\-parameterized form used by
vitalBayes, the model is:

\\M(t) = \frac{k\_{VB} \cdot L\_\infty}{L(t)}\\

where \\k\_{VB}\\ is the von Bertalanffy growth coefficient serving as
the asymptotic mortality rate (\\M\_\infty\\), \\L\_\infty\\ is
asymptotic length, and \\L(t)\\ is predicted length at age \\t\\.
Mortality is inversely proportional to body size: small, young
individuals experience the highest mortality, which declines as the
organism grows toward \\L\_\infty\\.

This equation involves \\k\\ in two distinct roles. The \\k\_{VB}\\ in
the numerator is a *theoretical constant* derived from the CW
formulation under von Bertalanffy assumptions — it must always be the VB
growth coefficient, regardless of which growth model was actually
fitted. The \\L(t)\\ in the denominator is a *growth trajectory
prediction* — it should come from whatever model best describes the
organism’s actual growth dynamics. When a Gompertz or Logistic model is
fitted, vitalBayes uses the native growth equation with native \\k\\ for
\\L(t)\\, while computing a separate VB-equivalent \\k\\ for
\\M\_\infty\\. The VB-equivalent is derived from the shared biological
milestones \\(L\_\infty, L_0, L\_{mat}, t\_{mat})\\:

\\k\_{VB}^{equiv} = \frac{1}{t\_{mat}} \ln\\\left(\frac{L\_\infty -
L_0}{L\_\infty - L\_{mat}}\right)\\

For k-based fits without maturity milestones, a derivative-matching
fallback is used (see *Derivative-Based Conversion* below).

The Chen-Watanabe model also supports a two-phase extension for
late-life senescence, in which mortality transitions from the declining
CW trajectory to an increasing senescence function (Gompertz or
Logistic) at a fraction of the age at maturity. For the complete
mathematical development of the CW reparameterization, the normalized
growth coefficient \\G(t)\\, and the senescence extensions, see
`vignette("chen_watanabe_reparameterization")`.

### Peterson-Wroblewski (1984)

The Peterson-Wroblewski model estimates mortality from an allometric
relationship between body weight and mortality rate, derived from
empirical data on pelagic fish populations:

\\M(W) = 1.92 \cdot W(t)^{-0.25}\\

where \\W(t)\\ is body weight in grams at age \\t\\, computed from
predicted length via a user-supplied length-weight function \\W =
f(L)\\. The \\-0.25\\ exponent reflects metabolic scaling theory:
metabolic rate scales approximately as \\W^{0.75}\\ (Kleiber’s law), and
mortality is assumed proportional to mass-specific metabolic rate,
yielding a \\W^{-0.25}\\ dependence.

Unlike the CW model, the PW model requires a length-weight relationship
function (provided via the `lw_fun` argument). Because the model
operates on predicted body weight rather than growth model parameters
directly, it is inherently growth-model-agnostic. Any growth model that
predicts \\L(t)\\ can be combined with any \\L\\-\\W\\ conversion to
produce mortality estimates. The length-at-age predictions use the
native growth model’s \\k\\ (not a VB-equivalent), since the PW model
was not derived under VB-specific assumptions.

### Lorenzen (1996, 2022)

The Lorenzen model relates natural mortality to body size through
empirical calibration across a broad taxonomic range. vitalBayes
supports two formulations.

**Weight-based** (Lorenzen, 1996):

\\M(W) = \alpha \cdot W(t)^{\beta}\\

where \\\alpha\\ and \\\beta\\ are allometric parameters estimated from
a meta-analysis of marine and freshwater fishes (\\\alpha \sim
\mathcal{N}(3.69, \\ 0.502)\\, \\\beta \sim \mathcal{N}(-0.305, \\
0.029)\\). Like the PW model, this is purely weight-based: \\k\\ does
not appear in the mortality equation and serves only to predict \\L(t)\\
via the native growth model. The native \\k\\ and native growth equation
should always be used.

**Growth-based** (Lorenzen, 2022):

\\\ln M = 0.28 - 1.30 \\ \ln\\\left(\frac{L(t)}{L\_\infty}\right) + 1.08
\\ \ln(k\_{VB})\\

This formulation has the same dual-\\k\\ structure as Chen-Watanabe. The
\\\ln(k\_{VB})\\ term was calibrated against von Bertalanffy parameters
(Lorenzen, 2022), so it must be the VB growth coefficient. The
\\L(t)/L\_\infty\\ ratio represents relative body size at age and should
use the native growth trajectory for the most accurate prediction. When
inputs come from a Gompertz or Logistic fit, vitalBayes passes the
native \\k\\ for \\L(t)\\ prediction and a separate \\k\_{VB}^{equiv}\\
for the regression term. The growth-based formulation has the advantage
of not requiring a length-weight relationship, making it useful when
\\L\\-\\W\\ data are unavailable.

### Derivative-Based \\k\_{VB}\\ Conversion

When a k-based Gompertz or Logistic fit is used without maturity
milestones, the VB-equivalent \\k\\ cannot be derived from \\(L\_{mat},
t\_{mat})\\. In this case, vitalBayes falls back to derivative matching
at birth — computing the VB growth coefficient that produces the same
instantaneous growth rate at age 0:

\\k\_{VB}^{equiv} = k_g \cdot \frac{L_0 \cdot \ln(L\_\infty /
L_0)}{L\_\infty - L_0} \quad \text{(Gompertz)}\\

\\k\_{VB}^{equiv} = k_l \cdot \frac{L_0}{L\_\infty} \quad
\text{(Logistic)}\\

These closed-form expressions follow from equating the VB growth rate at
birth, \\dL/dt\|\_{t=0} = k\_{VB}(L\_\infty - L_0)\\, with the
corresponding Gompertz or Logistic growth rate. The milestone-based
derivation is preferred when maturity data are available, since it
anchors the conversion at a biologically meaningful interior point of
the growth trajectory.

## Scaling Mortality Schedules

### Why Scale?

Theoretical mortality models (CW, PW, Lorenzen) produce age-specific
*shapes* that are biologically informative, but their absolute levels
often do not match empirical observations. This is expected: the models
were derived from general relationships across broad taxonomic groups,
not calibrated to the specific population under study. Scaling adjusts
the overall mortality level while preserving the age-specific pattern,
anchoring the curve to an independent empirical estimate.

### The Scaling Transformation

Given a raw mortality schedule \\M\_{raw}(t)\\ and a target mean
mortality \\\bar{M}\_{target}\\, the scaled schedule is:

\\M\_{scaled}(t) = M\_{raw}(t) \times
\frac{\bar{M}\_{target}}{\bar{M}\_{raw}}\\

where \\\bar{M}\_{raw} = \frac{1}{n} \sum_t M\_{raw}(t)\\ is the mean of
the unscaled schedule. This multiplicative adjustment preserves the
relative differences between ages — if neonatal mortality is 5\\\times\\
adult mortality in the raw schedule, it remains 5\\\times\\ in the
scaled schedule.

### Target Specification

The scaling target \\\bar{M}\_{target}\\ can be specified in three ways.

**Survival probability** (`M_target = NULL, p = ...`): Specify the
proportion of a cohort expected to survive to maximum age \\t\_{max}\\.
The target is derived as:

\\\bar{M}\_{target} = -\frac{\ln(p)}{t\_{max}}\\

A common choice is \\p = 0.001\\ (0.1% survival to \\t\_{max}\\), which
produces moderate scaling for long-lived species.

**Empirical relationship** (`M_target = function(tmax) ...`): Supply a
function of \\t\_{max}\\, such as the Then et al. (2015) estimator:

\\\bar{M}\_{target} = 4.899 \cdot t\_{max}^{-0.916}\\

This approach leverages the strong empirical correlation between maximum
age and mean mortality across fish taxa.

**Fixed value** (`M_target = 0.15`): Directly specify a target mean
mortality from literature estimates or independent analyses.

``` r
library(vitalBayes)

# Survival probability scaling (0.1% survive to tmax)
mort_p <- get_stochastic_mortality(
  method = "CW", growth_fit = growth_fit, sex = 1,
  scaled = TRUE, M_target = NULL, p = 0.001
)

# Then et al. (2015) empirical relationship
then_2015 <- function(tmax) 4.899 * tmax^(-0.916)
mort_then <- get_stochastic_mortality(
  method = "CW", growth_fit = growth_fit, sex = 1,
  scaled = TRUE, M_target = then_2015
)

# Fixed target from independent analysis
mort_fixed <- get_stochastic_mortality(
  method = "CW", growth_fit = growth_fit, sex = 1,
  scaled = TRUE, M_target = 0.15
)
```

## Basic Usage

### Standard Workflow with vitalBayes Fits

The typical workflow begins with growth (and optionally maturity) model
fitting, then feeds posteriors into mortality estimation:

``` r
library(vitalBayes)
library(data.table)

# Step 1: Fit growth model (any of the three)
# growth_fit <- fit_bayesian_growth(...)
# maturity_fit <- fit_bayesian_maturity(...)

# Step 2: Generate Chen-Watanabe mortality schedules
mort <- get_stochastic_mortality(
  method       = "CW",
  growth_fit   = growth_fit,
  maturity_fit = maturity_fit,
  sex          = 2,        # Males
  iter         = 2000,     # Number of posterior draws
  scaled       = TRUE,
  p            = 0.001,    # 0.1% survival to tmax
  two_phase    = TRUE,     # Include senescence
  late_model   = "gompertz",
  print_plot   = TRUE
)
```

The function automatically detects the growth model type from the fit
object and extracts \\(L\_\infty, L_0, L\_{mat}, t\_{mat})\\ from the
joint posterior. When the growth model is Gompertz or Logistic, two
things happen: the native \\k\\ and growth equation are used for
\\L(t)\\ prediction, and a separate VB-equivalent \\k\\ is computed for
mortality models that require it (CW and growth-based Lorenzen). For
maturity-based fits, \\k\_{VB}^{equiv}\\ is derived from the biological
milestones; for k-based fits, a derivative-matching conversion is
applied. Joint posterior sampling preserves parameter correlations —
critically, the strong negative correlation between \\L\_\infty\\ and
\\k\\ that is typical of growth model fits.

### Output Structure

``` r
names(mort)
#> [1] "Schedules"  "Parameters" "Summary"    "Plot"

# Full age-specific schedules for every posterior draw
head(mort$Schedules)
#>    set_id  age     M_raw  M_scaled age_round
#>     <int> <num>     <num>    <num>     <num>
#> 1:      1   0.5  0.42318  0.31204       0.5
#> 2:      1   1.0  0.33721  0.24870       1.0

# Parameter draws used (includes both k values and growth model type)
head(mort$Parameters)
#>    set_id   Linf     L0   Lmat  tmat k_native k_vb_equiv growth_model  tmax
#>     <int>  <num>  <num>  <num> <num>    <num>      <num>       <char> <num>
#> 1:      1  98.42  24.31  68.92  11.2   0.1521     0.0945     gompertz  34.2

# Summary statistics across posterior draws
head(mort$Summary)
#>    age_round M_median M_mean M_lower M_upper
#>        <num>    <num>  <num>   <num>   <num>
#> 1:       0.5   0.3012 0.3124  0.2018  0.4532
```

### Manual Parameter Specification

For literature-based analyses where posterior draws are unavailable,
parameters can be specified as mean-SD pairs. Two sampling modes are
supported, mirroring the k-based and maturity-based parameterizations in
[`fit_bayesian_growth()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md).

#### k-Based (Independent Sampling)

When only `k` is provided without `Lmat`, all parameters are sampled
independently from their specified normal distributions. The `k` value
should be the native growth coefficient for the model specified by
`growth_model`. For VB fits (the default), the native \\k\\ is already
the VB growth coefficient. For Gompertz or Logistic fits, the function
automatically derives \\k\_{VB}^{equiv}\\ via derivative matching at
birth and uses the native growth equation for \\L(t)\\ prediction.

``` r
# VB parameters (backward-compatible: k IS the VB k)
mort_vb <- get_stochastic_mortality(
  method = "CW",
  Linf   = c(100, 8),     # mean, sd
  L0     = c(25, 2),
  k      = c(0.08, 0.01), # VB k
  tmat   = c(12, 1.5),
  iter   = 2000,
  scaled = TRUE,
  p      = 0.001
)

# Gompertz parameters (derivative matching for k_vb)
mort_gomp <- get_stochastic_mortality(
  method       = "CW",
  Linf         = c(100, 8),
  L0           = c(25, 2),
  k            = c(0.15, 0.02),  # Gompertz k (NOT VB k)
  tmat         = c(12, 1.5),
  growth_model = "gompertz",     # Must match source of k
  iter         = 2000,
  scaled       = TRUE,
  p            = 0.001
)
```

Independent sampling destroys any correlations between parameters. This
typically produces wider (less precise) uncertainty bands than joint
posterior sampling from a fitted model, and can occasionally generate
parameter combinations that are individually plausible but jointly
inconsistent.

#### Maturity-Based (Bivariate Sampling)

When both `Lmat` and `tmat` are provided, they are sampled jointly from
a bivariate normal distribution via Cholesky factorization. This
preserves the positive biological correlation between length- and
age-at-maturity: larger individuals tend to mature at older ages. The
VB-equivalent \\k\\ is then *derived* from the sampled milestones
\\(L\_\infty, L_0, L\_{mat}, t\_{mat})\\ rather than taken from the user
directly.

``` r
mort_mat <- get_stochastic_mortality(
  method        = "CW",
  Linf          = c(100, 8),
  L0            = c(25, 2),
  Lmat          = c(72, 4),       # length-at-maturity: mean, sd
  tmat          = c(12, 1.5),     # age-at-maturity: mean, sd
  k             = c(0.15, 0.02),  # Native Gompertz k (for L(t) trajectory)
  growth_model  = "gompertz",
  rho_Lmat_tmat = 0.5,            # Lmat-tmat correlation (default)
  iter          = 2000,
  scaled        = TRUE,
  p             = 0.001
)
```

The `rho_Lmat_tmat` argument controls the strength of the assumed
correlation. The default of 0.5 reflects the broad biological
expectation that larger individuals mature later. Alternative values can
be obtained from paired individual-level data, from computing the
empirical correlation between aligned posteriors of separate maturity
model fits, from published species-specific estimates, or by running
sensitivity analyses across a range of \\\rho\\ values.

When both `Lmat`/`tmat` and `k` are provided, the milestone-derived
\\k\_{VB}^{equiv}\\ is used for CW and growth-based Lorenzen (as the
\\M\_\infty\\ constant), while the user-supplied `k` is retained as the
native growth coefficient for \\L(t)\\ prediction under the specified
`growth_model`. This is the preferred configuration for non-VB growth
models: the maturity milestones anchor the VB-equivalent \\k\\ at a
biologically meaningful point, and the native \\k\\ produces the most
accurate growth trajectory.

If `k` is omitted in this mode, the milestone-derived
\\k\_{VB}^{equiv}\\ is used for *both* roles, and \\L(t)\\ defaults to
the VB equation regardless of `growth_model`. This is appropriate when
only maturity data are available, but the function will issue a warning
for non-VB `growth_model` settings.

#### Sensitivity to \\\rho\\

When the correlation between maturity parameters is unknown, running the
analysis across a range of \\\rho\\ values quantifies how sensitive the
mortality estimates are to this assumption:

``` r
rho_values <- c(0.3, 0.5, 0.7)
mort_rho <- lapply(rho_values, function(r) {
  get_stochastic_mortality(
    method        = "CW",
    Linf          = c(100, 8),
    L0            = c(25, 2),
    Lmat          = c(72, 4),
    tmat          = c(12, 1.5),
    rho_Lmat_tmat = r,
    iter          = 2000,
    print_plot    = FALSE
  )
})
```

In practice, sensitivity to \\\rho\\ is usually modest when the marginal
standard deviations of \\L\_{mat}\\ and \\t\_{mat}\\ are small relative
to their means. The primary effect of higher \\\rho\\ is to tighten the
uncertainty bands on derived \\k\\ (and therefore \\M\\), because
correlated sampling generates fewer extreme parameter combinations than
independent sampling.

## Comparing Mortality Models

### Across Mortality Methods

Comparing mortality estimates from different theoretical frameworks
provides a measure of structural uncertainty:

``` r
# Chen-Watanabe
mort_cw <- get_stochastic_mortality(
  method = "CW", growth_fit = growth_fit, sex = 1,
  scaled = TRUE, p = 0.001, print_plot = FALSE
)

# Peterson-Wroblewski (requires length-weight function)
lw_fun <- function(L) 0.0001 * L^3.1
mort_pw <- get_stochastic_mortality(
  method = "PW", growth_fit = growth_fit, sex = 1,
  lw_fun = lw_fun, scaled = TRUE, p = 0.001,
  print_plot = FALSE
)

# Lorenzen growth-based
mort_lor <- get_stochastic_mortality(
  method = "L", growth_fit = growth_fit, sex = 1,
  weight_based = FALSE, scaled = TRUE, p = 0.001,
  print_plot = FALSE
)

# Combine for comparison plot
library(ggplot2)
combined <- rbind(
  mort_cw$Summary[, model := "Chen-Watanabe"],
  mort_pw$Summary[, model := "Peterson-Wroblewski"],
  mort_lor$Summary[, model := "Lorenzen"]
)

ggplot(combined, aes(x = age_round, color = model, fill = model)) +
  geom_ribbon(aes(ymin = M_lower, ymax = M_upper), alpha = 0.2, color = NA) +
  geom_line(aes(y = M_median), linewidth = 1) +
  scale_color_manual(values = vital_palette(3)) +
  scale_fill_manual(values = vital_palette(3)) +
  labs(x = "Age (years)", y = "Natural Mortality (M)",
       title = "Mortality Model Comparison") +
  theme_bw() +
  theme(legend.position = "top")
```

All three models predict declining mortality with age, but they differ
in their curvature, early-age magnitude, and late-life behavior.
Comparing across models gives a sense of how sensitive management
conclusions are to the choice of theoretical framework.

### Across Growth Models

A complementary diagnostic tests whether mortality estimates are robust
to the choice of growth model:

``` r
# Fit all three growth models
# vb_fit   <- fit_bayesian_growth(..., model = "v")
# gomp_fit <- fit_bayesian_growth(..., model = "g")
# logis_fit <- fit_bayesian_growth(..., model = "l")

# CW mortality from each growth model (auto-detects model type)
mort_from_vb <- get_stochastic_mortality(
  method = "CW", growth_fit = vb_fit, sex = 1, print_plot = FALSE
)
mort_from_gomp <- get_stochastic_mortality(
  method = "CW", growth_fit = gomp_fit, sex = 1, print_plot = FALSE
)
mort_from_logis <- get_stochastic_mortality(
  method = "CW", growth_fit = logis_fit, sex = 1, print_plot = FALSE
)

combined_growth <- rbind(
  mort_from_vb$Summary[, growth_model := "von Bertalanffy"],
  mort_from_gomp$Summary[, growth_model := "Gompertz"],
  mort_from_logis$Summary[, growth_model := "Logistic"]
)

ggplot(combined_growth, aes(x = age_round, color = growth_model, fill = growth_model)) +
  geom_ribbon(aes(ymin = M_lower, ymax = M_upper), alpha = 0.2, color = NA) +
  geom_line(aes(y = M_median), linewidth = 1) +
  labs(x = "Age (years)", y = "Natural Mortality (M)",
       title = "CW Mortality by Growth Model",
       subtitle = "Native L(t) trajectories with VB-equivalent M_inf") +
  theme_bw() +
  theme(legend.position = "top")
```

Differences between growth models in CW mortality now arise from two
sources: the VB-equivalent \\k\\ (which may differ across models if the
milestone-based inversions produce slightly different values), and —
more importantly — the native \\L(t)\\ trajectories, which have
different curvatures between anchor points. Models with faster early
growth (like Gompertz) predict larger \\L(t)\\ at young ages and
therefore lower juvenile mortality, while models that approach
\\L\_\infty\\ less tightly (like Logistic) may produce slightly
different late-age behavior. If the three curves are similar, mortality
estimates are robust to growth model choice. Divergence indicates
genuine model uncertainty that should be reported.

## Modular Mortality Functions

For applications outside the Monte Carlo framework — quick calculations,
integration into custom analyses, or pedagogical exploration —
vitalBayes exports the individual mortality model functions:

``` r
ages <- seq(0.5, 30, by = 0.5)

# Chen-Watanabe with VB (k serves both roles)
M_cw_vb <- M_chen_watanabe_L0(
  age = ages, Linf = 100, L0 = 25, k = 0.08,
  two_phase = TRUE, tmat = 12, late_model = "gompertz"
)

# Chen-Watanabe with Gompertz (separated k roles)
M_cw_gomp <- M_chen_watanabe_L0(
  age = ages, Linf = 100, L0 = 25,
  k = 0.08,              # VB-equivalent k (for M_inf)
  k_native = 0.15,       # Gompertz k (for L(t))
  growth_model = "gompertz",
  two_phase = TRUE, tmat = 12
)

# Peterson-Wroblewski (native k only, no VB theory)
lw_func <- function(L) 0.0001 * L^3.1
M_pw <- M_peterson_wroblewski(
  age = ages, Linf = 100, L0 = 25, k = 0.15,
  lw_fun = lw_func, growth_model = "gompertz"
)

# Lorenzen growth-based (both k values)
M_lor <- M_lorenzen(
  age = ages, Linf = 100, L0 = 25,
  k = 0.15,          # Gompertz k (for L(t)/Linf ratio)
  k_vb = 0.08,       # VB-equivalent (for ln(k) regression term)
  weight_based = FALSE,
  growth_model = "gompertz"
)

# Lorenzen weight-based (native k only)
M_lor_wt <- M_lorenzen(
  age = ages, Linf = 100, L0 = 25, k = 0.15,
  lw_fun = lw_func, weight_based = TRUE,
  growth_model = "gompertz"
)

# Scale any of them
M_scaled <- scale_mortality(M_cw_gomp, M_target = NULL, tmax = 35, p = 0.001)
```

## Downstream Integration: Survival Simulation

Mortality schedules from
[`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md)
feed directly into
[`simulate_survivorship()`](https://brian-j-moe.github.io/vitalBayes/reference/simulate_survivorship.md)
for cohort survival analysis:

``` r
surv <- simulate_survivorship(
  mc_object = mort,
  n         = 50000,    # Starting cohort size
  n_iter    = 2000,     # Simulation iterations
  mode      = "random"  # Sample parameter sets randomly
)

# Key survival metrics
surv$Aggregate$Age_of_Death
#>     mean    sd  lower  upper
#>    <num> <num>  <num>  <num>
#> 1:  8.42  2.31   4.89  13.21

surv$Aggregate$Survival_to_tmat
#>     mean      sd   lower   upper
#>    <num>   <num>   <num>   <num>
#> 1: 0.312  0.0891   0.168   0.512
```

See
[`vignette("survivorship_simulation")`](https://brian-j-moe.github.io/vitalBayes/articles/survivorship_simulation.md)
for comprehensive coverage of survival analysis options, visualization,
and interpretation.

## Troubleshooting

| Issue                       | Possible Cause                               | Solution                                                                                |
|-----------------------------|----------------------------------------------|-----------------------------------------------------------------------------------------|
| Negative mortality values   | CW two-phase Taylor expansion breakdown      | Function automatically caps at 0; consider single-phase or different `late_model`       |
| Very wide uncertainty bands | Weak growth model posterior                  | Check growth model convergence and parameter correlations                               |
| Unrealistic tmax estimates  | `Linf_factor` too extreme                    | Reduce from 0.99 to 0.95 for exploratory analysis                                       |
| Missing tmat error          | CW needs age-at-maturity                     | Provide `maturity_fit`, manual `tmat`, or both `Lmat` and `tmat` for bivariate sampling |
| `lw_fun` error              | PW/Lorenzen weight-based needs length-weight | Provide function: `lw_fun = function(L) a * L^b`                                        |
| CW and PW disagree strongly | Different theoretical assumptions            | Report both; consider which assumptions best match your species                         |

## Reporting Results

When reporting mortality estimates in publications, include:

1.  **Model specification**: Which mortality model(s), growth model
    used, scaling approach, and parameter sources.
2.  **Growth parameter source**: vitalBayes fit (specify model: VB,
    Gompertz, or Logistic) or literature values.
3.  **Key metrics**: Mean \\M\\, age-specific \\M\\ at key life stages
    (e.g., neonate, sub-adult, adult), \\t\_{max}\\ estimate.
4.  **Uncertainty**: 95% credible intervals throughout, noting whether
    these reflect joint posterior sampling or independent parameter
    variation.

Example methods text:

> “Natural mortality was estimated using the Chen-Watanabe model (Chen &
> Watanabe, 1989) with Gompertz senescence, implemented via the
> vitalBayes package (v1.0.0; Author, Year). Growth parameters were
> drawn from the joint posterior of a Gompertz growth model (selected by
> LOO-CV). Predicted length-at-age used the native Gompertz trajectory,
> while the asymptotic mortality rate (\\M\_\infty = k\_{VB}\\) was
> computed from a VB-equivalent growth coefficient derived from
> biological milestones \\(L\_\infty, L_0, L\_{mat}, t\_{mat})\\.
> Mortality schedules were scaled to 0.1% survival at estimated
> \\t\_{max}\\ and propagated through 2,000 Monte Carlo simulations.
> Median age-specific mortality with 95% credible intervals is
> reported.”

## References

Chen, S., & Watanabe, S. (1989). Age dependence of natural mortality
coefficient in fish population dynamics. *Nippon Suisan Gakkaishi*,
55(2), 205–208.

Lorenzen, K. (1996). The relationship between body weight and natural
mortality in juvenile and adult fish. *Journal of Fish Biology*, 49(4),
627–642.

Lorenzen, K. (2022). Size- and age-dependent natural mortality in fish
populations. *Fisheries Research*, 255, 106454.

Peterson, I., & Wroblewski, J. S. (1984). Mortality rate of fishes in
the pelagic ecosystem. *Canadian Journal of Fisheries and Aquatic
Sciences*, 41(7), 1117–1120.

Then, A. Y., Hoenig, J. M., Hall, N. G., & Hewitt, D. A. (2015).
Evaluating the predictive performance of empirical estimators of natural
mortality rate using information on over 200 fish species. *ICES Journal
of Marine Science*, 72(1), 82–92.

## See Also

- `vignette("chen_watanabe_reparameterization")` — Mathematical
  derivation of the CW reparameterization and \\G(t)\\ framework
- [`vignette("survivorship_simulation")`](https://brian-j-moe.github.io/vitalBayes/articles/survivorship_simulation.md)
  — Cohort survival analysis using mortality schedules
- [`vignette("fit_bayesian_growth")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_growth.md)
  — Growth model fitting that feeds the mortality workflow
- [`vignette("fit_bayesian_maturity")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_maturity.md)
  — Maturity estimation for \\t\_{mat}\\ and \\L\_{mat}\\

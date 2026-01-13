# Stochastic Estimation of Age-Specific Natural Mortality

Monte Carlo simulation of age-specific natural mortality schedules with
full uncertainty propagation from growth model posteriors. Supports
Chen-Watanabe, Peterson-Wroblewski, and Lorenzen models with automatic
derivation of VB-equivalent \\k\\ from any growth model fit.

## Usage

``` r
get_stochastic_mortality(
  method = c("CW", "PW", "L"),
  growth_fit = NULL,
  maturity_fit = NULL,
  sex = NULL,
  Linf = NULL,
  L0 = NULL,
  k = NULL,
  tmat = NULL,
  Linf_factor = 0.99,
  age_seq = function(tmax) seq(0.1, ceiling(tmax), length.out = 500),
  iter = 2000,
  scaled = TRUE,
  M_target = NULL,
  p = 0.001,
  two_phase = TRUE,
  late_model = c("gompertz", "logistic"),
  tm_factor = 2/3,
  M_mult = 2,
  smooth_factor = 1/3,
  lw_fun = NULL,
  weight_based = FALSE,
  growth_model = c("vb", "gompertz", "logistic"),
  seed = 1234,
  palette = c("synthwave", "viridis", "okabe", "plasma", "inferno"),
  print_plot = TRUE,
  show_progress = TRUE
)
```

## Arguments

- method:

  Character. Mortality model: `"CW"`, `"PW"`, or `"L"`.

- growth_fit:

  Optional `CmdStanMCMC` object from
  [`fit_bayesian_growth`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md).
  If provided, parameters are extracted from the joint posterior.

- maturity_fit:

  Optional `CmdStanMCMC` object from
  [`fit_bayesian_maturity`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md)
  providing age-at-maturity for k-based growth fits or two-phase CW.

- sex:

  Integer. Sex code (1 = female, 2 = male) for hierarchical models.

- Linf, L0, k, tmat:

  Alternative to `growth_fit`: specify parameters directly as
  `c(mean, sd)` vectors for independent normal sampling.

- Linf_factor:

  Numeric in (0, 1). Fraction of \\L\_\infty\\ for \\t\_{max}\\
  estimation. Default 0.99.

- age_seq:

  Function or numeric vector defining ages for mortality calculation.
  Default `function(tmax) seq(0, ceiling(tmax), length.out = 500)`.

- iter:

  Number of Monte Carlo iterations. Default 2000.

- scaled:

  Logical. If `TRUE` (default), scales mortality to `M_target` or
  survival probability `p`.

- M_target:

  Target mean mortality. Can be numeric scalar, function of tmax, or
  `NULL` for survival-probability-based scaling.

- p:

  Survival probability to \\t\_{max}\\ for scaling. Default 0.001.

- two_phase:

  Logical. For CW model, use two-phase senescence? Default `TRUE`.

- late_model:

  Character. Senescence model: `"gompertz"` or `"logistic"`. Default
  `"gompertz"`.

- tm_factor, M_mult, smooth_factor:

  Two-phase model parameters.

- lw_fun:

  Length-weight function for PW and weight-based Lorenzen.

- weight_based:

  Logical. For Lorenzen, use weight-based formulation? Default `FALSE`.

- growth_model:

  Character. Growth model type when using manual parameters: `"vb"`,
  `"gompertz"`, or `"logistic"`.

- seed:

  Random seed for reproducibility. Default 1234.

- palette:

  Color palette for plot: `"synthwave"`, `"viridis"`, `"okabe"`,
  `"plasma"`, or `"inferno"`.

- print_plot:

  Logical. Print plot on completion? Default `TRUE`.

- show_progress:

  Logical. Show progress messages? Default `TRUE`.

## Value

A list with components: Schedules (data.table of all mortality schedules
with columns set_id, age, M, M_scaled), Parameters (data.table of
sampled life history parameters), Summary (data.table with median and 95
percent CI by age), and Plot (ggplot2 object).

## Details

A key feature of this function is growth-model-agnostic mortality
estimation. When a growth fit from
[`fit_bayesian_growth`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md)
is provided, the function extracts biological milestones (Linf, L0,
Lmat, tmat) and computes the VB-equivalent k needed for Chen-Watanabe
and growth-based Lorenzen models.

This allows users to fit whichever growth model (von Bertalanffy,
Gompertz, or Logistic) best describes their data, then estimate
mortality without theoretical compromise.

When `growth_fit` is provided, parameters are drawn from the joint
posterior distribution, preserving correlations. This yields mortality
estimates with appropriate (often narrower) uncertainty bounds compared
to independent sampling of each parameter.

Three mortality models are available: CW (Chen-Watanabe 1989 with L0
parameterization and optional two-phase senescence), PW
(Peterson-Wroblewski 1984 weight-based allometric model), and L
(Lorenzen 1996/2022 in weight-based or growth-based form).

## Examples

``` r
if (FALSE) { # \dontrun{
# From a Gompertz growth fit (maturity-based parameterization)
mort <- get_stochastic_mortality(
  method     = "CW",
  growth_fit = gomp_fit,  # Any growth model works!
  sex        = 1,
  iter       = 2000,
  scaled     = TRUE,
  p          = 0.001
)

# View plot
mort$Plot

# Check VB-equivalent k distribution
hist(mort$Parameters$k_vb_equiv)

# Manual specification for sensitivity analysis
mort <- get_stochastic_mortality(
  method = "CW",
  Linf = c(100, 5),
  L0   = c(25, 2),
  k    = c(0.1, 0.02),  # VB k or VB-equivalent k
  tmat = c(10, 1)
)
} # }
```

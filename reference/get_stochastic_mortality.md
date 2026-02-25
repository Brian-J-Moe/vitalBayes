# Stochastic Estimation of Age-Specific Natural Mortality

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
  Lmat = NULL,
  tmat = NULL,
  rho_Lmat_tmat = 0.5,
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
  If provided, growth model type is auto-detected and parameters are
  extracted from the joint posterior.

- maturity_fit:

  Optional `CmdStanMCMC` object from
  [`fit_bayesian_maturity`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md)
  providing age-at-maturity.

- sex:

  Integer. Sex code (1 = female, 2 = male) for hierarchical models.

- Linf, L0:

  Required for manual mode: `c(mean, sd)` vectors.

- k:

  Optional `c(mean, sd)` for native growth coefficient (manual mode).
  Required for PW and weight-based Lorenzen. Optional for CW and
  growth-based Lorenzen when `Lmat` and `tmat` are provided (the
  VB-equivalent will be derived from milestones).

- Lmat, tmat:

  Optional `c(mean, sd)` vectors for maturity milestones. When both are
  provided, enables bivariate sampling and milestone-based
  \\k\_{VB}^{equiv}\\ derivation.

- rho_Lmat_tmat:

  Numeric in (-1, 1). Correlation between \\L\_{mat}\\ and \\t\_{mat}\\
  for bivariate sampling. Default 0.5.

- Linf_factor:

  Numeric in (0, 1). Fraction of \\L\_\infty\\ for \\t\_{max}\\
  estimation. Default 0.99.

- age_seq:

  Function or numeric vector defining age grid. Default
  `function(tmax) seq(0.1, ceiling(tmax), length.out = 500)`.

- iter:

  Number of Monte Carlo iterations. Default 2000.

- scaled:

  Logical. Scale mortality to target? Default `TRUE`.

- M_target:

  Target mean mortality (numeric, function of tmax, or `NULL` for
  survival-based scaling).

- p:

  Survival probability to \\t\_{max}\\. Default 0.001.

- two_phase:

  Logical. For CW, use two-phase senescence? Default `TRUE`.

- late_model:

  Character. Senescence model: `"gompertz"` or `"logistic"`.

- tm_factor, M_mult, smooth_factor:

  Two-phase model parameters.

- lw_fun:

  Length-weight function for PW and weight-based Lorenzen.

- weight_based:

  Logical. For Lorenzen, use weight-based? Default `FALSE`.

- growth_model:

  Character. Growth model type for manual specification: `"vb"`,
  `"gompertz"`, or `"logistic"`. Ignored when `growth_fit` is provided
  (auto-detected).

- seed:

  Random seed. Default 1234.

- palette:

  Color palette for plot.

- print_plot:

  Logical. Print plot on completion? Default `TRUE`.

- show_progress:

  Logical. Show messages? Default `TRUE`.

## Value

A list with components:

- Schedules:

  data.table of all mortality schedules (set_id, age, M_raw, M_scaled)

- Parameters:

  data.table of sampled life history parameters including k_original
  (native), k_vb_equiv (VB-equivalent), and growth_model

- Summary:

  data.table with median and 95\\ Plotggplot2 object

Monte Carlo simulation of age-specific natural mortality schedules with
full uncertainty propagation from growth model posteriors. Supports
Chen-Watanabe, Peterson-Wroblewski, and Lorenzen models with automatic
derivation of VB-equivalent \\k\\ from any growth model fit. A key
feature of this function is growth-model-agnostic mortality estimation.
When a growth fit from
[`fit_bayesian_growth`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md)
is provided, the function automatically detects the growth model type
and extracts posterior draws of the native growth coefficient and
biological milestones. All mortality models use the native growth
trajectory for body-size predictions, while CW and growth-based Lorenzen
additionally receive the VB-equivalent \\k\\ for their model-specific
parameters that were derived or calibrated under VB assumptions. Growth
Coefficient RolesThe three mortality models use \\k\\ in fundamentally
different ways:

- CW:

  \\M(t) = k\_{VB} \cdot L\_\infty / L(t)\\: the \\k\_{VB}\\ appears as
  the asymptotic mortality rate constant (derived under VB assumptions),
  while \\L(t)\\ is predicted by the native growth model.

- Lorenzen growth-based:

  \\\ln M = a_0 + a_1 \ln(L/L\_\infty) + a_2 \ln(k\_{VB})\\: the
  \\k\_{VB}\\ enters as a calibration coefficient (fitted against VB
  parameters), while \\L(t)/L\_\infty\\ uses the native trajectory.

- PW and Lorenzen weight-based:

  Operate entirely on predicted body weight: \\M(W(t))\\. No \\k\\
  appears in the mortality equation itself; only the native growth
  trajectory matters.

Manual Parameter SpecificationWhen a growth fit is not available,
parameters can be specified manually. Two sampling modes are supported:

- Bivariate (preferred):

  When both `Lmat` and `tmat` are provided, they are sampled jointly
  from a bivariate normal with correlation `rho_Lmat_tmat`. The
  VB-equivalent \\k\\ is then derived from the sampled milestones.

- Independent (fallback):

  When only `k` is provided, it is sampled independently from a normal
  distribution and assumed to be VB-equivalent for CW and growth-based
  Lorenzen. For non-VB growth models, the derivative-matching approach
  provides a VB-equivalent.

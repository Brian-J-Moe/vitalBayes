# Stochastic Estimation of Age-Specific Natural Mortality

## Usage

    get_stochastic_mortality(
      method = c("CW", "PW", "L"),
      growth_fit = NULL,
      birth_fit = NULL,
      length_maturity_fit = NULL,
      age_maturity_fit = NULL,
      sex = 1L,
      Linf = NULL,
      L0 = NULL,
      Lmat = NULL,
      tmat = NULL,
      maturity_cor = 0.5,
    <<<<<<< HEAD
      growth_model = c("vb", "gompertz", "logistic"),
      Linf_factor = 0.99,
      age_seq = function(tmax) seq(0.1, ceiling(tmax), length.out = 500),
      iter = 2000,
      scaled = TRUE,
      M_target = NULL,
      p = 0.001,
      two_phase = FALSE,
      late_model = c("gompertz", "logistic"),
      tm_factor = 2/3,
      M_mult = 2,
      smooth_factor = 1/3,
      weight_based = FALSE,
      lw_fun = NULL,
      seed = 1234,
      palette = c("synthwave", "viridis", "okabe", "plasma", "inferno"),
      print_plot = TRUE,
      show_progress = TRUE
    )

    get_stochastic_mortality(
      method = c("CW", "PW", "L"),
      growth_fit = NULL,
      birth_fit = NULL,
      length_maturity_fit = NULL,
      age_maturity_fit = NULL,
      sex = 1L,
      Linf = NULL,
      L0 = NULL,
      Lmat = NULL,
      tmat = NULL,
      maturity_cor = 0.5,
    =======
    >>>>>>> 89361dd64d8e03d181b8dff5161237d2b73d660d
      growth_model = c("vb", "gompertz", "logistic"),
      Linf_factor = 0.99,
      age_seq = function(tmax) seq(0.1, ceiling(tmax), length.out = 500),
      iter = 2000,
      scaled = TRUE,
      M_target = NULL,
      p = 0.001,
      two_phase = FALSE,
      late_model = c("gompertz", "logistic"),
      tm_factor = 2/3,
      M_mult = 2,
      smooth_factor = 1/3,
      weight_based = FALSE,
      lw_fun = NULL,
      seed = 1234,
      palette = c("synthwave", "viridis", "okabe", "plasma", "inferno"),
      print_plot = TRUE,
      show_progress = TRUE
    )

## Arguments

- method:

  Character. Mortality model: `"CW"`, `"PW"`, or `"L"`.

- growth_fit:

  CmdStanMCMC object from `fit_bayesian_growth`. If provided, extracts
  Linf, L0, and (for maturity-based fits) Lmat, tmat with preserved
  correlations. Default NULL.

- birth_fit:

  CmdStanMCMC object from `fit_bayesian_birth`. Provides L0 (birth size)
  if not available from growth_fit. Default NULL.

- length_maturity_fit:

  CmdStanMCMC object from `fit_bayesian_maturity` for
  length-at-maturity. Provides Lmat if not in growth_fit. Default NULL.

- age_maturity_fit:

  CmdStanMCMC object from `fit_bayesian_maturity` for age-at-maturity.
  Provides tmat if not in growth_fit. Default NULL.

- sex:

  Integer. Sex code (1 = female, 2 = male) for extracting from two-sex
  model fits. Default 1.

\<\<\<\<\<\<\< HEAD

- Linf, L0:

  Required for manual mode: `c(mean, sd)` vectors.

- Lmat, tmat:

  Optional `c(mean, sd)` vectors for maturity milestones. When both are
  provided, enables bivariate sampling and milestone-based
  \\k\_{VB}^{equiv}\\ derivation.

=======

- Linf, L0, Lmat, tmat:

  Numeric vectors of length 2: `c(mean, sd)`. Used when corresponding
  stanfit is not provided.

\>\>\>\>\>\>\> 89361dd64d8e03d181b8dff5161237d2b73d660d

- maturity_cor:

  Numeric. Correlation between Lmat and tmat when both come from manual
  specification or separate fits. Default 0.5 (positive correlation is
  biologically reasonable). Set to 0 for independent sampling.

- growth_model:

  Character. Growth model for G(t) and L(t): `"vb"`, `"gompertz"`, or
  `"logistic"`. Default `"vb"`.

- Linf_factor:

  Fraction of Linf used to estimate tmax. Default 0.99.

- age_seq:

  Function or numeric vector for age grid.

- iter:

  Number of Monte Carlo iterations. Default 2000.

- scaled:

  Logical. Scale mortality to M_target? Default TRUE.

- M_target:

  Target mean mortality. Can be scalar, function of tmax, or NULL.

- p:

  Survival probability for scaling when M_target is NULL. Default 0.001.

- two_phase:

  Logical. Use CW two-phase senescence? Default FALSE.

- late_model:

  Character. Senescence model. Default `"gompertz"`.

\<\<\<\<\<\<\< HEAD

- tm_factor, M_mult, smooth_factor:

  Two-phase model parameters.

=======

- tm_factor:

  Numeric. Transition age as fraction of tmat. Default 2/3.

\>\>\>\>\>\>\> 89361dd64d8e03d181b8dff5161237d2b73d660d

- M_mult:

  Numeric. Mortality multiplier for senescence. Default 2.

- smooth_factor:

  Numeric. Transition smoothness. Default 1/3.

- weight_based:

  Logical. For Lorenzen: use weight-based? Default FALSE.

- lw_fun:

  Function. Length-weight relationship for PW and weight-based Lorenzen.

- seed:

  Integer. Random seed. Default 1234.

- palette:

  Character. Color palette for plot.

- print_plot:

  Logical. Print the plot? Default TRUE.

- show_progress:

  Logical. Show progress messages? Default TRUE.

\<\<\<\<\<\<\< HEAD

- maturity_fit:

  Optional `CmdStanMCMC` object from
  [`fit_bayesian_maturity`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md)
  providing age-at-maturity.

- k:

  Optional `c(mean, sd)` for native growth coefficient (manual mode).
  Required for PW and weight-based Lorenzen. Optional for CW and
  growth-based Lorenzen when `Lmat` and `tmat` are provided (the
  VB-equivalent will be derived from milestones).

- rho_Lmat_tmat:

  Numeric in (-1, 1). Correlation between \\L\_{mat}\\ and \\t\_{mat}\\
  for bivariate sampling. Default 0.5.

- Linf, L0, Lmat, tmat:

  Numeric vectors of length 2: `c(mean, sd)`. Used when corresponding
  stanfit is not provided.

- tm_factor:

  Numeric. Transition age as fraction of tmat. Default 2/3.

## Value

A list with components:

- Schedules:

  data.table of all mortality schedules (set_id, age, M_raw, M_scaled)

- Parameters:

  data.table of sampled life history parameters including k_original
  (native), k_vb_equiv (VB-equivalent), and growth_model

- Summary:

  data.table with median and 95\\ Plotggplot2 object

A list with components: \>\>\>\>\>\>\>
89361dd64d8e03d181b8dff5161237d2b73d660d

- Schedules:

  data.table with columns: set_id, age, M

- Parameters:

  data.table with life history parameters and derived quantities

- Summary:

  data.table with age-wise summary statistics

- Plot:

  ggplot2 object

Monte Carlo simulation of age-specific natural mortality schedules with
\<\<\<\<\<\<\< HEAD full uncertainty propagation from growth model
posteriors. Supports Chen-Watanabe, Peterson-Wroblewski, and Lorenzen
models with automatic derivation of VB-equivalent \\k\\ from any growth
model fit.Monte Carlo simulation of age-specific natural mortality
schedules with ======= \>\>\>\>\>\>\>
89361dd64d8e03d181b8dff5161237d2b73d660d uncertainty propagation from
growth and maturity parameters. Accepts either vitalBayes stanfit
objects (preserving posterior correlations) or manual `c(mean, sd)`
specifications. \<\<\<\<\<\<\< HEAD A key feature of this function is
growth-model-agnostic mortality estimation. When a growth fit from
[`fit_bayesian_growth`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md)
is provided, the function automatically detects the growth model type
and extracts posterior draws of the native growth coefficient and
biological milestones. All mortality models use the native growth
trajectory for body-size predictions, while CW and growth-based Lorenzen
additionally receive the VB-equivalent \\k\\ for their model-specific
parameters that were derived or calibrated under VB assumptions.This
function samples life history parameters and computes mortality
schedules using the chosen method. When vitalBayes fit objects are
provided, posterior correlations between parameters are preserved,
yielding more realistic uncertainty bounds than independent
sampling.======= This function samples life history parameters and
computes mortality schedules using the chosen method. When vitalBayes
fit objects are provided, posterior correlations between parameters are
preserved, yielding more realistic uncertainty bounds than independent
sampling.\>\>\>\>\>\>\> 89361dd64d8e03d181b8dff5161237d2b73d660d When
maturity parameters (Lmat, tmat) come from separate model fits, the
`maturity_cor` argument allows specification of assumed correlation
between these parameters. Larger individuals typically mature at older
ages, so a positive correlation (default 0.5) is biologically
reasonable. \<\<\<\<\<\<\< HEAD Growth Coefficient RolesThe three
mortality models use \\k\\ in fundamentally different ways:

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

Parameter SourcesParameters can be supplied from multiple sources with
the following priority:

- Linf:

  growth_fit \> manual

- L0:

  growth_fit \> birth_fit \> manual

- Lmat:

  growth_fit \> length_maturity_fit \> manual

- tmat:

  growth_fit \> age_maturity_fit \> manual

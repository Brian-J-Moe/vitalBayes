# Stochastic Estimation of Age-Specific Natural Mortality

Monte Carlo simulation of age-specific natural mortality schedules with
uncertainty propagation from growth and maturity parameters. Accepts
either vitalBayes stanfit objects (preserving posterior correlations) or
manual `c(mean, sd)` specifications.

## Usage

``` r
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
```

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

- Linf, L0, Lmat, tmat:

  Numeric vectors of length 2: `c(mean, sd)`. Used when corresponding
  stanfit is not provided.

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

- tm_factor:

  Numeric. Transition age as fraction of tmat. Default 2/3.

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

## Value

A list with components:

- Schedules:

  data.table with columns: set_id, age, M

- Parameters:

  data.table with life history parameters and derived quantities

- Summary:

  data.table with age-wise summary statistics

- Plot:

  ggplot2 object

## Details

This function samples life history parameters and computes mortality
schedules using the chosen method. When vitalBayes fit objects are
provided, posterior correlations between parameters are preserved,
yielding more realistic uncertainty bounds than independent sampling.

When maturity parameters (Lmat, tmat) come from separate model fits, the
`maturity_cor` argument allows specification of assumed correlation
between these parameters. Larger individuals typically mature at older
ages, so a positive correlation (default 0.5) is biologically
reasonable.

## Parameter Sources

Parameters can be supplied from multiple sources with the following
priority:

- Linf:

  growth_fit \> manual

- L0:

  growth_fit \> birth_fit \> manual

- Lmat:

  growth_fit \> length_maturity_fit \> manual

- tmat:

  growth_fit \> age_maturity_fit \> manual

## Examples

``` r
if (FALSE) { # \dontrun{
# Using vitalBayes fits (full correlation preservation)
mort <- get_stochastic_mortality(
  method = "CW",
  growth_fit = growth_fit,  # From maturity-based fit_bayesian_growth
  sex = 1,
  growth_model = "vb"
)

# Using separate maturity fits with specified correlation
mort <- get_stochastic_mortality(
  method = "CW",
  growth_fit = growth_fit,           # For Linf, L0
  length_maturity_fit = L50_fit,     # For Lmat
  age_maturity_fit = t50_fit,        # For tmat
  maturity_cor = 0.6,                # Assumed Lmat-tmat correlation
  sex = 1
)

# Using manual parameters with correlation
mort <- get_stochastic_mortality(
  method = "CW",
  Linf = c(126, 10),
  L0 = c(35, 3),
  Lmat = c(83, 5),
  tmat = c(47, 4),
  maturity_cor = 0.5,
  growth_model = "gompertz"
)
} # }
```

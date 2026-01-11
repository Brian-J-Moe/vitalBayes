# Stochastic Simulation of Age-Specific Natural Mortality

Monte Carlo simulation of age-specific natural mortality schedules under
parameter uncertainty. Life-history parameters can be drawn from
vitalBayes model posteriors (preserving correlations) or specified
manually.

## Usage

``` r
get_stochastic_mortality(
  method = c("CW", "PW", "L"),
  growth_fit = NULL,
  maturity_fit = NULL,
  sex = NULL,
  Linf = NULL,
  k = NULL,
  t0 = NULL,
  tmat = NULL,
  Linf_factor = 0.999,
  age_seq = function(tmax) seq(0, ceiling(tmax), length.out = 2000),
  iter = 2000,
  scaled = TRUE,
  M_target = NULL,
  p = 0.001,
  two_phase = TRUE,
  late_model = "gompertz",
  tm_factor = 2/3,
  M_mult = 2,
  M_cap_factor = 1/4,
  smooth_factor = 1/3,
  mode = "K_mult",
  alpha = 3,
  r_given = NULL,
  weight_based = FALSE,
  lw_fun = NULL,
  seed = 1234,
  palette = c("synthwave", "viridis", "okabe", "plasma", "inferno"),
  fill_color = NULL,
  line_color = NULL,
  ribbon_alpha = 0.6,
  linewidth = 1,
  plot_round = 2,
  xlab = "Age (yrs)",
  ylab = "Instantaneous Mortality (M)",
  title = "Natural Mortality Rates",
  subtitle = "Median with 95% credible interval",
  base_size = 11,
  font_family = "serif",
  print_plot = TRUE,
  show_progress = TRUE
)
```

## Arguments

- method:

  Character. One of `"CW"`, `"PW"`, `"L"`.

- growth_fit:

  Optional CmdStanMCMC object from
  [`fit_bayesian_growth()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md).
  If provided, life-history parameters are drawn from the joint
  posterior.

- maturity_fit:

  Optional CmdStanMCMC object from
  [`fit_bayesian_maturity()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md).
  Provides age-at-maturity (tmat) for CW model.

- sex:

  Integer. Sex code (1 = female, 2 = male) for hierarchical models.
  Required if fits are hierarchical.

- Linf, k, t0, tmat:

  Manual parameter specification as `c(mean, sd)`. Ignored if
  corresponding fit objects are provided.

- Linf_factor:

  Proportion of asymptotic length for tmax estimation.
  `tmax = -(1/k) * log(1 - Linf_factor) + t0`. Default 0.999.

- age_seq:

  Function returning evaluation ages given tmax. Default:
  `function(tmax) seq(0, ceiling(tmax), length.out = 2000)`.

- iter:

  Integer. Number of Monte Carlo iterations. Default 2000.

- scaled:

  Logical. Scale schedules to target mortality? Default TRUE.

- M_target:

  Target mean mortality. Can be NULL (derive from p), numeric, or
  function of tmax (e.g., Hoenig). Default NULL.

- p:

  Probability of survival to tmax (if M_target = NULL). Default 0.001.

- two_phase:

  Logical. Use two-phase CW model? Default TRUE.

- late_model:

  CW late-phase: `"gompertz"` or `"logistic"`.

- tm_factor, M_mult, M_cap_factor, smooth_factor:

  CW late-phase parameters.

- mode, alpha, r_given:

  CW logistic late-phase parameters.

- weight_based:

  For Lorenzen: TRUE = weight-based, FALSE = growth-based.

- lw_fun:

  Length-weight function `lw_fun(length)` returning grams. Required for
  PW and Lorenzen weight-based.

- seed:

  Random seed. Default 1234.

- palette:

  Color palette: `"synthwave"`, `"viridis"`, `"okabe"`.

- fill_color, line_color:

  Override specific colors. If NULL, from palette.

- ribbon_alpha:

  Ribbon transparency. Default 0.6.

- linewidth:

  Line width. Default 1.

- plot_round:

  Decimal places for age rounding in plot summary. Default 2.

- xlab, ylab, title, subtitle:

  Plot labels.

- base_size:

  Base font size. Default 11.

- font_family:

  Font family. Default "serif".

- print_plot:

  Print plot on completion? Default TRUE.

- show_progress:

  Show progress messages? Default TRUE.

## Value

A list with:

- Schedules:

  data.table of all mortality schedules (set_id, age, M)

- Parameters:

  data.table of sampled life-history parameters

- Summary:

  data.table with mean, 2.5 percent, 97.5 percent by age

- Plot:

  ggplot2 object

## Details

Three mortality models are available:

- **CW** - Chen & Watanabe (1989), with optional late-life senescence

- **PW** - Peterson & Wroblewski (1984), weight-based

- **L** - Lorenzen (1996/2022), weight- or growth-based

Schedules can be scaled to match a target mean mortality derived from:
empirical relationships (e.g., Hoenig 1983, Then et al. 2015), a fixed
value, or probability of survival to maximum age.

## Examples

``` r
if (FALSE) { # \dontrun{
# Using vitalBayes fits
mort <- get_stochastic_mortality(
  method = "CW",
  growth_fit = growth_fit,
  maturity_fit = maturity_fit,
  sex = 2,
  scaled = TRUE,
  M_target = function(tmax) 4.899 * tmax^(-0.916)
)

# Manual specification
mort <- get_stochastic_mortality(
  method = "L",
  Linf = c(484, 10),
  k = c(0.17, 0.02),
  t0 = c(-0.97, 0.05),
  weight_based = FALSE,
  palette = "okabe"
)
} # }
```

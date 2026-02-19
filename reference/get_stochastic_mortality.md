# Stochastic Estimation of Age-Specific Natural Mortality

Monte Carlo simulation of age-specific natural mortality schedules with
uncertainty propagation from growth and maturity parameters.

## Usage

``` r
get_stochastic_mortality(
  method = c("CW", "PW", "L"),
  Linf,
  L0,
  Lmat,
  tmat,
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

- Linf, L0, Lmat, tmat:

  Numeric vectors of length 2: `c(mean, sd)`.

- growth_model:

  Character. Growth model for G(t) and L(t): `"vb"`, `"gompertz"`, or
  `"logistic"`.

- Linf_factor:

  Fraction of Linf used to estimate tmax. Default 0.99.

- age_seq:

  Function or numeric vector for age grid. Default creates 500 points
  from 0.1 to tmax.

- iter:

  Number of Monte Carlo iterations. Default 2000.

- scaled:

  Scale mortality to M_target? Default TRUE.

- M_target:

  Target mean mortality. Can be scalar, function of tmax, or NULL.

- p:

  Survival probability for scaling if M_target NULL. Default 0.001.

- two_phase:

  Use CW two-phase senescence? Default FALSE.

- late_model:

  Senescence model: `"gompertz"` or `"logistic"`.

- tm_factor:

  Transition age factor. Default 2/3.

- M_mult:

  Mortality multiplier for senescence. Default 2.

- smooth_factor:

  Transition smoothness. Default 1/3.

- weight_based:

  For Lorenzen: use weight-based formulation? Default FALSE.

- lw_fun:

  Length-weight function (required for PW and weight-based Lorenzen).

- seed:

  Random seed. Default 1234.

- palette:

  Color palette: `"synthwave"`, `"viridis"`, `"okabe"`, `"plasma"`, or
  `"inferno"`.

- print_plot:

  Print plot? Default TRUE.

- show_progress:

  Show progress messages? Default TRUE.

## Value

A list with:

- Schedules:

  data.table with columns: set_id, age, M

- Parameters:

  data.table with columns: set_id, Linf, L0, Lmat, tmat, Minf, k_native,
  tmax

- Summary:

  data.table with age-wise median and 95% CI

- Plot:

  ggplot2 object

## Details

This function samples life history parameters from specified
distributions and computes mortality schedules using the chosen method.
The output format is compatible with
[`simulate_survivorship`](https://brian-j-moe.github.io/vitalBayes/reference/simulate_survivorship.md).

Three mortality models are available:

- **CW**: Chen-Watanabe (1989) with model-dependent G(t)

- **PW**: Peterson-Wroblewski (1984) weight-based

- **L**: Lorenzen (1996/2022) weight- or growth-based

## Examples

``` r
if (FALSE) { # \dontrun{
mort <- get_stochastic_mortality(
  method = "CW",
  Linf = c(108, 10),
  L0 = c(35, 2),
  Lmat = c(83, 5),
  tmat = c(47, 3),
  growth_model = "gompertz",
  iter = 2000
)
mort$Plot
} # }
```

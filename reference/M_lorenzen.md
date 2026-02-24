# Lorenzen Natural Mortality Model

Computes size-dependent natural mortality following Lorenzen (1996,
2022). Supports both weight-based and growth-based formulations.

## Usage

``` r
M_lorenzen(
  age,
  Linf,
  L0,
  Lmat,
  tmat,
  lw_fun = NULL,
  weight_based = FALSE,
  growth_model = c("vb", "gompertz", "logistic"),
  sample_params = TRUE
)

M_lorenzen(
  age,
  Linf,
  L0,
  Lmat,
  tmat,
  lw_fun = NULL,
  weight_based = FALSE,
  growth_model = c("vb", "gompertz", "logistic"),
  sample_params = TRUE
)
```

## Arguments

- age:

  Numeric vector of ages.

- Linf:

  Asymptotic length.

- L0:

  Length at birth.

- Lmat:

  Length at maturity.

- tmat:

  Age at maturity.

- lw_fun:

  Length-weight function (required if weight_based = TRUE).

- weight_based:

  Use weight-based formulation?

- growth_model:

  Growth model for L(t).

- sample_params:

  Sample parameters from their distributions?

- k:

  Growth coefficient native to the fitted model (used for \\L(t)\\
  prediction via `growth_model`).

- k_vb:

  VB-equivalent growth coefficient for the growth-based formulation's
  \\\ln(k)\\ term. Defaults to `k`, which is correct for VB fits. For
  non-VB fits, use `compute_k_vb_equivalent` or the derivative-matching
  approach. Ignored for weight-based formulation.

## Value

Numeric vector of instantaneous mortality rates.

Numeric vector of instantaneous mortality rates.

## Details

Two formulations are available:

**Weight-based (Lorenzen 1996):** \$\$M(W) = \alpha \cdot W^{\beta}\$\$
where \\\alpha \sim N(3.69, 0.502)\\ and \\\beta \sim N(-0.305,
0.029)\\. Body weight \\W(t)\\ is computed from predicted length via a
user-supplied length-weight function using the native growth model
trajectory.

**Growth-based (Lorenzen 2022):** \$\$\ln M = 0.28 - 1.30
\ln(L(t)/L\_\infty) + 1.08 \ln(k\_{VB})\$\$ This formulation was
calibrated using von Bertalanffy parameters. The \\L(t)/L\_\infty\\
ratio is computed from the native growth model for accurate
relative-size estimates, while \\k\_{VB}\\ must be the VB-equivalent
growth coefficient (or native VB \\k\\) because the regression
coefficients were calibrated against VB parameters.

## Two Roles of k (Growth-Based Formulation)

Like Chen-Watanabe, the growth-based Lorenzen involves \\k\\ in two
distinct roles:

- Calibration coefficient:

  \\\ln(k\_{VB})\\ enters the regression equation as a predictor. Since
  Lorenzen (2022) calibrated this term against VB parameters, it must be
  the VB growth coefficient or VB-equivalent. Passed via `k_vb`.

- Growth trajectory:

  \\L(t)/L\_\infty\\ represents relative body size. The native growth
  model produces more accurate predictions than a VB approximation.
  Controlled by `k` and `growth_model`.

## References

Lorenzen, K. (1996). The relationship between body weight and natural
mortality in juvenile and adult fish. *Journal of Fish Biology*, 49(4),
627-642.

Lorenzen, K. (2022). Size- and age-dependent natural mortality in fish
populations. *Fisheries Research*, 255, 106454.

## Examples

``` r
if (FALSE) { # \dontrun{
ages <- seq(0.5, 30, by = 0.5)
lw_fun <- function(L) 0.0001 * L^3.1

# Weight-based with Gompertz trajectory
M_wt <- M_lorenzen(
  age = ages, Linf = 100, L0 = 25, k = 0.15,
  lw_fun = lw_fun, weight_based = TRUE,
  growth_model = "gompertz"
)

# Growth-based with Gompertz trajectory + VB-equivalent k
M_gr <- M_lorenzen(
  age = ages, Linf = 100, L0 = 25,
  k = 0.15,                         # Native Gompertz k for L(t)
  k_vb = 0.08,                      # VB-equivalent k for coefficient
  growth_model = "gompertz"
)
} # }
```

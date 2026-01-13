# Lorenzen Natural Mortality Model

Computes size-dependent natural mortality following Lorenzen (1996,
2022). Supports both weight-based and growth-based formulations.

## Usage

``` r
M_lorenzen(
  age,
  Linf,
  L0,
  k,
  lw_fun = NULL,
  weight_based = FALSE,
  growth_model = c("vb", "gompertz", "logistic"),
  sample_params = FALSE
)
```

## Arguments

- age:

  Numeric vector of ages at which to compute mortality.

- Linf:

  Asymptotic length.

- L0:

  Length at birth.

- k:

  Growth coefficient. For `weight_based = FALSE`, should be
  VB-equivalent \\k\\ (use
  [`compute_k_vb_equivalent`](https://brian-j-moe.github.io/vitalBayes/reference/compute_k_vb_equivalent.md)).

- lw_fun:

  Function mapping length to weight (required if `weight_based = TRUE`).

- weight_based:

  Logical. If `TRUE`, uses weight-based formulation. If `FALSE`
  (default), uses growth-based formulation.

- growth_model:

  Character. Growth model for length prediction (only used for
  weight-based formulation): `"vb"`, `"gompertz"`, or `"logistic"`.

- sample_params:

  Logical. If `TRUE`, samples allometric parameters from their
  distributions. If `FALSE`, uses mean values.

## Value

Numeric vector of instantaneous mortality rates.

## Details

Two formulations are available:

Weight-based (Lorenzen 1996): \$\$M(W) = \alpha \cdot W^{\beta}\$\$
where \\\alpha \sim N(3.69, 0.502)\\ and \\\beta \sim N(-0.305,
0.029)\\.

Growth-based (Lorenzen 2022): \$\$\ln M = 0.28 - 1.30 \ln(L/L\_\infty) +
1.08 \ln(k)\$\$ This formulation was calibrated using von Bertalanffy
parameters, so \\k\\ should be the VB-equivalent \\k\\ when using fits
from other growth models.

## References

Lorenzen, K. (1996). The relationship between body weight and natural
mortality in juvenile and adult fish. *Journal of Fish Biology*, 49(4),
627-642.

Lorenzen, K. (2022). Size- and age-dependent natural mortality in fish
populations. *Fisheries Research*, 255, 106454.

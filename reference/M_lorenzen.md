# Lorenzen Natural Mortality Model

Computes size-dependent natural mortality following Lorenzen (1996,
2022).

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

## Value

Numeric vector of instantaneous mortality rates.

## Details

Two formulations are available:

**Weight-based** (Lorenzen 1996): \$\$M(W) = \alpha \cdot W^{\beta}\$\$

**Growth-based** (Lorenzen 2022): \$\$\ln M = 0.28 - 1.30
\ln(L/L\_\infty) + 1.08 \ln(k)\$\$

For the growth-based formulation, M_inf (VB-derived) is used as k for
consistency with the Chen-Watanabe framework.

## References

Lorenzen, K. (1996). The relationship between body weight and natural
mortality in juvenile and adult fish. *Journal of Fish Biology*, 49(4),
627-642.

Lorenzen, K. (2022). Size- and age-dependent natural mortality in fish
populations. *Fisheries Research*, 255, 106454.

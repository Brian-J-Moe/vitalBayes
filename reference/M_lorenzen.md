# Lorenzen Natural Mortality Model

Computes size-dependent natural mortality following Lorenzen (1996,
2022). Supports both weight-based and von Bertalanffy growth-based
formulations.

## Usage

``` r
M_lorenzen(
  age,
  Linf,
  k,
  t0,
  weight_based = FALSE,
  lw_fun = NULL,
  sample_params = TRUE
)
```

## Arguments

- age:

  Numeric vector of ages at which to compute mortality.

- Linf:

  Asymptotic length.

- k:

  von Bertalanffy growth coefficient.

- t0:

  Theoretical age at length zero.

- weight_based:

  Logical. If TRUE, uses weight-based formulation requiring `lw_fun`. If
  FALSE, uses growth-based formulation.

- lw_fun:

  Function mapping length to weight in grams (required if weight_based =
  TRUE).

- sample_params:

  Logical. If TRUE, samples allometric parameters from their
  distributions (for Monte Carlo simulation). If FALSE, uses mean
  values.

## Value

Numeric vector of instantaneous mortality rates.

## Details

**Weight-based formulation:** \$\$M(W) = \alpha W^{\beta}\$\$ where
\\\alpha \sim N(3.69, 0.502)\\ and \\\beta \sim N(-0.305, 0.029)\\.

**Growth-based formulation (Lorenzen 2022):** \$\$\ln M = 0.28 - 1.30
\ln(L/L\_\infty) + 1.08 \ln(k)\$\$ with uncertainty incorporated via
parameter distributions.

## References

Lorenzen, K. (1996). The relationship between body weight and natural
mortality in juvenile and adult fish: a comparison of natural ecosystems
and aquaculture. Journal of Fish Biology, 49(4), 627-642.

Lorenzen, K. (2022). Size- and age-dependent natural mortality in fish
populations: Biology, models, implications, and a generalized
length-inverse model. Fisheries Research, 255, 106454.

## Examples

``` r
if (FALSE) { # \dontrun{
# Growth-based formulation
ages <- seq(0.1, 30, length.out = 100)
M <- M_lorenzen(ages, Linf = 200, k = 0.15, t0 = -1, weight_based = FALSE)
plot(ages, M, type = "l")
} # }
```

# Compute Asymptotic Mortality Rate (M_inf)

Computes the asymptotic mortality rate \\M\_\infty\\ from biological
milestones using the VB formula. This serves as the unified mortality
scaling factor in the Chen-Watanabe framework.

## Usage

``` r
compute_Minf(Linf, L0, Lmat, tmat, warn = TRUE)
```

## Arguments

- Linf:

  Numeric vector. Asymptotic length.

- L0:

  Numeric vector. Length at birth.

- Lmat:

  Numeric vector. Length at maturity.

- tmat:

  Numeric vector. Age at maturity.

- warn:

  Logical. Warn on invalid values?

## Value

Numeric vector of asymptotic mortality rates.

## Details

\$\$M\_\infty = \frac{1}{t\_{mat}} \ln\left(\frac{L\_\infty -
L_0}{L\_\infty - L\_{mat}}\right)\$\$

This formula is used regardless of growth model, providing a consistent
mortality anchor that prevents the ~10-fold survival differences arising
from model-native k formulas.

## Examples

``` r
Minf <- compute_Minf(Linf = 126, L0 = 35, Lmat = 83, tmat = 47)
```

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

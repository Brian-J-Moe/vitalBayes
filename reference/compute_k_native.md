# Compute Native Growth Coefficient for a Specific Model

Derives the growth coefficient \\k\\ for a specified growth model from
biological milestones \\(L\_\infty, L_0, L\_{mat}, t\_{mat})\\.

## Usage

``` r
compute_k_native(
  Linf,
  L0,
  Lmat,
  tmat,
  growth_model = c("vb", "gompertz", "logistic"),
  warn = TRUE
)
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

- growth_model:

  Character. Growth model: `"vb"`, `"gompertz"`, or `"logistic"`.

- warn:

  Logical. Warn on invalid values?

## Value

Numeric vector of native k values. Invalid values returned as NA.

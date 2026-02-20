# Compute Maximum Age from Growth Parameters

Compute Maximum Age from Growth Parameters

## Usage

``` r
compute_tmax(
  Linf,
  L0,
  k,
  growth_model = c("vb", "gompertz", "logistic"),
  Linf_factor = 0.99
)
```

## Arguments

- Linf:

  Numeric. Asymptotic length.

- L0:

  Numeric. Length at birth.

- k:

  Numeric. Native growth coefficient.

- growth_model:

  Character. Growth model.

- Linf_factor:

  Numeric in (0,1). Fraction of Linf (default 0.99).

## Value

Numeric tmax value.

# Compute Maximum Age from Growth Parameters

Estimates tmax as age when L(t) reaches a specified fraction of Linf.

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

## Examples

``` r
tmax <- compute_tmax(Linf = 126, L0 = 35, k = 0.016,
                     growth_model = "vb", Linf_factor = 0.99)
```

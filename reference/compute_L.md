# Compute Length at Age L(t)

Computes predicted length at age using a specified growth model.

## Usage

``` r
compute_L(age, Linf, L0, k, growth_model = c("vb", "gompertz", "logistic"))
```

## Arguments

- age:

  Numeric vector. Ages.

- Linf:

  Numeric. Asymptotic length.

- L0:

  Numeric. Length at birth.

- k:

  Numeric. Native growth coefficient.

- growth_model:

  Character. Growth model.

## Value

Numeric vector of predicted lengths.

## Examples

``` r
L_vb <- compute_L(0:50, Linf = 126, L0 = 35, k = 0.016, growth_model = "vb")
```

# Compute Length at Age L(t)

Compute Length at Age L(t)

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

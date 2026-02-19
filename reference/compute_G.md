# Compute Relative Size G(t) = L(t)/Linf

Computes the relative size at specified ages using a given growth model.
This is the denominator in the mortality formula M(t) = M_inf / G(t).

## Usage

``` r
compute_G(age, Linf, L0, k, growth_model = c("vb", "gompertz", "logistic"))
```

## Arguments

- age:

  Numeric vector. Ages at which to compute G(t).

- Linf:

  Numeric. Asymptotic length.

- L0:

  Numeric. Length at birth.

- k:

  Numeric. Native growth coefficient.

- growth_model:

  Character. Growth model.

## Value

Numeric vector of G(t) values bounded in (0, 1\].

## Examples

``` r
ages <- seq(0, 100, by = 1)
G_vb <- compute_G(ages, Linf = 126, L0 = 35, k = 0.016, growth_model = "vb")
```

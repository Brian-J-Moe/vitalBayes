# Scale Mortality Schedule to Target Mean

Scale Mortality Schedule to Target Mean

## Usage

``` r
scale_mortality(M, M_target = NULL, tmax = NULL, p = 0.001)
```

## Arguments

- M:

  Numeric vector of mortality rates.

- M_target:

  Target mean mortality. Can be scalar, function of tmax, or NULL.

- tmax:

  Maximum age (required if M_target is function or NULL).

- p:

  Survival probability (used if M_target NULL).

## Value

Scaled mortality vector.

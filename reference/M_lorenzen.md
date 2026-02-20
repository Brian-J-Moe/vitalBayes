# Lorenzen Natural Mortality Model

Lorenzen Natural Mortality Model

## Usage

``` r
M_lorenzen(
  age,
  Linf,
  L0,
  Lmat,
  tmat,
  lw_fun = NULL,
  weight_based = FALSE,
  growth_model = c("vb", "gompertz", "logistic"),
  sample_params = TRUE
)
```

## Arguments

- age:

  Numeric vector of ages.

- Linf:

  Asymptotic length.

- L0:

  Length at birth.

- Lmat:

  Length at maturity.

- tmat:

  Age at maturity.

- lw_fun:

  Length-weight function (required if weight_based = TRUE).

- weight_based:

  Use weight-based formulation?

- growth_model:

  Growth model for L(t).

- sample_params:

  Sample parameters from their distributions?

## Value

Numeric vector of instantaneous mortality rates.

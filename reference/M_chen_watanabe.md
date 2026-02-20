# Chen-Watanabe Natural Mortality (Model-Dependent)

Computes age-specific natural mortality using a generalized
Chen-Watanabe framework where G(t) is derived from the native growth
model trajectory.

## Usage

``` r
M_chen_watanabe(
  age,
  Linf,
  L0,
  Lmat,
  tmat,
  growth_model = c("vb", "gompertz", "logistic"),
  tmax = NULL,
  Linf_factor = 0.99,
  two_phase = FALSE,
  late_model = c("gompertz", "logistic"),
  tm_factor = 2/3,
  M_mult = 2,
  smooth_factor = 1/3
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

- growth_model:

  Character. Growth model for G(t).

- tmax:

  Maximum age (computed if NULL).

- Linf_factor:

  Fraction of Linf for tmax estimation.

- two_phase:

  Use two-phase senescence model?

- late_model:

  Senescence model: `"gompertz"` or `"logistic"`.

- tm_factor:

  Transition age as fraction of tmat.

- M_mult:

  Mortality multiplier for senescence.

- smooth_factor:

  Transition smoothness.

## Value

Numeric vector of instantaneous mortality rates.

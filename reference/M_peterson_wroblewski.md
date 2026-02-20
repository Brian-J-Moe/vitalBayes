# Peterson-Wroblewski Natural Mortality Model

Peterson-Wroblewski Natural Mortality Model

## Usage

``` r
M_peterson_wroblewski(
  age,
  Linf,
  L0,
  Lmat,
  tmat,
  lw_fun,
  growth_model = c("vb", "gompertz", "logistic")
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

  Length-weight function.

- growth_model:

  Character. Growth model for L(t).

## Value

Numeric vector of instantaneous mortality rates.

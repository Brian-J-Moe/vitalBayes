# Peterson-Wroblewski Natural Mortality Model

Computes weight-based natural mortality following Peterson & Wroblewski
(1984).

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

  Length-weight function: `lw_fun(L)` returns weight in grams.

- growth_model:

  Character. Growth model for L(t).

## Value

Numeric vector of instantaneous mortality rates.

## Details

\$\$M(W) = 1.92 \cdot W^{-0.25}\$\$

**Warning**: This model was calibrated on teleost fishes and may produce
biologically implausible mortality rates for elasmobranchs.

## References

Peterson, I., & Wroblewski, J. S. (1984). Mortality rate of fishes in
the pelagic ecosystem. *Canadian Journal of Fisheries and Aquatic
Sciences*, 41(7), 1117-1120.

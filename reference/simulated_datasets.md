# Simulated von Bertalanffy datasets

Three example datasets generated with
[`simulate_vb_growth_data()`](https://brian-j-moe.github.io/vitalBayes/reference/simulate_vb_growth_data.md)
to illustrate common sample-size scenarios (imbalanced, larger balanced,
and limited sample sizes) for fitting and testing von Bertalanffy growth
and maturity relationships.

## Usage

``` r
imbalanced_data

growth_data

limited_data
```

## Format

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with 5 columns:

- sex:

  Character. Sex for non-embryo observations (`"female"` or `"male"`);
  `NA` for embryo rows.

- mat:

  Integer. Maturity indicator for non-embryo observations (`0` =
  immature, `1` = mature); `NA` for embryo rows.

- fl:

  Numeric. Fork length (cm). For embryos, this is the embryo length
  measurement.

- age:

  Numeric. Age (years) for non-embryo observations; `NA` for embryo
  rows.

- embryo:

  Logical. `TRUE` for embryo observations; `FALSE` otherwise.

An object of class `data.table` (inherits from `data.frame`) with 197
rows and 5 columns.

An object of class `data.table` (inherits from `data.frame`) with 391
rows and 5 columns.

An object of class `data.table` (inherits from `data.frame`) with 47
rows and 5 columns.

## Details

All three datasets share the same structure and are intended as
lightweight, reproducible examples for package vignettes, unit tests,
and demonstrations. Non-embryo observations are simulated from a von
Bertalanffy growth model with sex-specific parameter distributions,
observation error on length, a truncated-exponential age distribution
(older ages increasingly rare), and a sex-specific logistic maturity
ogive on observed length. Embryo rows include lengths only and are
flagged with `embryo = TRUE`.

Dataset-specific sizes:

- `imbalanced_data`:

  Default call to
  [`simulate_vb_growth_data()`](https://brian-j-moe.github.io/vitalBayes/reference/simulate_vb_growth_data.md)
  (defaults: 150 females, 34 males, 13 embryos).

- `growth_data`:

  Larger dataset (189 females, 176 males, 26 embryos; `seed = 1234`).

- `limited_data`:

  Small dataset (24 females, 18 males, 5 embryos; `seed = 12345`).

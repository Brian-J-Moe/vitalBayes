# Create LOO Comparison Table

Creates a publication-ready table from LOO comparison results.

## Usage

``` r
create_loo_table(comp_dt, format = "data.table", digits = 2)
```

## Arguments

- comp_dt:

  A data.table from
  [`compare_loo`](https://brian-j-moe.github.io/vitalBayes/reference/compare_loo.md).

- format:

  Character. Output format: `"data.table"` (default), `"kable"`, or
  `"gt"`.

- digits:

  Integer. Number of decimal places. Default 2.

## Value

Formatted table in specified format.

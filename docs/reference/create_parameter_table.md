# Create Parameter Summary Table

Creates a publication-ready parameter summary table from a vitalBayes
fit.

## Usage

``` r
create_parameter_table(
  fit,
  params = NULL,
  probs = c(0.025, 0.5, 0.975),
  digits = 3,
  format = "data.table"
)
```

## Arguments

- fit:

  A CmdStanMCMC object from any vitalBayes fitting function.

- params:

  Character vector of parameter names. If `NULL`, auto-selects key
  parameters based on model type.

- probs:

  Numeric vector of quantiles. Default c(0.025, 0.5, 0.975).

- digits:

  Integer. Decimal places. Default 3.

- format:

  Character. Output format: `"data.table"`, `"kable"`, or `"gt"`.

## Value

Formatted parameter table.

## Examples

``` r
if (FALSE) { # \dontrun{
create_parameter_table(growth_fit)
create_parameter_table(growth_fit, format = "kable")
} # }
```

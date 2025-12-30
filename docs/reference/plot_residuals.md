# Plot Residual Diagnostics

Creates diagnostic plots for growth model residuals.

## Usage

``` r
plot_residuals(
  fit,
  data = NULL,
  age_col = "age",
  length_col = "length",
  sex_col = "sex",
  type = c("all", "fitted", "qq", "histogram", "age"),
  colors = NULL,
  colorblind = FALSE,
  base_size = 12
)
```

## Arguments

- fit:

  A CmdStanMCMC object from
  [`fit_bayesian_growth`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md).

- data:

  Optional data.table with observed data.

- age_col:

  Character. Column name for age.

- length_col:

  Character. Column name for length.

- sex_col:

  Character. Column name for sex.

- type:

  Character. Type of plot: `"all"` (default), `"fitted"`, `"qq"`,
  `"histogram"`, or `"age"`.

- colors:

  Color palette.

- colorblind:

  Logical. Use colorblind-safe palette? Default `FALSE`.

- base_size:

  Base font size.

## Value

A ggplot object or list of ggplot objects.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_residuals(growth_fit, data = shark_data)
plot_residuals(growth_fit, type = "qq")
} # }
```

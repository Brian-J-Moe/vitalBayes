# Compare Multiple Growth Models

Creates a faceted comparison of multiple growth model fits.

## Usage

``` r
compare_growth_models(
  ...,
  data = NULL,
  age_col = "age",
  length_col = "length",
  sex_col = "sex",
  sex_labels = c(`1` = "Female", `2` = "Male"),
  age_range = NULL,
  k_based_list = NULL,
  which_model_list = NULL,
  colors = NULL,
  colorblind = FALSE,
  base_size = 12
)
```

## Arguments

- ...:

  Named CmdStanMCMC objects from
  [`fit_bayesian_growth`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md).
  Names become facet labels.

- data:

  Optional data.table/data.frame with observed data.

- age_col, length_col, sex_col:

  Column names in data.

- sex_labels:

  Named character vector. Labels for sex in plot legend and facets.
  Default `c("1" = "Female", "2" = "Male")`.

- age_range:

  Numeric vector of length 2. If `NULL`, auto-determined.

- k_based_list:

  Logical vector indicating which models are k-based. If `NULL`, assumes
  all are maturity-based.

- which_model_list:

  Integer vector of model types. If `NULL`, assumes VBGM.

- colors:

  Color palette.

- colorblind:

  Logical. Use colorblind-safe palette? Default `FALSE`.

- base_size:

  Base font size.

## Value

A ggplot object with faceted model comparison.

## See also

[`vignette("visualization")`](https://brian-j-moe.github.io/vitalBayes/articles/visualization.md)
for comprehensive plotting examples.

[`vignette("model_diagnostics")`](https://brian-j-moe.github.io/vitalBayes/articles/model_diagnostics.md)
for model comparison via LOO-CV.

[Statistical Methods: Model
Assessment](https://brian-j-moe.github.io/vitalBayes/doc/Understanding_vitalBayes.html#assessment)
for LOO-CV theory.

[`plot_growth_curve`](https://brian-j-moe.github.io/vitalBayes/reference/plot_growth_curve.md),
[`compute_loo`](https://brian-j-moe.github.io/vitalBayes/reference/compute_loo.md),
[`compare_loo`](https://brian-j-moe.github.io/vitalBayes/reference/compare_loo.md)

## Examples

``` r
if (FALSE) { # \dontrun{
compare_growth_models(
  "von Bertalanffy" = vb_fit,
  "Gompertz" = gomp_fit,
  "Logistic" = log_fit,
  data = shark_data
)

# With French labels
compare_growth_models(
  "von Bertalanffy" = vb_fit,
  "Gompertz" = gomp_fit,
  data = shark_data,
  sex_labels = c("1" = "Femelle", "2" = "Mâle")
)
} # }
```

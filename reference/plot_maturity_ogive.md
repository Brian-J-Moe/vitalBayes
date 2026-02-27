# Plot Maturity Ogive

Creates maturity ogive plots with credible intervals. Supports both
length-at-maturity and age-at-maturity models.

## Usage

``` r
plot_maturity_ogive(
  fit,
  type = c("length", "age"),
  data = NULL,
  x_col = NULL,
  maturity_col = "maturity",
  sex_col = "sex",
  sex_labels = c(`1` = "Female", `2` = "Male"),
  x_range = NULL,
  n_points = 500,
  show_data = TRUE,
  show_rug = TRUE,
  show_x50_line = TRUE,
  data_alpha = 0.5,
  ribbon_alpha = 0.25,
  line_width = 1.2,
  colors = NULL,
  colorblind = FALSE,
  facet_sex = TRUE,
  x_lab = NULL,
  y_lab = "Probability of Maturity",
  title = NULL,
  base_size = 12,
  theme_args = list(),
  additional_layers = list()
)
```

## Arguments

- fit:

  A CmdStanMCMC object from
  [`fit_bayesian_maturity`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md).

- type:

  Character. Type of maturity model: `"length"` or `"age"`.

- data:

  Optional data.table/data.frame with observed data.

- x_col:

  Character. Column name for predictor (length or age).

- maturity_col:

  Character. Column name for maturity status.

- sex_col:

  Character. Column name for sex.

- sex_labels:

  Named character vector. Labels for sex in plot legend and facets.
  Default `c("1" = "Female", "2" = "Male")`.

- x_range:

  Numeric vector of length 2. Range for predictions.

- n_points:

  Integer. Number of points for ogive. Default 500.

- show_data:

  Logical. Show observed data? Default `TRUE`.

- show_rug:

  Logical. Show rug plot at bottom? Default `TRUE`.

- show_x50_line:

  Logical. Show vertical line at x50? Default `TRUE`.

- data_alpha:

  Numeric. Transparency for data. Default 0.5.

- ribbon_alpha:

  Numeric. Transparency for ribbon. Default 0.25.

- line_width:

  Numeric. Width of ogive line. Default 1.2.

- colors:

  Character vector of colors.

- colorblind:

  Logical. Use colorblind-safe palette? Default `FALSE`.

- facet_sex:

  Logical. Facet by sex for two-sex models? Default `TRUE`.

- x_lab, y_lab:

  Character. Axis labels. Auto-generated if `NULL`.

- title:

  Character. Plot title.

- base_size:

  Numeric. Base font size. Default 12.

- theme_args, additional_layers:

  Lists for customization.

## Value

A ggplot object.

## See also

[`vignette("visualization")`](https://brian-j-moe.github.io/vitalBayes/articles/visualization.md)
for comprehensive plotting examples.

[`vignette("fit_bayesian_maturity")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_maturity.md)
for maturity model fitting.

[`fit_bayesian_maturity`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md),
[`vital_colors`](https://brian-j-moe.github.io/vitalBayes/reference/vital_colors.md),
[`theme_vital`](https://brian-j-moe.github.io/vitalBayes/reference/theme_vital.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Length-at-maturity
plot_maturity_ogive(L50_fit, type = "length", data = shark_data)

# French labels
plot_maturity_ogive(
  L50_fit,
  type = "length",
  data = shark_data,
  sex_labels = c("1" = "Femelle", "2" = "Mâle")
)

# Age-at-maturity with custom styling
plot_maturity_ogive(
  t50_fit,
  type = "age",
  data = shark_data,
  colors = c("#E63946", "#1D3557"),
  show_rug = TRUE
)
} # }
```

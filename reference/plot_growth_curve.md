# Plot Growth Curves

Creates publication-quality growth curve plots with credible intervals.
Automatically handles single-sex and two-sex models.

## Usage

``` r
plot_growth_curve(
  fit,
  data = NULL,
  age_col = "age",
  length_col = "length",
  sex_col = "sex",
  sex_labels = c(`1` = "Female", `2` = "Male"),
  age_range = NULL,
  n_points = 100,
  show_data = TRUE,
  show_50_ci = TRUE,
  data_alpha = 0.4,
  data_size = 1.5,
  ribbon_alpha = 0.2,
  ribbon_alpha_50 = 0.3,
  line_width = 1.2,
  colors = NULL,
  colorblind = FALSE,
  facet_sex = TRUE,
  x_lab = "Age (years)",
  y_lab = "Length (cm)",
  title = NULL,
  k_based = FALSE,
  which_model = 1L,
  base_size = 12,
  theme_args = list(),
  additional_layers = list()
)
```

## Arguments

- fit:

  A CmdStanMCMC object from
  [`fit_bayesian_growth`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md).

- data:

  Optional data.table/data.frame with observed data for overlay.

- age_col:

  Character. Column name for age in `data`. Default "age".

- length_col:

  Character. Column name for length in `data`. Default "length".

- sex_col:

  Character. Column name for sex in `data`. Default "sex".

- sex_labels:

  Named character vector. Labels for sex in plot legend and facets.
  Default `c("1" = "Female", "2" = "Male")`. Can be customized for other
  languages, e.g., `c("1" = "Femelle", "2" = "Mâle")` (French),
  `c("1" = "Hembra", "2" = "Macho")` (Spanish), or
  `c("1" = "メス", "2" = "オス")` (Japanese).

- age_range:

  Numeric vector of length 2. Age range for predictions. If `NULL`,
  auto-determined from data or defaults to c(0, 30).

- n_points:

  Integer. Number of points for prediction grid. Default 100.

- show_data:

  Logical. Overlay observed data points? Default `TRUE`.

- show_50_ci:

  Logical. Show 50% credible interval ribbon? Default `TRUE`.

- data_alpha:

  Numeric. Transparency for data points. Default 0.4.

- data_size:

  Numeric. Size for data points. Default 1.5.

- ribbon_alpha:

  Numeric. Transparency for 95% CI ribbon. Default 0.2.

- ribbon_alpha_50:

  Numeric. Transparency for 50% CI ribbon. Default 0.3.

- line_width:

  Numeric. Width of median line. Default 1.2.

- colors:

  Character vector of colors. If `NULL`, uses
  [`vital_colors()`](https://brian-j-moe.github.io/vitalBayes/reference/vital_colors.md).

- colorblind:

  Logical. Use colorblind-safe palette? Default `FALSE`. When `TRUE`,
  overrides `colors` with `vital_colors(, "sex_cb")` for two-sex models
  or `vital_colors(, "okabe_ito")` otherwise.

- facet_sex:

  Logical. For two-sex models, facet by sex? Default `TRUE`.

- x_lab, y_lab:

  Character. Axis labels.

- title:

  Character. Plot title. Default `NULL`.

- k_based:

  Logical. Was the model k-based? Default `FALSE`.

- which_model:

  Integer. 1=VBGM, 2=Gompertz, 3=Logistic. Default 1.

- base_size:

  Numeric. Base font size for theme. Default 12.

- theme_args:

  List of additional arguments for
  [`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html).

- additional_layers:

  List of additional ggplot2 layers to add.

## Value

A ggplot object.

## See also

[`vignette("visualization")`](https://brian-j-moe.github.io/vitalBayes/articles/visualization.md)
for comprehensive plotting examples.

[`vignette("fit_bayesian_growth")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_growth.md)
for growth model fitting.

[`compare_growth_models`](https://brian-j-moe.github.io/vitalBayes/reference/compare_growth_models.md),
[`plot_residuals`](https://brian-j-moe.github.io/vitalBayes/reference/plot_residuals.md),
[`vital_colors`](https://brian-j-moe.github.io/vitalBayes/reference/vital_colors.md),
[`theme_vital`](https://brian-j-moe.github.io/vitalBayes/reference/theme_vital.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
plot_growth_curve(growth_fit, data = shark_data)

# French labels
plot_growth_curve(
  growth_fit,
  data = shark_data,
  sex_labels = c("1" = "Femelle", "2" = "Mâle")
)

# Spanish labels
plot_growth_curve(
  growth_fit,
  data = shark_data,
  sex_labels = c("1" = "Hembra", "2" = "Macho")
)

# Customize appearance
plot_growth_curve(
  growth_fit,
  data = shark_data,
  colors = c("#E63946", "#1D3557"),
  theme_args = list(legend.position = "bottom"),
  additional_layers = list(
    geom_hline(yintercept = 100, linetype = "dashed")
  )
)
} # }
```

# Plot Birth Ogive

Creates publication-quality birth transition plots.

## Usage

``` r
plot_birth_ogive(
  fit,
  data = NULL,
  length_col = "length",
  status_col = "status",
  length_range = NULL,
  n_points = 100,
  show_data = TRUE,
  show_rug = TRUE,
  show_b50_line = TRUE,
  color = NULL,
  ribbon_alpha = 0.25,
  x_lab = "Length (cm)",
  y_lab = "Probability of Free-Swimming",
  title = NULL,
  base_size = 12,
  theme_args = list(),
  additional_layers = list()
)
```

## Arguments

- fit:

  A CmdStanMCMC object from
  [`fit_bayesian_birth`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_birth.md).

- data:

  Optional data.table/data.frame with observed data.

- length_col:

  Character. Column name for length.

- status_col:

  Character. Column name for birth status.

- length_range:

  Numeric vector of length 2.

- n_points:

  Integer. Number of points. Default 100.

- show_data:

  Logical. Show data? Default `TRUE`.

- show_rug:

  Logical. Show rug? Default `TRUE`.

- show_b50_line:

  Logical. Show b50 line? Default `TRUE`.

- color:

  Character. Line color. Default uses primary vitalBayes magenta.

- ribbon_alpha:

  Numeric. Ribbon transparency.

- x_lab, y_lab:

  Character. Axis labels.

- title:

  Character. Plot title.

- base_size:

  Numeric. Base font size.

- theme_args, additional_layers:

  Lists for customization.

## Value

A ggplot object.

# Monte Carlo Simulation of Survivorship

Simulates cohort survivorship trajectories using age-specific mortality
schedules. Computes cumulative survival, mean age at death, and
probabilities of reaching maturity and maximum age.

## Usage

``` r
simulate_survivorship(
  mc_object,
  n = 50000,
  n_iter = 5000,
  mode = c("random", "per_set"),
  seed = 1234,
  palette = c("synthwave", "viridis", "okabe", "plasma", "inferno"),
  bar_fill = NULL,
  bar_color = NULL,
  bar_alpha = 0.7,
  error_color = "black",
  error_width = 0.4,
  vline_color = NULL,
  vline_width = 1,
  vline_type = "solid",
  rect_fill = NULL,
  rect_alpha = 0.4,
  x_breaks = function(tmax) seq(0, ceiling(tmax), by = 5),
  title = NULL,
  subtitle = "Aggregated across simulations",
  base_size = 11,
  font_family = "serif",
  print_plot = TRUE,
  show_progress = TRUE
)
```

## Arguments

- mc_object:

  Output from
  [`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md),
  containing `$Schedules` and `$Parameters`.

- n:

  Starting cohort size. Default 50000.

- n_iter:

  Number of simulation iterations. Default 5000.

- mode:

  Character. Either:

  - `"random"` - Each iteration uses a randomly sampled parameter set

  - `"per_set"` - Run n_iter simulations for each parameter set

- seed:

  Random seed. Default 1234.

- palette:

  Color palette: `"synthwave"`, `"viridis"`, `"okabe"`.

- bar_fill, bar_color:

  Colors for survivorship bars.

- bar_alpha:

  Transparency for bars. Default 0.7.

- error_color:

  Color for error bars. Default "black".

- error_width:

  Width of error bars. Default 0.4.

- vline_color, vline_width, vline_type:

  Mean age-at-death line aesthetics.

- rect_fill, rect_alpha:

  Uncertainty band for mean age-at-death.

- x_breaks:

  Function returning x-axis breaks given max tmax.

- title, subtitle:

  Plot title and subtitle.

- base_size:

  Base font size. Default 11.

- font_family:

  Font family. Default "serif".

- print_plot:

  Print plot? Default TRUE.

- show_progress:

  Show progress messages? Default TRUE.

## Value

A list with:

- Aggregate:

  List with overall survival curve, age-at-death, and survival
  probabilities to tmat/tmax

- Per_Set:

  List with per-parameter-set summaries

- Raw:

  data.table of all trajectory data

- Plot:

  ggplot2 object

## Details

For each iteration, the function:

1.  Integrates the mortality schedule to get cumulative survival: \\S(t)
    = \exp(-\int_0^t M(a) da)\\

2.  Simulates discrete cohort transitions using binomial survival

3.  Tracks deaths to compute age-at-death distribution

## Examples

``` r
if (FALSE) { # \dontrun{
mort <- get_stochastic_mortality(method = "CW", ...)
surv <- simulate_survivorship(mort, n = 1e5, n_iter = 2000)
surv$Plot
} # }
```

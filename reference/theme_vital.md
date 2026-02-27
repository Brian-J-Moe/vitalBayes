# vitalBayes ggplot2 Theme

A ggplot2 theme inspired by the vitalBayes retro-wave aesthetic.
Provides a clean, modern look suitable for scientific publications while
maintaining visual consistency with the package branding.

## Usage

``` r
theme_vital(base_size = 12, base_family = "serif", grid = FALSE, dark = FALSE)
```

## Arguments

- base_size:

  Numeric. Base font size. Default 12.

- base_family:

  Character. Base font family. Default "".

- grid:

  Logical. Show grid lines? Default `TRUE` for major only.

- dark:

  Logical. Use dark theme variant? Default `FALSE`.

## Value

A ggplot2 theme object.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) +
  geom_point(color = vital_colors(1)) +
  theme_vital()
} # }
```

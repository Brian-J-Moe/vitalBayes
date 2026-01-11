# vitalBayes Extended Color Palette

Returns color palettes for vitalBayes visualizations with support for
synthwave (default), viridis, and colorblind-friendly (Okabe-Ito)
schemes.

## Usage

``` r
vital_palette(
  n = NULL,
  type = c("synthwave", "viridis", "okabe", "plasma", "inferno"),
  alpha = 1
)
```

## Arguments

- n:

  Number of colors needed. If NULL, returns the full palette.

- type:

  Character. One of:

  - `"synthwave"` - Retro-wave palette from vitalBayes hex sticker
    (default)

  - `"viridis"` - Viridis colorblind-friendly sequential palette

  - `"okabe"` - Okabe-Ito colorblind-friendly qualitative palette

  - `"plasma"` - Plasma sequential palette

  - `"inferno"` - Inferno sequential palette

- alpha:

  Numeric. Transparency value (0-1). Default 1 (opaque).

## Value

A character vector of hex color codes.

## Examples

``` r
if (FALSE) { # \dontrun{
vital_palette(4, "synthwave")
vital_palette(6, "okabe")
vital_palette(10, "viridis", alpha = 0.8)
} # }
```

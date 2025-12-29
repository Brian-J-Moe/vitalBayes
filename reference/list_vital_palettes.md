# List Available Color Palettes

Prints information about all available color palettes in vitalBayes,
including their colorblind accessibility status.

## Usage

``` r
list_vital_palettes(show_colors = FALSE)
```

## Arguments

- show_colors:

  Logical. Display color swatches? Default `FALSE`. Requires a graphics
  device if `TRUE`.

## Value

Invisibly returns a data.frame with palette information.

## Examples

``` r
list_vital_palettes()
#> 
#> vitalBayes Color Palettes
#> ==================================================
#>   default      (5 colors)
#>     vitalBayes logo retro-wave
#>   sex          (2 colors)
#>     Two-sex (magenta/cyan)
#>   gradient     (6 colors)
#>     Continuous gradient (pink to cyan)
#>   sunset       (5 colors)
#>     Warm sunset tones
#>   okabe_ito    (8 colors) [CB-safe]
#>     Okabe-Ito universal design
#>   viridis      (5 colors) [CB-safe]
#>     Viridis perceptually uniform
#>   tol_bright   (7 colors) [CB-safe]
#>     Paul Tol bright scheme
#>   sex_cb       (2 colors) [CB-safe]
#>     Two-sex colorblind-safe (orange/blue)
```

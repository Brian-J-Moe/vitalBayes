# Check if Palette is Colorblind-Friendly

Utility function to check if a palette name is colorblind-accessible.

## Usage

``` r
is_colorblind_safe(palette)
```

## Arguments

- palette:

  Character. Palette name to check.

## Value

Logical. TRUE if colorblind-friendly.

## Examples

``` r
is_colorblind_safe("okabe_ito")  # TRUE
#> [1] TRUE
is_colorblind_safe("default")    # FALSE
#> [1] FALSE
```

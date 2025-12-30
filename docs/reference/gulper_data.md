# Gulper Shark Life History Data

Length, age, and maturity data for the little gulper shark,
*Centrophorus uyato*. Contains `NA`s to illustrate that
[`fit_bayesian_growth()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md)
and
[`fit_bayesian_maturity()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md)
can handle missing data.

## Usage

``` r
data(gulper_data)
```

## Format

A `data.table` with 668 rows and 6 variables:

- sex:

  Integer. Sex indicator: `1` = female, `2` = male.

- mat:

  Integer. Maturity status: `0` = immature, `1` = mature.

- fl:

  Numeric. Fork length in cm.

- age1:

  Numeric. Age in years estimated from the dorsal fin spine "inner"
  growth band region.

- age2:

  Numeric. Age in years estimated from the sum of dorsal fin spine
  "inner" and "outer" growth band regions.

- embryo:

  Logical. If `TRUE`, individual was an embryo rather than
  free-swimming.

## Details

The sex and maturity variables use integer coding for compatibility with
Stan models:

- `sex`: 1 = female, 2 = male (auto-detected by vitalBayes functions)

- `mat`: 0 = immature, 1 = mature (standard binary coding)

Age estimates from two band-counting methodologies are provided to allow
comparison of growth models under different aging assumptions.

## Data Subsets

Common data subsets for analysis:

- Embryos: `gulper_data[embryo == TRUE]`

- Free-swimming: `gulper_data[embryo == FALSE]`

- Females: `gulper_data[sex == 1]`

- Males: `gulper_data[sex == 2]`

- Mature females: `gulper_data[sex == 1 & mat == 1]`

## References

Moe, B.J. (unpublished). Life history of the little gulper shark
(*Centrophorus uyato*) in the Gulf of Mexico.

## Examples

``` r
data("gulper_data")

# Data overview
str(gulper_data)
#> Classes 'data.table' and 'data.frame':   668 obs. of  6 variables:
#>  $ sex   : chr  NA NA NA NA ...
#>  $ mat   : num  NA NA NA NA NA NA NA NA NA NA ...
#>  $ fl    : num  8.5 16.5 16.5 19.2 25.4 26.3 28.6 29.6 31.7 36 ...
#>  $ age1  : num  NA NA NA NA NA NA NA NA NA NA ...
#>  $ age2  : num  NA NA NA NA NA NA NA NA NA NA ...
#>  $ embryo: logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  - attr(*, ".internal.selfref")=<externalptr> 
#>  - attr(*, "index")= int(0) 
#>   ..- attr(*, "__embryo")= int [1:668] 13 14 15 16 17 18 19 20 21 22 ...

# Sample sizes by sex
gulper_data[embryo == FALSE, .N, by = sex]
#>       sex     N
#>    <char> <int>
#> 1:      1   252
#> 2:      2   404

# Maturity by sex (free-swimming only)
gulper_data[embryo == FALSE, .(
  n = .N,
  n_mature = sum(mat, na.rm = TRUE),
  prop_mature = mean(mat, na.rm = TRUE)
), by = sex]
#>       sex     n n_mature prop_mature
#>    <char> <int>    <num>       <num>
#> 1:      1   252       85   0.3632479
#> 2:      2   404      280   0.7161125

# Extract embryo and free-swimming lengths for birth model
embryo_fl <- gulper_data[embryo == TRUE, fl]
freeswim_fl <- gulper_data[embryo == FALSE, fl]

# Subset females for single-sex analysis
female_data <- gulper_data[sex == 1]
```

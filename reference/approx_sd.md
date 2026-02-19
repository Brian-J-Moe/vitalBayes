# Approximate Standard Deviation from Confidence Interval

Estimates standard deviation from reported confidence interval bounds
assuming a normal distribution.

## Usage

``` r
approx_sd(lower, upper, level = 0.95)
```

## Arguments

- lower:

  Lower bound of confidence interval.

- upper:

  Upper bound of confidence interval.

- level:

  Confidence level (default 0.95).

## Value

Estimated standard deviation.

## Examples

``` r
# If 95% CI is (10, 20), approximate SD
approx_sd(10, 20, 0.95)
#> [1] 2.551067
```

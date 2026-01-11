# Approximate a Standard Deviation from a Confidence Interval

Uses a normal approximation to convert confidence interval bounds to a
standard deviation: SD = (upper - lower) / (2 \* z), where z = qnorm(1 -
(1-CI)/2). Useful for extracting uncertainty from published confidence
intervals.

## Usage

``` r
approx_sd(lower, upper, CI = 0.95)
```

## Arguments

- lower:

  Numeric vector of lower CI bounds.

- upper:

  Numeric vector of upper CI bounds (same length as lower).

- CI:

  Scalar or numeric vector in (0, 1); default 0.95.

## Value

A numeric vector of approximate standard deviations.

## Examples

``` r
if (FALSE) { # \dontrun{
approx_sd(8, 12)
approx_sd(c(8, 10), c(12, 15), 0.90)
} # }
```

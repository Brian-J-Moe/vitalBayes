# Sample from Bivariate Normal with Specified Correlation

Generates correlated samples from two normal distributions.

## Usage

``` r
sample_bivariate_normal(n, mu1, sd1, mu2, sd2, rho, seed = NULL)
```

## Arguments

- n:

  Number of samples.

- mu1, mu2:

  Means.

- sd1, sd2:

  Standard deviations.

- rho:

  Correlation coefficient.

- seed:

  Random seed.

## Value

List with components x1 and x2.

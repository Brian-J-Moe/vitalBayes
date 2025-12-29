# Compare Pooling Effects

Compares parameter estimates between pooled and unpooled two-sex models
to assess the effect of partial pooling.

## Usage

``` r
compare_pooling(fit_pooled, fit_unpooled, params = NULL)
```

## Arguments

- fit_pooled:

  CmdStanMCMC object from a pooled two-sex model.

- fit_unpooled:

  CmdStanMCMC object from an unpooled two-sex model.

- params:

  Character vector of parameter names to compare. If `NULL`, uses
  c("Linf", "L0", "k").

## Value

A data.table with parameter comparisons.

## Examples

``` r
if (FALSE) { # \dontrun{
comp <- compare_pooling(growth_pooled, growth_unpooled)
print(comp)
} # }
```

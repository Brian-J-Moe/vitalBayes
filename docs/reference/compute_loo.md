# LOO Cross-Validation for vitalBayes Models

Computes approximate leave-one-out cross-validation using PSIS-LOO.
Requires the `loo` package.

## Usage

``` r
compute_loo(fit, ...)
```

## Arguments

- fit:

  A CmdStanMCMC object with log_lik generated quantities.

- ...:

  Additional arguments passed to
  [`loo`](https://mc-stan.org/loo/reference/loo.html).

## Value

A `loo` object.

## See also

[`vignette("model_diagnostics")`](https://brian-j-moe.github.io/vitalBayes/articles/model_diagnostics.md)
for comprehensive model comparison workflow.

[Statistical Methods: Model
Assessment](https://brian-j-moe.github.io/vitalBayes/doc/vitalBayes_stats_explained.html#assessment)
for LOO-CV theory and interpretation.

[`compare_loo`](https://brian-j-moe.github.io/vitalBayes/reference/compare_loo.md),
[`create_loo_table`](https://brian-j-moe.github.io/vitalBayes/reference/create_loo_table.md),
[`ppc_summary`](https://brian-j-moe.github.io/vitalBayes/reference/ppc_summary.md)

## Examples

``` r
if (FALSE) { # \dontrun{
loo_growth <- compute_loo(growth_fit)
print(loo_growth)
} # }
```

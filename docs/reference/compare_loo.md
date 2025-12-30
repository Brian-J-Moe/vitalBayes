# Compare Models Using LOO

Compares multiple models using LOO-CV. Returns a comparison table with
expected log predictive density (elpd) differences.

## Usage

``` r
compare_loo(..., criterion = c("loo", "waic"))
```

## Arguments

- ...:

  Named CmdStanMCMC objects or loo objects to compare.

- criterion:

  Character. Comparison criterion: `"loo"` (default) or `"waic"`.

## Value

A data.table with model comparison statistics.

## See also

[`vignette("model_diagnostics")`](https://brian-j-moe.github.io/vitalBayes/articles/model_diagnostics.md)
for comprehensive model comparison workflow.

[Statistical Methods: Model
Assessment](https://brian-j-moe.github.io/vitalBayes/doc/vitalBayes_stats_explained.html#assessment)
for interpreting elpd differences.

[`compute_loo`](https://brian-j-moe.github.io/vitalBayes/reference/compute_loo.md),
[`create_loo_table`](https://brian-j-moe.github.io/vitalBayes/reference/create_loo_table.md),
[`compare_growth_models`](https://brian-j-moe.github.io/vitalBayes/reference/compare_growth_models.md)

## Examples

``` r
if (FALSE) { # \dontrun{
comp <- compare_loo(
  "von Bertalanffy" = vb_fit,
  "Gompertz" = gomp_fit,
  "Logistic" = log_fit
)
print(comp)
} # }
```

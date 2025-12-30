# Posterior Predictive Summary

Generates a comprehensive summary of posterior predictive checks for any
vitalBayes model fit.

## Usage

``` r
ppc_summary(fit, model_type = NULL, data = NULL, ...)
```

## Arguments

- fit:

  A CmdStanMCMC object from any vitalBayes fitting function.

- model_type:

  Character. Type of model: `"growth"`, `"maturity"`, or `"birth"`. If
  `NULL`, auto-detected.

- data:

  Optional data.table with observed data for additional checks.

- ...:

  Additional arguments (currently unused).

## Value

A list with class `"vitalBayes_ppc"` containing diagnostic summaries.

## See also

[`vignette("model_diagnostics")`](https://brian-j-moe.github.io/vitalBayes/articles/model_diagnostics.md)
for comprehensive diagnostic workflow.

[Statistical Methods: Model
Assessment](https://brian-j-moe.github.io/vitalBayes/doc/vitalBayes_stats_explained.html#assessment)
for posterior predictive check theory.

[`compute_loo`](https://brian-j-moe.github.io/vitalBayes/reference/compute_loo.md),
[`compare_loo`](https://brian-j-moe.github.io/vitalBayes/reference/compare_loo.md),
[`plot_residuals`](https://brian-j-moe.github.io/vitalBayes/reference/plot_residuals.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ppc <- ppc_summary(growth_fit, model_type = "growth")
print(ppc)
} # }
```

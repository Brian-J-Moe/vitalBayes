# Extract Life History Parameters from vitalBayes Fits

Extracts von Bertalanffy growth parameters and maturity information from
vitalBayes model fits. Can return full posterior draws (preserving
correlations) or summary statistics.

## Usage

``` r
extract_lh_params(
  growth_fit = NULL,
  maturity_fit = NULL,
  sex = NULL,
  format = c("draws", "summary"),
  n_draws = NULL
)
```

## Arguments

- growth_fit:

  Optional CmdStanMCMC object from
  [`fit_bayesian_growth()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md).

- maturity_fit:

  Optional CmdStanMCMC object from
  [`fit_bayesian_maturity()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md)
  for age-at-maturity (t50).

- sex:

  Integer. Sex code (1 = female, 2 = male) for hierarchical models. If
  NULL and model is hierarchical, extracts both sexes.

- format:

  Character. Either `"draws"` for full posterior draws as data.table
  (default), or `"summary"` for mean and SD for each parameter.

- n_draws:

  Integer. Number of posterior draws to sample (if format = "draws"). If
  NULL, uses all draws. Default NULL.

## Value

If format = "draws": data.table with columns for each parameter and
optional `.draw` index and `sex` indicator. If format = "summary":
data.table with columns `parameter`, `mean`, `sd`, and optional `sex`.

## Details

This function preserves the joint posterior distribution including
correlations between parameters. This is statistically superior to
independent resampling from marginal `c(mean, sd)` specifications,
particularly when correlations are substantial (as is often the case for
Linf, k, and t0).

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract full draws for males
pars <- extract_lh_params(growth_fit, maturity_fit, sex = 2, format = "draws")

# Extract summaries for both sexes
pars <- extract_lh_params(growth_fit, maturity_fit, sex = NULL, format = "summary")
} # }
```

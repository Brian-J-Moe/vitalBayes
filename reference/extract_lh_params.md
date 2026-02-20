# Extract Life History Parameters from vitalBayes Fits or Manual Specification

Extracts life history parameters from a combination of vitalBayes
stanfit objects and/or manual `c(mean, sd)` specifications. When stanfit
objects are provided, posterior correlations are preserved. When
maturity parameters come from separate fits, optional correlation can be
specified.

## Usage

``` r
extract_lh_params(
  growth_fit = NULL,
  birth_fit = NULL,
  length_maturity_fit = NULL,
  age_maturity_fit = NULL,
  Linf = NULL,
  L0 = NULL,
  Lmat = NULL,
  tmat = NULL,
  maturity_cor = 0.5,
  sex = 1L,
  iter = 2000,
  seed = 1234,
  show_progress = TRUE
)
```

## Arguments

- growth_fit:

  CmdStanMCMC object from `fit_bayesian_growth`, or NULL.

- birth_fit:

  CmdStanMCMC object from `fit_bayesian_birth`, or NULL.

- length_maturity_fit:

  CmdStanMCMC object from `fit_bayesian_maturity` for
  length-at-maturity, or NULL.

- age_maturity_fit:

  CmdStanMCMC object from `fit_bayesian_maturity` for age-at-maturity,
  or NULL.

- Linf, L0, Lmat, tmat:

  Manual specification as `c(mean, sd)`, or NULL.

- maturity_cor:

  Correlation between Lmat and tmat when both are from manual
  specification or separate fits. Default 0.5 (positive correlation
  reflecting that larger individuals tend to mature later). Set to 0 for
  independent sampling, or NA to attempt estimation from data.

- sex:

  Integer. Sex code (1 = female, 2 = male) for extracting from two-sex
  models. Default 1.

- iter:

  Number of parameter sets to generate. Default 2000.

- seed:

  Random seed. Default 1234.

- show_progress:

  Logical. Show messages? Default TRUE.

## Value

A data.table with columns: set_id, Linf, L0, Lmat, tmat, plus attribute
"sources" indicating parameter origins.

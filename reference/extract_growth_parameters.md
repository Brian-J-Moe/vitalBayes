# Extract Life History Parameters from Growth Model Posterior

Extracts posterior draws of \\(L\_\infty, L_0, L\_{mat}, t\_{mat})\\
from a vitalBayes growth model fit. Works identically for von
Bertalanffy, Gompertz, and Logistic models fitted via
[`fit_bayesian_growth`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md).

## Usage

``` r
extract_growth_parameters(
  growth_fit,
  maturity_fit = NULL,
  sex = NULL,
  n_draws = NULL,
  seed = 1234
)
```

## Arguments

- growth_fit:

  A `CmdStanMCMC` object from
  [`fit_bayesian_growth`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md).

- maturity_fit:

  Optional `CmdStanMCMC` object from
  [`fit_bayesian_maturity`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md)
  providing age-at-maturity. Required for k-based growth fits if `tmat`
  is needed.

- sex:

  Integer. Sex code (1 = female, 2 = male) for hierarchical models. If
  `NULL`, extracts from single-sex model or uses column 1.

- n_draws:

  Integer. Number of posterior draws to return. If `NULL`, returns all
  available draws. If specified, draws are subsampled randomly.

- seed:

  Integer. Random seed for reproducible subsampling.

## Value

A `data.table` with columns:

- draw:

  Integer draw index

- Linf:

  Asymptotic length

- L0:

  Length at birth

- Lmat:

  Length at maturity (if available)

- tmat:

  Age at maturity (if available)

- k:

  Growth coefficient (native to the fitted model)

- k_vb_equiv:

  VB-equivalent k (always computed; see Details)

- growth_model:

  Growth model type: `"vb"`, `"gompertz"`, or `"logistic"`

The `growth_model` column is constant across all rows. The `k_vb_equiv`
derivation method is stored as an attribute `"k_vb_method"`: one of
`"milestone"`, `"derivative"`, or `"identity"`.

## Details

This function provides a unified interface for extracting the biological
parameters common to all growth models. These parameters — asymptotic
length, birth size, and maturity milestones — represent real biological
quantities that exist independently of the mathematical model used to
describe growth.

For maturity-based growth fits (`k_based = FALSE`), \\L\_{mat}\\ and
\\t\_{mat}\\ are directly estimated parameters. For k-based fits, these
must be supplied separately via `maturity_fit`.

The function also computes the VB-equivalent \\k\\ needed for mortality
models that were derived under VB assumptions (Chen-Watanabe and
growth-based Lorenzen). Three derivation strategies are available,
applied in order of preference:

1.  **Milestone-based**: If \\L\_{mat}\\ and \\t\_{mat}\\ are available
    (from maturity-based fits or a separate `maturity_fit`),
    \\k\_{VB}^{equiv}\\ is derived by inverting the VB equation through
    the biological milestones \\(0, L_0)\\ and \\(t\_{mat}, L\_{mat})\\.

2.  **Derivative matching**: For k-based Gompertz or Logistic fits
    without maturity data, \\k\_{VB}^{equiv}\\ is derived by matching
    the instantaneous growth rate at birth (\\t = 0\\).

3.  **Identity**: For VB fits, the native \\k\\ is itself the VB growth
    coefficient.

## Examples

``` r
if (FALSE) { # \dontrun{
# From a Gompertz fit with maturity-based parameterization
params <- extract_growth_parameters(gomp_fit, sex = 1, n_draws = 2000)

# Check VB-equivalent k distribution
hist(params$k_vb_equiv, main = "VB-Equivalent k from Gompertz Fit")

# Verify the growth model was auto-detected
unique(params$growth_model)
#> [1] "gompertz"

# Check which derivation method was used
attr(params, "k_vb_method")
#> [1] "milestone"
} # }
```

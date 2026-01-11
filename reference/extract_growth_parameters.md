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

  Growth coefficient (original model's k, not VB-equivalent)

- k_vb_equiv:

  VB-equivalent k (if Lmat and tmat available)

## Details

This function provides a unified interface for extracting the biological
parameters common to all growth models. These parameters — asymptotic
length, birth size, and maturity milestones — represent real biological
quantities that exist independently of the mathematical model used to
describe growth.

For maturity-based growth fits (`k_based = FALSE`), \\L\_{mat}\\ and
\\t\_{mat}\\ are directly estimated parameters. For k-based fits, these
must be supplied separately via `maturity_fit`.

## Examples

``` r
if (FALSE) { # \dontrun{
# From a Gompertz fit with maturity-based parameterization
params <- extract_growth_parameters(gomp_fit, sex = 1, n_draws = 2000)

# Check VB-equivalent k distribution
hist(params$k_vb_equiv, main = "VB-Equivalent k from Gompertz Fit")

# Compare to original Gompertz k
plot(params$k, params$k_vb_equiv,
     xlab = "Gompertz k", ylab = "VB-equivalent k")
} # }
```

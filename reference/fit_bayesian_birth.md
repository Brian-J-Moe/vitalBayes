# Fit Bayesian Length-at-Birth Model

Fits a binomial regression with probit link function for estimating the
length at which 50% of individuals are expected to transition from
embryo to free-swimming status (b50). The model is implemented in Stan
via precompiled models using the instantiate package.

The probit link is chosen based on a threshold-crossing interpretation
where latent developmental readiness (influenced by maternal size,
condition, temperature, etc.) is assumed normally distributed. This
contrasts with the logit link which assumes log-odds, which are less
interpretable in this biological context.

## Usage

``` r
fit_bayesian_birth(
  embryo_lts,
  free_swimming_lts,
  mean_b50 = NULL,
  cv_b50 = 0.3,
  mean_slope = 0,
  sd_slope = 1,
  parallel = TRUE,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 500,
  seed = 1234,
  ...
)
```

## Arguments

- embryo_lts:

  Numeric vector of embryo lengths (internally coded as status = 0).

- free_swimming_lts:

  Numeric vector of free-swimming lengths (internally coded as status =
  1).

- mean_b50:

  Numeric. Prior mean for b50. If `NULL` (default), computed as the
  midpoint between min(free_swimming_lts) and max(embryo_lts).

- cv_b50:

  Numeric. Coefficient of variation for b50 prior. Default 0.3.

- mean_slope:

  Numeric. Prior mean for slope on log scale. Default 0.

- sd_slope:

  Numeric. Prior SD for slope on log scale. Default 1.

- parallel:

  Logical. Run chains in parallel? Default `TRUE`.

- chains:

  Integer. Number of MCMC chains. Default 4.

- iter_warmup:

  Integer. Warmup iterations per chain. Default 1000.

- iter_sampling:

  Integer. Sampling iterations per chain. Default 1000.

- refresh:

  Integer. Progress update frequency. Set to 0 for no output. Default
  500.

- seed:

  Integer. Random seed for reproducibility. Default 1234.

- ...:

  Additional arguments passed to `$sample()`.

## Value

A `CmdStanMCMC` object containing:

- b50:

  Length at 50% birth probability

- slope:

  Transition steepness (probit scale)

- b05, b95:

  Lengths at 5% and 95% birth probability

- transition_width:

  Width of transition zone (b95 - b05)

- log_lik:

  Log-likelihood values for LOO-CV

- p_pred:

  Predicted probabilities

- status_rep:

  Posterior predictive replications

- mean_p_embryo, mean_p_freeswim:

  Mean predicted probabilities by group

- prop_correct_rep:

  Classification accuracy in replications

## Details

The linear predictor is parameterized directly in terms of b50:
\$\$\eta_i = \text{slope} \times (\text{length}\_i - b\_{50})\$\$

Parameters are estimated using lognormal priors to ensure positivity.

## Prior Specification

Priors are constructed using a coefficient of variation approach. The
prior mean for b50 is automatically computed as the midpoint between the
maximum embryo length and minimum free-swimming length (the transition
zone). The prior SD is then `mean_b50 * cv_b50`.

## Initialization

Initial values are set to the data-derived midpoint estimate for b50 and
a reasonable default for slope, ensuring chains start near the
high-probability region of the posterior.

## See also

[`vignette("fit_bayesian_birth")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_birth.md)
for usage examples with gulper shark data.

`vignette("complete_workflow")` for the full three-stage analysis
pipeline.

[Statistical Methods: Birth Size
Estimation](https://brian-j-moe.github.io/vitalBayes/doc/vitalBayes_stats_explained.html#birth)
for the mathematical derivation.

[Statistical Methods: Probit
Link](https://brian-j-moe.github.io/vitalBayes/doc/vitalBayes_stats_explained.html#probit)
for justification of the probit over logit link.

[`fit_bayesian_maturity`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md),
[`fit_bayesian_growth`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md)
for downstream models that use birth estimates as priors.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit birth model with data-derived priors
birth_fit <- fit_bayesian_birth(
  embryo_lts        = sharks[embryo == TRUE, length],
  free_swimming_lts = sharks[embryo == FALSE, length]
)

# View summary
birth_fit$summary(c("b50", "slope", "transition_width"))

# Use b50 posterior as L0 prior for growth model
growth_fit <- fit_bayesian_growth(
  ...,
  birth_stanfit = birth_fit
)
} # }
```

# Extract Posterior Draws from a Stanfit Object

Extracts posterior draws for a parameter, handling both single-sex and
two-sex models.

## Usage

``` r
extract_param_draws(fit, param, sex = NULL)
```

## Arguments

- fit:

  CmdStanMCMC object.

- param:

  Character. Parameter name.

- sex:

  Integer or NULL. Sex index for two-sex models.

## Value

Numeric vector of draws.

# Scale Mortality Schedule to Target

Rescales an age-specific mortality schedule so the mean mortality equals
a target value.

## Usage

``` r
scale_mortality(M, M_target = NULL, tmax = NULL, p = 0.001)
```

## Arguments

- M:

  Numeric vector of instantaneous mortality rates.

- M_target:

  Target mean mortality. Can be a numeric scalar, a function of tmax
  (e.g., `function(tmax) ...`), or NULL to derive from survival
  probability `p`.

- tmax:

  Maximum age (required if M_target is a function or NULL).

- p:

  Probability of surviving to tmax. Used only if M_target is NULL.
  Default 0.001 (0.1 percent survival to tmax).

## Value

Numeric vector of scaled mortality rates.

## Details

The scaling applies: \\M\_{scaled} = M\_{raw} / \bar{M}\_{raw} \times
M\_{target}\\

If `M_target = NULL`, it is derived as \\-\ln(p) / t\_{max}\\.

## Examples

``` r
if (FALSE) { # \dontrun{
M_raw <- M_chen_watanabe(0:30, Linf = 200, k = 0.15, t0 = -1, tmax = 30)
M_scaled <- scale_mortality(M_raw, M_target = 0.2)

# Using Hoenig-style target
hoenig_target <- function(tmax) 4.899 * tmax^(-0.916)
M_scaled <- scale_mortality(M_raw, M_target = hoenig_target, tmax = 30)
} # }
```

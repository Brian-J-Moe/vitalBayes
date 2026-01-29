# Scale Mortality Schedule to Target Mean

Rescales an age-specific mortality schedule so its mean equals a target
value derived from empirical relationships (e.g., Hoenig, Then et al.)
or survival probability constraints.

## Usage

``` r
scale_mortality(M, M_target = NULL, tmax = NULL, p = 0.001)
```

## Arguments

- M:

  Numeric vector of instantaneous mortality rates.

- M_target:

  Target mean mortality. Can be a numeric scalar (fixed target), a
  function of tmax (e.g., `function(tmax) 4.899 * tmax^(-0.916)`), or
  `NULL` to derive from survival probability `p`.

- tmax:

  Maximum age (required if `M_target` is a function or `NULL`).

- p:

  Probability of surviving to `tmax`. Used only if `M_target = NULL`.
  Default 0.001 (0.1% survival).

## Value

Numeric vector of scaled mortality rates (same length as `M`).

## Details

The scaling applies: \$\$M\_{scaled}(t) = M\_{raw}(t) \times
\frac{M\_{target}}{\bar{M}\_{raw}}\$\$

This preserves the *shape* of the age-specific mortality curve while
adjusting its overall level. Scaling is useful because theoretical
mortality models often produce absolute levels that don't match
empirical observations, but the relative age pattern may still be
informative.

## Examples

``` r
if (FALSE) { # \dontrun{
M_raw <- M_chen_watanabe_L0(0:30, Linf = 100, L0 = 25, k = 0.1,
                             two_phase = FALSE)

# Scale to fixed target
M_scaled <- scale_mortality(M_raw, M_target = 0.2)

# Scale using Then et al. (2015) relationship
then_2015 <- function(tmax) 4.899 * tmax^(-0.916)
M_scaled <- scale_mortality(M_raw, M_target = then_2015, tmax = 30)

# Scale to survival probability
M_scaled <- scale_mortality(M_raw, M_target = NULL, tmax = 30, p = 0.01)
} # }
```

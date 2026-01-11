# Scale Mortality Schedule to Target Mean

## Usage

``` r
scale_mortality(M, M_target = NULL, tmax = NULL, p = 0.001)
```

## Arguments

- M:

  Numeric vector of instantaneous mortality rates.

- M_target:

  Target mean mortality. Can be:

  - Numeric scalar (fixed target)

  - Function of tmax: `function(tmax) ...`

  - `NULL` to derive from survival probability `p`

- tmax:

  Maximum age (required if `M_target` is a function or `NULL`).

- p:

  Probability of surviving to `tmax`. Used only if `M_target = NULL`.
  Default 0.001 (0.1\\

Numeric vector of scaled mortality rates (same length as `M`). Rescales
an age-specific mortality schedule so its mean equals a target value
derived from empirical relationships (e.g., Hoenig, Then et al.) or
survival probability constraints. The scaling applies:
\$\$M\_{scaled}(t) = M\_{raw}(t) \times
\frac{M\_{target}}{\bar{M}\_{raw}}\$\$This preserves the *shape* of the
age-specific mortality curve while adjusting its overall level. Scaling
is useful because theoretical mortality models often produce absolute
levels that don't match empirical observations, but the relative age
pattern may still be informative.

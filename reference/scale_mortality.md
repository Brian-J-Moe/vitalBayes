# Scale Mortality Schedule to Target via Cumulative Hazard

## Usage

``` r
scale_mortality(M, age, M_target = NULL, tmax = NULL, p = 0.001)
```

## Arguments

- M:

  Numeric vector of instantaneous mortality rates.

- age:

  Numeric vector of ages corresponding to `M`. Required for trapezoidal
  integration. Must be the same length as `M` and sorted in ascending
  order.

- M_target:

  Target mean mortality. Can be a numeric scalar (fixed target), a
  function of tmax (e.g., `function(tmax) 4.899 * tmax^(-0.916)` for
  Then et al. 2015), or `NULL` to derive from survival probability `p`.

- tmax:

  Maximum age (required if `M_target` is a function or `NULL`).

- p:

  Probability of surviving to `tmax`. Used only if `M_target = NULL`.
  Default 0.001 (0.1\\

Numeric vector of scaled mortality rates (same length as `M`). Rescales
an age-specific mortality schedule by matching the cumulative hazard to
an empirical target (from Hoenig, Then et al., or a specified survival
probability). The scaling finds a proportional constant \\c\\ such that:
\$\$M\_{scaled}(t) = c \times M\_{raw}(t)\$\$The constant \\c\\ is
determined by matching the cumulative hazard to the target. The
cumulative hazard is computed via trapezoidal integration: \$\$H\_{raw}
= \int_0^{t\_{max}} M\_{raw}(a) \\ da \approx \sum_i \frac{\Delta
a_i}{2} \left\[ M(a_i) + M(a\_{i+1}) \right\]\$\$For the
**survival-probability target** (`M_target = NULL`): \$\$c =
\frac{-\ln(p)}{H\_{raw}}\$\$ This ensures \\S(t\_{max}) = \exp\\\left(-c
\cdot H\_{raw}\right) = p\\ exactly.For a **fixed mean-mortality
target** (numeric or function of `tmax`): \$\$c =
\frac{\bar{M}\_{target} \times t\_{max}}{H\_{raw}}\$\$ This interprets
the target as the average mortality over the lifespan, so the cumulative
hazard of the scaled schedule equals \\\bar{M}\_{target} \times
t\_{max}\\.This approach is preferred over arithmetic-mean-based scaling
because it directly constrains the biologically relevant quantity
(cumulative survival) rather than a proxy, and is invariant to the age
grid spacing.

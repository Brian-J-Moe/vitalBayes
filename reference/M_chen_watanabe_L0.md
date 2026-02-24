# Chen-Watanabe Natural Mortality (\\L_0\\ Parameterization)

Computes age-specific natural mortality using the Chen & Watanabe (1989)
model with an \\L_0\\ parameterization that eliminates dependence on the
theoretical parameter \\t_0\\.

## Usage

``` r
M_chen_watanabe_L0(
  age,
  Linf,
  L0,
  k,
  k_native = k,
  growth_model = c("vb", "gompertz", "logistic"),
  tmax = NULL,
  Linf_factor = 0.99,
  two_phase = TRUE,
  tmat = NULL,
  late_model = c("gompertz", "logistic"),
  tm_factor = 2/3,
  M_mult = 2,
  smooth_factor = 1/3
)
```

## Arguments

- age:

  Numeric vector of ages at which to compute mortality.

- Linf:

  Asymptotic length.

- L0:

  Length at birth.

- k:

  VB-equivalent growth coefficient (used as \\M\_\infty\\ in the CW
  formula). For VB fits, this is the native \\k\\. For Gompertz or
  Logistic fits, use `compute_k_vb_equivalent` or
  `.derive_k_vb_derivative`.

- k_native:

  Growth coefficient native to the fitted model (used for \\L(t)\\
  prediction). Defaults to `k`, which is correct for VB fits. For non-VB
  fits, pass the model-specific growth coefficient here.

- growth_model:

  Character. Growth model for \\L(t)\\ prediction: `"vb"` (default),
  `"gompertz"`, or `"logistic"`.

- tmax:

  Maximum age. If `NULL`, estimated from growth parameters as age when
  \\L(t) = \\ `Linf_factor` \\\times L\_\infty\\ (using the VB equation
  with `k` for consistency with the CW theoretical framework).

- Linf_factor:

  Numeric in (0, 1). Fraction of \\L\_\infty\\ used to estimate
  \\t\_{max}\\. Default 0.99.

- two_phase:

  Logical. If `TRUE`, applies two-phase model with late-life senescence.
  Default `TRUE`.

- tmat:

  Age at maturity. Required if `two_phase = TRUE`.

- late_model:

  Character. Senescence model: `"gompertz"` (default) or `"logistic"`.

- tm_factor:

  Numeric. Fraction of \\t\_{mat}\\ at which transition to senescence
  begins. Default 2/3.

- M_mult:

  Numeric. Multiplier for senescence mortality plateau relative to
  mortality at \\t_m\\. Default 2.

- smooth_factor:

  Numeric. Controls smoothness of transition between phases. Default
  1/3.

## Value

Numeric vector of instantaneous mortality rates (same length as `age`).

## Details

The standard Chen-Watanabe formulation expresses mortality as: \$\$M(t)
= \frac{k}{1 - e^{-k(t - t_0)}}\$\$

where \\t_0\\ is the theoretical age at length zero — a parameter with
no direct biological interpretation that can take implausible values,
particularly when growth data are sparse.

We reparameterize using the relationship between \\t_0\\ and \\L_0\\
(birth length) under von Bertalanffy dynamics: \$\$L_0 = L\_\infty(1 -
e^{kt_0})\$\$

After algebraic manipulation (see vignette), the \\L_0\\-parameterized
form becomes: \$\$M(t) = \frac{k\_{VB} \cdot L\_\infty}{L(t)}\$\$

where \\k\_{VB}\\ is the von Bertalanffy growth coefficient (or
VB-equivalent for non-VB fits) and \\L(t)\\ is the predicted length at
age \\t\\ from the native growth model.

## Two Roles of k

The Chen-Watanabe equation involves \\k\\ in two conceptually distinct
roles:

- Asymptotic mortality rate:

  The VB-specific constant \\M\_\infty = k\_{VB}\\ that governs the
  theoretical minimum mortality as \\t \to \infty\\. This enters the CW
  formula as the numerator. Because the derivation assumes VB dynamics,
  this must always be the VB growth coefficient (or VB-equivalent for
  non-VB fits). Passed via the `k` argument.

- Growth trajectory:

  The length-at-age prediction \\L(t)\\ that determines how far an
  individual is from asymptotic size. For non-VB growth models, the
  native trajectory with native \\k\\ produces more accurate body-size
  estimates than the VB approximation. Controlled by `k_native` and
  `growth_model`.

For VB fits, both roles use the same \\k\\, and the default
`k_native = k` preserves full backward compatibility.

## Two-Phase Extension

The original CW model produces unrealistic mortality trajectories at old
ages (approaching zero asymptotically). When `two_phase = TRUE`, the
model adds a senescence component where mortality increases after
maturity, more realistically capturing late-life dynamics. Mortality
follows the standard CW model until age \\t_m\\ (a fraction of
\\t\_{mat}\\), then transitions to a senescence model (Gompertz or
logistic) that increases mortality toward \\t\_{max}\\.

## References

Chen, S., & Watanabe, S. (1989). Age dependence of natural mortality
coefficient in fish population dynamics. *Nippon Suisan Gakkaishi*,
55(2), 205-208.

## See also

`compute_k_vb_equivalent` for deriving \\k\_{VB}\\ from maturity
milestones,
[`get_stochastic_mortality`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md)
for Monte Carlo mortality estimation with uncertainty.

## Examples

``` r
if (FALSE) { # \dontrun{
ages <- seq(0.1, 30, by = 0.5)

# VB fit (k serves both roles)
M_vb <- M_chen_watanabe_L0(
  age = ages, Linf = 100, L0 = 25, k = 0.1,
  two_phase = TRUE, tmat = 10
)

# Gompertz fit (separate k for each role)
M_gomp <- M_chen_watanabe_L0(
  age = ages, Linf = 100, L0 = 25,
  k = 0.08,                       # VB-equivalent k for M_inf
  k_native = 0.12,                # Native Gompertz k for L(t)
  growth_model = "gompertz",
  two_phase = TRUE, tmat = 10
)
} # }
```

# Chen-Watanabe Natural Mortality (L₀ Parameterization)

## Usage

``` r
M_chen_watanabe(
  age,
  Linf,
  L0,
  k,
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

  VB-equivalent growth coefficient. Can be computed from any growth
  model using
  [`compute_k_vb_equivalent`](https://brian-j-moe.github.io/vitalBayes/reference/compute_k_vb_equivalent.md).

- tmax:

  Maximum age. If `NULL`, estimated from growth parameters as age when
  \\L(t) = \\ `Linf_factor` \\\times L\_\infty\\.

- Linf_factor:

  Numeric in (0, 1). Fraction of \\L\_\infty\\ used to estimate
  \\t\_{max}\\. Default 0.99 (age at 99\\

  two_phaseLogical. If `TRUE`, applies two-phase model with late-life
  senescence. Default `TRUE`.

  tmatAge at maturity. Required if `two_phase = TRUE`.

  late_modelCharacter. Senescence model: `"gompertz"` (default) or
  `"logistic"`.

  tm_factorNumeric. Fraction of \\t\_{mat}\\ at which transition to
  senescence begins. Default 2/3.

  M_multNumeric. Multiplier for senescence mortality plateau relative to
  mortality at \\t_m\\. Default 2.

  smooth_factorNumeric. Controls smoothness of transition between
  phases. Default 1/3.

Numeric vector of instantaneous mortality rates (same length as `age`).
Computes age-specific natural mortality using the Chen & Watanabe (1989)
model with an \\L_0\\ parameterization that eliminates dependence on the
theoretical parameter \\t_0\\. The standard Chen-Watanabe formulation
expresses mortality as: \$\$M(t) = \frac{k}{1 - e^{-k(t -
t_0)}}\$\$where \\t_0\\ is the theoretical age at length zero — a
parameter with no direct biological interpretation that can take
implausible values, particularly when growth data are sparse.We
reparameterize using the relationship between \\t_0\\ and \\L_0\\ (birth
length) under von Bertalanffy dynamics: \$\$L_0 = L\_\infty(1 -
e^{kt_0})\$\$After algebraic manipulation (see vignette), the
\\L_0\\-parameterized form becomes: \$\$M(t) = \frac{k \cdot
L\_\infty}{L(t)}\$\$where \\L(t) = L\_\infty - (L\_\infty -
L_0)e^{-kt}\\ is the predicted length at age \\t\\. Mathematical
FoundationThis reformulation reveals that Chen-Watanabe mortality is
inversely proportional to body size — smaller (younger) individuals
experience higher mortality. The ratio \\L\_\infty / L(t)\\ represents
how far an individual is from asymptotic size, with mortality declining
as this ratio approaches 1.

Two-Phase ExtensionThe original CW model produces unrealistic mortality
trajectories at old ages (approaching zero asymptotically). The
two-phase extension adds a senescence component where mortality
increases after maturity, more realistically capturing late-life
dynamics.When `two_phase = TRUE`, mortality follows the standard CW
model until age \\t_m\\ (a fraction of \\t\_{mat}\\), then transitions
to a senescence model (Gompertz or logistic) that increases mortality
toward \\t\_{max}\\.

Chen, S., & Watanabe, S. (1989). Age dependence of natural mortality
coefficient in fish population dynamics. *Nippon Suisan Gakkaishi*,
55(2), 205-208.
[`compute_k_vb_equivalent`](https://brian-j-moe.github.io/vitalBayes/reference/compute_k_vb_equivalent.md)
for deriving \\k\\ from any growth model,
[`get_stochastic_mortality`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md)
for Monte Carlo mortality estimation with uncertainty.

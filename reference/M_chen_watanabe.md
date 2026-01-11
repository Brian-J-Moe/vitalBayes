# Chen-Watanabe Natural Mortality Model

Computes age-specific natural mortality using the Chen & Watanabe (1989)
model, with optional late-life mortality phase following Gompertz or
logistic senescence.

## Usage

``` r
M_chen_watanabe(
  age,
  Linf,
  k,
  t0,
  tmax,
  tmat = NULL,
  two_phase = TRUE,
  late_model = c("gompertz", "logistic"),
  tm_factor = 2/3,
  M_mult = 2,
  M_cap_factor = 1/4,
  smooth_factor = 1/3,
  mode = c("K_mult", "r_target"),
  alpha = 3,
  r_given = NULL
)
```

## Arguments

- age:

  Numeric vector of ages at which to compute mortality.

- Linf:

  Asymptotic length.

- k:

  von Bertalanffy growth coefficient.

- t0:

  Theoretical age at length zero.

- tmax:

  Maximum age (used for late-phase parameterization).

- tmat:

  Optional age at maturity (used to estimate transition age).

- two_phase:

  Logical. If TRUE, includes late-life senescence phase.

- late_model:

  Character. Senescence model: `"gompertz"` or `"logistic"`.

- tm_factor:

  Proportion of adult lifespan where late-phase transition occurs.
  Default 2/3.

- M_mult:

  Multiplier for mortality at tmax relative to M at transition age.

- M_cap_factor:

  Cap for M_tmax as proportion of maximum early-phase M.

- smooth_factor:

  Controls width of transition zone between phases.

- mode:

  Logistic mode: `"K_mult"` or `"r_target"`.

- alpha:

  Logistic K multiplier (if mode = "K_mult").

- r_given:

  Fixed logistic rate (if mode = "r_target").

## Value

Numeric vector of instantaneous mortality rates.

## Details

The Chen-Watanabe model expresses mortality as a function of the von
Bertalanffy growth parameters: \$\$M(t) = k / (1 - e^{-k(t - t_0)})\$\$

For the two-phase model, mortality transitions from the stable CW
formulation to a senescence phase (Gompertz or logistic) after a
transition age \\t_m\\.

## References

Chen, S., & Watanabe, S. (1989). Age dependence of natural mortality
coefficient in fish population dynamics. Nippon Suisan Gakkaishi, 55(2),
205-208.

## Examples

``` r
if (FALSE) { # \dontrun{
ages <- seq(0.1, 30, length.out = 100)
M <- M_chen_watanabe(ages, Linf = 200, k = 0.15, t0 = -1, tmax = 30)
plot(ages, M, type = "l")
} # }
```

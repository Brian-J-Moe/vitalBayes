# The Chen-Watanabe Reparameterization

## Overview

This vignette documents the mathematical framework underlying
vitalBayes’s implementation of the Chen & Watanabe (1989) natural
mortality model. Two methodological design choices distinguish this
implementation from standard practice:

1.  **\\L_0\\ parameterization**: Expressing the CW model directly in
    terms of predicted body length \\L(t)\\ and asymptotic length
    \\L\_\infty\\, eliminating dependence on the theoretical parameter
    \\t_0\\.

2.  **Growth-model-agnostic estimation**: Defining mortality through a
    normalized growth coefficient \\G(t) = L(t) / L\_\infty\\ and an
    asymptotic mortality rate \\M\_\infty\\ that is equivalent to the
    von Bertalanffy \\k\\, regardless of which growth model is used to
    describe the data.

Together, these choices mean that users can fit whichever growth model
best describes their data — von Bertalanffy, Gompertz, or Logistic — and
still use the Chen-Watanabe mortality framework without theoretical
compromise. For the mathematical details of the \\L_0\\-based growth
model equations themselves, see
[`vignette("fit_bayesian_growth")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_growth.md).

## The Chen-Watanabe Model

### \\L_0\\-Parameterized Form

The Chen & Watanabe (1989) model relates instantaneous natural mortality
\\M(t)\\ to the von Bertalanffy growth coefficient \\k\\. After
reparameterizing from \\t_0\\ to \\L_0\\ (the algebraic steps are
detailed in
[`vignette("fit_bayesian_growth")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_growth.md)),
the CW model simplifies to:

\\M(t) = \frac{k \cdot L\_\infty}{L(t)}\\

where \\L(t) = L\_\infty - (L\_\infty - L_0)e^{-kt}\\ is the von
Bertalanffy predicted length at age \\t\\, \\L\_\infty\\ is asymptotic
length, and \\k\\ is the VB growth coefficient.

This compact form reveals the core biological mechanism: **mortality is
inversely proportional to body size**. Small, young individuals are more
vulnerable to predation, starvation, and environmental stochasticity; as
an organism grows toward its asymptotic size, its instantaneous
mortality declines toward a lower bound.

### The Normalized Growth Coefficient \\G(t)\\

To generalize the CW model beyond von Bertalanffy growth, it is useful
to define the **normalized growth coefficient**:

\\G(t) = \frac{L(t)}{L\_\infty}\\

This is the fraction of asymptotic size achieved at age \\t\\. At birth
(\\t = 0\\), \\G(0) = L_0 / L\_\infty\\, which is some small fraction.
As \\t \to \infty\\, \\G(t) \to 1\\ for all three growth models — the
organism approaches its asymptotic size.

Using \\G(t)\\, the CW model becomes:

\\M(t) = \frac{M\_\infty}{G(t)}\\

where \\M\_\infty = k\\ is the **asymptotic mortality rate** — the
mortality rate an individual would experience at asymptotic length. The
entirety of the age-specific mortality pattern is captured by \\G(t)\\:
how quickly an organism reaches its full size determines how quickly its
mortality declines.

### Why \\M\_\infty\\ Equals the Von Bertalanffy \\k\\

The identification \\M\_\infty \equiv k\_{VB}\\ emerges directly from
the Chen-Watanabe derivation, which was built on von Bertalanffy growth
dynamics. As age increases, the VB growth equation approaches
\\L\_\infty\\ and the CW mortality equation approaches \\k\\:

\\\lim\_{t \to \infty} M(t) = \lim\_{t \to \infty} \frac{k \cdot
L\_\infty}{L\_\infty - (L\_\infty - L_0)e^{-kt}} = k\\

This makes biological sense. The von Bertalanffy \\k\\ governs the rate
at which an organism approaches its asymptotic size; organisms with
higher \\k\\ grow faster but have shorter lifespans and higher baseline
mortality (Beverton & Holt, 1959; Charnov, 1993). The CW model
formalizes this life history trade-off: \\k\\ simultaneously determines
the growth rate (how fast \\G(t)\\ approaches 1) and the floor of the
mortality curve.

Critically, **this relationship holds regardless of which growth model
was used to fit the data**. The Gompertz and Logistic models have their
own growth coefficients \\k_g\\ and \\k_l\\, but these are *not*
\\M\_\infty\\. The CW mortality model was derived under VB assumptions,
and the parameter \\k\\ in the CW equation is specifically the VB growth
coefficient. When we use Gompertz or Logistic fits, we do not substitute
their native \\k\\ into the CW equation — instead, we compute the
VB-equivalent \\k\\ that encodes the same biological growth information
(described below).

## Growth-Specific \\G(t)\\ Functions

Although \\M\_\infty = k\_{VB}\\ is common to all growth models, the
shape of the mortality curve is governed by \\G(t)\\, which differs
across models. Each growth model implies a different trajectory toward
asymptotic size, and therefore a different pattern of declining
mortality. Below are the normalized growth coefficients for each
\\L_0\\-based growth model implemented in vitalBayes.

### Von Bertalanffy

\\L\_{VB}(t) = L\_\infty - (L\_\infty - L_0) \\ e^{-k\_{VB} \\ t}\\

\\G\_{VB}(t) = 1 - \left(1 - \frac{L_0}{L\_\infty}\right) e^{-k\_{VB} \\
t}\\

The VB model describes growth as a constant exponential approach to
\\L\_\infty\\. The absolute growth rate \\dL/dt\\ is highest at birth
and declines monotonically. This produces a mortality curve that is
steepest in early life and decelerates smoothly — a pattern consistent
with strong predation pressure on neonates that diminishes as
individuals grow through vulnerable size classes.

### Gompertz

\\L\_{Gomp}(t) = L\_\infty
\exp\\\left\[-\ln\\\left(\frac{L\_\infty}{L_0}\right) e^{-k_g \\
t}\right\]\\

\\G\_{Gomp}(t) = \exp\\\left\[-\ln\\\left(\frac{L\_\infty}{L_0}\right)
e^{-k_g \\ t}\right\]\\

The Gompertz model describes growth where the *rate of deceleration*
itself decelerates — an exponentially declining growth rate rather than
a linearly declining one. Growth is initially rapid but slows earlier
and more abruptly than under VB dynamics. The inflection point occurs at
\\L\_\infty / e \approx 0.368 \\ L\_\infty\\, much earlier than the VB
inflection at \\0.632 \\ L\_\infty\\. In terms of mortality, the
Gompertz \\G(t)\\ rises faster in early life than \\G\_{VB}(t)\\,
producing a steeper initial decline in mortality but a slower final
approach to \\M\_\infty\\.

### Logistic

\\L\_{Log}(t) = \frac{L\_\infty}{1 + \left(\frac{L\_\infty}{L_0} -
1\right) e^{-k_l \\ t}}\\

\\G\_{Log}(t) = \frac{1}{1 + \left(\frac{L\_\infty}{L_0} - 1\right)
e^{-k_l \\ t}}\\

The Logistic model is symmetric around its inflection point at
\\L\_\infty / 2\\, with growth accelerating before the midpoint and
decelerating after. This produces the most distinct mortality pattern:
mortality declines slowly in very early life (when growth is still
accelerating), then drops steeply as the organism passes through the
inflection point, before leveling off toward \\M\_\infty\\. This pattern
may suit species where size-dependent mortality relief is concentrated
in a relatively narrow size window.

### Comparing \\G(t)\\ Across Models

The three \\G(t)\\ functions all share the same boundary conditions
(\\G(0) = L_0/L\_\infty\\ and \\G(\infty) = 1\\), but they differ in how
they traverse the interval between these bounds. Because \\M(t) =
M\_\infty / G(t)\\, the shape of \\G(t)\\ directly determines the
mortality schedule:

- **VB**: Steady exponential approach. Mortality declines fastest at
  young ages and progressively slows.
- **Gompertz**: Rapid early approach, slower later convergence.
  Mortality drops quickly in early life but the final decline toward
  \\M\_\infty\\ is drawn out.
- **Logistic**: Sigmoidal approach. Mortality is relatively flat in very
  early life, then drops sharply around the growth inflection, and
  converges toward \\M\_\infty\\.

In practice, for many elasmobranch species and typical ratios of
\\L_0/L\_\infty\\ (roughly 0.15–0.40), the differences in \\G(t)\\
between models are modest provided the models fit the data comparably
well. The choice of growth model has a larger effect on the early-age
mortality curve (where \\G(t)\\ values are small and thus \\M(t)\\
values are large) than on mid- to late-life mortality.

``` r
library(ggplot2)
library(data.table)

# Example parameters
Linf <- 100; L0 <- 25
k_vb <- 0.08; k_g <- 0.12; k_l <- 0.15
ages <- seq(0, 40, by = 0.5)

# Compute G(t) for each model
G_vb   <- 1 - (1 - L0/Linf) * exp(-k_vb * ages)
G_gomp <- exp(-log(Linf/L0) * exp(-k_g * ages))
G_log  <- 1 / (1 + (Linf/L0 - 1) * exp(-k_l * ages))

dt <- data.table(
  age  = rep(ages, 3),
  G    = c(G_vb, G_gomp, G_log),
  model = rep(c("von Bertalanffy", "Gompertz", "Logistic"), each = length(ages))
)

ggplot(dt, aes(x = age, y = G, color = model)) +
  geom_line(linewidth = 1) +
  labs(x = "Age (years)", y = expression(G(t) == L(t)/L[infinity]),
       title = "Normalized Growth Coefficient by Model") +
  theme_bw() +
  theme(legend.position = "top")
```

## Computing VB-Equivalent \\k\\ from Any Growth Model

### The Biological Milestones Approach

Since \\M\_\infty = k\_{VB}\\, users who fit a Gompertz or Logistic
model need a way to obtain the VB-equivalent \\k\\ without fitting a
separate VB model. The solution exploits the fact that all three growth
models estimate the same biological quantities — \\L\_\infty\\, \\L_0\\,
\\L\_{mat}\\, and \\t\_{mat}\\ — regardless of their functional form. A
shark has a true birth size, a true size at maturity, and a true
asymptotic size; these do not change based on which equation we use to
model growth.

Given these shared biological milestones, the VB-equivalent \\k\\ is the
growth coefficient that would produce a von Bertalanffy curve passing
through \\(0, L_0)\\ and \\(t\_{mat}, L\_{mat})\\ with asymptote
\\L\_\infty\\. Substituting \\L(t\_{mat}) = L\_{mat}\\ into the VB
equation and solving:

\\k\_{VB}^{equiv} = \frac{1}{t\_{mat}} \ln\\\left(\frac{L\_\infty -
L_0}{L\_\infty - L\_{mat}}\right)\\

This formula is applied draw-by-draw to the posterior of any growth
model. When the growth model *is* von Bertalanffy with maturity-based
parameterization, \\k\_{VB}^{equiv}\\ exactly reproduces the fitted
\\k\\ by construction. For Gompertz and Logistic fits,
\\k\_{VB}^{equiv}\\ will differ from their native \\k\\ values, since
the three models define \\k\\ differently. But \\k\_{VB}^{equiv}\\ is
the correct input for the CW mortality model regardless.

In vitalBayes, `compute_k_vb_equivalent()` performs this computation and
[`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md)
calls it automatically when processing posteriors from non-VB growth
fits.

### Verification: Consistency When Using VB

When the underlying growth model is von Bertalanffy with maturity-based
parameterization, \\k\_{VB}^{equiv}\\ is algebraically identical to the
estimated \\k\\. This serves as an internal consistency check:

``` r
library(vitalBayes)

# Extract VB posterior draws
params <- extract_growth_parameters(vb_fit, sex = 1, n_draws = 2000)

# Compare native k to VB-equivalent k
all.equal(params$k, params$k_vb_equiv, tolerance = 1e-10)
#> [1] TRUE
```

For Gompertz and Logistic fits, the VB-equivalent \\k\\ and native \\k\\
will be correlated but not identical:

``` r
params_gomp <- extract_growth_parameters(gomp_fit, sex = 1, n_draws = 2000)

# Native Gompertz k vs VB-equivalent k
cor(params_gomp$k, params_gomp$k_vb_equiv, use = "complete.obs")
#> [1] 0.89  # Correlated but distinct
```

## Late-Life Senescence Modifications

### The Problem with Asymptotic Mortality

The single-phase CW model predicts that mortality monotonically
decreases toward \\M\_\infty\\ as \\G(t) \to 1\\. This is biologically
unrealistic for long-lived species: senescence, accumulated
physiological damage, and declining immune function all contribute to
increasing mortality at advanced ages. The monotonically declining CW
curve captures juvenile and sub-adult mortality patterns well, but fails
to represent late-life dynamics.

### Two-Phase Extension

vitalBayes implements a two-phase CW model that blends the standard
declining-mortality trajectory with a senescence component. The
transition occurs at a user-defined age \\t_m\\, typically set as a
fraction of the age at maturity (\\t_m = \frac{2}{3} t\_{mat}\\ by
default). Before \\t_m\\, mortality follows the standard CW model. After
\\t_m\\, a senescence function gradually increases mortality toward
\\t\_{max}\\.

Two senescence models are available:

**Gompertz senescence** (`late_model = "gompertz"`):

\\M\_{late}(t) = M_s \exp\\\left\[r \cdot (t - t\_{max})\right\]\\

where \\M_s = c \cdot M(t_m)\\ is the mortality plateau at \\t\_{max}\\
(with \\c\\ = `M_mult`, default 2), and \\r\\ is solved from continuity
at \\t_m\\:

\\r = \frac{\ln\left(M(t_m) / M_s\right)}{t_m - t\_{max}}\\

This produces an accelerating increase in mortality — modest at first,
then more pronounced as the organism approaches maximum age. This
pattern aligns with actuarial models of biological senescence in which
the hazard function increases exponentially (Gompertz, 1825).

**Logistic senescence** (`late_model = "logistic"`):

\\M\_{late}(t) = \frac{K}{1 + \exp\\\left\[-r \cdot (t -
t\_{max})\right\]}\\

where \\K = 2 \\ c \cdot M(t_m)\\ and \\r\\ is again solved from
continuity at \\t_m\\.

Logistic senescence produces a bounded mortality increase that plateaus
at \\K\\ rather than growing without limit. This may be more appropriate
when there is reason to believe that late-life mortality stabilizes
rather than continuing to accelerate — for example, in populations
subject to a physiological upper bound on hazard rates.

### Smooth Blending

The transition between the early-phase CW model and the late-phase
senescence model is smoothed using a logistic weighting function:

\\w(t) = \frac{1}{1 + \exp\\\left\[-\frac{2(t - t_m)}{s \cdot
\|t\_{max} - t_m\|}\right\]}\\

where \\s\\ is the `smooth_factor` (default \\\frac{1}{3}\\). The
blended mortality at any age is:

\\M(t) = (1 - w) \cdot M\_{CW}(t) + w \cdot M\_{late}(t)\\

This ensures a smooth, differentiable transition rather than a sharp
kink at \\t_m\\.

``` r
ages <- seq(0.1, 35, by = 0.25)

# Single-phase: monotonic decline
M_single <- M_chen_watanabe_L0(
  age = ages, Linf = 100, L0 = 25, k = 0.08,
  two_phase = FALSE
)

# Two-phase: Gompertz senescence
M_gomp <- M_chen_watanabe_L0(
  age = ages, Linf = 100, L0 = 25, k = 0.08,
  two_phase = TRUE, tmat = 12, late_model = "gompertz"
)

# Two-phase: Logistic senescence
M_logis <- M_chen_watanabe_L0(
  age = ages, Linf = 100, L0 = 25, k = 0.08,
  two_phase = TRUE, tmat = 12, late_model = "logistic"
)

dt <- data.table(
  age = rep(ages, 3),
  M   = c(M_single, M_gomp, M_logis),
  Model = rep(c("Single-phase", "Gompertz senescence", "Logistic senescence"),
              each = length(ages))
)

ggplot(dt, aes(x = age, y = M, color = Model)) +
  geom_line(linewidth = 1) +
  labs(x = "Age (years)", y = "Instantaneous Mortality (M)",
       title = "Chen-Watanabe Mortality with Senescence Modifications") +
  theme_bw() +
  theme(legend.position = "top")
```

## Practical Implications

### Using the Best Growth Model for Mortality

The growth-model-agnostic framework means researchers no longer face a
choice between statistical fit and mortality estimation capability. In
elasmobranch datasets with sparse adult observations, the VB model
frequently produces biologically implausible \\L\_\infty\\ estimates
(sometimes below observed maximum lengths) and unstable \\k\\ values
driven by the strong \\L\_\infty\\-\\k\\ correlation (Katsanevakis &
Maravelias, 2008). Gompertz and Logistic models often provide more
reliable fits under these conditions.

With the reparameterization, the workflow is straightforward: fit all
candidate growth models, select the best model by LOO-CV or other
criteria, and pass that model’s posterior directly to
[`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md).
The function handles the VB-equivalent \\k\\ computation automatically.

### When Growth Models Agree and Disagree

If all three growth models fit the data comparably, their CW mortality
estimates will be similar — the biological milestones \\(L\_\infty, L_0,
L\_{mat}, t\_{mat})\\ are well-constrained and produce nearly identical
\\k\_{VB}^{equiv}\\ values. Divergence between mortality estimates from
different growth models indicates genuine model uncertainty: the data do
not strongly favor one growth trajectory, and the implied mortality
schedules differ as a result. In such cases, reporting sensitivity to
growth model choice is informative.

### Relationship to Other Mortality Models

The \\M(t) = M\_\infty / G(t)\\ formulation is specific to the
Chen-Watanabe model. The Peterson-Wroblewski and Lorenzen mortality
models in vitalBayes use different theoretical foundations (allometric
scaling of metabolic rate with body weight) and do not share this
\\G(t)\\ structure. See
[`vignette("mortality_estimation")`](https://brian-j-moe.github.io/vitalBayes/articles/mortality_estimation.md)
for equations and practical guidance on all three mortality models.

## References

Beverton, R. J. H., & Holt, S. J. (1959). A review of the lifespans and
mortality rates of fish in nature, and their relation to growth and
other physiological characteristics. In G. E. W. Wolstenholme & M.
O’Connor (Eds.), *CIBA Foundation Colloquia on Ageing* (Vol. 5,
pp. 142–180). Churchill.

Charnov, E. L. (1993). *Life History Invariants: Some Explorations of
Symmetry in Evolutionary Ecology*. Oxford University Press.

Chen, S., & Watanabe, S. (1989). Age dependence of natural mortality
coefficient in fish population dynamics. *Nippon Suisan Gakkaishi*,
55(2), 205–208.

Gompertz, B. (1825). On the nature of the function expressive of the law
of human mortality, and on a new mode of determining the value of life
contingencies. *Philosophical Transactions of the Royal Society of
London*, 115, 513–583.

Katsanevakis, S., & Maravelias, C. D. (2008). Modelling fish growth:
multi-model inference as a better alternative to a priori using von
Bertalanffy equation. *Fish and Fisheries*, 9(2), 178–187.

## See Also

- [`vignette("fit_bayesian_growth")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_growth.md)
  — Growth model equations and the \\L_0\\ parameterization
- [`vignette("mortality_estimation")`](https://brian-j-moe.github.io/vitalBayes/articles/mortality_estimation.md)
  — Stochastic mortality estimation with
  [`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md)
- [`vignette("survivorship_simulation")`](https://brian-j-moe.github.io/vitalBayes/articles/survivorship_simulation.md)
  — Cohort survival analysis using mortality schedules

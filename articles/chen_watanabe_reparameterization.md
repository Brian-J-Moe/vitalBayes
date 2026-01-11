# Chen-Watanabe Reparameterization: From t₀ to L₀

## Overview

This vignette documents two methodological innovations in vitalBayes’s
approach to natural mortality estimation:

1.  **L₀ parameterization of Chen-Watanabe**: Replacing the theoretical
    parameter \\t_0\\ with the observable quantity \\L_0\\ (birth
    length)

2.  **Growth-model-agnostic k derivation**: Computing a VB-equivalent
    growth coefficient from biological milestones, enabling
    Chen-Watanabe mortality estimation from *any* growth model fit (von
    Bertalanffy, Gompertz, or Logistic)

These innovations align mortality estimation with vitalBayes’s core
philosophy: parameterize models using biologically meaningful,
observable quantities rather than abstract mathematical constructs.

## Motivation

### The Problem with \\t_0\\

The Chen-Watanabe (1989) natural mortality model was derived under von
Bertalanffy growth dynamics and is traditionally expressed as:

\\M(t) = \frac{k}{1 - e^{-k(t - t_0)}}\\

where \\k\\ is the VB growth coefficient and \\t_0\\ is the “theoretical
age at length zero” — the age at which an organism would have zero
length if VB growth extended backward in time.

The parameter \\t_0\\ presents several challenges:

**Biological ambiguity.** No fish has zero length at any age. The
parameter is a mathematical artifact of extrapolating the VB curve
backward, not a biological reality. Values can be positive, negative, or
zero, and interpretation varies.

**Estimation instability.** When growth data are sparse at young ages
(common in elasmobranch research), \\t_0\\ estimates become unreliable.
Strong correlations with other parameters (\\k\\, \\L\_\infty\\) make it
poorly identified.

**Awkward uncertainty propagation.** Specifying priors on \\t_0\\ is
unintuitive. What prior belief does one have about the “age at length
zero”? This contrasts sharply with birth length (\\L_0\\), where
researchers often have strong prior information from embryo
measurements.

### The Problem with Growth Model Dependence

The standard Chen-Watanabe formulation requires von Bertalanffy
parameters specifically. This creates a practical problem: **the VB
model often performs poorly on elasmobranch data**.

In datasets with sparse adult observations — typical for slow-growing,
long-lived species — the VB model frequently produces:

- Biologically implausible \\L\_\infty\\ estimates (sometimes below
  observed maximum lengths)
- Unstable \\k\\ estimates driven by the \\L\_\infty\\-\\k\\ correlation
- Poor predictive performance compared to Gompertz or Logistic
  alternatives

Researchers face an uncomfortable choice: use the best-fitting growth
model for biological inference but abandon Chen-Watanabe mortality
estimation, or force a VB fit specifically for mortality calculation
despite knowing it describes the data poorly.

## Mathematical Development

### Reparameterization: \\t_0 \rightarrow L_0\\

We begin with the relationship between \\t_0\\ and \\L_0\\ under von
Bertalanffy dynamics. At age \\t = 0\\:

\\L_0 = L\_\infty \left(1 - e^{kt_0}\right)\\

Solving for the exponential term:

\\e^{kt_0} = 1 - \frac{L_0}{L\_\infty} = \frac{L\_\infty -
L_0}{L\_\infty}\\

Now consider the denominator of the Chen-Watanabe equation at age \\t\\:

\\1 - e^{-k(t - t_0)} = 1 - e^{-kt} \cdot e^{kt_0}\\

Substituting our expression for \\e^{kt_0}\\:

\\1 - e^{-k(t-t_0)} = 1 - e^{-kt} \cdot \frac{L\_\infty -
L_0}{L\_\infty}\\

\\= \frac{L\_\infty - (L\_\infty - L_0)e^{-kt}}{L\_\infty}\\

The numerator is precisely the von Bertalanffy growth equation:

\\L(t) = L\_\infty - (L\_\infty - L_0)e^{-kt}\\

Therefore:

\\1 - e^{-k(t-t_0)} = \frac{L(t)}{L\_\infty}\\

Substituting back into the Chen-Watanabe equation:

\\\boxed{M(t) = \frac{k \cdot L\_\infty}{L(t)}}\\

This is the **L₀-parameterized Chen-Watanabe model**. It expresses
mortality as inversely proportional to body size, with no dependence on
\\t_0\\.

### Biological Interpretation

The reparameterized form reveals the biological mechanism underlying
Chen-Watanabe mortality: smaller individuals experience higher
mortality, and this relationship scales with how close an organism is to
its asymptotic size.

The ratio \\L\_\infty / L(t)\\ can be interpreted as a “size deficit” —
how far the organism is from its maximum potential size. At birth (\\t =
0\\), this ratio is \\L\_\infty / L_0\\, its maximum value. As the
organism grows, the ratio approaches 1, and mortality approaches \\k\\.

This interpretation aligns with life history theory: juvenile mortality
is dominated by predation (size-dependent), while adult mortality
approaches a baseline rate related to senescence and physiological
constraints.

### Growth-Model-Agnostic k: The Key Insight

The L₀-parameterized CW model still requires a von Bertalanffy \\k\\.
But here’s the crucial observation: **all three growth models estimate
the same biological quantities**, just with different functional forms.

Regardless of whether you fit von Bertalanffy, Gompertz, or Logistic,
your growth model produces posterior distributions for:

- \\L\_\infty\\ — asymptotic length
- \\L_0\\ — length at birth  
- \\L\_{mat}\\ — length at maturity (in maturity-based parameterization)
- \\t\_{mat}\\ — age at maturity

These are real biological quantities that exist independently of the
mathematical model. A shark has some true birth size and some true
size/age at maturity — these don’t change based on which equation we use
to describe growth.

Given these biological milestones, we can ask: **what von Bertalanffy
\\k\\ would produce a growth curve passing through these points?**

The VB model with \\L_0\\ parameterization is:

\\L(t) = L\_\infty - (L\_\infty - L_0)e^{-kt}\\

At maturity (\\t = t\_{mat}\\, \\L = L\_{mat}\\):

\\L\_{mat} = L\_\infty - (L\_\infty - L_0)e^{-k \cdot t\_{mat}}\\

Solving for \\k\\:

\\e^{-k \cdot t\_{mat}} = \frac{L\_\infty - L\_{mat}}{L\_\infty - L_0}\\

\\-k \cdot t\_{mat} = \ln\left(\frac{L\_\infty - L\_{mat}}{L\_\infty -
L_0}\right)\\

\\\boxed{k\_{VB}^{equiv} = \frac{1}{t\_{mat}} \ln\left(\frac{L\_\infty -
L_0}{L\_\infty - L\_{mat}}\right)}\\

This **VB-equivalent \\k\\** can be computed from the posterior of *any*
growth model, as long as that model provides \\(L\_\infty, L_0,
L\_{mat}, t\_{mat})\\.

### Verification: VB Consistency

When the growth model *is* von Bertalanffy with maturity-based
parameterization, the VB-equivalent \\k\\ computed from this formula
exactly equals the \\k\\ estimated by Stan. This must be true by
construction — the maturity-based VB model derives \\k\\ using this same
relationship.

For Gompertz and Logistic fits, the VB-equivalent \\k\\ will differ from
the native \\k\\ of those models. This is expected: the three models use
different functional forms, so their \\k\\ parameters have different
meanings. But the VB-equivalent \\k\\ captures the “growth rate
information” in a form that Chen-Watanabe can use.

## Implementation in vitalBayes

### Computing VB-Equivalent k

``` r
library(vitalBayes)
library(data.table)

# Fit a Gompertz growth model (assume it fits your data better than VB)
gomp_fit <- fit_bayesian_growth(
  lt        = "fl",
  age       = "age", 
  sex       = "sex",
  data      = growth_data[embryo == FALSE & !is.na(age)],
  model     = "gompertz",
  k_based   = FALSE,  # Maturity-based parameterization
  birth_fit = birth_fit,
  L50_fit   = L50_fit,
  t50_fit   = t50_fit
)

# Extract biological parameters from Gompertz posterior
params <- extract_growth_parameters(gomp_fit, sex = 1, n_draws = 2000)

# Examine the relationship between native Gompertz k and VB-equivalent k
head(params[, .(k_gompertz = k, k_vb_equiv)])
#>    k_gompertz k_vb_equiv
#> 1:     0.0823     0.0945
#> 2:     0.0791     0.0912
#> 3:     0.0856     0.0987
#> ...

# They're correlated but not identical
cor(params$k, params$k_vb_equiv, use = "complete.obs")
#> [1] 0.89
```

### Mortality Estimation from Any Growth Model

``` r
# Chen-Watanabe mortality from a Gompertz fit — this now works!
mort <- get_stochastic_mortality(
  method     = "CW",
  growth_fit = gomp_fit,  # Gompertz, not VB
  sex        = 1,
  iter       = 2000,
  scaled     = TRUE,
  p          = 0.001
)

# The function automatically:
# 1. Extracts (Linf, L0, Lmat, tmat) from the Gompertz posterior
# 2. Computes VB-equivalent k for each posterior draw
# 3. Uses the L0-parameterized CW model

# Check what k values were used
summary(mort$Parameters$k_vb_equiv)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.0712  0.0891  0.0943  0.0951  0.1008  0.1234
```

### Comparing Mortality from Different Growth Models

A useful diagnostic is comparing mortality estimates derived from
different growth models. If the models are all reasonable fits to the
data, the mortality estimates should be similar:

``` r
# Fit all three growth models
vb_fit <- fit_bayesian_growth(
  lt = "fl", age = "age", sex = "sex",
  data = growth_data[embryo == FALSE & !is.na(age)],
  model = "vb", k_based = FALSE,
  birth_fit = birth_fit, L50_fit = L50_fit, t50_fit = t50_fit
)

gomp_fit <- fit_bayesian_growth(
  lt = "fl", age = "age", sex = "sex",
  data = growth_data[embryo == FALSE & !is.na(age)],
  model = "gompertz", k_based = FALSE,
  birth_fit = birth_fit, L50_fit = L50_fit, t50_fit = t50_fit
)

logis_fit <- fit_bayesian_growth(
  lt = "fl", age = "age", sex = "sex",
  data = growth_data[embryo == FALSE & !is.na(age)],
  model = "logistic", k_based = FALSE,
  birth_fit = birth_fit, L50_fit = L50_fit, t50_fit = t50_fit
)

# Mortality from each
mort_vb <- get_stochastic_mortality(
  method = "CW", growth_fit = vb_fit, sex = 1, 
  print_plot = FALSE
)
mort_gomp <- get_stochastic_mortality(
  method = "CW", growth_fit = gomp_fit, sex = 1,
  print_plot = FALSE
)
mort_logis <- get_stochastic_mortality(
  method = "CW", growth_fit = logis_fit, sex = 1,
  print_plot = FALSE
)

# Combine for comparison
combined <- rbind(
  mort_vb$Summary[, model := "von Bertalanffy"],
  mort_gomp$Summary[, model := "Gompertz"],
  mort_logis$Summary[, model := "Logistic"]
)

# Plot
library(ggplot2)
ggplot(combined, aes(x = age_round, color = model, fill = model)) +
  geom_ribbon(aes(ymin = M_lower, ymax = M_upper), alpha = 0.2, color = NA) +
  geom_line(aes(y = M_median), linewidth = 1) +
  labs(
    x = "Age (years)",
    y = "Natural Mortality (M)",
    title = "CW Mortality Estimates by Growth Model",
    subtitle = "All derived from same biological milestones"
  ) +
  theme_bw() +
  theme(legend.position = "top")
```

If the three curves are similar, this suggests the mortality estimates
are robust to growth model choice — the biological information
constrains the result regardless of functional form. If they diverge
substantially, this reveals genuine model uncertainty that should be
acknowledged.

## Why This Matters

### Practical Benefits

**Use the best model.** Researchers can select the growth model that
best fits their data (by LOO-CV or other criteria) without sacrificing
mortality estimation capability.

**Biologically meaningful priors.** Specifying priors on \\L_0\\ (birth
length) is straightforward — researchers often have direct measurements
from embryos or neonates. Priors on \\t_0\\ were always awkward.

**Interpretable parameters.** The L₀-parameterized CW model makes the
size-mortality relationship explicit. The mortality at any age can be
understood as a function of current body size relative to asymptotic
size.

**Consistent workflow.** The maturity-based parameterization philosophy
extends naturally from growth to mortality. All models in the vitalBayes
workflow use the same biological milestones.

### Theoretical Coherence

A legitimate concern is whether using a non-VB growth model with
Chen-Watanabe mortality is theoretically sound. After all, CW was
*derived* under VB assumptions.

Our position is that the VB-equivalent \\k\\ approach maintains
theoretical coherence because:

1.  **CW is fundamentally about size-dependent mortality.** The
    biological mechanism — smaller individuals die faster — doesn’t
    depend on the VB functional form. It depends on the growth
    trajectory.

2.  **The biological milestones are model-agnostic.** An organism’s
    birth size, maturity size/age, and asymptotic size exist
    independently of how we mathematically describe growth. We’re using
    these real quantities, not model artifacts.

3.  **VB-equivalent \\k\\ captures the relevant information.** What CW
    needs is “how fast does this organism approach its asymptotic size?”
    The VB-equivalent \\k\\ answers this question using information that
    any growth model provides.

4.  **Consistency check available.** If VB, Gompertz, and Logistic all
    produce similar mortality estimates (as they should when all fit the
    data reasonably), this validates the approach. If they diverge
    wildly, that reveals genuine uncertainty that shouldn’t be hidden.

## Limitations and Caveats

### When Models Disagree

If different growth models produce substantially different mortality
estimates, this indicates:

- The growth models encode meaningfully different information about
  growth rate
- Genuine model uncertainty that should be reported
- Possibly that some models fit poorly and shouldn’t be used

In such cases, we recommend: (1) comparing growth models via LOO-CV, (2)
using the best-fitting model for primary inference, and (3) reporting
sensitivity to model choice.

### Very Different Growth Trajectories

The approach works best when all three growth models produce similar
predicted trajectories across the observed age range. If a dataset
strongly favors one model — e.g., Gompertz fits vastly better than VB —
the VB-equivalent \\k\\ from Gompertz is still valid, but it’s
essentially computing “what would VB k be if we forced this biology into
VB form?”

This is statistically defensible but philosophically imperfect. In
extreme cases, researchers might consider whether a non-CW mortality
model (Peterson-Wroblewski, Lorenzen) is more appropriate.

### Extrapolation Beyond Data

All mortality estimation involves extrapolation — we predict mortality
at ages beyond our data range. The VB-equivalent \\k\\ approach doesn’t
change this fundamental limitation. If biological milestones are poorly
estimated (wide posteriors on \\L\_{mat}\\ or \\t\_{mat}\\), mortality
estimates will have appropriately wide uncertainty.

## Conclusions

The L₀-parameterized, growth-model-agnostic Chen-Watanabe implementation
in vitalBayes addresses long-standing practical limitations:

1.  **Observable parameters**: Birth length replaces theoretical “age at
    length zero”
2.  **Model flexibility**: Any growth model can feed mortality
    estimation  
3.  **Biological grounding**: All inference flows from real, measurable
    quantities
4.  **Uncertainty propagation**: Full posterior information preserved

This approach embodies the vitalBayes philosophy: parameterize with
biology, not mathematics.

## References

Chen, S., & Watanabe, S. (1989). Age dependence of natural mortality
coefficient in fish population dynamics. *Nippon Suisan Gakkaishi*,
55(2), 205-208.

Katsanevakis, S., & Maravelias, C. D. (2008). Modelling fish growth:
multi-model inference as a better alternative to a priori using von
Bertalanffy equation. *Fish and Fisheries*, 9(2), 178-187.

Lorenzen, K. (1996). The relationship between body weight and natural
mortality in juvenile and adult fish. *Journal of Fish Biology*, 49(4),
627-642.

Thorson, J. T., & Prager, M. H. (2011). Better catch curves:
Incorporating age-specific natural mortality and logistic selectivity.
*Transactions of the American Fisheries Society*, 140(2), 356-366.

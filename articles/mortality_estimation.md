# Natural Mortality Estimation with get_stochastic_mortality()

## Overview

Natural mortality (\\M\\) — death from causes other than fishing — is
among the most difficult parameters to estimate in fisheries biology.
Direct observation is rarely possible, so we rely on indirect methods
that derive mortality from life history parameters. These approaches
introduce uncertainty that must be properly propagated to downstream
analyses.

The
[`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md)
function implements Monte Carlo simulation of age-specific mortality
schedules with full uncertainty propagation from growth model
posteriors. A key innovation is **growth-model-agnostic estimation**:
mortality can be derived from von Bertalanffy, Gompertz, or Logistic
growth fits without theoretical compromise.

This vignette covers the three available mortality models, the
growth-model-agnostic architecture, scaling approaches, and integration
with the broader vitalBayes workflow.

## Mortality Models

### Chen-Watanabe (1989)

The Chen-Watanabe model assumes that natural mortality is inversely
related to body size under von Bertalanffy growth dynamics. vitalBayes
implements an **L₀-parameterized** version that eliminates dependence on
the theoretical parameter \\t_0\\:

\\M(t) = \frac{k \cdot L\_\infty}{L(t)}\\

where \\L(t) = L\_\infty - (L\_\infty - L_0)e^{-kt}\\ is the predicted
length at age \\t\\.

This reformulation uses observable quantities (\\L_0\\, birth length)
rather than mathematical abstractions (\\t_0\\, age at length zero). See
[`vignette("chen_watanabe_reparameterization")`](https://brian-j-moe.github.io/vitalBayes/articles/chen_watanabe_reparameterization.md)
for the full mathematical derivation.

**Two-phase extension.** The original CW model produces mortality that
declines asymptotically to \\k\\, which is biologically unrealistic —
mortality should increase at old ages due to senescence. The two-phase
extension models late-life mortality increase via Gompertz or logistic
functions, producing more realistic bathtub-shaped mortality curves.

**When to use CW.** Species where juvenile mortality is high and
declines with size; species with documented senescence patterns;
situations where the theoretical grounding of CW (von Bertalanffy
derivation) is appropriate.

### Peterson-Wroblewski (1984)

The Peterson-Wroblewski model takes a purely allometric approach:

\\M(W) = 1.92 \cdot W^{-0.25}\\

where \\W\\ is body weight in grams. This captures the well-documented
pattern that smaller individuals experience higher mortality, with the
-0.25 exponent reflecting metabolic scaling.

**Requirements.** A length-weight function is needed to convert
predicted lengths to body mass.

**When to use PW.** Species with well-characterized length-weight
relationships; comparative studies across taxa; situations where you
want a simple, empirically-grounded model without growth model
assumptions.

### Lorenzen (1996, 2022)

Lorenzen offers two formulations:

**Weight-based (1996):** \\M(W) = \alpha \cdot W^{\beta}\\ where
\\\alpha \sim N(3.69, 0.502)\\ and \\\beta \sim N(-0.305, 0.029)\\.

**Growth-based (2022):** \\\ln M = 0.28 - 1.30 \ln(L/L\_\infty) + 1.08
\ln(k)\\

The growth-based formulation was calibrated using von Bertalanffy
parameters, so it should use VB-equivalent \\k\\ when the underlying
growth model is Gompertz or Logistic.

**When to use Lorenzen.** General applications; meta-analytic contexts;
when you want built-in parameter uncertainty from the original
calibration.

## Growth-Model-Agnostic Estimation

A central innovation in vitalBayes is the ability to estimate
Chen-Watanabe mortality from *any* growth model fit — not just von
Bertalanffy.

### The Challenge

Chen-Watanabe was derived under von Bertalanffy assumptions, so its
\\k\\ parameter specifically means the VB growth coefficient. But the VB
model often performs poorly on elasmobranch data, especially when adult
observations are sparse. Gompertz and Logistic models frequently provide
better fits.

Historically, researchers faced an uncomfortable choice: use the
best-fitting growth model for biological inference but abandon CW
mortality estimation, or force a VB fit despite knowing it describes the
data poorly.

### The Solution: VB-Equivalent k

All three growth models estimate the same biological quantities:

- \\L\_\infty\\ — asymptotic length
- \\L_0\\ — length at birth
- \\L\_{mat}\\ — length at maturity
- \\t\_{mat}\\ — age at maturity

From these milestones, we can compute what the von Bertalanffy \\k\\
would need to be:

\\k\_{VB}^{equiv} = \frac{1}{t\_{mat}} \ln\left(\frac{L\_\infty -
L_0}{L\_\infty - L\_{mat}}\right)\\

This VB-equivalent \\k\\ captures the “growth rate information” in a
form Chen-Watanabe can use, regardless of which growth model generated
the biological parameters.

``` r
library(vitalBayes)

# Fit Gompertz (assume it fits your data better than VB)
gomp_fit <- fit_bayesian_growth(
  lt    = "fl",
  age   = "age",
  sex   = "sex",
  data  = growth_data[embryo == FALSE & !is.na(age)],
  model = "gompertz",
  k_based = FALSE,  # Maturity-based parameterization required
  birth_fit = birth_fit,
  L50_fit = L50_fit,
  t50_fit = t50_fit
)

# CW mortality from Gompertz fit — this works automatically
mort <- get_stochastic_mortality(
  method     = "CW",
  growth_fit = gomp_fit,  # Any growth model
  sex        = 1,
  iter       = 2000
)

# The function extracts (Linf, L0, Lmat, tmat) from Gompertz posterior
# and computes VB-equivalent k internally
summary(mort$Parameters$k_vb_equiv)
```

### Requirements for Growth-Model-Agnostic Estimation

This approach requires **maturity-based growth parameterization**
(`k_based = FALSE` in
[`fit_bayesian_growth()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md)).
The maturity parameters \\L\_{mat}\\ and \\t\_{mat}\\ must be available
in the growth model posterior.

For k-based growth fits, you can still supply a separate maturity fit
via the `maturity_fit` argument, but the integration is less seamless.

## Basic Usage

### From Growth Model Posterior

The recommended workflow extracts parameters from vitalBayes growth
fits, preserving correlations:

``` r
# Complete workflow
birth_fit <- fit_bayesian_birth(
  embryo_lts = growth_data[embryo == TRUE, fl],
  free_swimming_lts = growth_data[embryo == FALSE, fl]
)

L50_fit <- fit_bayesian_maturity(
  maturity = "mat", lt = "fl", sex = "sex",
  data = growth_data[embryo == FALSE & !is.na(mat)],
  use_pooling = TRUE
)

t50_fit <- fit_bayesian_maturity(
  maturity = "mat", age = "age", sex = "sex",
  data = growth_data[embryo == FALSE & !is.na(mat) & !is.na(age)],
  use_pooling = TRUE
)

growth_fit <- fit_bayesian_growth(
  lt = "fl", age = "age", sex = "sex",
  data = growth_data[embryo == FALSE & !is.na(age)],
  model = "gompertz",  # Use whichever model fits best
  k_based = FALSE,
  birth_fit = birth_fit,
  L50_fit = L50_fit,
  t50_fit = t50_fit
)

# Mortality with full uncertainty propagation
mort <- get_stochastic_mortality(
  method     = "CW",
  growth_fit = growth_fit,
  sex        = 1,
  iter       = 2000,
  scaled     = TRUE,
  p          = 0.001
)

# View results
mort$Plot
mort$Summary
```

### Manual Parameter Specification

For sensitivity analysis or when using literature values:

``` r
mort <- get_stochastic_mortality(
  method = "CW",
  Linf = c(100, 5),    # mean, sd
  L0   = c(25, 2),
  k    = c(0.10, 0.02),  # Should be VB or VB-equivalent k
  tmat = c(10, 1),
  iter = 2000,
  scaled = TRUE,
  p = 0.001
)
```

Note that manual specification samples parameters independently,
ignoring correlations. This typically produces wider uncertainty bounds
than posterior-based estimation.

## Scaling Mortality Schedules

Theoretical mortality models often produce absolute levels that don’t
match empirical observations. Scaling adjusts the overall level while
preserving the age-specific pattern.

### Survival Probability Scaling

Scale so that a specified fraction survives to maximum age:

``` r
# 0.1% survive to tmax (common assumption)
mort <- get_stochastic_mortality(
  method = "CW",
  growth_fit = growth_fit,
  sex = 1,
  scaled = TRUE,
  M_target = NULL,  # Derive from p
  p = 0.001
)
```

### Empirical Relationship Scaling

Scale to match empirical M estimators like Then et al. (2015):

``` r
# Then et al. (2015) tmax-based estimator
then_2015 <- function(tmax) 4.899 * tmax^(-0.916)

mort <- get_stochastic_mortality(
  method = "CW",
  growth_fit = growth_fit,
  sex = 1,
  scaled = TRUE,
  M_target = then_2015
)
```

### Fixed Target Scaling

Scale to a specific mean mortality value:

``` r
mort <- get_stochastic_mortality(
  method = "CW",
  growth_fit = growth_fit,
  sex = 1,
  scaled = TRUE,
  M_target = 0.15  # Fixed target
)
```

## Comparing Models

### Across Mortality Models

``` r
# Chen-Watanabe
mort_cw <- get_stochastic_mortality(
  method = "CW", growth_fit = growth_fit, sex = 1,
  print_plot = FALSE
)

# Peterson-Wroblewski
lw_fun <- function(L) 0.0001 * L^3.1
mort_pw <- get_stochastic_mortality(
  method = "PW", growth_fit = growth_fit, sex = 1,
  lw_fun = lw_fun, print_plot = FALSE
)

# Lorenzen (growth-based)
mort_lor <- get_stochastic_mortality(
  method = "L", growth_fit = growth_fit, sex = 1,
  weight_based = FALSE, print_plot = FALSE
)

# Combine
combined <- rbind(
  mort_cw$Summary[, model := "Chen-Watanabe"],
  mort_pw$Summary[, model := "Peterson-Wroblewski"],
  mort_lor$Summary[, model := "Lorenzen"]
)

# Plot comparison
library(ggplot2)
ggplot(combined, aes(x = age_round, color = model, fill = model)) +
  geom_ribbon(aes(ymin = M_lower, ymax = M_upper), alpha = 0.2, color = NA) +
  geom_line(aes(y = M_median), linewidth = 1) +
  labs(x = "Age (years)", y = "Natural Mortality (M)",
       title = "Mortality Model Comparison") +
  theme_bw() +
  theme(legend.position = "top")
```

### Across Growth Models

Test sensitivity of mortality to growth model choice:

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

# CW mortality from each
mort_from_vb <- get_stochastic_mortality(
  method = "CW", growth_fit = vb_fit, sex = 1, print_plot = FALSE
)
mort_from_gomp <- get_stochastic_mortality(
  method = "CW", growth_fit = gomp_fit, sex = 1, print_plot = FALSE
)
mort_from_logis <- get_stochastic_mortality(
  method = "CW", growth_fit = logis_fit, sex = 1, print_plot = FALSE
)

# If estimates are similar, mortality is robust to growth model choice
# If they differ substantially, report this sensitivity
```

## Two-Phase Chen-Watanabe Options

### Gompertz Senescence (Default)

``` r
mort <- get_stochastic_mortality(
  method = "CW",
  growth_fit = growth_fit,
  sex = 1,
  two_phase = TRUE,
  late_model = "gompertz",
  tm_factor = 2/3,   # Transition starts at 2/3 of tmat
  M_mult = 2         # Senescence mortality = 2x early-phase M
)
```

### Logistic Senescence

``` r
mort <- get_stochastic_mortality(
  method = "CW",
  growth_fit = growth_fit,
  sex = 1,
  two_phase = TRUE,
  late_model = "logistic",
  tm_factor = 2/3,
  M_mult = 2
)
```

### Single-Phase (Original CW)

``` r
mort <- get_stochastic_mortality(
  method = "CW",
  growth_fit = growth_fit,
  sex = 1,
  two_phase = FALSE  # Original CW, no senescence
)
```

## Output Structure

[`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md)
returns a list with four components:

``` r
names(mort)
#> [1] "Schedules"  "Parameters" "Summary"    "Plot"

# Schedules: All mortality curves
head(mort$Schedules)
#>    set_id  age     M_raw  M_scaled
#> 1:      1  0.1 0.4521032 0.3892145
#> 2:      1  0.5 0.3201893 0.2756123
#> ...

# Parameters: Life history values used
head(mort$Parameters)
#>    set_id     Linf       L0     Lmat     tmat k_original k_vb_equiv k_for_mort     tmax
#> 1:      1 98.23451 24.12893 68.92341 9.823145  0.0891234  0.0934521  0.0934521 32.14521
#> ...

# Summary: Median and 95% CI by age
head(mort$Summary)
#>    age_round  M_median    M_mean   M_lower   M_upper
#> 1:       0.1 0.3892145 0.3912893 0.2891234 0.5123456
#> ...

# Plot: ggplot2 object
class(mort$Plot)
#> [1] "gg"     "ggplot"
```

## Downstream Integration

Mortality schedules feed directly into survival simulation:

``` r
surv <- simulate_survivorship(
  mc_object = mort,
  n = 50000,
  n_iter = 2000
)

# Key metrics
surv$Aggregate$Age_of_Death
surv$Aggregate$Survival_to_tmat
```

See
[`vignette("survivorship_simulation")`](https://brian-j-moe.github.io/vitalBayes/articles/survivorship_simulation.md)
for comprehensive survival analysis.

## References

Chen, S., & Watanabe, S. (1989). Age dependence of natural mortality
coefficient in fish population dynamics. *Nippon Suisan Gakkaishi*,
55(2), 205-208.

Lorenzen, K. (1996). The relationship between body weight and natural
mortality in juvenile and adult fish. *Journal of Fish Biology*, 49(4),
627-642.

Lorenzen, K. (2022). Size- and age-dependent natural mortality in fish
populations. *Fisheries Research*, 255, 106454.

Peterson, I., & Wroblewski, J. S. (1984). Mortality rate of fishes in
the pelagic ecosystem. *Canadian Journal of Fisheries and Aquatic
Sciences*, 41(7), 1117-1120.

Then, A. Y., Hoenig, J. M., Hall, N. G., & Hewitt, D. A. (2015).
Evaluating the predictive performance of empirical estimators of natural
mortality rate. *ICES Journal of Marine Science*, 72(1), 82-92.

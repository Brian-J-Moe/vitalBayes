# vitalBayes <img src="man/figures/logo.png" align="right" height="139" />
## From Birth to Death: Bayesian Estimation of Vital Rates for Elasmobranchs

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen.svg)](https://github.com/Brian-J-Moe/vitalBayes)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->


`vitalBayes` provides a coherent Bayesian framework for estimating population 
vital rates in elasmobranchs — from birth and maturity through growth to natural 
mortality and survival. The package implements hierarchical models via 
[Stan](https://mc-stan.org/) with precompiled models for fast, reliable inference.

The package addresses critical methodological challenges in elasmobranch research: 
sparse datasets, imbalanced sex ratios, strong parameter correlations, and the 
need for biologically coherent estimates that propagate uncertainty across life 
stages.

## Why vitalBayes?

Elasmobranch vital rate estimation presents unique challenges:

- **Sparse embryo samples** make birth size estimation uncertain
- **Imbalanced sex ratios** lead to unreliable parameter estimates for the 
undersampled sex
- **Strong parameter correlations** (especially $L_\infty$ and $k$) complicate 
growth model inference
- **Mortality is rarely measured directly** — it must be derived from life 
history parameters with appropriate uncertainty propagation

vitalBayes addresses these challenges through:

| Challenge | Solution |
|-----------|----------|
| Sparse data | **Partial pooling** borrows strength across sexes without forcing identical parameters |
| Parameter correlation | **Maturity-based parameterization** derives $k$ from observable maturity milestones |
| Threshold estimation | **Probit link** provides biologically intuitive interpretation for birth and maturity |
| Prior specification | **CV-based priors** offer scale-invariant, intuitive parameter uncertainty |
| Mortality uncertainty | **Joint posterior sampling** preserves correlations when deriving mortality schedules |
| Reproducibility | **Precompiled Stan models** ensure consistent, fast inference |

<!--
## Interactive Demo

Explore how the three growth models (von Bertalanffy, Gompertz, Logistic) respond to changes in life history parameters. 
All curves pass through the same maturity point ($t_{mat}$, $L_{mat}$) — the growth coefficient $k$ is derived from these observable milestones rather than estimated directly.

<div align="center">
<iframe 
  src="growth_explorer.html" 
  width="100%" 
  height="650" 
  style="border: 2px solid #3D4466; border-radius: 12px; max-width: 1000px;"
  loading="lazy">
</iframe>
</div>
<p align="center"><em>First load may take 10-20 seconds while WebR initializes</em></p>
-->

## Installation

vitalBayes requires a working C++ toolchain for Stan. See the 
[RStan Getting Started Guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) 
for setup instructions.

```r
# Install dependencies
install.packages(c("cmdstanr", "instantiate", "data.table", "ggplot2", "loo"))

# Install cmdstan (if not already installed)
cmdstanr::install_cmdstan()

# Install vitalBayes from GitHub
# install.packages("pak")
pak::pak("Brian-J-Moe/vitalBayes")
```

## The Four-Stage Workflow

vitalBayes implements an integrated workflow where information flows forward 
through life stages — from individual-level parameters to population-level 
vital rates:

```
  Birth      ───▶      Maturity      ───▶      Growth      ───▶    Mortality/Survival    
  (b₅₀)               (L₅₀, t₅₀)              (L∞, k, L₀)              (M, S(t))

    │                     │                       │                       │
    │                     │                       │                       │
    ▼                     ▼                       ▼                       ▼
 Informs L₀            Informs               Derived from           Derived from
  prior              Lmat, tmat             maturity params        growth posterior
```

Each stage produces posterior distributions that inform the next, creating a 
biologically coherent analysis with fully propagated uncertainty.

## Quick Start

```r
library(vitalBayes)
library(data.table)

# Load example data (simulated von Bertalanffy growth data)
data(growth_data)

# ─────────────────────────────────────────────────────────────
# Stage 1: Birth Size
# ─────────────────────────────────────────────────────────────
birth_fit <- fit_bayesian_birth(
  embryo_lts        = growth_data[embryo == TRUE, fl],
  free_swimming_lts = growth_data[embryo == FALSE, fl]
)

birth_fit$summary(c("b50", "transition_width"))

# ─────────────────────────────────────────────────────────────
# Stage 2: Maturity
# ─────────────────────────────────────────────────────────────
mat_data <- growth_data[embryo == FALSE & !is.na(mat)]

# Length-at-maturity with partial pooling for imbalanced sexes
L50_fit <- fit_bayesian_maturity(
  maturity    = "mat",
  lt          = "fl",
  sex         = "sex",
  data        = mat_data,
  use_pooling = TRUE
)

L50_fit$summary(c("L50", "L50_diff"))

# Age-at-maturity
t50_fit <- fit_bayesian_maturity(
  maturity    = "mat",
  age         = "age",
  sex         = "sex",
  data        = mat_data[!is.na(age)],
  use_pooling = TRUE
)

# ─────────────────────────────────────────────────────────────
# Stage 3: Growth (maturity-based parameterization)
# ─────────────────────────────────────────────────────────────
growth_fit <- fit_bayesian_growth(
  lt        = "fl",
  age       = "age",
  sex       = "sex",
  data      = growth_data[embryo == FALSE & !is.na(age)],
  model     = "vb",
  k_based   = FALSE,        # Use maturity-based parameterization
  birth_fit = birth_fit,    # L₀ prior from birth model
  L50_fit   = L50_fit,      # Lmat prior
  t50_fit   = t50_fit,      # tmat prior
  use_pooling = TRUE
)

growth_fit$summary(c("Linf", "L0", "k", "Lmat", "tmat"))

# ─────────────────────────────────────────────────────────────
# Stage 4: Mortality & Survival
# ─────────────────────────────────────────────────────────────
# Generate age-specific mortality schedules with uncertainty
mort <- get_stochastic_mortality(
  method       = "CW",           # Chen-Watanabe model
  growth_fit   = growth_fit,     # Joint posterior sampling
  maturity_fit = t50_fit,        # Age-at-maturity for two-phase CW
  sex          = 2,              # Males
  iter         = 2000,
  scaled       = TRUE,           # Scale to survival probability
  p            = 0.001           # 0.1% survive to tmax
)

# Simulate cohort survival
surv <- simulate_survivorship(
  mc_object = mort,
  n         = 50000,             # Cohort size
  n_iter    = 2000               # Simulation iterations
)

# Key demographic metrics
surv$Aggregate$Age_of_Death      # Mean lifespan with 95% CI
surv$Aggregate$Survival_to_tmat  # Probability of reaching maturity
```

## Example Datasets

vitalBayes includes three simulated datasets for demonstrations and testing:

| Dataset | Description | Sample Sizes | Best For |
|---------|-------------|--------------|----------|
| `growth_data` | Balanced, comprehensive | 189F, 176M, 26 embryos | Main workflow, documentation |
| `imbalanced_data` | Imbalanced sex ratio | 150F, 34M, 13 embryos | Partial pooling demonstrations |
| `limited_data` | Small sample sizes | 24F, 18M, 5 embryos | Prior sensitivity analysis |

All datasets share the same structure with columns: `sex`, `mat`, `fl`, `age`, `embryo`.

## Visualization

vitalBayes includes publication-ready plotting functions with its signature 
synthwave color palette:

```r
# Growth curves with credible intervals
plot_growth_curve(
  fit        = growth_fit,
  data       = growth_data[embryo == FALSE & !is.na(age)],
  age_col    = "age",
  length_col = "fl",
  sex_col    = "sex"
)

# Maturity ogives
plot_maturity_ogive(
  fit  = L50_fit,
  type = "length",
  data = mat_data,
  x_col = "fl",
  maturity_col = "mat",
  sex_col = "sex"
)

# Mortality schedules
mort$Plot  # Auto-generated by get_stochastic_mortality()

# Survivorship curves
surv$Plot  # Auto-generated by simulate_survivorship()

# Compare growth models
compare_growth_models(
  "von Bertalanffy" = vb_fit,
  "Gompertz"        = gomp_fit,
  "Logistic"        = logis_fit,
  data = growth_data[embryo == FALSE & !is.na(age)]
)
```

### Colorblind-Safe Palettes

```r
# Check available palettes
list_vital_palettes()

# Use colorblind-safe alternatives
plot_growth_curve(growth_fit, data = growth_data, palette = "okabe")

# Available palettes: synthwave (default), viridis, okabe, plasma, inferno
vital_palette(n = 4, type = "okabe")
```

## Model Comparison

Compare alternative models using leave-one-out cross-validation:

```r
# Compute LOO-CV
loo_vb   <- compute_loo(vb_fit)
loo_gomp <- compute_loo(gomp_fit)

# Compare models
compare_loo(
  "von Bertalanffy" = loo_vb,
  "Gompertz"        = loo_gomp
)

# Publication-ready table
create_loo_table(
  "von Bertalanffy" = loo_vb,
  "Gompertz"        = loo_gomp
)
```

## Key Innovations

### 1. Maturity-Based Growth Parameterization

Traditional growth models estimate $k$ directly, but $k$ and $L_\infty$ are 
strongly negatively correlated (often *r* < −0.9), making both poorly identified 
when data don't extend to near-asymptotic sizes.

**Our solution:** Derive $k$ from observable maturity milestones:

$$k = \frac{1}{t_{mat}} \ln\left(\frac{L_\infty - L_0}{L_\infty - L_{mat}}\right)$$

This parameterization uses quantities that fall within the data range ($L_{mat}$, 
$t_{mat}$), ensures biological consistency by forcing the growth curve through 
the maturity point, and propagates uncertainty from upstream maturity models.

### 2. Probit Link for Threshold Models

Both birth and maturity models use a probit link function rather than the more 
commonly used logit. The probit link assumes that underlying developmental 
readiness is normally distributed across individuals. An individual transitions 
(births or matures) when this latent readiness crosses a threshold.

**Why probit?** Normal variation in developmental readiness is biologically 
plausible, the transition width ($x_{95}$ − $x_{05}$) has units of the predictor 
(cm or years), and derived quantities ($L_{05}$, $L_{95}$) have direct biological 
interpretation.

### 3. Partial Pooling for Imbalanced Sex Ratios

**The problem:** Elasmobranch sampling frequently yields imbalanced sex ratios. 
A dataset might contain 150 females but only 34 males. Fitting separate models 
produces wide, unreliable credible intervals for the sparse sex, estimates 
driven by a few influential observations, and inefficient use of information 
(both sexes are the same species).

Complete pooling (ignoring sex) is biologically unrealistic. No pooling (fully 
separate) wastes information.

**Our solution:** Partial pooling models sex-specific parameters as draws from 
a common distribution:

$$\log(\theta_s) = \mu + \tau \cdot \eta_s, \quad \eta_s \sim \mathcal{N}(0, 1)$$

where μ is the population mean, τ is between-sex standard deviation, and 
η<sub>s</sub> are standardized sex deviations. The model **learns τ from the 
data**: if sexes appear similar, estimates shrink together; if sexes appear 
different, estimates stay separated; sparse sex borrows strength from data-rich sex.

**Non-centered parameterization:** We use the non-centered form (estimating η 
rather than θ directly) because it dramatically improves MCMC sampling when τ 
is small.

**Half-normal prior on τ:** We place τ ~ Half-Normal(0, σ<sub>τ</sub>) rather 
than the commonly-used half-Cauchy. The half-Cauchy's heavy tails can pull τ 
toward implausibly large values, especially with only two groups.

```r
# Enable partial pooling (especially useful for imbalanced data)
data(imbalanced_data)  # 150 females, 34 males

fit <- fit_bayesian_maturity(
  maturity    = "mat",
  lt          = "fl", 
  sex         = "sex",
  data        = imbalanced_data[embryo == FALSE],
  use_pooling = TRUE,    # Hierarchical structure
  prior_tau   = 0.5      # Half-normal scale (on log scale)
)

# Sex-specific estimates with appropriate shrinkage
fit$summary("L50")

# Between-sex difference (still estimable!)
fit$summary("L50_diff")
```

See `vignette("partial_pooling")` for a comprehensive treatment.

### 4. Constrained Asymptotic Length ($L_\infty$ > $L_{max}$)

**The problem:** Unconstrained growth models frequently converge on $L_\infty$ 
values *below* the largest observed individuals — biologically impossible, yet 
common when data don't extend to ages near asymptotic size.

**Our solution:** vitalBayes enforces $L_\infty$ > $L_{max}$ through a shifted 
lognormal prior:

$$L_\infty = L_{max} + \exp(\nu), \quad \nu \sim \mathcal{N}(\mu_\nu, \sigma_\nu)$$

This ensures $L_\infty$ always exceeds the maximum observed length while still 
allowing uncertainty about *how much* it exceeds it.

```r
# Lmax is auto-detected from data
growth_fit <- fit_bayesian_growth(
  lt   = "fl",
  age  = "age",
  data = growth_data[embryo == FALSE]
)
# Message: "Lmax from data: 98.5 cm (female), 87.2 cm (male)"

# Or specify manually
growth_fit <- fit_bayesian_growth(
  lt   = "fl",
  age  = "age",
  data = growth_data[embryo == FALSE],
  Lmax = c(100, 90)  # Female, Male
)
```

### 5. Joint Posterior Sampling for Mortality

**The problem:** Mortality models (Chen-Watanabe, Peterson-Wroblewski, Lorenzen) 
require growth parameters as inputs. Traditional approaches specify these as 
independent point estimates with uncertainty — e.g., $L_\infty \sim N(100, 5)$ 
and $k \sim N(0.1, 0.02)$ — ignoring that these parameters are typically strongly 
correlated in the posterior.

Independent sampling occasionally generates biologically implausible combinations 
(high $L_\infty$ with high $k$) that never appeared in the original growth model 
posterior, artificially inflating mortality uncertainty.

**Our solution:** `get_stochastic_mortality()` accepts vitalBayes growth fits 
directly, drawing from the **joint posterior distribution**. Parameter correlations 
are preserved, yielding mortality schedules with narrower but more honest 
uncertainty intervals.

```r
# Traditional approach (ignores correlations)
mort_naive <- get_stochastic_mortality(
  method = "CW",
  Linf   = c(100, 5),   # Specified as mean, sd
  k      = c(0.1, 0.02),
  t0     = c(-1, 0.2)
)

# vitalBayes approach (preserves correlations)
mort_joint <- get_stochastic_mortality(
  method     = "CW",
  growth_fit = growth_fit,  # Full posterior
  sex        = 2
)
```

### 6. Multiple Mortality Models

vitalBayes implements three natural mortality models, each with distinct assumptions:

| Model | Best For | Key Feature |
|-------|----------|-------------|
| **Chen-Watanabe** | Species with distinct life phases | Two-phase: declining juvenile M, late-life senescence |
| **Peterson-Wroblewski** | Weight-mortality relationships | Allometric: $M \propto W^{-0.25}$ |
| **Lorenzen** | General applications | Growth-based or weight-based formulations |

The Chen-Watanabe model supports multiple late-life senescence options (Gompertz, 
logistic) and can be scaled to empirical mortality estimators (Hoenig, Then et al.).

```r
# Chen-Watanabe with Gompertz senescence
mort_cw <- get_stochastic_mortality(
  method     = "CW",
  growth_fit = growth_fit,
  sex        = 1,
  two_phase  = TRUE,
  late_model = "gompertz"
)

# Peterson-Wroblewski (requires length-weight function)
lw_fun <- function(L) 0.0001 * L^3.1
mort_pw <- get_stochastic_mortality(
  method     = "PW",
  growth_fit = growth_fit,
  sex        = 1,
  lw_fun     = lw_fun
)

# Lorenzen growth-based (no L-W function needed)
mort_lor <- get_stochastic_mortality(
  method       = "L",
  growth_fit   = growth_fit,
  sex          = 1,
  weight_based = FALSE
)
```

### 7. CV-Based Prior Elicitation

**The problem:** Specifying priors on log-transformed parameters is unintuitive. 
What does σ = 0.5 on log($L_{50}$) mean in practical terms?

**Our solution:** vitalBayes uses **coefficient of variation (CV)** for prior 
specification. The CV expresses uncertainty as a proportion of the mean — a 
scale-invariant, intuitive quantity.

**Example:** Setting `CV_Linf = 0.2` for a species with expected $L_\infty$ ≈ 100 cm 
means prior SD ≈ 20 cm and 95% prior interval ≈ 60–140 cm.

```r
growth_fit <- fit_bayesian_growth(
  lt   = "fl",
  age  = "age",
  data = growth_data[embryo == FALSE],
  
  # Intuitive CV-based priors
  CV_Linf = 0.20,    # 20% uncertainty on asymptotic length
  CV_L0   = 0.30,    # 30% uncertainty on birth length
  CV_Lmat = 0.20,    # 20% uncertainty on length-at-maturity
  CV_tmat = 0.30     # 30% uncertainty on age-at-maturity
)
```

### 8. Automatic Multilingual Sex Coding

vitalBayes automatically detects sex coding in multiple languages:

```r
# All of these work automatically:
# English:    "F"/"M", "Female"/"Male"
# Spanish:    "H"/"M", "Hembra"/"Macho"
# Portuguese: "F"/"M", "Fêmea"/"Macho"
# French:     "F"/"M", "Femelle"/"Mâle"
# German:     "W"/"M", "Weibchen"/"Männchen"
# Japanese:   "メス"/"オス", "雌"/"雄"
# Chinese:    "母"/"公", "雌"/"雄"
# Russian:    "Ж"/"М", "Самка"/"Самец"
# Symbols:    "♀"/"♂"
# Numeric:    1/2 (female/male)
```

## Vignettes

| Vignette | Description |
|----------|-------------|
| `vignette("fit_bayesian_birth")` | Birth size estimation |
| `vignette("fit_bayesian_maturity")` | L₅₀ and t₅₀ estimation |
| `vignette("fit_bayesian_growth")` | Growth models with maturity-based parameterization |
| `vignette("partial_pooling")` | Hierarchical modeling for imbalanced data |
| `vignette("mortality_estimation")` | Natural mortality models and uncertainty propagation |
| `vignette("survivorship_simulation")` | Cohort survival analysis |
| `vignette("visualization")` | Plotting functions and color palettes |
| `vignette("model_diagnostics")` | Convergence, PPC, and LOO-CV |

For a comprehensive statistical background, see the 
[Statistical Methods Guide](articles/Understanding_vitalBayes.html).

## Citation

If you use vitalBayes in your research, please cite:

```
@manual{,
  title = {vitalBayes: Bayesian Vital Rate Estimation for Elasmobranchs},
  author = {Brian J Moe},
  year = {2025},
  note = {R package version 0.4.0},
  url = {https://github.com/Brian-J-Moe/vitalBayes}
}
```

## Dependencies

**Core:**

- [cmdstanr](https://mc-stan.org/cmdstanr/) — CmdStan interface
- [instantiate](https://github.com/wlandau/instantiate) — Precompiled Stan models
- [data.table](https://rdatatable.gitlab.io/data.table/) — Fast data manipulation

**Visualization:**

- [ggplot2](https://ggplot2.tidyverse.org/) — Publication-quality graphics

**Model Comparison:**

- [loo](https://mc-stan.org/loo/) — Leave-one-out cross-validation

**Numerical Integration:**

- [pracma](https://CRAN.R-project.org/package=pracma) — Trapezoidal integration for survival curves

## Contributing

Contributions are welcome! Please open an issue to discuss proposed changes or 
submit a pull request.

## License

[MIT License](LICENSE.md)

---

*vitalBayes was developed for the elasmobranch research community to support 
robust, reproducible vital rate analyses — from birth to mortality.*

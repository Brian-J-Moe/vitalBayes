# vitalBayes

## Bayesian Life History Parameter Estimation for Elasmobranchs

`vitalBayes` provides a coherent statistical framework for estimating
birth size, maturity, and growth parameters in sharks and rays using
Bayesian methods. The package implements hierarchical models via
[Stan](https://mc-stan.org/) with precompiled models for fast, reliable
inference.

The package addresses critical methodological challenges in shark and
ray research, including sparse datasets, imbalanced sex ratios,
uncertain aging methodologies, and the need for biologically coherent
parameter estimation across life stages.

## Why vitalBayes?

Elasmobranch life history estimation presents unique challenges:

- **Sparse embryo samples** make birth size estimation uncertain
- **Imbalanced sex ratios** lead to unreliable parameter estimates for
  the undersampled sex
- **Strong parameter correlations** (especially *L*_(∞) and *k*)
  complicate growth model inference
- **Multiple aging methods** may yield different biological conclusions

vitalBayes addresses these challenges through:

| Challenge             | Solution                                                                               |
|-----------------------|----------------------------------------------------------------------------------------|
| Sparse data           | **Partial pooling** borrows strength across sexes without forcing identical parameters |
| Parameter correlation | **Maturity-based parameterization** derives *k* from observable maturity milestones    |
| Threshold estimation  | **Probit link** provides biologically intuitive interpretation for birth and maturity  |
| Prior specification   | **CV-based priors** offer scale-invariant, intuitive parameter uncertainty             |
| Reproducibility       | **Precompiled Stan models** ensure consistent, fast inference                          |

## Installation

vitalBayes requires a working C++ toolchain for Stan. See the [RStan
Getting Started
Guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for
setup instructions.

``` r
# Install dependencies
install.packages(c("cmdstanr", "instantiate", "data.table", "ggplot2", "loo"))

# Install cmdstan (if not already installed)
cmdstanr::install_cmdstan()

# Install vitalBayes from GitHub
# install.packages("pak")
pak::pak("Brian-J-Moe/vitalBayes")
```

## The Three-Stage Workflow

vitalBayes implements an integrated workflow where information flows
forward through life history stages:

            Birth      ───▶      Maturity      ───▶      Growth
            (b₅₀)               (L₅₀, t₅₀)              (L∞, k, L₀)

              │                     │                       │
              │                     │                       │
              ▼                     ▼                       ▼
          Informs L₀             Informs               Derived from
            prior              Lmat, tmat             maturity params

Each stage produces posterior distributions that become informative
priors for the next, creating a biologically coherent analysis.

## Quick Start

``` r
library(vitalBayes)
library(data.table)

# Load example data (little gulper shark)
data(gulper_data)

# ─────────────────────────────────────────────────────────────
# Stage 1: Birth Size
# ─────────────────────────────────────────────────────────────
birth_fit <- fit_bayesian_birth(
  embryo_lts        = gulper_data[embryo == TRUE, fl],
  free_swimming_lts = gulper_data[embryo == FALSE, fl]
)

birth_fit$summary(c("b50", "transition_width"))

# ─────────────────────────────────────────────────────────────
# Stage 2: Maturity
# ─────────────────────────────────────────────────────────────
mat_data <- gulper_data[embryo == FALSE & !is.na(mat)]

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
  age         = "age1",
  sex         = "sex",
  data        = mat_data[!is.na(age1)],
  use_pooling = TRUE
)

# ─────────────────────────────────────────────────────────────
# Stage 3: Growth (maturity-based parameterization)
# ─────────────────────────────────────────────────────────────
growth_fit <- fit_bayesian_growth(
  lt        = "fl",
  age       = "age1",
  sex       = "sex",
  data      = gulper_data[embryo == FALSE & !is.na(age1)],
  model     = "vb",
  k_based   = FALSE,        # Use maturity-based parameterization
  birth_fit = birth_fit,    # L₀ prior from birth model
  L50_fit   = L50_fit,      # Lmat prior
  t50_fit   = t50_fit,      # tmat prior
  use_pooling = TRUE
)

growth_fit$summary(c("Linf", "L0", "k", "Lmat", "tmat"))
```

## Visualization

vitalBayes includes publication-ready plotting functions with a
distinctive retro-wave color palette:

``` r
# Growth curves with credible intervals
plot_growth_curve(
  fit        = growth_fit,
  data       = gulper_data[embryo == FALSE & !is.na(age1)],
  age_col    = "age1",
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

# Compare growth models
compare_growth_models(
  "von Bertalanffy" = vb_fit,
  "Gompertz"        = gomp_fit,
  "Logistic"        = logis_fit,
  data = growth_data
)
```

### Colorblind-Safe Palettes

``` r
# Check available palettes
list_vital_palettes()

# Use colorblind-safe alternatives
plot_growth_curve(growth_fit, data = growth_data, colorblind = TRUE)

# Or specify directly
vital_colors(2, "sex_cb")     # Orange/blue for sex
vital_colors(5, "okabe_ito")  # Okabe-Ito palette
```

### Multilingual Support

Sex labels in plots support multiple languages:

``` r
# French
plot_maturity_ogive(L50_fit, sex_labels = c("1" = "Femelle", "2" = "Mâle"))

# Spanish
plot_growth_curve(growth_fit, sex_labels = c("1" = "Hembra", "2" = "Macho"))

# Japanese
plot_maturity_ogive(L50_fit, sex_labels = c("1" = "メス", "2" = "オス"))
```

## Model Comparison

Compare alternative models using leave-one-out cross-validation:

``` r
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

**The problem:** Traditional growth models estimate the growth
coefficient *k* directly, but *k* is notoriously difficult to estimate.
It exhibits strong negative correlation with *L*_(∞) (often *r* \<
-0.9), and in datasets that don’t extend to ages approaching asymptotic
size, both parameters become poorly identified. Small changes in one can
be compensated by changes in the other, inflating posterior uncertainty.

**Our solution:** Rather than estimating *k* directly, vitalBayes
derives it from observable maturity milestones. The key insight is that
any growth curve must pass through the point (*t*_(mat), *L*_(mat))—the
age and length at which an individual reaches maturity. We solve each
growth equation for *k*:

#### von Bertalanffy

At maturity, \\L(t\_{mat}) = L\_{mat}\\: \\L\_{mat} = L\_\infty -
(L\_\infty - L_0) e^{-k \cdot t\_{mat}}\\

Solving for *k*: \\k = \frac{1}{t\_{mat}} \ln\left(\frac{L\_\infty -
L_0}{L\_\infty - L\_{mat}}\right)\\

#### Gompertz

At maturity: \\L\_{mat} = L\_\infty
\exp\left\[-\ln\left(\frac{L\_\infty}{L_0}\right) e^{-k \cdot
t\_{mat}}\right\]\\

Solving for *k*: \\k = \frac{1}{t\_{mat}} \ln\left(\frac{\ln(L\_\infty /
L_0)}{\ln(L\_\infty / L\_{mat})}\right)\\

#### Logistic

At maturity: \\L\_{mat} = \frac{L\_\infty}{1 +
\left(\frac{L\_\infty}{L_0} - 1\right) e^{-k \cdot t\_{mat}}}\\

Solving for *k*: \\k = \frac{1}{t\_{mat}}
\ln\left(\frac{L\_{mat}(L\_\infty - L_0)}{L_0(L\_\infty -
L\_{mat})}\right)\\

**Why this matters:**

| Aspect                  | Traditional (*k*-based)      | Maturity-based                       |
|-------------------------|------------------------------|--------------------------------------|
| Parameter correlation   | Strong (*L*_(∞), *k* ~ -0.9) | Reduced (maturity params observable) |
| Data requirements       | Ages near asymptote          | Maturity data within observed range  |
| Prior information       | Difficult to specify         | Informed by maturity analysis        |
| Uncertainty propagation | Independent                  | Flows from upstream models           |
| Biological coherence    | Not enforced                 | Curve passes through maturity        |

The maturity parameters (*L*_(mat), *t*_(mat)) typically fall well
within the observed data range, making them far more identifiable than
asymptotic parameters. By anchoring the growth curve at this observable
milestone, we transform a poorly-identified parameter (*k*) into a
derived quantity constrained by well-informed inputs.

``` r
# Maturity-based parameterization (recommended)
growth_fit <- fit_bayesian_growth(
  lt      = "fl",
  age     = "age1",
  sex     = "sex",
  data    = growth_data,
  k_based = FALSE,      # Derive k from maturity parameters
  L50_fit = L50_fit,    # Posterior for Lmat
  t50_fit = t50_fit     # Posterior for tmat
)

# k is now a derived quantity with appropriate uncertainty
growth_fit$summary("k")
```

### 2. Probit Link for Threshold Models

**The problem:** Maturity and birth-status models are typically fitted
using logistic regression (logit link). While computationally
convenient, the logit link implies that the log-odds of maturity change
linearly with size—a mathematical abstraction without clear biological
meaning.

**Our solution:** vitalBayes uses the probit link, which has a natural
**threshold-crossing interpretation**. The probit model assumes there
exists a latent continuous variable (developmental readiness,
reproductive investment) that is normally distributed across
individuals. An individual transitions when this latent variable crosses
a threshold:

\\P(\text{mature}\_i) = \Phi\[\beta \cdot (L_i - L\_{50})\]\\

where Φ(·) is the standard normal CDF.

**Biological justification:** Developmental readiness is influenced by
many factors—maternal condition, temperature, food availability, genetic
variation, sibling competition (for embryos). The central limit theorem
suggests that the aggregate effect of many small, independent factors
tends toward normality. The probit link formalizes this intuition.

**Practical implications:** - The slope parameter β has units of inverse
length (or inverse time), representing how quickly transition
probability changes - The transition width (5% to 95% probability)
equals 3.29/β, providing an intuitive measure of individual variation -
Derived quantities (*L*₀₅, *L*₉₅) have direct biological interpretation

### 3. Partial Pooling for Imbalanced Sex Ratios

**The problem:** Elasmobranch sampling frequently yields imbalanced sex
ratios. A dataset might contain 150 females but only 23 males. Fitting
separate models produces: - Wide, unreliable credible intervals for the
sparse sex - Estimates driven by a few influential observations  
- Inefficient use of information (both sexes are the same species)

Complete pooling (ignoring sex) is biologically unrealistic. No pooling
(fully separate) wastes information.

**Our solution:** Partial pooling models sex-specific parameters as
draws from a common distribution:

\\\log(\theta_s) = \mu + \tau \cdot \eta_s, \quad \eta_s \sim
\mathcal{N}(0, 1)\\

where μ is the population mean, τ is between-sex standard deviation, and
η_(s) are standardized sex deviations. The model **learns τ from the
data**:

- If sexes appear similar → small τ → estimates shrink together
- If sexes appear different → large τ → estimates stay separated
- Sparse sex → borrows strength from data-rich sex

**Non-centered parameterization:** We use the non-centered form
(estimating η rather than θ directly) because it dramatically improves
MCMC sampling when τ is small. When between-sex variation is minimal,
centered parameterizations create a “funnel” geometry that causes
divergent transitions.

**Half-normal prior on τ:** We place τ ~ Half-Normal(0, σ_(τ)) rather
than the commonly-used half-Cauchy. The half-Cauchy’s heavy tails can
pull τ toward implausibly large values, especially with only two groups
(sexes). The half-normal provides gentle regularization while still
allowing substantial between-sex variation when warranted.

``` r
# Enable partial pooling
fit <- fit_bayesian_maturity(
  maturity    = "mat",
  lt          = "fl", 
  sex         = "sex",
  data        = mat_data,
  use_pooling = TRUE,    # Hierarchical structure
  prior_tau   = 0.5      # Half-normal scale (on log scale)
)

# Sex-specific estimates with appropriate shrinkage
fit$summary("L50")

# Between-sex difference (still estimable!)
fit$summary("L50_diff")
```

See
[`vignette("partial_pooling")`](https://brian-j-moe.github.io/vitalBayes/articles/partial_pooling.md)
for a comprehensive treatment.

### 4. Constrained Asymptotic Length (*L*_(∞) \> *L*_(max))

**The problem:** Unconstrained growth models frequently converge on
*L*_(∞) values *below* the largest observed individuals. This is
biologically impossible—asymptotic length must exceed any realized
length—yet it occurs regularly when: - Data don’t extend to ages near
asymptotic size - Older individuals are undersampled (fishing
selectivity, natural mortality) - Strong *L*_(∞)–*k* correlation allows
trade-offs

**Our solution:** vitalBayes enforces *L*_(∞) \> *L*_(max) through a
shifted lognormal prior:

\\L\_\infty = L\_{max} + \exp(\nu), \quad \nu \sim \mathcal{N}(\mu\_\nu,
\sigma\_\nu)\\

This ensures *L*_(∞) always exceeds the maximum observed length while
still allowing uncertainty about *how much* it exceeds it.

**Biological rationale:** The largest individual we’ve observed
represents a lower bound on the species’ potential size. Natural
mortality, fishing pressure, and sampling limitations virtually
guarantee we haven’t captured the true maximum. Our prior should encode
this biological reality rather than allowing impossible values.

``` r
# Lmax is auto-detected from data
growth_fit <- fit_bayesian_growth(
  lt   = "fl",
  age  = "age1",
  data = growth_data
)
# Message: "Lmax from data: 98.5 cm (female), 87.2 cm (male)"

# Or specify manually
growth_fit <- fit_bayesian_growth(
  lt   = "fl",
  age  = "age1",
  data = growth_data,
  Lmax = c(100, 90)  # Female, Male
)
```

### 5. CV-Based Prior Elicitation

**The problem:** Bayesian models require prior distributions, but
specifying priors on log-transformed parameters is unintuitive. What
does σ = 0.5 on log(*L*₅₀) mean in practical terms? Researchers often
resort to “weakly informative” priors without considering whether they
encode reasonable biological information.

**Our solution:** vitalBayes uses **coefficient of variation (CV)** for
prior specification. The CV expresses uncertainty as a proportion of the
mean—a scale-invariant, intuitive quantity.

For a parameter θ with prior mean μ and CV = c: - The prior SD on the
original scale is σ_(θ) = μ · c - Approximately 95% of prior mass falls
within μ ± 2μc

**Example:** Setting `CV_Linf = 0.2` for a species with expected *L*_(∞)
≈ 100 cm means: - Prior SD ≈ 20 cm - 95% prior interval ≈ 60–140 cm

This is far more interpretable than specifying σ = 0.2 on log(*L*_(∞)).

``` r
growth_fit <- fit_bayesian_growth(
  lt   = "fl",
  age  = "age1",
  data = growth_data,
  
  # Intuitive CV-based priors
  CV_Linf = 0.20,    # 20% uncertainty on asymptotic length
  CV_L0   = 0.30,    # 30% uncertainty on birth length
  CV_Lmat = 0.20,    # 20% uncertainty on length-at-maturity
  CV_tmat = 0.30     # 30% uncertainty on age-at-maturity
)
```

### 6. Integrated Three-Stage Workflow

**The problem:** Life history parameters are typically estimated
independently, ignoring biological connections. Birth size is estimated
from embryo data, maturity from reproductive assessments, and growth
from age-length pairs—each analysis producing point estimates and
confidence intervals that don’t propagate into downstream analyses.

**Our solution:** vitalBayes implements an integrated workflow where
**posterior distributions flow forward** through life history stages:

    Birth (b₅₀)  ──▶  prior on L₀ for growth model
                     ╱
    Maturity (L₅₀) ──  prior on Lmat for growth model
                     ╲
    Maturity (t₅₀)  ──▶  prior on tmat for growth model

When you pass `birth_fit`, `L50_fit`, and `t50_fit` to
[`fit_bayesian_growth()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md),
the function extracts posterior summaries and constructs informative
priors. Uncertainty from upstream models propagates naturally into
growth parameter estimates.

**Why this matters:** - **Biological coherence:** Growth curves are
constrained to pass through maturity milestones - **Proper
uncertainty:** Imprecise maturity estimates → wider growth parameter
CIs - **Information efficiency:** All available data inform the final
estimates - **Reproducibility:** The entire workflow is encapsulated in
a few function calls

``` r
# Complete integrated workflow
birth_fit <- fit_bayesian_birth(embryo_lts, freeswim_lts)
L50_fit   <- fit_bayesian_maturity(maturity = "mat", lt = "fl", ...)
t50_fit   <- fit_bayesian_maturity(maturity = "mat", age = "age1", ...)

growth_fit <- fit_bayesian_growth(
  lt        = "fl",
  age       = "age1",
  data      = growth_data,
  birth_fit = birth_fit,  # Informs L₀
  L50_fit   = L50_fit,    # Informs Lmat
  t50_fit   = t50_fit,    # Informs tmat
  k_based   = FALSE       # Derive k from maturity
)
```

### 7. Automatic Multilingual Sex Coding

vitalBayes automatically detects sex coding in 11+ languages:

``` r
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

| Vignette                                                                                                          | Description                                        |
|-------------------------------------------------------------------------------------------------------------------|----------------------------------------------------|
| [`vignette("fit_bayesian_birth")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_birth.md)       | Birth size estimation                              |
| [`vignette("fit_bayesian_maturity")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_maturity.md) | L₅₀ and t₅₀ estimation                             |
| [`vignette("fit_bayesian_growth")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_growth.md)     | Growth models with maturity-based parameterization |
| [`vignette("partial_pooling")`](https://brian-j-moe.github.io/vitalBayes/articles/partial_pooling.md)             | Hierarchical modeling for imbalanced data          |
| [`vignette("visualization")`](https://brian-j-moe.github.io/vitalBayes/articles/visualization.md)                 | Plotting functions and color palettes              |
| [`vignette("model_diagnostics")`](https://brian-j-moe.github.io/vitalBayes/articles/model_diagnostics.md)         | Convergence, PPC, and LOO-CV                       |

For comprehensive statistical background, see the [Statistical Methods
Guide](https://brian-j-moe.github.io/vitalBayes/articles/vitalBayes_stats_explained.md).

## Citation

If you use vitalBayes in your research, please cite:

    @software{vitalBayes,
      author = {Moe, Brian J.},
      title = {vitalBayes: Bayesian Life History Parameter Estimation for Elasmobranchs},
      year = {2025},
      url = {https://github.com/Brian-J-Moe/vitalBayes},
      version = {0.2.1}
    }

## Dependencies

**Core:** - [cmdstanr](https://mc-stan.org/cmdstanr/) — CmdStan
interface - [instantiate](https://github.com/wlandau/instantiate) —
Precompiled Stan models -
[data.table](https://rdatatable.gitlab.io/data.table/) — Fast data
manipulation

**Visualization:** - [ggplot2](https://ggplot2.tidyverse.org/) —
Publication-quality graphics

**Model Comparison:** - [loo](https://mc-stan.org/loo/) — Leave-one-out
cross-validation

## Contributing

Contributions are welcome! Please open an issue to discuss proposed
changes or submit a pull request.

## License

[MIT License](https://brian-j-moe.github.io/vitalBayes/LICENSE.md)

------------------------------------------------------------------------

*vitalBayes is developed for the elasmobranch research community to
support robust, reproducible life history analyses.*

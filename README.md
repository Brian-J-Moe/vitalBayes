# vitalBayes <img src="man/figures/logo.png" align="right" height="139" />

## Bayesian Vital Rates for Elasmobranchs

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

`vitalBayes` provides a coherent Bayesian framework for estimating vital rates in elasmobranchs—from birth and maturity through growth to natural mortality and survival. The package implements hierarchical models using [Stan software](https://mc-stan.org/) via the [CmdStan R package](https://mc-stan.org/cmdstanr/reference/cmdstanr-package.html). For fast and reliable inference, the user will precompile several Stan models upon instillation.

## Installation

```r
# Install dependencies
install.packages(c("cmdstanr", "data.table", "ggplot2", "loo", "pracma"))

# Install cmdstan (if not already installed)
cmdstanr::install_cmdstan()

# Install vitalBayes from GitHub
pak::pak("Brian-J-Moe/vitalBayes")

# compile Stan models
vitalBayes::precompile_models()
```

## The Workflow

vitalBayes implements an integrated four-stage workflow where posterior distributions flow forward through life stages:

```
Birth  ──▶  Maturity  ──▶  Growth  ──▶  Mortality/Survival
(L₀)       (L₅₀, t₅₀)      (L∞, k)         (M, S(t))
```

```r
library(vitalBayes)
library(data.table)

data(growth_data)

# Stage 1: Birth size ─────────────────────────────────────────
birth_fit <- fit_bayesian_birth(
  embryo_lts        = growth_data[embryo == TRUE, fl],
  free_swimming_lts = growth_data[embryo == FALSE, fl]
)

# Stage 2: Maturity ───────────────────────────────────────────
mat_data <- growth_data[embryo == FALSE & !is.na(mat)]

L50_fit <- fit_bayesian_maturity(
  maturity = "mat", lt = "fl", sex = "sex",
  data = mat_data, use_pooling = TRUE
)

t50_fit <- fit_bayesian_maturity(
  maturity = "mat", age = "age", sex = "sex",
  data = mat_data[!is.na(age)], use_pooling = TRUE
)

# Stage 3: Growth ─────────────────────────────────────────────
growth_fit <- fit_bayesian_growth(
  lt = "fl", age = "age", sex = "sex",
  data = growth_data[embryo == FALSE & !is.na(age)],
  model   = "vb",
  k_based = FALSE,      # Maturity-based parameterization
  birth_stanfit         = birth_fit,
  length.mature_stanfit = L50_fit,
  age.mature_stanfit    = t50_fit,
  use_pooling = TRUE
)

# Stage 4: Mortality & Survival ───────────────────────────────
mort <- get_stochastic_mortality(
  method       = "CW",
  Linf = c(100, 8), L0 = c(30, 3), Lmat = c(70, 5), tmat = c(12, 2),
  growth_model = "vb",
  iter = 2000, scaled = TRUE, p = 0.001
)

surv <- simulate_survivorship(mc_object = mort, n = 50000, n_iter = 2000)
surv$Aggregate$Age_of_Death
surv$Aggregate$Survival_to_tmat
```

## Key Methodological Features

### Maturity-Based Growth Parameterization

The growth coefficient *k* is derived from observable maturity milestones (L₅₀, t₅₀) rather than estimated directly. This breaks the notorious L∞–k correlation and anchors the growth curve to data within the observed range. Supports von Bertalanffy, Gompertz, and Logistic models.

### Partial Pooling for Imbalanced Data

Hierarchical models borrow strength across sexes when sample sizes are unequal—common in elasmobranch research. The sparse sex gets regularized estimates without assuming identical parameters. Selective pooling prevents "double-pooling" when maturity priors are themselves pooled.

### Unified Mortality Framework

Natural mortality uses a consistent formulation across growth models: M(t) = M∞ / G(t), where M∞ is computed from biological milestones and G(t) uses the native growth trajectory. This ensures comparable survival predictions regardless of which growth model best fits your data.

### Probit Link for Thresholds

Birth size and maturity ogives use probit rather than logit links, providing direct interpretation as cumulative normal distributions of developmental thresholds.

### CV-Based Priors

Prior uncertainty is specified via coefficient of variation—intuitive, scale-invariant, and easily communicated (e.g., "20% uncertainty on L∞").

## Supported Models

| Component | Options |
|-----------|---------|
| **Growth** | von Bertalanffy, Gompertz, Logistic |
| **Mortality** | Chen-Watanabe (with optional two-phase senescence), Peterson-Wroblewski, Lorenzen |
| **Pooling** | None, partial (hierarchical), or selective |

## Example Datasets
 
| Dataset | Sex Ratio | Use Case |
|---------|-----------|----------|
| `growth_data` | 189F, 176M | Main workflow examples |
| `imbalanced_data` | 150F, 34M | Partial pooling demonstrations |
| `limited_data` | 24F, 18M | Prior sensitivity analysis |


## Vignettes

| Topic | Vignette |
|-------|----------|
| Birth size | `vignette("fit_bayesian_birth")` |
| Maturity ogives | `vignette("fit_bayesian_maturity")` |
| Growth models | `vignette("fit_bayesian_growth")` |
| Partial pooling | `vignette("partial_pooling")` |
| Mortality estimation | `vignette("mortality_estimation")` |
| Mortality framework | `vignette("chen_watanabe_reparameterization")` |
| Survival simulation | `vignette("survivorship_simulation")` |
| Diagnostics | `vignette("model_diagnostics")` |

## Citation

```bibtex
@manual{vitalBayes,
  title = {{vitalBayes}: Bayesian Vital Rate Estimation for Elasmobranchs},
  author = {Brian J. Moe},
  year = {2025},
  note = {R package version 0.7.0},
  url = {https://github.com/Brian-J-Moe/vitalBayes}
}
```

## License

[MIT](LICENSE.md)



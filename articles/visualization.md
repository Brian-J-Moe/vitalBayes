# Visualizing vitalBayes Results

## Overview

vitalBayes provides publication-ready visualization functions for all
model types, with built-in support for colorblind-safe palettes and
multilingual labels.

## Color Palettes

``` r
library(vitalBayes)

# List available palettes
list_vital_palettes()

# Default retro-wave palette
vital_colors(5, "default")

# Two-sex palette (pink/cyan)
vital_colors(2, "sex")

# Colorblind-safe alternatives
vital_colors(5, "okabe_ito")
vital_colors(2, "sex_cb")        # Orange/blue for sex
vital_colors(5, "viridis")
vital_colors(5, "tol_bright")

# Check if a palette is colorblind-safe
is_colorblind_safe("okabe_ito")  # TRUE
is_colorblind_safe("sex")        # FALSE
```

## The vitalBayes Theme

All plots use
[`theme_vital()`](https://brian-j-moe.github.io/vitalBayes/reference/theme_vital.md),
a clean ggplot2 theme optimized for publications:

``` r
library(ggplot2)

# Use standalone
ggplot(mtcars, aes(mpg, wt)) +
 geom_point() +
 theme_vital()

# Customize
ggplot(mtcars, aes(mpg, wt)) +
 geom_point() +
 theme_vital(base_size = 14) +
 theme(legend.position = "bottom")
```

------------------------------------------------------------------------

## plot_birth_ogive()

Visualize birth size transition probability:

``` r
library(data.table)
data(gulper_data)

# Fit birth model
birth_fit <- fit_bayesian_birth(
 embryo_lts = gulper_data[embryo == TRUE, fl],
 free_swimming_lts = gulper_data[embryo == FALSE, fl]
)

# Basic plot
plot_birth_ogive(
 fit = birth_fit,
 embryo_lengths = gulper_data[embryo == TRUE, fl],
 freeswim_lengths = gulper_data[embryo == FALSE, fl]
)

# Customized
plot_birth_ogive(
 fit = birth_fit,
 embryo_lengths = gulper_data[embryo == TRUE, fl],
 freeswim_lengths = gulper_data[embryo == FALSE, fl],
 show_data = TRUE,
 show_b50_line = TRUE,
 colors = vital_colors(1, "sunset"),
 x_lab = "Fork Length (cm)",
 title = "Little Gulper Shark Birth Size"
)
```

------------------------------------------------------------------------

## plot_maturity_ogive()

Visualize maturity probability curves:

``` r
# Fit two-sex maturity model
mat_data <- gulper_data[embryo == FALSE & !is.na(mat)]

L50_fit <- fit_bayesian_maturity(
 maturity = "mat",
 lt = "fl",
 sex = "sex",
 data = mat_data,
 use_pooling = TRUE
)

# Basic usage
plot_maturity_ogive(
 fit = L50_fit,
 type = "length",
 data = mat_data,
 x_col = "fl",
 maturity_col = "mat",
 sex_col = "sex"
)
```

### Key Arguments

``` r
plot_maturity_ogive(
 fit = L50_fit,
 type = "length",           # or "age"
 data = mat_data,
 
 # Data columns
 x_col = "fl",
 maturity_col = "mat",
 sex_col = "sex",
 
 # Display options
 show_data = TRUE,          # Overlay observed points
 show_rug = TRUE,           # Rug plot at top/bottom
 show_x50_line = TRUE,      # Vertical line at x50
 
 # Appearance
 colors = NULL,             # Auto from vital_colors("sex")
 colorblind = FALSE,        # Use colorblind-safe palette
 facet_sex = TRUE,          # Separate panels by sex
 
 # Labels
 x_lab = "Fork Length (cm)",
 y_lab = "Probability of Maturity",
 title = "Length-at-Maturity Ogive"
)
```

### Multilingual Labels

``` r
# French
plot_maturity_ogive(
 fit = L50_fit,
 type = "length",
 data = mat_data,
 sex_labels = c("1" = "Femelle", "2" = "Mâle"),
 x_lab = "Longueur à la fourche (cm)",
 y_lab = "Probabilité de maturité"
)

# Spanish
plot_maturity_ogive(
 fit = L50_fit,
 type = "length",
 data = mat_data,
 sex_labels = c("1" = "Hembra", "2" = "Macho"),
 x_lab = "Longitud furcal (cm)",
 y_lab = "Probabilidad de madurez"
)

# Japanese
plot_maturity_ogive(
 fit = L50_fit,
 type = "length",
 data = mat_data,
 sex_labels = c("1" = "メス", "2" = "オス"),
 x_lab = "尾叉長 (cm)",
 y_lab = "成熟確率"
)
```

------------------------------------------------------------------------

## plot_growth_curve()

Visualize fitted growth curves with uncertainty:

``` r
growth_data <- gulper_data[embryo == FALSE & !is.na(age1)]

# Fit growth model
growth_fit <- fit_bayesian_growth(
 lt = "fl",
 age = "age1",
 sex = "sex",
 data = growth_data,
 model = "vb",
 k_based = FALSE,
 L50_fit = L50_fit
)

# Basic plot
plot_growth_curve(
 fit = growth_fit,
 data = growth_data,
 age_col = "age1",
 length_col = "fl",
 sex_col = "sex"
)
```

### Key Arguments

``` r
plot_growth_curve(
 fit = growth_fit,
 data = growth_data,
 
 # Data columns
 age_col = "age1",
 length_col = "fl",
 sex_col = "sex",
 
 # Prediction range
 age_range = c(0, 50),      # Override auto-detection
 n_points = 100,            # Points for smooth curve
 
 # Display options
 show_data = TRUE,          # Overlay observed points
 show_50_ci = TRUE,         # Show 50% credible interval
 facet_sex = TRUE,          # Separate panels by sex
 
 # Appearance
 data_alpha = 0.4,
 ribbon_alpha = 0.2,
 ribbon_alpha_50 = 0.3,
 line_width = 1.2,
 colors = NULL,             # Auto from vital_colors("sex")
 colorblind = FALSE,
 
 # Model specification (must match fitting call)
 k_based = FALSE,
 which_model = 1,           # 1=VB, 2=Gompertz, 3=Logistic
 
 # Labels
 x_lab = "Age (years)",
 y_lab = "Fork Length (cm)",
 title = NULL
)
```

### Additional Layers

Add custom ggplot2 layers:

``` r
library(ggplot2)

plot_growth_curve(
 fit = growth_fit,
 data = growth_data,
 additional_layers = list(
   # Add horizontal line at L50
   geom_hline(yintercept = 75, linetype = "dashed", color = "gray50"),
   # Add annotation
   annotate("text", x = 5, y = 77, label = "L50", size = 3)
 ),
 theme_args = list(
   legend.position = "bottom",
   plot.title = element_text(hjust = 0.5)
 )
)
```

------------------------------------------------------------------------

## compare_growth_models()

Side-by-side comparison of multiple growth models:

``` r
# Fit multiple models
vb_fit <- fit_bayesian_growth(
 lt = "fl", age = "age1", sex = "sex", data = growth_data,
 model = "vb", k_based = FALSE, L50_fit = L50_fit
)

gomp_fit <- fit_bayesian_growth(
 lt = "fl", age = "age1", sex = "sex", data = growth_data,
 model = "gompertz", k_based = FALSE, L50_fit = L50_fit
)

logis_fit <- fit_bayesian_growth(
 lt = "fl", age = "age1", sex = "sex", data = growth_data,
 model = "logistic", k_based = FALSE, L50_fit = L50_fit
)

# Compare visually
compare_growth_models(
 "von Bertalanffy" = vb_fit,
 "Gompertz" = gomp_fit,
 "Logistic" = logis_fit,
 data = growth_data,
 age_col = "age1",
 length_col = "fl",
 sex_col = "sex",
 k_based_list = c(FALSE, FALSE, FALSE),
 which_model_list = c(1, 2, 3)
)

# Colorblind-safe
compare_growth_models(
 "von Bertalanffy" = vb_fit,
 "Gompertz" = gomp_fit,
 data = growth_data,
 colorblind = TRUE
)
```

------------------------------------------------------------------------

## plot_residuals()

Diagnostic residual plots for growth models:

``` r
# All diagnostic plots
plot_residuals(
 fit = growth_fit,
 data = growth_data,
 age_col = "age1",
 length_col = "fl",
 sex_col = "sex",
 type = "all"              # Combined panel
)

# Individual plot types
plot_residuals(fit = growth_fit, data = growth_data, type = "fitted")
plot_residuals(fit = growth_fit, data = growth_data, type = "qq")
plot_residuals(fit = growth_fit, data = growth_data, type = "histogram")
plot_residuals(fit = growth_fit, data = growth_data, type = "age")
```

------------------------------------------------------------------------

## Saving Publication Figures

``` r
library(patchwork)

# Create multi-panel figure
p1 <- plot_maturity_ogive(L50_fit, type = "length", data = mat_data)
p2 <- plot_growth_curve(growth_fit, data = growth_data)

combined <- p1 / p2 + plot_annotation(tag_levels = "A")

# Save at publication resolution
ggsave("Figure1.pdf", combined, width = 10, height = 12)
ggsave("Figure1.png", combined, width = 10, height = 12, dpi = 300)
```

## See Also

- [Statistical Methods
  Guide](https://brian-j-moe.github.io/vitalBayes/articles/vitalBayes_stats_explained.md)
  — Mathematical background for all models
- [`vital_colors()`](https://brian-j-moe.github.io/vitalBayes/reference/vital_colors.md),
  [`theme_vital()`](https://brian-j-moe.github.io/vitalBayes/reference/theme_vital.md)
  — Styling utilities
- [`list_vital_palettes()`](https://brian-j-moe.github.io/vitalBayes/reference/list_vital_palettes.md),
  [`is_colorblind_safe()`](https://brian-j-moe.github.io/vitalBayes/reference/is_colorblind_safe.md)
  — Palette information
- [`vignette("fit_bayesian_maturity")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_maturity.md)
  — Maturity model fitting
- [`vignette("fit_bayesian_growth")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_growth.md)
  — Growth model fitting

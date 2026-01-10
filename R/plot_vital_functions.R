# =============================================================================
# vitalBayes Visualization Functions
# =============================================================================
# Publication-quality plotting functions for growth curves, maturity ogives,
# and birth models. Uses ggplot2 with customizable theming.
# =============================================================================

#' Vital Package Color Palette
#'
#' @description
#' Returns color palettes for vitalBayes plots. The default palette is derived
#' directly from the package hex sticker's retro-wave/synthwave aesthetic.
#' Colorblind-friendly alternatives are provided based on established research
#' (Okabe-Ito, viridis, and Paul Tol's schemes).
#'
#' @param n Integer. Number of colors to return. Default returns all available
#'   colors in the selected palette.
#' @param palette Character. Which palette to use:
#'   \describe{
#'     \item{\code{"default"}}{Full retro-wave palette from logo (5 colors)}
#'     \item{\code{"sex"}}{Female (magenta) / Male (cyan) for two-sex models}
#'     \item{\code{"gradient"}}{Smooth gradient from pink to cyan (6 colors)}
#'     \item{\code{"sunset"}}{Warm sunset tones from the logo's sun element}
#'     \item{\code{"okabe_ito"}}{Colorblind-safe Okabe-Ito palette (8 colors)}
#'     \item{\code{"viridis"}}{Colorblind-safe viridis-inspired (5 colors)}
#'     \item{\code{"tol_bright"}}{Paul Tol's colorblind-safe bright scheme (7 colors)}
#'     \item{\code{"sex_cb"}}{Colorblind-safe sex palette (orange/blue)}
#'   }
#'
#' @return Character vector of hex color codes.
#'
#' @section Colorblind Accessibility:
#' The \code{"okabe_ito"}, \code{"viridis"}, \code{"tol_bright"}, and \code{"sex_cb"}
#' palettes are designed to be distinguishable by individuals with the most common
#' forms of color vision deficiency (deuteranopia, protanopia, tritanopia).
#'
#' For two-sex models in publications requiring accessibility, use \code{"sex_cb"}
#' which employs orange (female) and blue (male) - the most robust color combination
#' for all forms of color blindness.
#'
#' @references
#' Okabe M, Ito K (2008). "Color Universal Design (CUD): How to make figures and
#' presentations that are friendly to Colorblind people."
#' \url{https://jfly.uni-koeln.de/color/}
#'
#' Tol P (2021). "Colour Schemes." SRON Technical Note.
#' \url{https://personal.sron.nl/~pault/}
#'
#' @examples
#' # Default vitalBayes palette
#' vital_colors()
#'
#' # Two-sex models
#' vital_colors(2, "sex")       # Default pink/cyan
#' vital_colors(2, "sex_cb")    # Colorblind-safe orange/blue
#'
#' # Colorblind-safe palettes for multi-group comparisons
#' vital_colors(4, "okabe_ito")
#' vital_colors(5, "viridis")
#'
#' # Preview palettes
#' barplot(rep(1, 5), col = vital_colors(5), border = NA, main = "Default")
#' barplot(rep(1, 8), col = vital_colors(8, "okabe_ito"), border = NA, main = "Okabe-Ito")
#'
#' @export
vital_colors <- function(n = 5, palette = c("default", "sex", "gradient", "sunset",
                                             "okabe_ito", "viridis", "tol_bright",
                                             "sex_cb")) {

  palette <- match.arg(palette)

  colors <- switch(palette,

    # =========================================================================
    # vitalBayes Branding Palettes (from hex sticker)
    # =========================================================================

    # Primary palette - full retro-wave aesthetic
    default = c(
      "#E930A7",  # Vivid magenta (shark/primary accent)
      "#00CFFF",  # Electric cyan (growth curve/secondary)
      "#7B2D8E",  # Deep purple (background gradient)
      "#FF6B35",  # Sunset orange (sun element)
      "#F72585"   # Hot pink (gradient highlight)
    ),

    # Two-sex models: Female = magenta, Male = cyan
    sex = c(
      "#E930A7",  # Female - vivid magenta
      "#00CFFF"   # Male - electric cyan
    ),

    # Smooth gradient for continuous scales
    gradient = c(
      "#F72585",  # Hot pink
      "#B5179E",  # Purple-pink
      "#7209B7",  # Purple
      "#3A0CA3",  # Deep purple
      "#4361EE",  # Blue
      "#4CC9F0"   # Cyan
    ),

    # Warm sunset tones (for age/maturity visualizations)
    sunset = c(
      "#FF6B35",  # Bright orange
      "#F7931E",  # Golden orange
      "#FFB627",  # Amber
      "#FF006E",  # Magenta
      "#8338EC"   # Purple
    ),

    # =========================================================================
    # Colorblind-Friendly Palettes
    # =========================================================================

    # Okabe-Ito palette - excellent colorblind accessibility
    # Reference: https://jfly.uni-koeln.de/color/
    okabe_ito = c(
      "#E69F00",  # Orange
      "#56B4E9",  # Sky blue
      "#009E73",  # Bluish green
      "#F0E442",  # Yellow
      "#0072B2",  # Blue
      "#D55E00",  # Vermillion
      "#CC79A7",  # Reddish purple
      "#000000"   # Black
    ),

    # Viridis-inspired discrete palette
    # Perceptually uniform, colorblind-safe
    viridis = c(
      "#440154",  # Dark purple
      "#3B528B",  # Blue-purple
      "#21918C",  # Teal
      "#5DC863",  # Green
      "#FDE725"   # Yellow
    ),

    # Paul Tol's bright scheme - high contrast, colorblind-safe
    # Reference: https://personal.sron.nl/~pault/
    tol_bright = c(
      "#4477AA",  # Blue
      "#EE6677",  # Red
      "#228833",  # Green
      "#CCBB44",  # Yellow
      "#66CCEE",  # Cyan
      "#AA3377",  # Purple
      "#BBBBBB"   # Grey
    ),

    # Colorblind-safe sex palette
    # Orange and blue are the most robust pair for all color vision types
    sex_cb = c(
      "#E69F00",  # Female - orange (warm)
      "#0072B2"   # Male - blue (cool)
    )
  )

  # Expand palette if more colors needed
  if (n > length(colors)) {
    colors <- grDevices::colorRampPalette(colors)(n)
  }

  return(colors[1:n])
}

#' vitalBayes Extended Color Palette
#'
#' @description
#' Returns color palettes for vitalBayes visualizations with support for
#' synthwave (default), viridis, and colorblind-friendly (Okabe-Ito) schemes.
#'
#' @param n Number of colors needed. If NULL, returns the full palette.
#' @param type Character. One of:
#'   \itemize{
#'     \item \code{"synthwave"} - Retro-wave palette from vitalBayes hex sticker (default)
#'     \item \code{"viridis"} - Viridis colorblind-friendly sequential palette
#'     \item \code{"okabe"} - Okabe-Ito colorblind-friendly qualitative palette
#'     \item \code{"plasma"} - Plasma sequential palette
#'     \item \code{"inferno"} - Inferno sequential palette
#'   }
#' @param alpha Numeric. Transparency value (0-1). Default 1 (opaque).
#'
#' @return A character vector of hex color codes.
#'
#' @examples
#' \dontrun{
#' vital_palette(4, "synthwave")
#' vital_palette(6, "okabe")
#' vital_palette(10, "viridis", alpha = 0.8)
#' }
#'
#' @importFrom grDevices colorRampPalette adjustcolor
#' @export
vital_palette <- function(n = NULL, type = c("synthwave", "viridis", "okabe", "plasma", "inferno"),
                          alpha = 1) {

  type <- match.arg(type)

  # Define base palettes
  palettes <- list(
    synthwave = c(
      vital_pink    = "#FF1D8E",
      bayes_blue    = "#00D9FF",
      sunset_purple = "#5B3A8F",
      deep_ocean    = "#1A1F3A",
      coral_orange  = "#FF6B35",
      teal_accent   = "#4ECDC4",
      hot_magenta   = "#E930A7",
      royal_purple  = "#7B2D8E"
    ),
    okabe = c(
      orange    = "#E69F00",
      sky_blue  = "#56B4E9",
      green     = "#009E73",
      yellow    = "#F0E442",
      blue      = "#0072B2",
      vermilion = "#D55E00",
      purple    = "#CC79A7",
      black     = "#000000"
    )
  )

  # Handle viridis family
  if (type %in% c("viridis", "plasma", "inferno")) {
    if (!requireNamespace("viridis", quietly = TRUE)) {
      stop("Package 'viridis' required for this palette type. Install with install.packages('viridis').")
    }
    if (is.null(n)) n <- 6

    pal <- switch(type,
                  viridis = viridis::viridis(n, alpha = alpha),
                  plasma  = viridis::plasma(n, alpha = alpha),
                  inferno = viridis::inferno(n, alpha = alpha))
    return(pal)
  }

  # Named palettes
  base <- unname(palettes[[type]])

  if (is.null(n)) {
    if (alpha < 1) return(grDevices::adjustcolor(base, alpha.f = alpha))
    return(base)
  }

  # Interpolate or subset
  if (n <= length(base)) {
    pal <- base[1:n]
  } else {
    pal <- grDevices::colorRampPalette(base)(n)
  }

  if (alpha < 1) pal <- grDevices::adjustcolor(pal, alpha.f = alpha)
  pal
}


#' Check if Palette is Colorblind-Friendly
#'
#' @description
#' Utility function to check if a palette name is colorblind-accessible.
#'
#' @param palette Character. Palette name to check.
#'
#' @return Logical. TRUE if colorblind-friendly.
#'
#' @examples
#' is_colorblind_safe("okabe_ito")  # TRUE
#' is_colorblind_safe("default")    # FALSE
#'
#' @export
is_colorblind_safe <- function(palette) {
  cb_palettes <- c("okabe_ito", "viridis", "tol_bright", "sex_cb")
  tolower(palette) %in% cb_palettes
}


#' List Available Color Palettes
#'
#' @description
#' Prints information about all available color palettes in vitalBayes,
#' including their colorblind accessibility status.
#'
#' @param show_colors Logical. Display color swatches? Default \code{FALSE}.
#'   Requires a graphics device if \code{TRUE}.
#'
#' @return Invisibly returns a data.frame with palette information.
#'
#' @examples
#' list_vital_palettes()
#'
#' @export
list_vital_palettes <- function(show_colors = FALSE) {

  palettes <- data.frame(
    palette = c("default", "sex", "gradient", "sunset",
                "okabe_ito", "viridis", "tol_bright", "sex_cb"),
    n_colors = c(5, 2, 6, 5, 8, 5, 7, 2),
    colorblind_safe = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
    description = c(
      "vitalBayes logo retro-wave",
      "Two-sex (magenta/cyan)",
      "Continuous gradient (pink to cyan)",
      "Warm sunset tones",
      "Okabe-Ito universal design",
      "Viridis perceptually uniform",
      "Paul Tol bright scheme",
      "Two-sex colorblind-safe (orange/blue)"
    ),
    stringsAsFactors = FALSE
  )

  message("\nvitalBayes Color Palettes")
  message(paste(rep("=", 50), collapse = ""))

  for (i in 1:nrow(palettes)) {
    cb_marker <- if (palettes$colorblind_safe[i]) " [CB-safe]" else ""
    message(sprintf("  %-12s (%d colors)%s",
                    palettes$palette[i],
                    palettes$n_colors[i],
                    cb_marker))
    message(sprintf("    %s\n", palettes$description[i]))
  }

  if (show_colors) {
    grDevices::dev.new(width = 10, height = 6)
    graphics::par(mfrow = c(4, 2), mar = c(1, 1, 2, 1))

    for (pal in palettes$palette) {
      n <- palettes$n_colors[palettes$palette == pal]
      cols <- vital_colors(n, pal)
      graphics::barplot(rep(1, n), col = cols, border = NA,
                        main = pal, axes = FALSE)
    }
  }

  invisible(palettes)
}


#' vitalBayes ggplot2 Theme
#'
#' @description
#' A publication-ready ggplot2 theme inspired by the vitalBayes retro-wave
#' aesthetic. Provides a clean, modern look suitable for scientific publications
#' while maintaining visual consistency with the package branding.
#'
#' @param base_size Numeric. Base font size. Default 12.
#' @param base_family Character. Base font family. Default "".
#' @param grid Logical. Show grid lines? Default \code{TRUE} for major only.
#' @param dark Logical. Use dark theme variant? Default \code{FALSE}.
#'
#' @return A ggplot2 theme object.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point(color = vital_colors(1)) +
#'   theme_vital()
#' }
#'
#' @import ggplot2
#' @export
theme_vital <- function(base_size = 12, base_family = "", grid = TRUE, dark = FALSE) {

  if (dark) {
    # Dark variant - closer to logo background
    bg_color <- "#1A1A2E"
    text_color <- "#E0E0E0"
    grid_color <- "#3D3D5C"
    axis_color <- "#6B6B8D"
  } else {
    # Light variant - publication-friendly
    bg_color <- "#FFFFFF"
    text_color <- "#2D2D2D"
    grid_color <- "#E0E0E0"
    axis_color <- "#4A4A4A"
  }

  theme_base <- ggplot2::theme_bw(base_size = base_size, base_family = base_family)

  theme_base +
    ggplot2::theme(
      # Background
      plot.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.background = ggplot2::element_rect(fill = bg_color, color = NA),

      # Grid
      panel.grid.major = if (grid) {
        ggplot2::element_line(color = grid_color, linewidth = 0.3)
      } else {
        ggplot2::element_blank()
      },
      panel.grid.minor = ggplot2::element_blank(),

      # Axes
      axis.line = ggplot2::element_line(color = axis_color, linewidth = 0.5),
      axis.ticks = ggplot2::element_line(color = axis_color, linewidth = 0.3),
      axis.text = ggplot2::element_text(color = text_color),
      axis.title = ggplot2::element_text(color = text_color, face = "bold"),

      # Panel border
      panel.border = ggplot2::element_rect(color = axis_color, fill = NA, linewidth = 0.5),

      # Legend
      legend.background = ggplot2::element_rect(fill = bg_color, color = NA),
      legend.key = ggplot2::element_rect(fill = bg_color, color = NA),
      legend.text = ggplot2::element_text(color = text_color),
      legend.title = ggplot2::element_text(color = text_color, face = "bold"),

      # Title/subtitle
      plot.title = ggplot2::element_text(
        color = if (dark) "#E930A7" else "#7B2D8E",
        face = "bold",
        size = base_size * 1.2,
        hjust = 0
      ),
      plot.subtitle = ggplot2::element_text(
        color = text_color,
        size = base_size * 0.9,
        hjust = 0
      ),

      # Strip for facets
      strip.background = ggplot2::element_rect(
        fill = if (dark) "#3D3D5C" else "#F0F0F0",
        color = axis_color
      ),
      strip.text = ggplot2::element_text(
        color = text_color,
        face = "bold"
      )
    )
}


# -----------------------------------------------------------------------------
# Growth Curve Visualization
# -----------------------------------------------------------------------------

#' Generate Posterior Predictions for Growth Curves
#'
#' @description
#' Internal function to generate growth curve predictions from posterior draws.
#'
#' @param fit CmdStanMCMC object from fit_bayesian_growth().
#' @param age_seq Numeric vector of ages for predictions.
#' @param sex_code Integer. Sex code (1=female, 2=male) for two-sex models.
#' @param k_based Logical. Was the model k-based?
#' @param which_model Integer. 1=VBGM, 2=Gompertz, 3=Logistic.
#' @param n_draws Integer. Number of posterior draws to use. Default 500.
#'
#' @return data.table with columns: age, median, lower, upper.
#'
#' @noRd
.posterior_predict_growth <- function(
    fit,
    age_seq,
    sex_code = 1L,
    k_based = FALSE,
    which_model = 1L,
    n_draws = 500
) {

  # Extract draws
  Linf_draws <- fit$draws("Linf", format = "matrix")
  L0_draws <- fit$draws("L0", format = "matrix")

  is_hierarchical <- ncol(Linf_draws) > 1

  # Subset draws for efficiency
  total_draws <- nrow(Linf_draws)
  draw_idx <- sample.int(total_draws, min(n_draws, total_draws))

  # Get appropriate column for sex
  s <- if (is_hierarchical) sex_code else 1L

  Linf <- Linf_draws[draw_idx, s]
  L0 <- L0_draws[draw_idx, s]

  if (k_based) {
    k_draws <- fit$draws("k", format = "matrix")
    k <- k_draws[draw_idx, s]
  } else {
    Lmat_draws <- fit$draws("Lmat", format = "matrix")
    tmat_draws <- fit$draws("tmat", format = "matrix")
    Lmat <- Lmat_draws[draw_idx, s]
    tmat <- tmat_draws[draw_idx, s]

    # Derive k
    k <- .derive_k_from_maturity(Linf, L0, Lmat, tmat, which_model)
  }

  # Compute predictions
  n_age <- length(age_seq)
  pred_mat <- matrix(NA_real_, nrow = length(draw_idx), ncol = n_age)

  for (i in seq_along(age_seq)) {
    pred_mat[, i] <- .growth_curve(Linf, L0, k, age_seq[i], which_model)
  }

  # Summarize
  data.table::data.table(
    age = age_seq,
    median = apply(pred_mat, 2, stats::median),
    lower = apply(pred_mat, 2, stats::quantile, 0.025),
    upper = apply(pred_mat, 2, stats::quantile, 0.975),
    lower_50 = apply(pred_mat, 2, stats::quantile, 0.25),
    upper_50 = apply(pred_mat, 2, stats::quantile, 0.75)
  )
}


#' Plot Growth Curves
#'
#' @description
#' Creates publication-quality growth curve plots with credible intervals.
#' Automatically handles single-sex and two-sex models.
#'
#' @param fit A CmdStanMCMC object from \code{\link{fit_bayesian_growth}}.
#' @param data Optional data.table/data.frame with observed data for overlay.
#' @param age_col Character. Column name for age in \code{data}. Default "age".
#' @param length_col Character. Column name for length in \code{data}. Default "length".
#' @param sex_col Character. Column name for sex in \code{data}. Default "sex".
#' @param sex_labels Named character vector. Labels for sex in plot legend and facets.
#'   Default \code{c("1" = "Female", "2" = "Male")}. Can be customized for other
#'   languages, e.g., \code{c("1" = "Femelle", "2" = "Mâle")} (French),
#'   \code{c("1" = "Hembra", "2" = "Macho")} (Spanish), or
#'   \code{c("1" = "メス", "2" = "オス")} (Japanese).
#' @param age_range Numeric vector of length 2. Age range for predictions.
#'   If \code{NULL}, auto-determined from data or defaults to c(0, 30).
#' @param n_points Integer. Number of points for prediction grid. Default 100.
#' @param show_data Logical. Overlay observed data points? Default \code{TRUE}.
#' @param show_50_ci Logical. Show 50% credible interval ribbon? Default \code{TRUE}.
#' @param data_alpha Numeric. Transparency for data points. Default 0.4.
#' @param data_size Numeric. Size for data points. Default 1.5.
#' @param ribbon_alpha Numeric. Transparency for 95% CI ribbon. Default 0.2.
#' @param ribbon_alpha_50 Numeric. Transparency for 50% CI ribbon. Default 0.3.
#' @param line_width Numeric. Width of median line. Default 1.2.
#' @param colors Character vector of colors. If \code{NULL}, uses \code{vital_colors()}.
#' @param colorblind Logical. Use colorblind-safe palette? Default \code{FALSE}.
#'   When \code{TRUE}, overrides \code{colors} with \code{vital_colors(, "sex_cb")}
#'   for two-sex models or \code{vital_colors(, "okabe_ito")} otherwise.
#' @param facet_sex Logical. For two-sex models, facet by sex? Default \code{TRUE}.
#' @param x_lab,y_lab Character. Axis labels.
#' @param title Character. Plot title. Default \code{NULL}.
#' @param k_based Logical. Was the model k-based? Default \code{FALSE}.
#' @param which_model Integer. 1=VBGM, 2=Gompertz, 3=Logistic. Default 1.
#' @param base_size Numeric. Base font size for theme. Default 12.
#' @param theme_args List of additional arguments for \code{ggplot2::theme()}.
#' @param additional_layers List of additional ggplot2 layers to add.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' plot_growth_curve(growth_fit, data = shark_data)
#'
#' # French labels
#' plot_growth_curve(
#'   growth_fit,
#'   data = shark_data,
#'   sex_labels = c("1" = "Femelle", "2" = "Mâle")
#' )
#'
#' # Spanish labels
#' plot_growth_curve(
#'   growth_fit,
#'   data = shark_data,
#'   sex_labels = c("1" = "Hembra", "2" = "Macho")
#' )
#'
#' # Customize appearance
#' plot_growth_curve(
#'   growth_fit,
#'   data = shark_data,
#'   colors = c("#E63946", "#1D3557"),
#'   theme_args = list(legend.position = "bottom"),
#'   additional_layers = list(
#'     geom_hline(yintercept = 100, linetype = "dashed")
#'   )
#' )
#' }
#'
#' @seealso
#' \code{vignette("visualization")} for comprehensive plotting examples.
#'
#' \code{vignette("fit_bayesian_growth")} for growth model fitting.
#'
#' \code{\link{compare_growth_models}}, \code{\link{plot_residuals}},
#' \code{\link{vital_colors}}, \code{\link{theme_vital}}
#'
#' @import ggplot2
#' @import data.table
#' @export
plot_growth_curve <- function(
    fit,
    data          = NULL,
    age_col       = "age",
    length_col    = "length",
    sex_col       = "sex",
    sex_labels    = c("1" = "Female", "2" = "Male"),
    age_range     = NULL,
    n_points      = 100,
    show_data     = TRUE,
    show_50_ci    = TRUE,
    data_alpha    = 0.4,
    data_size     = 1.5,
    ribbon_alpha  = 0.2,
    ribbon_alpha_50 = 0.3,
    line_width    = 1.2,
    colors        = NULL,
    colorblind    = FALSE,
    facet_sex     = TRUE,
    x_lab         = "Age (years)",
    y_lab         = "Length (cm)",
    title         = NULL,
    k_based       = FALSE,
    which_model   = 1L,
    base_size     = 12,
    theme_args    = list(),
    additional_layers = list()
) {

  # Check if hierarchical
  Linf_test <- fit$draws("Linf", format = "matrix")
  is_hierarchical <- ncol(Linf_test) > 1

  # Determine age range
  if (is.null(age_range)) {
    if (!is.null(data) && age_col %in% names(data)) {
      age_range <- range(data[[age_col]], na.rm = TRUE)
      age_range[1] <- max(0, age_range[1])
      age_range[2] <- age_range[2] * 1.05
    } else {
      age_range <- c(0, 30)
    }
  }

  age_seq <- seq(age_range[1], age_range[2], length.out = n_points)

  # Generate predictions
  if (is_hierarchical) {
    pred_f <- .posterior_predict_growth(fit, age_seq, 1L, k_based, which_model)
    pred_f[, sex := sex_labels["1"]]

    pred_m <- .posterior_predict_growth(fit, age_seq, 2L, k_based, which_model)
    pred_m[, sex := sex_labels["2"]]

    pred_dt <- data.table::rbindlist(list(pred_f, pred_m))

    # Set factor levels to maintain order
    pred_dt[, sex := factor(sex, levels = c(sex_labels["1"], sex_labels["2"]))]
  } else {
    pred_dt <- .posterior_predict_growth(fit, age_seq, 1L, k_based, which_model)
    pred_dt[, sex := "Combined"]
  }

  # Set up colors - handle colorblind option
  if (colorblind) {
    colors <- vital_colors(2, "sex_cb")
  } else if (is.null(colors)) {
    colors <- vital_colors(2, "sex")
  }

  # Build plot
  p <- ggplot2::ggplot()

  # Add additional layers first (behind ribbons)
  for (layer in additional_layers) {
    if (inherits(layer, "Layer") || inherits(layer, "gg")) {
      p <- p + layer
    }
  }

  # Add 95% CI ribbon
  p <- p + ggplot2::geom_ribbon(
    data = pred_dt,
    ggplot2::aes(x = age, ymin = lower, ymax = upper, fill = sex),
    alpha = ribbon_alpha
  )

  # Add 50% CI ribbon
  if (show_50_ci) {
    p <- p + ggplot2::geom_ribbon(
      data = pred_dt,
      ggplot2::aes(x = age, ymin = lower_50, ymax = upper_50, fill = sex),
      alpha = ribbon_alpha_50
    )
  }

  # Add median line
  p <- p + ggplot2::geom_line(
    data = pred_dt,
    ggplot2::aes(x = age, y = median, color = sex),
    linewidth = line_width
  )

  # Add observed data
  if (show_data && !is.null(data)) {
    data_dt <- data.table::as.data.table(data)

    if (is_hierarchical && sex_col %in% names(data_dt)) {
      # Map sex values to labels (handles both numeric and character)
      sex_vals <- data_dt[[sex_col]]
      if (is.numeric(sex_vals)) {
        data_dt[, sex := factor(
          ifelse(get(sex_col) == 1, sex_labels["1"], sex_labels["2"]),
          levels = c(sex_labels["1"], sex_labels["2"])
        )]
      } else {
        data_dt[, sex := factor(
          ifelse(grepl("^[Ff]", get(sex_col)), sex_labels["1"], sex_labels["2"]),
          levels = c(sex_labels["1"], sex_labels["2"])
        )]
      }
    } else {
      data_dt[, sex := "Combined"]
    }

    p <- p + ggplot2::geom_point(
      data = data_dt,
      ggplot2::aes(x = .data[[age_col]], y = .data[[length_col]], color = sex),
      alpha = data_alpha,
      size = data_size
    )
  }

  # Apply colors
  if (is_hierarchical) {
    p <- p + ggplot2::scale_color_manual(
      values = colors,
      name = "Sex",
      labels = c(sex_labels["1"], sex_labels["2"])
    )
    p <- p + ggplot2::scale_fill_manual(
      values = colors,
      name = "Sex",
      labels = c(sex_labels["1"], sex_labels["2"])
    )
  } else {
    p <- p + ggplot2::scale_color_manual(values = colors[1], guide = "none")
    p <- p + ggplot2::scale_fill_manual(values = colors[1], guide = "none")
  }

  # Faceting
  if (is_hierarchical && facet_sex) {
    p <- p + ggplot2::facet_wrap(~sex)
  }

  # Labels and vitalBayes theme
  p <- p + ggplot2::labs(x = x_lab, y = y_lab, title = title)
  p <- p + theme_vital(base_size = base_size)

  # Apply custom theme arguments
  if (length(theme_args) > 0) {
    p <- p + do.call(ggplot2::theme, theme_args)
  }

  return(p)
}


#' Compare Multiple Growth Models
#'
#' @description
#' Creates a faceted comparison of multiple growth model fits.
#'
#' @param ... Named CmdStanMCMC objects from \code{\link{fit_bayesian_growth}}.
#'   Names become facet labels.
#' @param data Optional data.table/data.frame with observed data.
#' @param age_col,length_col,sex_col Column names in data.
#' @param sex_labels Named character vector. Labels for sex in plot legend and facets.
#'   Default \code{c("1" = "Female", "2" = "Male")}. Can be customized for other
#'   languages, e.g., \code{c("1" = "Femelle", "2" = "Mâle")} (French).
#' @param age_range Numeric vector of length 2. If \code{NULL}, auto-determined.
#' @param k_based_list Logical vector indicating which models are k-based.
#'   If \code{NULL}, assumes all are maturity-based.
#' @param which_model_list Integer vector of model types. If \code{NULL}, assumes VBGM.
#' @param colors Color palette.
#' @param colorblind Logical. Use colorblind-safe palette? Default \code{FALSE}.
#' @param base_size Base font size.
#'
#' @return A ggplot object with faceted model comparison.
#'
#' @examples
#' \dontrun{
#' compare_growth_models(
#'   "von Bertalanffy" = vb_fit,
#'   "Gompertz" = gomp_fit,
#'   "Logistic" = log_fit,
#'   data = shark_data
#' )
#'
#' # With French labels
#' compare_growth_models(
#'   "von Bertalanffy" = vb_fit,
#'   "Gompertz" = gomp_fit,
#'   data = shark_data,
#'   sex_labels = c("1" = "Femelle", "2" = "Mâle")
#' )
#' }
#'
#' @seealso
#' \code{vignette("visualization")} for comprehensive plotting examples.
#'
#' \code{vignette("model_diagnostics")} for model comparison via LOO-CV.
#'
#' \href{../doc/vitalBayes_stats_explained.html#assessment}{Statistical Methods: Model Assessment}
#' for LOO-CV theory.
#'
#' \code{\link{plot_growth_curve}}, \code{\link{compute_loo}}, \code{\link{compare_loo}}
#'
#' @import ggplot2
#' @import data.table
#' @export
compare_growth_models <- function(
    ...,
    data            = NULL,
    age_col         = "age",
    length_col      = "length",
    sex_col         = "sex",
    sex_labels      = c("1" = "Female", "2" = "Male"),
    age_range       = NULL,
    k_based_list    = NULL,
    which_model_list = NULL,
    colors          = NULL,
    colorblind      = FALSE,
    base_size       = 12
) {

  fits <- list(...)
  n_models <- length(fits)
  model_names <- names(fits)

  if (is.null(model_names) || any(model_names == "")) {
    model_names <- paste0("Model ", seq_len(n_models))
  }

  if (is.null(k_based_list)) {
    k_based_list <- rep(FALSE, n_models)
  }

  if (is.null(which_model_list)) {
    which_model_list <- rep(1L, n_models)
  }

  # Determine age range
  if (is.null(age_range)) {
    if (!is.null(data) && age_col %in% names(data)) {
      age_range <- range(data[[age_col]], na.rm = TRUE)
      age_range[2] <- age_range[2] * 1.05
    } else {
      age_range <- c(0, 30)
    }
  }

  age_seq <- seq(age_range[1], age_range[2], length.out = 100)

  # Generate predictions for each model
  all_preds <- list()

  for (i in seq_len(n_models)) {
    fit <- fits[[i]]
    is_hier <- ncol(fit$draws("Linf", format = "matrix")) > 1

    if (is_hier) {
      pred_f <- .posterior_predict_growth(
        fit, age_seq, 1L, k_based_list[i], which_model_list[i]
      )
      pred_f[, `:=`(sex = sex_labels["1"], model = model_names[i])]

      pred_m <- .posterior_predict_growth(
        fit, age_seq, 2L, k_based_list[i], which_model_list[i]
      )
      pred_m[, `:=`(sex = sex_labels["2"], model = model_names[i])]

      all_preds[[i]] <- data.table::rbindlist(list(pred_f, pred_m))
    } else {
      pred <- .posterior_predict_growth(
        fit, age_seq, 1L, k_based_list[i], which_model_list[i]
      )
      pred[, `:=`(sex = "Combined", model = model_names[i])]
      all_preds[[i]] <- pred
    }
  }

  pred_dt <- data.table::rbindlist(all_preds)
  pred_dt[, model := factor(model, levels = model_names)]

  # Set factor levels for sex to maintain order
  if (any(pred_dt$sex != "Combined")) {
    pred_dt[, sex := factor(sex, levels = c(sex_labels["1"], sex_labels["2"], "Combined"))]
  }

  # Set colors - handle colorblind option
  if (colorblind) {
    colors <- vital_colors(2, "sex_cb")
  } else if (is.null(colors)) {
    colors <- vital_colors(2, "sex")
  }

  # Build plot
  p <- ggplot2::ggplot()

  p <- p + ggplot2::geom_ribbon(
    data = pred_dt,
    ggplot2::aes(x = age, ymin = lower, ymax = upper, fill = sex),
    alpha = 0.2
  )

  p <- p + ggplot2::geom_line(
    data = pred_dt,
    ggplot2::aes(x = age, y = median, color = sex),
    linewidth = 1
  )

  # Add data if provided
  if (!is.null(data)) {
    data_dt <- data.table::as.data.table(data)

    if (sex_col %in% names(data_dt)) {
      # Map sex values to labels (handles both numeric and character)
      sex_vals <- data_dt[[sex_col]]
      if (is.numeric(sex_vals)) {
        data_dt[, sex := factor(
          ifelse(get(sex_col) == 1, sex_labels["1"], sex_labels["2"]),
          levels = c(sex_labels["1"], sex_labels["2"])
        )]
      } else {
        data_dt[, sex := factor(
          ifelse(grepl("^[Ff]", get(sex_col)), sex_labels["1"], sex_labels["2"]),
          levels = c(sex_labels["1"], sex_labels["2"])
        )]
      }
    } else {
      data_dt[, sex := "Combined"]
    }

    p <- p + ggplot2::geom_point(
      data = data_dt,
      ggplot2::aes(x = .data[[age_col]], y = .data[[length_col]], color = sex),
      alpha = 0.3,
      size = 1
    )
  }

  p <- p + ggplot2::scale_color_manual(
    values = colors,
    name = "Sex",
    labels = c(sex_labels["1"], sex_labels["2"])
  )
  p <- p + ggplot2::scale_fill_manual(
    values = colors,
    name = "Sex",
    labels = c(sex_labels["1"], sex_labels["2"])
  )
  p <- p + ggplot2::facet_wrap(~model, ncol = 2)
  p <- p + ggplot2::labs(x = "Age (years)", y = "Length (cm)")
  p <- p + theme_vital(base_size = base_size)
  p <- p + ggplot2::theme(legend.position = "bottom")

  return(p)
}


# -----------------------------------------------------------------------------
# Maturity Ogive Visualization
# -----------------------------------------------------------------------------

#' Generate Posterior Predictions for Maturity Ogive
#'
#' @param fit CmdStanMCMC object from fit_bayesian_maturity().
#' @param x_seq Numeric vector of x values (length or age).
#' @param sex_code Integer. 1=female, 2=male.
#' @param type Character. "length" or "age".
#' @param n_draws Number of posterior draws to use.
#'
#' @return data.table with columns: x, median, lower, upper.
#'
#' @noRd
.posterior_predict_maturity <- function(
    fit,
    x_seq,
    sex_code = 1L,
    type = "length",
    n_draws = 500
) {

  param_name <- if (type == "length") "L50" else "t50"

  x50_draws <- fit$draws(param_name, format = "matrix")
  slope_draws <- fit$draws("slope", format = "matrix")

  is_hierarchical <- ncol(x50_draws) > 1

  total_draws <- nrow(x50_draws)
  draw_idx <- sample.int(total_draws, min(n_draws, total_draws))

  s <- if (is_hierarchical) sex_code else 1L

  x50 <- x50_draws[draw_idx, s]
  slope <- slope_draws[draw_idx, s]

  # Compute probabilities
  n_x <- length(x_seq)
  prob_mat <- matrix(NA_real_, nrow = length(draw_idx), ncol = n_x)

  for (i in seq_along(x_seq)) {
    eta <- slope * (x_seq[i] - x50)
    prob_mat[, i] <- stats::plogis(eta)  # logit link
  }

  data.table::data.table(
    x = x_seq,
    median = apply(prob_mat, 2, stats::median),
    lower = apply(prob_mat, 2, stats::quantile, 0.025),
    upper = apply(prob_mat, 2, stats::quantile, 0.975)
  )
}


#' Plot Maturity Ogive
#'
#' @description
#' Creates publication-quality maturity ogive plots with credible intervals.
#' Supports both length-at-maturity and age-at-maturity models.
#'
#' @param fit A CmdStanMCMC object from \code{\link{fit_bayesian_maturity}}.
#' @param type Character. Type of maturity model: \code{"length"} or \code{"age"}.
#' @param data Optional data.table/data.frame with observed data.
#' @param x_col Character. Column name for predictor (length or age).
#' @param maturity_col Character. Column name for maturity status.
#' @param sex_col Character. Column name for sex.
#' @param x_range Numeric vector of length 2. Range for predictions.
#' @param n_points Integer. Number of points for ogive. Default 100.
#' @param show_data Logical. Show observed data? Default \code{TRUE}.
#' @param show_rug Logical. Show rug plot at bottom? Default \code{TRUE}.
#' @param show_x50_line Logical. Show vertical line at x50? Default \code{TRUE}.
#' @param data_alpha Numeric. Transparency for data. Default 0.5.
#' @param ribbon_alpha Numeric. Transparency for ribbon. Default 0.25.
#' @param line_width Numeric. Width of ogive line. Default 1.2.
#' @param colors Character vector of colors.
#' @param colorblind Logical. Use colorblind-safe palette? Default \code{FALSE}.
#' @param facet_sex Logical. Facet by sex for two-sex models? Default \code{TRUE}.
#' @param x_lab,y_lab Character. Axis labels. Auto-generated if \code{NULL}.
#' @param title Character. Plot title.
#' @param sex_labels Named character vector. Labels for sex in plot legend and facets.
#'   Default \code{c("1" = "Female", "2" = "Male")}. Can be customized for other
#'   languages, e.g., \code{c("1" = "Femelle", "2" = "Mâle")} (French),
#'   \code{c("1" = "Hembra", "2" = "Macho")} (Spanish).
#' @param base_size Numeric. Base font size. Default 12.
#' @param theme_args,additional_layers Lists for customization.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # Length-at-maturity
#' plot_maturity_ogive(L50_fit, type = "length", data = shark_data)
#'
#' # French labels
#' plot_maturity_ogive(
#'   L50_fit,
#'   type = "length",
#'   data = shark_data,
#'   sex_labels = c("1" = "Femelle", "2" = "Mâle")
#' )
#'
#' # Age-at-maturity with custom styling
#' plot_maturity_ogive(
#'   t50_fit,
#'   type = "age",
#'   data = shark_data,
#'   colors = c("#E63946", "#1D3557"),
#'   show_rug = TRUE
#' )
#' }
#'
#' @seealso
#' \code{vignette("visualization")} for comprehensive plotting examples.
#'
#' \code{vignette("fit_bayesian_maturity")} for maturity model fitting.
#'
#' \code{\link{fit_bayesian_maturity}}, \code{\link{vital_colors}}, \code{\link{theme_vital}}
#'
#' @import ggplot2
#' @import data.table
#' @export
plot_maturity_ogive <- function(
    fit,
    type           = c("length", "age"),
    data           = NULL,
    x_col          = NULL,
    maturity_col   = "maturity",
    sex_col        = "sex",
    sex_labels     = c("1" = "Female", "2" = "Male"),
    x_range        = NULL,
    n_points       = 100,
    show_data      = TRUE,
    show_rug       = TRUE,
    show_x50_line  = TRUE,
    data_alpha     = 0.5,
    ribbon_alpha   = 0.25,
    line_width     = 1.2,
    colors         = NULL,
    colorblind     = FALSE,
    facet_sex      = TRUE,
    x_lab          = NULL,
    y_lab          = "Probability of Maturity",
    title          = NULL,
    base_size      = 12,
    theme_args     = list(),
    additional_layers = list()
) {

  type <- match.arg(type)

  # Set defaults based on type
  if (is.null(x_col)) {
    x_col <- if (type == "length") "length" else "age"
  }

  if (is.null(x_lab)) {
    x_lab <- if (type == "length") "Length (cm)" else "Age (years)"
  }

  # Check if hierarchical
  param_name <- if (type == "length") "L50" else "t50"
  x50_test <- fit$draws(param_name, format = "matrix")
  is_hierarchical <- ncol(x50_test) > 1

  # Determine x range
  if (is.null(x_range)) {
    if (!is.null(data) && x_col %in% names(data)) {
      x_range <- range(data[[x_col]], na.rm = TRUE)
      buffer <- diff(x_range) * 0.1
      x_range <- x_range + c(-buffer, buffer)
    } else {
      x50_median <- stats::median(x50_test)
      x_range <- x50_median * c(0.5, 1.5)
    }
  }

  x_seq <- seq(x_range[1], x_range[2], length.out = n_points)

  # Generate predictions
  if (is_hierarchical) {
    pred_f <- .posterior_predict_maturity(fit, x_seq, 1L, type)
    pred_f[, sex := sex_labels["1"]]

    pred_m <- .posterior_predict_maturity(fit, x_seq, 2L, type)
    pred_m[, sex := sex_labels["2"]]

    pred_dt <- data.table::rbindlist(list(pred_f, pred_m))

    # Set factor levels to maintain order
    pred_dt[, sex := factor(sex, levels = c(sex_labels["1"], sex_labels["2"]))]
  } else {
    pred_dt <- .posterior_predict_maturity(fit, x_seq, 1L, type)
    pred_dt[, sex := "Combined"]
  }

  # Set colors - handle colorblind option
  if (colorblind) {
    colors <- vital_colors(2, "sex_cb")
  } else if (is.null(colors)) {
    colors <- if (is_hierarchical) vital_colors(2, "sex") else vital_colors(1)
  }

  # Build plot
  p <- ggplot2::ggplot()

  # Additional layers first
  for (layer in additional_layers) {
    p <- p + layer
  }

  # Ribbon
  p <- p + ggplot2::geom_ribbon(
    data = pred_dt,
    ggplot2::aes(x = x, ymin = lower, ymax = upper, fill = sex),
    alpha = ribbon_alpha
  )

  # Ogive line
  p <- p + ggplot2::geom_line(
    data = pred_dt,
    ggplot2::aes(x = x, y = median, color = sex),
    linewidth = line_width
  )

  # x50 reference lines
  if (show_x50_line) {
    x50_summary <- fit$summary(param_name)

    if (is_hierarchical) {
      x50_dt <- data.table::data.table(
        sex = factor(c(sex_labels["1"], sex_labels["2"]),
                     levels = c(sex_labels["1"], sex_labels["2"])),
        x50 = x50_summary$median
      )
    } else {
      x50_dt <- data.table::data.table(
        sex = "Combined",
        x50 = x50_summary$median
      )
    }

    p <- p + ggplot2::geom_vline(
      data = x50_dt,
      ggplot2::aes(xintercept = x50, color = sex),
      linetype = "dashed",
      linewidth = 0.7,
      alpha = 0.7
    )

    p <- p + ggplot2::geom_hline(
      yintercept = 0.5,
      linetype = "dotted",
      color = "gray50"
    )
  }

  # Add observed data
  if (show_data && !is.null(data)) {
    data_dt <- data.table::as.data.table(data)

    if (is_hierarchical && sex_col %in% names(data_dt)) {
      # Map sex values to labels (handles both numeric and character)
      sex_vals <- data_dt[[sex_col]]
      if (is.numeric(sex_vals)) {
        data_dt[, sex := factor(
          ifelse(get(sex_col) == 1, sex_labels["1"], sex_labels["2"]),
          levels = c(sex_labels["1"], sex_labels["2"])
        )]
      } else {
        data_dt[, sex := factor(
          ifelse(grepl("^[Ff]", get(sex_col)), sex_labels["1"], sex_labels["2"]),
          levels = c(sex_labels["1"], sex_labels["2"])
        )]
      }
    } else {
      data_dt[, sex := "Combined"]
    }

    # Jitter maturity for visualization
    data_dt[, maturity_jitter := get(maturity_col) + stats::runif(.N, -0.03, 0.03)]

    p <- p + ggplot2::geom_point(
      data = data_dt,
      ggplot2::aes(x = .data[[x_col]], y = maturity_jitter, color = sex),
      alpha = data_alpha,
      size = 1.5
    )

    if (show_rug) {
      p <- p + ggplot2::geom_rug(
        data = data_dt[get(maturity_col) == 0],
        ggplot2::aes(x = .data[[x_col]], color = sex),
        sides = "b",
        alpha = 0.3
      )
      p <- p + ggplot2::geom_rug(
        data = data_dt[get(maturity_col) == 1],
        ggplot2::aes(x = .data[[x_col]], color = sex),
        sides = "t",
        alpha = 0.3
      )
    }
  }

  # Colors and scales
  p <- p + ggplot2::scale_color_manual(
    values = colors,
    name = "Sex",
    labels = c(sex_labels["1"], sex_labels["2"])
  )
  p <- p + ggplot2::scale_fill_manual(
    values = colors,
    name = "Sex",
    labels = c(sex_labels["1"], sex_labels["2"])
  )
  p <- p + ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25))

  # Faceting
  if (is_hierarchical && facet_sex) {
    p <- p + ggplot2::facet_wrap(~sex)
  }

  # Labels and vitalBayes theme
  p <- p + ggplot2::labs(x = x_lab, y = y_lab, title = title)
  p <- p + theme_vital(base_size = base_size)

  if (length(theme_args) > 0) {
    p <- p + do.call(ggplot2::theme, theme_args)
  }

  return(p)
}


# -----------------------------------------------------------------------------
# Birth Model Visualization
# -----------------------------------------------------------------------------

#' Plot Birth Ogive
#'
#' @description
#' Creates publication-quality birth transition plots.
#'
#' @param fit A CmdStanMCMC object from \code{\link{fit_bayesian_birth}}.
#' @param data Optional data.table/data.frame with observed data.
#' @param length_col Character. Column name for length.
#' @param status_col Character. Column name for birth status.
#' @param length_range Numeric vector of length 2.
#' @param n_points Integer. Number of points. Default 100.
#' @param show_data Logical. Show data? Default \code{TRUE}.
#' @param show_rug Logical. Show rug? Default \code{TRUE}.
#' @param show_b50_line Logical. Show b50 line? Default \code{TRUE}.
#' @param color Character. Line color. Default uses primary vitalBayes magenta.
#' @param ribbon_alpha Numeric. Ribbon transparency.
#' @param x_lab,y_lab Character. Axis labels.
#' @param title Character. Plot title.
#' @param base_size Numeric. Base font size.
#' @param theme_args,additional_layers Lists for customization.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @import data.table
#' @export
plot_birth_ogive <- function(
    fit,
    data          = NULL,
    length_col    = "length",
    status_col    = "status",
    length_range  = NULL,
    n_points      = 100,
    show_data     = TRUE,
    show_rug      = TRUE,
    show_b50_line = TRUE,
    color         = NULL,
    ribbon_alpha  = 0.25,
    x_lab         = "Length (cm)",
    y_lab         = "Probability of Free-Swimming",
    title         = NULL,
    base_size     = 12,
    theme_args    = list(),
    additional_layers = list()
) {

  # Default to primary vitalBayes magenta
  if (is.null(color)) {
    color <- vital_colors(1)
  }

  # Extract draws
  b50_draws <- as.vector(fit$draws("b50", format = "matrix"))
  slope_draws <- as.vector(fit$draws("slope", format = "matrix"))

  # Determine length range
  if (is.null(length_range)) {
    b50_median <- stats::median(b50_draws)
    length_range <- b50_median * c(0.7, 1.3)
  }

  length_seq <- seq(length_range[1], length_range[2], length.out = n_points)

  # Subsample draws
  n_draws <- min(500, length(b50_draws))
  idx <- sample.int(length(b50_draws), n_draws)

  # Compute probabilities (probit link)
  prob_mat <- matrix(NA_real_, nrow = n_draws, ncol = n_points)

  for (i in seq_along(length_seq)) {
    eta <- slope_draws[idx] * (length_seq[i] - b50_draws[idx])
    prob_mat[, i] <- stats::pnorm(eta)  # probit link
  }

  pred_dt <- data.table::data.table(
    length = length_seq,
    median = apply(prob_mat, 2, stats::median),
    lower = apply(prob_mat, 2, stats::quantile, 0.025),
    upper = apply(prob_mat, 2, stats::quantile, 0.975)
  )

  # Build plot
  p <- ggplot2::ggplot()

  for (layer in additional_layers) {
    p <- p + layer
  }

  p <- p + ggplot2::geom_ribbon(
    data = pred_dt,
    ggplot2::aes(x = length, ymin = lower, ymax = upper),
    fill = color,
    alpha = ribbon_alpha
  )

  p <- p + ggplot2::geom_line(
    data = pred_dt,
    ggplot2::aes(x = length, y = median),
    color = color,
    linewidth = 1.2
  )

  if (show_b50_line) {
    b50_median <- stats::median(b50_draws)

    p <- p + ggplot2::geom_vline(
      xintercept = b50_median,
      linetype = "dashed",
      color = color,
      alpha = 0.7
    )

    p <- p + ggplot2::geom_hline(
      yintercept = 0.5,
      linetype = "dotted",
      color = "gray50"
    )
  }

  if (show_data && !is.null(data)) {
    data_dt <- data.table::as.data.table(data)
    data_dt[, status_jitter := get(status_col) + stats::runif(.N, -0.03, 0.03)]

    p <- p + ggplot2::geom_point(
      data = data_dt,
      ggplot2::aes(x = .data[[length_col]], y = status_jitter),
      color = color,
      alpha = 0.5,
      size = 1.5
    )

    if (show_rug) {
      p <- p + ggplot2::geom_rug(
        data = data_dt[get(status_col) == 0],
        ggplot2::aes(x = .data[[length_col]]),
        sides = "b",
        color = color,
        alpha = 0.3
      )
      p <- p + ggplot2::geom_rug(
        data = data_dt[get(status_col) == 1],
        ggplot2::aes(x = .data[[length_col]]),
        sides = "t",
        color = color,
        alpha = 0.3
      )
    }
  }

  p <- p + ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25))
  p <- p + ggplot2::labs(x = x_lab, y = y_lab, title = title)
  p <- p + theme_vital(base_size = base_size)

  if (length(theme_args) > 0) {
    p <- p + do.call(ggplot2::theme, theme_args)
  }

  return(p)
}

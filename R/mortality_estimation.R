# =============================================================================
# vitalBayes Mortality Estimation Functions
# =============================================================================
#
# This module implements natural mortality estimation using biological
# milestones (Linf, L0, Lmat, tmat) with model-dependent G(t).
#
# Core mortality formula:
#   M(t) = M_inf / G(t)
#
# where:
#   M_inf = (1/tmat) * ln[(Linf - L0)/(Linf - Lmat)]  (VB-derived, unified anchor)
#   G(t) = L(t)/Linf  (computed using native growth model trajectory)
#
# Supported mortality models:
#   - CW: Chen-Watanabe (1989) with optional two-phase senescence
#   - PW: Peterson-Wroblewski (1984) weight-based
#   - L:  Lorenzen (1996/2022) weight- or growth-based
#
# =============================================================================


# =============================================================================
# SECTION 1: CORE HELPER FUNCTIONS
# =============================================================================


#' Compute Native Growth Coefficient for a Specific Model
#'
#' @description
#' Derives the growth coefficient \eqn{k} for a specified growth model from
#' biological milestones \eqn{(L_\infty, L_0, L_{mat}, t_{mat})}.
#'
#' @details
#' The native k formulas for each model are:
#'
#' \strong{Von Bertalanffy:}
#' \deqn{k = \frac{1}{t_{mat}} \ln\left(\frac{L_\infty - L_0}{L_\infty - L_{mat}}\right)}
#'
#' \strong{Gompertz:}
#' \deqn{k = -\frac{1}{t_{mat}} \ln\left(\frac{\ln(L_\infty / L_{mat})}{\ln(L_\infty / L_0)}\right)}
#'
#' \strong{Logistic:}
#' \deqn{k = -\frac{1}{t_{mat}} \ln\left(\frac{L_\infty/L_{mat} - 1}{L_\infty/L_0 - 1}\right)}
#'
#' @param Linf Numeric vector. Asymptotic length.
#' @param L0 Numeric vector. Length at birth.
#' @param Lmat Numeric vector. Length at maturity.
#' @param tmat Numeric vector. Age at maturity.
#' @param growth_model Character. Growth model: \code{"vb"}, \code{"gompertz"},
#'   or \code{"logistic"}.
#' @param warn Logical. Warn on invalid values?
#'
#' @return Numeric vector of native k values. Invalid values returned as NA.
#'
#' @examples
#' k_vb <- compute_k_native(Linf = 126, L0 = 35, Lmat = 83, tmat = 47,
#'                          growth_model = "vb")
#' k_gomp <- compute_k_native(Linf = 108, L0 = 35, Lmat = 83, tmat = 47,
#'                            growth_model = "gompertz")
#'
#' @export
compute_k_native <- function(Linf, L0, Lmat, tmat,
                             growth_model = c("vb", "gompertz", "logistic"),
                             warn = TRUE) {

  growth_model <- match.arg(growth_model)

  n <- max(length(Linf), length(L0), length(Lmat), length(tmat))
  Linf <- rep_len(Linf, n)
  L0   <- rep_len(L0, n)
  Lmat <- rep_len(Lmat, n)
  tmat <- rep_len(tmat, n)

  k <- switch(
    growth_model,
    "vb" = {
      numerator   <- Linf - L0
      denominator <- Linf - Lmat
      invalid <- numerator <= 0 | denominator <= 0 | tmat <= 0
      k_out <- (1 / tmat) * log(numerator / denominator)
      k_out[invalid] <- NA_real_
      k_out
    },
    "gompertz" = {
      ratio_mat <- Linf / Lmat
      ratio_0   <- Linf / L0
      invalid <- ratio_mat <= 1 | ratio_0 <= 1 | tmat <= 0
      k_out <- -(1 / tmat) * log(log(ratio_mat) / log(ratio_0))
      k_out[invalid] <- NA_real_
      k_out
    },
    "logistic" = {
      term_mat <- Linf / Lmat - 1
      term_0   <- Linf / L0 - 1
      invalid <- term_mat <= 0 | term_0 <= 0 | tmat <= 0
      k_out <- -(1 / tmat) * log(term_mat / term_0)
      k_out[invalid] <- NA_real_
      k_out
    }
  )

  if (warn && any(is.na(k))) {
    n_invalid <- sum(is.na(k))
    warning(sprintf("%d of %d values produced invalid k.", n_invalid, n), call. = FALSE)
  }

  k
}


#' Compute Asymptotic Mortality Rate (M_inf)
#'
#' @description
#' Computes the asymptotic mortality rate \eqn{M_\infty} from biological
#' milestones using the VB formula. This serves as the unified mortality
#' scaling factor in the Chen-Watanabe framework.
#'
#' @details
#' \deqn{M_\infty = \frac{1}{t_{mat}} \ln\left(\frac{L_\infty - L_0}{L_\infty - L_{mat}}\right)}
#'
#' This formula is used regardless of growth model, providing a consistent
#' mortality anchor that prevents the ~10-fold survival differences arising
#' from model-native k formulas.
#'
#' @param Linf Numeric vector. Asymptotic length.
#' @param L0 Numeric vector. Length at birth.
#' @param Lmat Numeric vector. Length at maturity.
#' @param tmat Numeric vector. Age at maturity.
#' @param warn Logical. Warn on invalid values?
#'
#' @return Numeric vector of asymptotic mortality rates.
#'
#' @examples
#' Minf <- compute_Minf(Linf = 126, L0 = 35, Lmat = 83, tmat = 47)
#'
#' @export
compute_Minf <- function(Linf, L0, Lmat, tmat, warn = TRUE) {

  n <- max(length(Linf), length(L0), length(Lmat), length(tmat))
  Linf <- rep_len(Linf, n)
  L0   <- rep_len(L0, n)
  Lmat <- rep_len(Lmat, n)
  tmat <- rep_len(tmat, n)

  numerator   <- Linf - L0
  denominator <- Linf - Lmat

  invalid <- numerator <= 0 | denominator <= 0 | tmat <= 0 |
    is.na(numerator) | is.na(denominator) | is.na(tmat)

  if (warn && any(invalid, na.rm = TRUE)) {
    n_invalid <- sum(invalid, na.rm = TRUE)
    warning(sprintf("%d of %d draws produced invalid M_inf.", n_invalid, n), call. = FALSE)
  }

  Minf <- (1 / tmat) * log(numerator / denominator)
  Minf[invalid] <- NA_real_

  Minf
}


#' Compute Relative Size G(t) = L(t)/Linf
#'
#' @description
#' Computes the relative size at specified ages using a given growth model.
#' This is the denominator in the mortality formula M(t) = M_inf / G(t).
#'
#' @param age Numeric vector. Ages at which to compute G(t).
#' @param Linf Numeric. Asymptotic length.
#' @param L0 Numeric. Length at birth.
#' @param k Numeric. Native growth coefficient.
#' @param growth_model Character. Growth model.
#'
#' @return Numeric vector of G(t) values bounded in (0, 1].
#'
#' @examples
#' ages <- seq(0, 100, by = 1)
#' G_vb <- compute_G(ages, Linf = 126, L0 = 35, k = 0.016, growth_model = "vb")
#'
#' @export
compute_G <- function(age, Linf, L0, k,
                      growth_model = c("vb", "gompertz", "logistic")) {

  growth_model <- match.arg(growth_model)

  G <- switch(
    growth_model,
    "vb" = 1 - ((Linf - L0) / Linf) * exp(-k * age),
    "gompertz" = {
      r0 <- log(Linf / L0)
      exp(-r0 * exp(-k * age))
    },
    "logistic" = {
      c <- Linf / L0 - 1
      1 / (1 + c * exp(-k * age))
    }
  )

  pmax(G, 1e-6)
}


#' Compute Length at Age L(t)
#'
#' @description
#' Computes predicted length at age using a specified growth model.
#'
#' @param age Numeric vector. Ages.
#' @param Linf Numeric. Asymptotic length.
#' @param L0 Numeric. Length at birth.
#' @param k Numeric. Native growth coefficient.
#' @param growth_model Character. Growth model.
#'
#' @return Numeric vector of predicted lengths.
#'
#' @examples
#' L_vb <- compute_L(0:50, Linf = 126, L0 = 35, k = 0.016, growth_model = "vb")
#'
#' @export
compute_L <- function(age, Linf, L0, k,
                      growth_model = c("vb", "gompertz", "logistic")) {

  growth_model <- match.arg(growth_model)

  switch(
    growth_model,
    "vb" = Linf - (Linf - L0) * exp(-k * age),
    "gompertz" = {
      r0 <- log(Linf / L0)
      Linf * exp(-r0 * exp(-k * age))
    },
    "logistic" = {
      c <- Linf / L0 - 1
      Linf / (1 + c * exp(-k * age))
    }
  )
}


#' Compute Maximum Age from Growth Parameters
#'
#' @description
#' Estimates tmax as age when L(t) reaches a specified fraction of Linf.
#'
#' @param Linf Numeric. Asymptotic length.
#' @param L0 Numeric. Length at birth.
#' @param k Numeric. Native growth coefficient.
#' @param growth_model Character. Growth model.
#' @param Linf_factor Numeric in (0,1). Fraction of Linf (default 0.99).
#'
#' @return Numeric tmax value.
#'
#' @examples
#' tmax <- compute_tmax(Linf = 126, L0 = 35, k = 0.016,
#'                      growth_model = "vb", Linf_factor = 0.99)
#'
#' @export
compute_tmax <- function(Linf, L0, k,
                         growth_model = c("vb", "gompertz", "logistic"),
                         Linf_factor = 0.99) {

  growth_model <- match.arg(growth_model)

  if (Linf_factor <= 0 || Linf_factor >= 1) {
    stop("Linf_factor must be in (0, 1).", call. = FALSE)
  }

  switch(
    growth_model,
    "vb" = -log(Linf * (1 - Linf_factor) / (Linf - L0)) / k,
    "gompertz" = {
      r0 <- log(Linf / L0)
      -log(-log(Linf_factor) / r0) / k
    },
    "logistic" = {
      c <- Linf / L0 - 1
      -log((1 / Linf_factor - 1) / c) / k
    }
  )
}


# =============================================================================
# SECTION 2: MORTALITY MODELS
# =============================================================================


#' Chen-Watanabe Natural Mortality (Model-Dependent)
#'
#' @description
#' Computes age-specific natural mortality using a generalized Chen-Watanabe
#' framework where G(t) is derived from the native growth model trajectory.
#'
#' @details
#' The mortality model is:
#' \deqn{M(t) = \frac{M_\infty}{G(t)}}
#'
#' where M_inf is VB-derived and G(t) uses the native growth model trajectory.
#' This formulation provides a unified mortality anchor while capturing
#' model-specific growth dynamics.
#'
#' @param age Numeric vector of ages.
#' @param Linf Asymptotic length.
#' @param L0 Length at birth.
#' @param Lmat Length at maturity.
#' @param tmat Age at maturity.
#' @param growth_model Character. Growth model for G(t): \code{"vb"},
#'   \code{"gompertz"}, or \code{"logistic"}.
#' @param tmax Maximum age (computed if NULL).
#' @param Linf_factor Fraction of Linf for tmax estimation.
#' @param two_phase Use two-phase senescence model?
#' @param late_model Senescence model: \code{"gompertz"} or \code{"logistic"}.
#' @param tm_factor Transition age as fraction of tmat.
#' @param M_mult Mortality multiplier for senescence.
#' @param smooth_factor Transition smoothness.
#'
#' @return Numeric vector of instantaneous mortality rates.
#'
#' @references
#' Chen, S., & Watanabe, S. (1989). Age dependence of natural mortality
#' coefficient in fish population dynamics. \emph{Nippon Suisan Gakkaishi},
#' 55(2), 205-208.
#'
#' @examples
#' ages <- seq(0.1, 150, by = 1)
#' M_vb <- M_chen_watanabe(ages, Linf = 126, L0 = 35, Lmat = 83, tmat = 47,
#'                         growth_model = "vb")
#'
#' @export
M_chen_watanabe <- function(
    age,
    Linf,
    L0,
    Lmat,
    tmat,
    growth_model = c("vb", "gompertz", "logistic"),
    tmax = NULL,
    Linf_factor = 0.99,
    two_phase = FALSE,
    late_model = c("gompertz", "logistic"),
    tm_factor = 2/3,
    M_mult = 2,
    smooth_factor = 1/3
) {

  growth_model <- match.arg(growth_model)
  late_model   <- match.arg(late_model)

  if (L0 >= Linf) stop("L0 must be less than Linf.", call. = FALSE)
  if (Lmat >= Linf) stop("Lmat must be less than Linf.", call. = FALSE)
  if (L0 >= Lmat) stop("L0 must be less than Lmat.", call. = FALSE)
  if (tmat <= 0) stop("tmat must be positive.", call. = FALSE)

  Minf <- compute_Minf(Linf, L0, Lmat, tmat, warn = FALSE)
  k_native <- compute_k_native(Linf, L0, Lmat, tmat, growth_model, warn = FALSE)

  if (is.na(Minf) || Minf <= 0) stop("Invalid M_inf.", call. = FALSE)
  if (is.na(k_native) || k_native <= 0) stop("Invalid native k.", call. = FALSE)

  if (is.null(tmax)) {
    tmax <- compute_tmax(Linf, L0, k_native, growth_model, Linf_factor)
  }

  G_t <- compute_G(age, Linf, L0, k_native, growth_model)
  M_cw <- Minf / G_t

  if (!two_phase) return(M_cw)

  # Two-phase extension
  tm <- tm_factor * tmat
  G_at_tm <- compute_G(tm, Linf, L0, k_native, growth_model)
  M_early_at_tm <- Minf / G_at_tm

  if (late_model == "gompertz") {
    M_s <- M_mult * M_early_at_tm
    r <- log(M_early_at_tm / M_s) / (tm - tmax)
    M_late <- function(t) M_s * exp(r * (t - tmax))
  } else {
    K <- 2 * M_mult * M_early_at_tm
    r <- -log(K / M_early_at_tm - 1) / (tm - tmax)
    M_late <- function(t) K / (1 + exp(-r * (t - tmax)))
  }

  smooth_width <- smooth_factor * abs(tmax - tm)
  weight <- 1 / (1 + exp(-2 * (age - tm) / smooth_width))

  M <- (1 - weight) * M_cw + weight * M_late(age)
  M[age < tm * 0.5] <- M_cw[age < tm * 0.5]
  pmax(M, 0)
}


#' Peterson-Wroblewski Natural Mortality Model
#'
#' @description
#' Computes weight-based natural mortality following Peterson & Wroblewski (1984).
#'
#' @details
#' \deqn{M(W) = 1.92 \cdot W^{-0.25}}
#'
#' \strong{Warning}: This model was calibrated on teleost fishes and may
#' produce biologically implausible mortality rates for elasmobranchs.
#'
#' @param age Numeric vector of ages.
#' @param Linf Asymptotic length.
#' @param L0 Length at birth.
#' @param Lmat Length at maturity.
#' @param tmat Age at maturity.
#' @param lw_fun Length-weight function: \code{lw_fun(L)} returns weight in grams.
#' @param growth_model Character. Growth model for L(t).
#'
#' @return Numeric vector of instantaneous mortality rates.
#'
#' @references
#' Peterson, I., & Wroblewski, J. S. (1984). Mortality rate of fishes in the
#' pelagic ecosystem. \emph{Canadian Journal of Fisheries and Aquatic Sciences},
#' 41(7), 1117-1120.
#'
#' @export
M_peterson_wroblewski <- function(
    age,
    Linf,
    L0,
    Lmat,
    tmat,
    lw_fun,
    growth_model = c("vb", "gompertz", "logistic")
) {

  growth_model <- match.arg(growth_model)

  if (is.null(lw_fun) || !is.function(lw_fun)) {
    stop("PW model requires a length-weight function 'lw_fun(L)'.", call. = FALSE)
  }

  k_native <- compute_k_native(Linf, L0, Lmat, tmat, growth_model, warn = FALSE)
  L_t <- compute_L(age, Linf, L0, k_native, growth_model)
  W_t <- lw_fun(L_t)
  W_t <- pmax(W_t, 1)

  1.92 * W_t^(-0.25)
}


#' Lorenzen Natural Mortality Model
#'
#' @description
#' Computes size-dependent natural mortality following Lorenzen (1996, 2022).
#'
#' @details
#' Two formulations are available:
#'
#' \strong{Weight-based} (Lorenzen 1996):
#' \deqn{M(W) = \alpha \cdot W^{\beta}}
#'
#' \strong{Growth-based} (Lorenzen 2022):
#' \deqn{\ln M = 0.28 - 1.30 \ln(L/L_\infty) + 1.08 \ln(k)}
#'
#' For the growth-based formulation, M_inf (VB-derived) is used as k for
#' consistency with the Chen-Watanabe framework.
#'
#' @param age Numeric vector of ages.
#' @param Linf Asymptotic length.
#' @param L0 Length at birth.
#' @param Lmat Length at maturity.
#' @param tmat Age at maturity.
#' @param lw_fun Length-weight function (required if weight_based = TRUE).
#' @param weight_based Use weight-based formulation?
#' @param growth_model Growth model for L(t).
#' @param sample_params Sample parameters from their distributions?
#'
#' @return Numeric vector of instantaneous mortality rates.
#'
#' @references
#' Lorenzen, K. (1996). The relationship between body weight and natural
#' mortality in juvenile and adult fish. \emph{Journal of Fish Biology},
#' 49(4), 627-642.
#'
#' Lorenzen, K. (2022). Size- and age-dependent natural mortality in fish
#' populations. \emph{Fisheries Research}, 255, 106454.
#'
#' @export
M_lorenzen <- function(
    age,
    Linf,
    L0,
    Lmat,
    tmat,
    lw_fun = NULL,
    weight_based = FALSE,
    growth_model = c("vb", "gompertz", "logistic"),
    sample_params = TRUE
) {

  growth_model <- match.arg(growth_model)

  k_native <- compute_k_native(Linf, L0, Lmat, tmat, growth_model, warn = FALSE)
  Minf <- compute_Minf(Linf, L0, Lmat, tmat, warn = FALSE)

  if (weight_based) {
    if (is.null(lw_fun) || !is.function(lw_fun)) {
      stop("Weight-based Lorenzen requires 'lw_fun'.", call. = FALSE)
    }

    L_t <- compute_L(age, Linf, L0, k_native, growth_model)
    W_t <- lw_fun(L_t)
    W_t <- pmax(W_t, 1)

    if (sample_params) {
      alpha <- stats::rnorm(1, 3.69, 0.502)
      beta  <- stats::rnorm(1, -0.305, 0.029)
    } else {
      alpha <- 3.69
      beta  <- -0.305
    }
    M <- alpha * W_t^beta

  } else {
    L_t <- compute_L(age, Linf, L0, k_native, growth_model)
    L_ratio <- L_t / Linf

    if (sample_params) {
      intercept <- stats::rnorm(1, 0.28, 0.105)
      coef_L    <- stats::rnorm(1, -1.30, 0.059)
      coef_k    <- stats::rnorm(1, 1.08, 0.082)
    } else {
      intercept <- 0.28
      coef_L    <- -1.30
      coef_k    <- 1.08
    }

    log_M <- intercept + coef_L * log(L_ratio) + coef_k * log(Minf)
    M <- exp(log_M)
  }

  M
}


#' Scale Mortality Schedule to Target Mean
#'
#' @description
#' Rescales an age-specific mortality schedule so its mean equals a target
#' value derived from empirical relationships or survival probability.
#'
#' @param M Numeric vector of mortality rates.
#' @param M_target Target mean mortality. Can be scalar, function of tmax, or NULL.
#' @param tmax Maximum age (required if M_target is function or NULL).
#' @param p Survival probability (used if M_target NULL).
#'
#' @return Scaled mortality vector.
#'
#' @export
scale_mortality <- function(M, M_target = NULL, tmax = NULL, p = 0.001) {

  if (is.null(M_target)) {
    if (is.null(tmax)) stop("tmax required when M_target is NULL.", call. = FALSE)
    M_target <- -log(p) / tmax
  } else if (is.function(M_target)) {
    if (is.null(tmax)) stop("tmax required when M_target is a function.", call. = FALSE)
    M_target <- M_target(tmax)
  }

  M_mean <- mean(M, na.rm = TRUE)
  if (M_mean <= 0) {
    warning("Mean mortality <= 0; scaling not applied.", call. = FALSE)
    return(M)
  }

  M * (M_target / M_mean)
}


# =============================================================================
# SECTION 3: STOCHASTIC MORTALITY ESTIMATION
# =============================================================================


#' Stochastic Estimation of Age-Specific Natural Mortality
#'
#' @description
#' Monte Carlo simulation of age-specific natural mortality schedules with
#' uncertainty propagation from growth and maturity parameters.
#'
#' @details
#' This function samples life history parameters from specified distributions
#' and computes mortality schedules using the chosen method. The output format
#' is compatible with \code{\link{simulate_survivorship}}.
#'
#' Three mortality models are available:
#' \itemize{
#'   \item \strong{CW}: Chen-Watanabe (1989) with model-dependent G(t)
#'   \item \strong{PW}: Peterson-Wroblewski (1984) weight-based
#'   \item \strong{L}: Lorenzen (1996/2022) weight- or growth-based
#' }
#'
#' @param method Character. Mortality model: \code{"CW"}, \code{"PW"}, or \code{"L"}.
#' @param Linf,L0,Lmat,tmat Numeric vectors of length 2: \code{c(mean, sd)}.
#' @param growth_model Character. Growth model for G(t) and L(t): \code{"vb"},
#'   \code{"gompertz"}, or \code{"logistic"}.
#' @param Linf_factor Fraction of Linf used to estimate tmax. Default 0.99.
#' @param age_seq Function or numeric vector for age grid. Default creates
#'   500 points from 0.1 to tmax.
#' @param iter Number of Monte Carlo iterations. Default 2000.
#' @param scaled Scale mortality to M_target? Default TRUE.
#' @param M_target Target mean mortality. Can be scalar, function of tmax, or NULL.
#' @param p Survival probability for scaling if M_target NULL. Default 0.001.
#' @param two_phase Use CW two-phase senescence? Default FALSE.
#' @param late_model Senescence model: \code{"gompertz"} or \code{"logistic"}.
#' @param tm_factor Transition age factor. Default 2/3.
#' @param M_mult Mortality multiplier for senescence. Default 2.
#' @param smooth_factor Transition smoothness. Default 1/3.
#' @param weight_based For Lorenzen: use weight-based formulation? Default FALSE.
#' @param lw_fun Length-weight function (required for PW and weight-based Lorenzen).
#' @param seed Random seed. Default 1234.
#' @param palette Color palette: \code{"synthwave"}, \code{"viridis"}, \code{"okabe"},
#'   \code{"plasma"}, or \code{"inferno"}.
#' @param print_plot Print plot? Default TRUE.
#' @param show_progress Show progress messages? Default TRUE.
#'
#' @return A list with:
#' \describe{
#'   \item{Schedules}{data.table with columns: set_id, age, M}
#'   \item{Parameters}{data.table with columns: set_id, Linf, L0, Lmat, tmat, Minf, k_native, tmax}
#'   \item{Summary}{data.table with age-wise median and 95% CI}
#'   \item{Plot}{ggplot2 object}
#' }
#'
#' @examples
#' \dontrun{
#' mort <- get_stochastic_mortality(
#'   method = "CW",
#'   Linf = c(108, 10),
#'   L0 = c(35, 2),
#'   Lmat = c(83, 5),
#'   tmat = c(47, 3),
#'   growth_model = "gompertz",
#'   iter = 2000
#' )
#' mort$Plot
#' }
#'
#' @import data.table
#' @importFrom stats rnorm quantile median
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line labs theme_bw theme
#' @importFrom ggplot2 element_text scale_x_continuous scale_y_continuous expansion
#' @export
get_stochastic_mortality <- function(
    method = c("CW", "PW", "L"),
    Linf,
    L0,
    Lmat,
    tmat,
    growth_model = c("vb", "gompertz", "logistic"),
    Linf_factor = 0.99,
    age_seq = function(tmax) seq(0.1, ceiling(tmax), length.out = 500),
    iter = 2000,
    scaled = TRUE,
    M_target = NULL,
    p = 0.001,
    # CW two-phase parameters
    two_phase = FALSE,
    late_model = c("gompertz", "logistic"),
    tm_factor = 2/3,
    M_mult = 2,
    smooth_factor = 1/3,
    # Lorenzen parameters
    weight_based = FALSE,
    lw_fun = NULL,
    seed = 1234,
    # Plot aesthetics
    palette = c("synthwave", "viridis", "okabe", "plasma", "inferno"),
    print_plot = TRUE,
    show_progress = TRUE
) {

  method <- match.arg(method)
  growth_model <- match.arg(growth_model)
  late_model <- match.arg(late_model)
  palette <- match.arg(palette)

  # Validate inputs
  if (length(Linf) != 2L || length(L0) != 2L || length(Lmat) != 2L || length(tmat) != 2L) {
    stop("Life-history parameters must each be c(mean, sd).", call. = FALSE)
  }

  set.seed(seed)

  # -------------------------------------------------------------------------
  # Parameter Sampling
  # -------------------------------------------------------------------------

  if (show_progress) message("Sampling life-history parameters...")

  Linf_draws <- stats::rnorm(iter, Linf[1], Linf[2])
  L0_draws   <- stats::rnorm(iter, L0[1], L0[2])
  Lmat_draws <- stats::rnorm(iter, Lmat[1], Lmat[2])
  tmat_draws <- stats::rnorm(iter, tmat[1], tmat[2])

  # Ensure biological constraints
  Linf_draws <- pmax(Linf_draws, Lmat_draws + 1)
  Lmat_draws <- pmax(Lmat_draws, L0_draws + 1)
  L0_draws   <- pmax(L0_draws, 0.1)
  tmat_draws <- pmax(tmat_draws, 0.1)

  # Compute M_inf for each draw (VB-derived)
  Minf_draws <- compute_Minf(Linf_draws, L0_draws, Lmat_draws, tmat_draws, warn = FALSE)

  # Compute native k for each draw
  k_native_draws <- compute_k_native(Linf_draws, L0_draws, Lmat_draws, tmat_draws,
                                     growth_model = growth_model, warn = FALSE)

  # Compute tmax for each draw
  tmax_draws <- mapply(
    compute_tmax,
    Linf = Linf_draws,
    L0 = L0_draws,
    k = k_native_draws,
    MoreArgs = list(growth_model = growth_model, Linf_factor = Linf_factor)
  )

  par_draws <- data.table::data.table(
    set_id   = seq_len(iter),
    Linf     = Linf_draws,
    L0       = L0_draws,
    Lmat     = Lmat_draws,
    tmat     = tmat_draws,
    Minf     = Minf_draws,
    k_native = k_native_draws,
    tmax     = tmax_draws
  )

  # Remove invalid draws
  valid_mask <- !is.na(par_draws$Minf) &
    par_draws$Minf > 0 &
    !is.na(par_draws$k_native) &
    par_draws$k_native > 0 &
    par_draws$tmax > 0 &
    is.finite(par_draws$tmax)

  n_invalid <- sum(!valid_mask)
  if (n_invalid > 0) {
    if (show_progress) {
      message(sprintf("Removing %d invalid parameter sets (%.1f%%)",
                      n_invalid, 100 * n_invalid / iter))
    }
    par_draws <- par_draws[valid_mask]
  }

  if (nrow(par_draws) < 100) {
    stop("Fewer than 100 valid parameter sets remain. Check input data.", call. = FALSE)
  }

  # -------------------------------------------------------------------------
  # Mortality Calculation
  # -------------------------------------------------------------------------

  if (show_progress) message(sprintf("Computing %s mortality schedules...", method))

  # Determine common age grid (based on median tmax)
  median_tmax <- stats::median(par_draws$tmax)
  if (is.function(age_seq)) {
    ages <- age_seq(median_tmax)
  } else {
    ages <- age_seq
  }
  ages <- ages[ages > 0]

  n_sets <- nrow(par_draws)
  step <- max(1L, as.integer(n_sets * 0.1))
  notify <- unique(c(seq(step, n_sets, by = step), n_sets))

  schedules_list <- vector("list", n_sets)

  for (i in seq_len(n_sets)) {

    p_i <- par_draws[i]

    # Compute mortality based on method
    M_raw <- switch(
      method,

      "CW" = {
        G_t <- compute_G(ages, p_i$Linf, p_i$L0, p_i$k_native, growth_model)
        M_cw <- p_i$Minf / G_t

        if (two_phase && !is.na(p_i$tmat)) {
          tm <- tm_factor * p_i$tmat
          G_at_tm <- compute_G(tm, p_i$Linf, p_i$L0, p_i$k_native, growth_model)
          M_early_at_tm <- p_i$Minf / G_at_tm

          if (late_model == "gompertz") {
            M_s <- M_mult * M_early_at_tm
            r <- log(M_early_at_tm / M_s) / (tm - p_i$tmax)
            M_late <- M_s * exp(r * (ages - p_i$tmax))
          } else {
            K <- 2 * M_mult * M_early_at_tm
            r <- -log(K / M_early_at_tm - 1) / (tm - p_i$tmax)
            M_late <- K / (1 + exp(-r * (ages - p_i$tmax)))
          }

          smooth_width <- smooth_factor * abs(p_i$tmax - tm)
          weight <- 1 / (1 + exp(-2 * (ages - tm) / smooth_width))

          M_cw <- (1 - weight) * M_cw + weight * M_late
          early_idx <- ages < tm * 0.5
          G_early <- compute_G(ages[early_idx], p_i$Linf, p_i$L0, p_i$k_native, growth_model)
          M_cw[early_idx] <- p_i$Minf / G_early
          M_cw <- pmax(M_cw, 0)
        }
        M_cw
      },

      "PW" = {
        if (is.null(lw_fun)) {
          stop("Peterson-Wroblewski requires 'lw_fun'.", call. = FALSE)
        }
        L_t <- compute_L(ages, p_i$Linf, p_i$L0, p_i$k_native, growth_model)
        W_t <- lw_fun(L_t)
        W_t <- pmax(W_t, 1)
        1.92 * W_t^(-0.25)
      },

      "L" = {
        if (weight_based) {
          if (is.null(lw_fun)) {
            stop("Weight-based Lorenzen requires 'lw_fun'.", call. = FALSE)
          }
          L_t <- compute_L(ages, p_i$Linf, p_i$L0, p_i$k_native, growth_model)
          W_t <- lw_fun(L_t)
          W_t <- pmax(W_t, 1)
          alpha <- stats::rnorm(1, 3.69, 0.502)
          beta  <- stats::rnorm(1, -0.305, 0.029)
          alpha * W_t^beta
        } else {
          L_t <- compute_L(ages, p_i$Linf, p_i$L0, p_i$k_native, growth_model)
          L_ratio <- L_t / p_i$Linf
          intercept <- stats::rnorm(1, 0.28, 0.105)
          coef_L    <- stats::rnorm(1, -1.30, 0.059)
          coef_k    <- stats::rnorm(1, 1.08, 0.082)
          log_M <- intercept + coef_L * log(L_ratio) + coef_k * log(p_i$Minf)
          exp(log_M)
        }
      }
    )

    # Scale if requested
    if (scaled) {
      if (is.null(M_target)) {
        M_target_i <- -log(p) / p_i$tmax
      } else if (is.function(M_target)) {
        M_target_i <- M_target(p_i$tmax)
      } else {
        M_target_i <- M_target
      }
      M_final <- M_raw / mean(M_raw, na.rm = TRUE) * M_target_i
    } else {
      M_final <- M_raw
    }

    # Store with column name "M" for compatibility with simulate_survivorship
    schedules_list[[i]] <- data.table::data.table(
      set_id = p_i$set_id,
      age    = ages,
      M      = M_final
    )

    if (show_progress && i %in% notify) {
      message(sprintf("  Progress: %d/%d (%.0f%%)", i, n_sets, 100 * i / n_sets))
    }
  }

  # Combine all schedules
  schedules <- data.table::rbindlist(schedules_list)

  # -------------------------------------------------------------------------
  # Summary Statistics
  # -------------------------------------------------------------------------

  if (show_progress) message("Computing summary statistics...")

  schedules[, age_round := round(age, 2)]

  summary_dt <- schedules[, .(
    M_median = stats::median(M, na.rm = TRUE),
    M_mean   = mean(M, na.rm = TRUE),
    M_lower  = stats::quantile(M, 0.025, na.rm = TRUE),
    M_upper  = stats::quantile(M, 0.975, na.rm = TRUE)
  ), by = age_round]

  tmax_summary <- par_draws[, .(
    mean  = mean(tmax),
    lower = stats::quantile(tmax, 0.025),
    upper = stats::quantile(tmax, 0.975)
  )]

  # -------------------------------------------------------------------------
  # Plotting
  # -------------------------------------------------------------------------

  if (show_progress) message("Generating plot...")

  pal <- switch(
    palette,
    "synthwave" = c("#FF6B9D", "#C490D1", "#9B6DFF", "#00D4AA"),
    "viridis"   = c("#440154", "#31688E", "#35B779", "#FDE725"),
    "okabe"     = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
    "plasma"    = c("#0D0887", "#7E03A8", "#CC4678", "#F89441"),
    "inferno"   = c("#000004", "#57106E", "#BC3754", "#F98E09")
  )

  fill_color <- pal[1]
  line_color <- pal[3]

  caption <- sprintf(
    "Estimated tmax: %.1f yrs (95%% CI: %.1f - %.1f) | Method: %s | Growth: %s",
    tmax_summary$mean, tmax_summary$lower, tmax_summary$upper,
    method, growth_model
  )

  mort_plot <- ggplot2::ggplot(summary_dt, ggplot2::aes(x = age_round)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = M_lower, ymax = M_upper),
      fill = fill_color, alpha = 0.4
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = M_median),
      color = line_color, linewidth = 1.2
    ) +
    ggplot2::labs(
      x       = "Age (years)",
      y       = "Instantaneous Mortality (M)",
      title   = "Age-Specific Natural Mortality",
      subtitle = "Median with 95% credible interval",
      caption = caption
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.05))) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(face = "italic"),
      plot.caption  = ggplot2::element_text(hjust = 0, size = 9),
      axis.title    = ggplot2::element_text(face = "bold")
    )

  if (print_plot) print(mort_plot)

  if (show_progress) message("Done.")

  # -------------------------------------------------------------------------
  # Return
  # -------------------------------------------------------------------------

  # Remove temporary age_round column from schedules
  schedules[, age_round := NULL]

  list(
    Schedules  = schedules,
    Parameters = par_draws,
    Summary    = summary_dt,
    Plot       = mort_plot
  )
}


# =============================================================================
# SECTION 4: HELPER FUNCTIONS
# =============================================================================


#' Approximate Standard Deviation from Confidence Interval
#'
#' @description
#' Estimates standard deviation from reported confidence interval bounds
#' assuming a normal distribution.
#'
#' @param lower Lower bound of confidence interval.
#' @param upper Upper bound of confidence interval.
#' @param level Confidence level (default 0.95).
#'
#' @return Estimated standard deviation.
#'
#' @examples
#' # If 95% CI is (10, 20), approximate SD
#' approx_sd(10, 20, 0.95)
#'
#' @export
approx_sd <- function(lower, upper, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  (upper - lower) / (2 * z)
}

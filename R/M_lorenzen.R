#' Lorenzen Natural Mortality Model
#'
#' @description
#' Computes size-dependent natural mortality following Lorenzen (1996, 2022).
#' Supports both weight-based and growth-based formulations.
#'
#' @details
#' Two formulations are available:
#'
#' \strong{Weight-based (Lorenzen 1996):}
#' \deqn{M(W) = \alpha \cdot W^{\beta}}
#' where \eqn{\alpha \sim N(3.69, 0.502)} and \eqn{\beta \sim N(-0.305, 0.029)}.
#' Body weight \eqn{W(t)} is computed from predicted length via a user-supplied
#' length-weight function using the native growth model trajectory.
#'
#' \strong{Growth-based (Lorenzen 2022):}
#' \deqn{\ln M = 0.28 - 1.30 \ln(L(t)/L_\infty) + 1.08 \ln(k_{VB})}
#' This formulation was calibrated using von Bertalanffy parameters. The
#' \eqn{L(t)/L_\infty} ratio is computed from the native growth model for
#' accurate relative-size estimates, while \eqn{k_{VB}} must be the
#' VB-equivalent growth coefficient (or native VB \eqn{k}) because the
#' regression coefficients were calibrated against VB parameters.
#'
#' @section Two Roles of k (Growth-Based Formulation):
#' Like Chen-Watanabe, the growth-based Lorenzen involves \eqn{k} in two
#' distinct roles:
#' \describe{
#'   \item{Calibration coefficient}{\eqn{\ln(k_{VB})} enters the regression
#'     equation as a predictor. Since Lorenzen (2022) calibrated this term
#'     against VB parameters, it must be the VB growth coefficient or
#'     VB-equivalent. Passed via \code{k_vb}.}
#'   \item{Growth trajectory}{\eqn{L(t)/L_\infty} represents relative body
#'     size. The native growth model produces more accurate predictions than
#'     a VB approximation. Controlled by \code{k} and \code{growth_model}.}
#' }
#'
#' @param age Numeric vector of ages at which to compute mortality.
#' @param Linf Asymptotic length.
#' @param L0 Length at birth.
#' @param k Growth coefficient native to the fitted model (used for \eqn{L(t)}
#'   prediction via \code{growth_model}).
#' @param k_vb VB-equivalent growth coefficient for the growth-based
#'   formulation's \eqn{\ln(k)} term. Defaults to \code{k}, which is correct
#'   for VB fits. For non-VB fits, use \code{\link{compute_k_vb_equivalent}} or
#'   the derivative-matching approach. Ignored for weight-based formulation.
#' @param lw_fun Function mapping length to weight in grams: \code{lw_fun(L)}.
#'   Required if \code{weight_based = TRUE}.
#' @param weight_based Logical. If \code{TRUE}, uses weight-based formulation.
#'   If \code{FALSE} (default), uses growth-based formulation.
#' @param growth_model Character. Growth model for length prediction:
#'   \code{"vb"}, \code{"gompertz"}, or \code{"logistic"}. Default \code{"vb"}.
#' @param sample_params Logical. If \code{TRUE}, samples allometric/regression
#'   parameters from their distributions for each call. If \code{FALSE}, uses
#'   mean values.
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
#' @examples
#' \dontrun{
#' ages <- seq(0.5, 30, by = 0.5)
#' lw_fun <- function(L) 0.0001 * L^3.1
#'
#' # Weight-based with Gompertz trajectory
#' M_wt <- M_lorenzen(
#'   age = ages, Linf = 100, L0 = 25, k = 0.15,
#'   lw_fun = lw_fun, weight_based = TRUE,
#'   growth_model = "gompertz"
#' )
#'
#' # Growth-based with Gompertz trajectory + VB-equivalent k
#' M_gr <- M_lorenzen(
#'   age = ages, Linf = 100, L0 = 25,
#'   k = 0.15,                         # Native Gompertz k for L(t)
#'   k_vb = 0.08,                      # VB-equivalent k for coefficient
#'   growth_model = "gompertz"
#' )
#' }
#'
#' @export
M_lorenzen <- function(
    age,
    Linf,
    L0,
    k,
    k_vb = k,
    lw_fun = NULL,
    weight_based = FALSE,
    growth_model = c("vb", "gompertz", "logistic"),
    sample_params = FALSE
) {

  growth_model <- match.arg(growth_model)

  # Predict length at age using NATIVE growth model
  L_t <- .predict_length_at_age(age, Linf, L0, k, growth_model)

  if (weight_based) {
    # ----- Weight-based formulation (Lorenzen 1996) -----
    # Pure weight-based: no k enters the mortality equation directly.
    # Only needs L(t) -> W(t).

    if (is.null(lw_fun) || !is.function(lw_fun)) {
      stop("Weight-based Lorenzen requires 'lw_fun'.", call. = FALSE)
    }

    W_t <- lw_fun(L_t)

    if (sample_params) {
      alpha <- stats::rnorm(1, 3.69, 0.502)
      beta  <- stats::rnorm(1, -0.305, 0.029)
    } else {
      alpha <- 3.69
      beta  <- -0.305
    }

    M <- alpha * W_t^beta

  } else {
    # ----- Growth-based formulation (Lorenzen 2022) -----
    # L(t)/Linf ratio: from native growth model (Role 2: body size)
    # ln(k_vb): VB-equivalent k (Role 1: calibration coefficient)

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

    log_M <- intercept + coef_L * log(L_ratio) + coef_k * log(k_vb)
    M <- exp(log_M)
  }

  M
}

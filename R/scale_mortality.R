# =============================================================================
# Mortality Scaling via Cumulative Hazard
# =============================================================================
#
# Replaces the arithmetic-mean-based scaling in mortality_estimation.R
# with trapezoidal numerical integration of the cumulative hazard.
# =============================================================================


# -----------------------------------------------------------------------------
# Trapezoidal Integration Helper
# -----------------------------------------------------------------------------

#' Trapezoidal numerical integration
#'
#' Computes \eqn{\int_a^b f(x) dx} via the trapezoidal rule given
#' \eqn{(x, y)} pairs. Assumes \code{x} is sorted in ascending order.
#'
#' @param x Numeric vector of grid points (must be sorted ascending).
#' @param y Numeric vector of function values at \code{x}.
#'
#' @return Scalar: the approximate integral.
#'
#' @noRd
.trapz <- function(x, y) {
  n <- length(x)
  if (n < 2) return(0)
  dx <- diff(x)
  sum(dx * (y[-n] + y[-1]) / 2)
}


# -----------------------------------------------------------------------------
# Revised Scale Mortality Function
# -----------------------------------------------------------------------------

#' Scale Mortality Schedule to Target via Cumulative Hazard
#'
#' @description
#' Rescales an age-specific mortality schedule by matching the cumulative
#' hazard to an empirical target (from Hoenig, Then et al., or a specified
#' survival probability).
#'
#' @details
#' The scaling finds a proportional constant \eqn{c} such that:
#' \deqn{M_{scaled}(t) = c \times M_{raw}(t)}
#'
#' The constant \eqn{c} is determined by matching the cumulative hazard to
#' the target. The cumulative hazard is computed via trapezoidal integration:
#' \deqn{H_{raw} = \int_0^{t_{max}} M_{raw}(a) \, da \approx \sum_i
#'   \frac{\Delta a_i}{2} \left[ M(a_i) + M(a_{i+1}) \right]}
#'
#' For the \strong{survival-probability target} (\code{M_target = NULL}):
#' \deqn{c = \frac{-\ln(p)}{H_{raw}}}
#' This ensures \eqn{S(t_{max}) = \exp\!\left(-c \cdot H_{raw}\right) = p}
#' exactly.
#'
#' For a \strong{fixed mean-mortality target} (numeric or function of
#' \code{tmax}):
#' \deqn{c = \frac{\bar{M}_{target} \times t_{max}}{H_{raw}}}
#' This interprets the target as the average mortality over the lifespan,
#' so the cumulative hazard of the scaled schedule equals
#' \eqn{\bar{M}_{target} \times t_{max}}.
#'
#' This approach is preferred over arithmetic-mean-based scaling because it
#' directly constrains the biologically relevant quantity (cumulative survival)
#' rather than a proxy, and is invariant to the age grid spacing.
#'
#' @param M Numeric vector of instantaneous mortality rates.
#' @param age Numeric vector of ages corresponding to \code{M}. Required for
#'   trapezoidal integration. Must be the same length as \code{M} and sorted
#'   in ascending order.
#' @param M_target Target mean mortality. Can be a numeric scalar (fixed
#'   target), a function of tmax (e.g.,
#'   \code{function(tmax) 4.899 * tmax^(-0.916)} for Then et al. 2015), or
#'   \code{NULL} to derive from survival probability \code{p}.
#' @param tmax Maximum age (required if \code{M_target} is a function or
#'   \code{NULL}).
#' @param p Probability of surviving to \code{tmax}. Used only if
#'   \code{M_target = NULL}. Default 0.001 (0.1% survival).
#'
#' @return Numeric vector of scaled mortality rates (same length as \code{M}).
#'
#' @examples
#' \dontrun{
#' ages  <- seq(0.1, 30, length.out = 500)
#' M_raw <- M_chen_watanabe_L0(ages, Linf = 100, L0 = 25, k = 0.1,
#'                              two_phase = FALSE)
#'
#' # Scale to survival probability
#' M_sp <- scale_mortality(M_raw, age = ages, tmax = 30, p = 0.01)
#'
#' # Scale using Then et al. (2015) relationship
#' then_2015 <- function(tmax) 4.899 * tmax^(-0.916)
#' M_then <- scale_mortality(M_raw, age = ages, M_target = then_2015, tmax = 30)
#'
#' # Scale to fixed target
#' M_fixed <- scale_mortality(M_raw, age = ages, M_target = 0.2, tmax = 30)
#' }
#'
#' @export
scale_mortality <- function(M, age, M_target = NULL, tmax = NULL, p = 0.001) {

  if (length(M) != length(age)) {
    stop("M and age must have the same length.", call. = FALSE)
  }

  # Compute cumulative hazard via trapezoidal integration
  H_raw <- .trapz(age, M)

  if (H_raw <= 0) {
    warning("Cumulative hazard is non-positive; scaling not applied.", call. = FALSE)
    return(M)
  }

  # Determine target cumulative hazard
  if (is.null(M_target)) {
    # Survival-probability target: S(tmax) = exp(-c * H_raw) = p
    # => c * H_raw = -log(p)
    if (is.null(tmax)) {
      stop("tmax required when M_target is NULL.", call. = FALSE)
    }
    H_target <- -log(p)

  } else {
    # Fixed or function-based mean mortality target
    if (is.function(M_target)) {
      if (is.null(tmax)) {
        stop("tmax required when M_target is a function.", call. = FALSE)
      }
      M_target <- M_target(tmax)
    }

    if (is.null(tmax)) {
      # If tmax not provided, use the span of the age grid
      tmax <- max(age) - min(age)
    }

    # Target cumulative hazard = M_bar_target * tmax
    H_target <- M_target * tmax
  }

  # Proportional scaling constant
  c_scale <- H_target / H_raw

  M * c_scale
}

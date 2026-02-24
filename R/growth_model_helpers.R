# =============================================================================
# Growth Model Helper Functions
# =============================================================================
#
# Add these functions to helper_functions.R
#
# Contents:
#   .predict_length_at_age()     - Native growth model L(t) prediction
#   .derive_k_vb_derivative()    - Derivative-matching VB-equivalent k
#   .sample_bivariate_normal()   - Cholesky-based bivariate normal sampling
# =============================================================================


# -----------------------------------------------------------------------------
# Predict Length-at-Age from Any Growth Model
# -----------------------------------------------------------------------------

#' Predict length at age using the native growth model equation
#'
#' Dispatches on growth model type to compute L(t) using the appropriate
#' functional form and native growth coefficient. Used internally by
#' mortality functions that need accurate body-size trajectories.
#'
#' @param age Numeric vector of ages.
#' @param Linf Asymptotic length.
#' @param L0 Length at birth.
#' @param k Growth coefficient (native to the specified model).
#' @param growth_model Character: \code{"vb"}, \code{"gompertz"}, or
#'   \code{"logistic"}.
#'
#' @return Numeric vector of predicted lengths (same length as \code{age}).
#'
#' @noRd
.predict_length_at_age <- function(age, Linf, L0, k,
                                   growth_model = c("vb", "gompertz", "logistic")) {

  growth_model <- match.arg(growth_model)

  L_t <- switch(
    growth_model,

    "vb" = Linf - (Linf - L0) * exp(-k * age),

    "gompertz" = {
      r0 <- log(Linf / L0)
      Linf * exp(-r0 * exp(-k * age))
    },

    "logistic" = {
      A <- Linf / L0 - 1
      Linf / (1 + A * exp(-k * age))
    }
  )

  # Numerical safety: ensure L_t is positive
  pmax(L_t, L0 * 0.01)
}


# -----------------------------------------------------------------------------
# Derive VB-Equivalent k via Derivative Matching at Birth
# -----------------------------------------------------------------------------

#' Derive VB-equivalent k from native growth coefficient via derivative matching
#'
#' @description
#' For k-based Gompertz or Logistic fits where maturity milestones are
#' unavailable, derives the VB-equivalent k by matching the instantaneous
#' growth rate at birth (age 0).
#'
#' @details
#' At \eqn{t = 0}, the VB growth rate is \eqn{dL/dt = k_{VB}(L_\infty - L_0)}.
#' Setting this equal to the native model's growth rate at birth and solving:
#'
#' Gompertz: \eqn{k_{VB} = k_g \cdot L_0 \cdot \ln(L_\infty / L_0) / (L_\infty - L_0)}
#'
#' Logistic: \eqn{k_{VB} = k_l \cdot L_0 / L_\infty}
#'
#' For VB models, returns \code{k_native} unchanged.
#'
#' This is a fallback when maturity milestones are unavailable. The
#' milestone-based derivation via \code{\link{compute_k_vb_equivalent}} is
#' preferred when \eqn{L_{mat}} and \eqn{t_{mat}} are known, as it anchors
#' at a biologically meaningful point within the observed data range.
#'
#' @param k_native Numeric vector. Native growth coefficient.
#' @param Linf Numeric vector. Asymptotic length.
#' @param L0 Numeric vector. Length at birth.
#' @param growth_model Character: \code{"vb"}, \code{"gompertz"}, or
#'   \code{"logistic"}.
#'
#' @return Numeric vector of VB-equivalent k values.
#'
#' @noRd
.derive_k_vb_derivative <- function(k_native, Linf, L0,
                                    growth_model = c("vb", "gompertz", "logistic")) {

  growth_model <- match.arg(growth_model)

  switch(
    growth_model,

    # VB: native k IS the VB k
    "vb" = k_native,

    # Gompertz: dL/dt|_0 = k_g * L0 * ln(Linf/L0)
    # VB:       dL/dt|_0 = k_vb * (Linf - L0)
    # => k_vb = k_g * L0 * ln(Linf/L0) / (Linf - L0)
    "gompertz" = k_native * L0 * log(Linf / L0) / (Linf - L0),

    # Logistic: dL/dt|_0 = k_l * L0 * (1 - L0/Linf)
    # VB:       dL/dt|_0 = k_vb * (Linf - L0)
    # => k_vb = k_l * L0 * (1 - L0/Linf) / (Linf - L0)
    #         = k_l * L0 * (Linf - L0) / (Linf * (Linf - L0))
    #         = k_l * L0 / Linf
    "logistic" = k_native * L0 / Linf
  )
}


# -----------------------------------------------------------------------------
# Bivariate Normal Sampling via Cholesky Factorization
# -----------------------------------------------------------------------------

#' Sample from a bivariate normal distribution
#'
#' @description
#' Generates correlated draws from a bivariate normal distribution using
#' Cholesky factorization. Used for joint sampling of maturity parameters
#' (Lmat, tmat) in manual mortality estimation.
#'
#' @param n Integer. Number of draws.
#' @param mu1,mu2 Numeric. Means of the two variables.
#' @param sd1,sd2 Numeric. Standard deviations of the two variables.
#' @param rho Numeric in (-1, 1). Correlation between the two variables.
#'
#' @return A matrix with \code{n} rows and 2 columns.
#'
#' @noRd
.sample_bivariate_normal <- function(n, mu1, sd1, mu2, sd2, rho = 0.5) {

  if (abs(rho) >= 1) {
    stop("rho must be strictly between -1 and 1.", call. = FALSE)
  }

  z1 <- stats::rnorm(n)
  z2 <- stats::rnorm(n)

  x1 <- mu1 + sd1 * z1
  x2 <- mu2 + sd2 * (rho * z1 + sqrt(1 - rho^2) * z2)

  cbind(x1, x2)
}

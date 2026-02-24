# =============================================================================
# Bivariate Normal Sampling via Cholesky Factorization
# =============================================================================
#
# Add this function to helper_functions.R
#


#' Sample from a Bivariate Normal Distribution
#'
#' @description
#' Generates correlated draws from a bivariate normal distribution using
#' Cholesky factorization. Used for joint sampling of maturity parameters
#' (Lmat, tmat) when specifying manual inputs to mortality estimation.
#'
#' @details
#' Correlated samples are generated via the Cholesky decomposition of the
#' correlation matrix. Given independent standard normals \eqn{Z_1, Z_2}:
#'
#' \deqn{x_1 = \mu_1 + \sigma_1 Z_1}
#' \deqn{x_2 = \mu_2 + \sigma_2 (\rho Z_1 + \sqrt{1 - \rho^2} Z_2)}
#'
#' This produces samples with the specified marginal distributions and
#' correlation \eqn{\rho}.
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

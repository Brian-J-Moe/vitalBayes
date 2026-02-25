#' Peterson-Wroblewski Natural Mortality Model
#'
#' @description
#' Computes weight-based natural mortality following Peterson & Wroblewski
#' (1984). Mortality scales allometrically with body weight.
#'
#' @details
#' The model expresses mortality as a power function of body weight:
#' \deqn{M(W) = 1.92 \cdot W^{-0.25}}
#' where \eqn{W} is body weight in grams.
#'
#' This model is growth-model-agnostic: it only requires predicted body weight
#' at age, which can be derived from any growth model via a length-weight
#' relationship. The \eqn{-0.25} exponent reflects metabolic scaling theory:
#' metabolic rate scales approximately as \eqn{W^{0.75}} (Kleiber's law), and
#' mortality is assumed proportional to mass-specific metabolic rate, yielding
#' a \eqn{W^{-0.25}} dependence.
#'
#' Because \eqn{k} does not appear in the mortality equation itself, the native
#' growth coefficient and native growth model should always be used for
#' \eqn{L(t)} prediction --- no VB-equivalent conversion is needed.
#'
#' @param age Numeric vector of ages at which to compute mortality.
#' @param Linf Asymptotic length.
#' @param L0 Length at birth.
#' @param k Growth coefficient (native to the fitted model).
#' @param lw_fun Function mapping length to weight in grams: \code{lw_fun(L)}.
#' @param growth_model Character. Growth model for length prediction:
#'   \code{"vb"}, \code{"gompertz"}, or \code{"logistic"}. Default \code{"vb"}.
#'
#' @return Numeric vector of instantaneous mortality rates.
#'
#' @references
#' Peterson, I., & Wroblewski, J. S. (1984). Mortality rate of fishes in the
#' pelagic ecosystem. \emph{Canadian Journal of Fisheries and Aquatic Sciences},
#' 41(7), 1117-1120.
#'
#' @examples
#' \dontrun{
#' lw_fun <- function(L) 0.0001 * L^3.1  # Length in cm, weight in g
#' ages <- seq(0.5, 30, by = 0.5)
#'
#' M <- M_peterson_wroblewski(
#'   age = ages, Linf = 100, L0 = 25, k = 0.1,
#'   lw_fun = lw_fun, growth_model = "gompertz"
#' )
#' }
#'
#' @export
M_peterson_wroblewski <- function(
    age,
    Linf,
    L0,
    k,
    lw_fun,
    growth_model = c("vb", "gompertz", "logistic")
) {

  growth_model <- match.arg(growth_model)

  if (is.null(lw_fun) || !is.function(lw_fun)) {
    stop("PW model requires a length-weight function 'lw_fun(L)'.", call. = FALSE)
  }

  # Predict length at age using native growth model
  L_t <- .predict_length_at_age(age, Linf, L0, k, growth_model)

  # Convert to weight and compute mortality
  W_t <- lw_fun(L_t)
  1.92 * W_t^(-0.25)
}

#' Chen-Watanabe Natural Mortality (\eqn{L_0} Parameterization)
#'
#' @description
#' Computes age-specific natural mortality using the Chen & Watanabe (1989)
#' model with an \eqn{L_0} parameterization that eliminates dependence on the
#' theoretical parameter \eqn{t_0}.
#'
#' @details
#' The standard Chen-Watanabe formulation expresses mortality as:
#' \deqn{M(t) = \frac{k}{1 - e^{-k(t - t_0)}}}
#'
#' where \eqn{t_0} is the theoretical age at length zero --- a parameter with no
#' direct biological interpretation that can take implausible values,
#' particularly when growth data are sparse.
#'
#' We reparameterize using the relationship between \eqn{t_0} and \eqn{L_0}
#' (birth length) under von Bertalanffy dynamics:
#' \deqn{L_0 = L_\infty(1 - e^{kt_0})}
#'
#' After algebraic manipulation (see vignette), the \eqn{L_0}-parameterized form
#' becomes:
#' \deqn{M(t) = \frac{k_{VB} \cdot L_\infty}{L(t)}}
#'
#' where \eqn{k_{VB}} is the von Bertalanffy growth coefficient (or
#' VB-equivalent for non-VB fits) and \eqn{L(t)} is the predicted length at
#' age \eqn{t} from the native growth model.
#'
#' @section Two Roles of k:
#' The Chen-Watanabe equation involves \eqn{k} in two conceptually distinct
#' roles:
#' \describe{
#'   \item{Asymptotic mortality rate}{The VB-specific constant
#'     \eqn{M_\infty = k_{VB}} that governs the theoretical minimum mortality
#'     as \eqn{t \to \infty}. This enters the CW formula as the numerator.
#'     Because the derivation assumes VB dynamics, this must always be the
#'     VB growth coefficient (or VB-equivalent for non-VB fits). Passed via
#'     the \code{k} argument.}
#'   \item{Growth trajectory}{The length-at-age prediction \eqn{L(t)} that
#'     determines how far an individual is from asymptotic size. For non-VB
#'     growth models, the native trajectory with native \eqn{k} produces
#'     more accurate body-size estimates than the VB approximation. Controlled
#'     by \code{k_native} and \code{growth_model}.}
#' }
#'
#' For VB fits, both roles use the same \eqn{k}, and the default
#' \code{k_native = k} preserves full backward compatibility.
#'
#' @section Two-Phase Extension:
#' The original CW model produces unrealistic mortality trajectories at old ages
#' (approaching zero asymptotically). When \code{two_phase = TRUE}, the model
#' adds a senescence component where mortality increases after maturity,
#' more realistically capturing late-life dynamics. Mortality follows the
#' standard CW model until age \eqn{t_m} (a fraction of \eqn{t_{mat}}), then
#' transitions to a senescence model (Gompertz or logistic) that increases
#' mortality toward \eqn{t_{max}}.
#'
#' @param age Numeric vector of ages at which to compute mortality.
#' @param Linf Asymptotic length.
#' @param L0 Length at birth.
#' @param k VB-equivalent growth coefficient (used as \eqn{M_\infty} in the CW
#'   formula). For VB fits, this is the native \eqn{k}. For Gompertz or
#'   Logistic fits, use \code{\link{compute_k_vb_equivalent}} or
#'   \code{.derive_k_vb_derivative}.
#' @param k_native Growth coefficient native to the fitted model (used for
#'   \eqn{L(t)} prediction). Defaults to \code{k}, which is correct for VB
#'   fits. For non-VB fits, pass the model-specific growth coefficient here.
#' @param growth_model Character. Growth model for \eqn{L(t)} prediction:
#'   \code{"vb"} (default), \code{"gompertz"}, or \code{"logistic"}.
#' @param tmax Maximum age. If \code{NULL}, estimated from growth parameters
#'   as age when \eqn{L(t) = } \code{Linf_factor} \eqn{\times L_\infty}
#'   (using the VB equation with \code{k} for consistency with the CW
#'   theoretical framework).
#' @param Linf_factor Numeric in (0, 1). Fraction of \eqn{L_\infty} used to
#'   estimate \eqn{t_{max}}. Default 0.99.
#' @param two_phase Logical. If \code{TRUE}, applies two-phase model with
#'   late-life senescence. Default \code{TRUE}.
#' @param tmat Age at maturity. Required if \code{two_phase = TRUE}.
#' @param late_model Character. Senescence model: \code{"gompertz"} (default)
#'   or \code{"logistic"}.
#' @param tm_factor Numeric. Fraction of \eqn{t_{mat}} at which transition to
#'   senescence begins. Default 2/3.
#' @param M_mult Numeric. Multiplier for senescence mortality plateau relative
#'   to mortality at \eqn{t_m}. Default 2.
#' @param smooth_factor Numeric. Controls smoothness of transition between
#'   phases. Default 1/3.
#'
#' @return Numeric vector of instantaneous mortality rates (same length as
#'   \code{age}).
#'
#' @references
#' Chen, S., & Watanabe, S. (1989). Age dependence of natural mortality
#' coefficient in fish population dynamics. \emph{Nippon Suisan Gakkaishi},
#' 55(2), 205-208.
#'
#' @examples
#' \dontrun{
#' ages <- seq(0.1, 30, by = 0.5)
#'
#' # VB fit (k serves both roles)
#' M_vb <- M_chen_watanabe_L0(
#'   age = ages, Linf = 100, L0 = 25, k = 0.1,
#'   two_phase = TRUE, tmat = 10
#' )
#'
#' # Gompertz fit (separate k for each role)
#' M_gomp <- M_chen_watanabe_L0(
#'   age = ages, Linf = 100, L0 = 25,
#'   k = 0.08,                       # VB-equivalent k for M_inf
#'   k_native = 0.12,                # Native Gompertz k for L(t)
#'   growth_model = "gompertz",
#'   two_phase = TRUE, tmat = 10
#' )
#' }
#'
#' @seealso \code{\link{compute_k_vb_equivalent}} for deriving \eqn{k_{VB}} from
#'   maturity milestones, \code{\link{get_stochastic_mortality}} for Monte Carlo
#'   mortality estimation with uncertainty.
#'
#' @export
M_chen_watanabe_L0 <- function(
    age,
    Linf,
    L0,
    k,
    k_native = k,
    growth_model = c("vb", "gompertz", "logistic"),
    tmax = NULL,
    Linf_factor = 0.99,
    two_phase = TRUE,
    tmat = NULL,
    late_model = c("gompertz", "logistic"),
    tm_factor = 2/3,
    M_mult = 2,
    smooth_factor = 1/3
) {

  late_model   <- match.arg(late_model)
  growth_model <- match.arg(growth_model)

  # Validate inputs
  if (L0 >= Linf) {
    stop("L0 must be less than Linf.", call. = FALSE)
  }
  if (k <= 0) {
    stop("k (VB-equivalent) must be positive.", call. = FALSE)
  }
  if (k_native <= 0) {
    stop("k_native must be positive.", call. = FALSE)
  }
  if (two_phase && is.null(tmat)) {
    stop("tmat required for two-phase model.", call. = FALSE)
  }

  # Estimate tmax if not provided.
  # Uses VB k for consistency with CW theoretical framework:
  # tmax represents the age at which the VB curve reaches Linf_factor * Linf,
  # the lifespan implied by the CW model's asymptotic structure.
  if (is.null(tmax)) {
    tmax <- -log(Linf * (1 - Linf_factor) / (Linf - L0)) / k
  }

  # Predicted length at each age (NATIVE growth model)
  L_t <- .predict_length_at_age(age, Linf, L0, k_native, growth_model)

  # Core CW mortality: M(t) = k_vb * Linf / L(t)
  #   k (= k_vb)  → asymptotic mortality rate constant (Role 1)
  #   L(t)         → native growth trajectory (Role 2)
  M_cw <- k * Linf / L_t

  if (!two_phase) {
    return(M_cw)
  }

  # ----- Two-Phase Extension -----

  # Transition age (fraction of maturity age)
  tm <- tm_factor * tmat

  # Early-phase CW mortality at transition point:
  # M_inf (from k_vb) evaluated at L(tm) from native growth model
  L_at_tm <- .predict_length_at_age(tm, Linf, L0, k_native, growth_model)
  M_early_at_tm <- k * Linf / L_at_tm

  # Late-phase mortality setup
  if (late_model == "gompertz") {
    # Gompertz senescence: M(t) = M_s * exp(r * (t - tmax))
    # Continuity: M(tm) = M_early_at_tm
    # Terminal: M(tmax) = M_mult * M_early_at_tm
    M_s <- M_mult * M_early_at_tm
    r <- log(M_early_at_tm / M_s) / (tm - tmax)

    M_late <- function(t) M_s * exp(r * (t - tmax))

  } else {
    # Logistic senescence: M(t) = K / (1 + exp(-r * (t - tmax)))
    K <- 2 * M_mult * M_early_at_tm
    r <- -log(K / M_early_at_tm - 1) / (tm - tmax)

    M_late <- function(t) K / (1 + exp(-r * (t - tmax)))
  }

  # Smooth transition between phases
  smooth_width <- smooth_factor * abs(tmax - tm)
  weight <- 1 / (1 + exp(-2 * (age - tm) / smooth_width))

  # Blend early and late phases
  M <- (1 - weight) * M_cw + weight * M_late(age)

  # For ages well below tm, use pure CW
  M[age < tm * 0.5] <- M_cw[age < tm * 0.5]

  # Cap any negative values (numerical artifacts)
  pmax(M, 0)
}

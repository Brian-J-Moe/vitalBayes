# =============================================================================
# vitalBayes Mortality Estimation Functions
# =============================================================================
#
# This module implements natural mortality estimation with key innovations:
#
#   1. L0-PARAMETERIZED CHEN-WATANABE: Eliminates dependence on t0, using
#      observable birth size instead. This aligns with the vitalBayes philosophy
#      of parameterizing models with biologically meaningful quantities.
#
#   2. GROWTH-MODEL-AGNOSTIC k DERIVATION: Computes VB-equivalent k from
#      (Linf, L0, Lmat, tmat) — parameters estimated by ALL growth models.
#      Users can fit whichever growth model best describes their data (VB,
#      Gompertz, or Logistic) and still use Chen-Watanabe mortality estimation.
#
#   3. JOINT POSTERIOR SAMPLING: Preserves correlations between growth
#      parameters when propagating uncertainty to mortality estimates.
#
# =============================================================================


# -----------------------------------------------------------------------------
# Core Helper: Compute VB-Equivalent k from Maturity Parameters
# -----------------------------------------------------------------------------

#' Compute Von Bertalanffy-Equivalent Growth Coefficient
#'
#' @description
#' Derives the von Bertalanffy growth coefficient \eqn{k} from biological
#' milestones \eqn{(L_\infty, L_0, L_{mat}, t_{mat})}. This enables Chen-Watanabe
#' mortality estimation using posteriors from \emph{any} growth model (von
#' Bertalanffy, Gompertz, or Logistic).
#'
#' @details
#' The key insight is that while the three growth models use different functional
#' forms and produce different numerical \eqn{k} values, they all estimate the
#' same underlying biological quantities: asymptotic length, birth size, and
#' maturity milestones. The VB-equivalent \eqn{k} is computed as:
#'
#' \deqn{k_{VB}^{equiv} = \frac{1}{t_{mat}} \ln\left(\frac{L_\infty - L_0}{L_\infty - L_{mat}}\right)}
#'
#' This is the growth coefficient that would produce a von Bertalanffy curve
#' passing through the biological milestones \eqn{(0, L_0)} and
#' \eqn{(t_{mat}, L_{mat})} with asymptote \eqn{L_\infty}.
#'
#' When the input growth fit is von Bertalanffy with maturity-based
#' parameterization, this computation exactly reproduces the fitted \eqn{k}.
#' When the input is Gompertz or Logistic, it produces the VB-equivalent
#' \eqn{k} that encodes the same biological growth information.
#'
#' The von Bertalanffy model tends to produce unstable \eqn{L_\infty} estimates
#' when data are sparse at older ages - a common situation in elasmobranch
#' research. Gompertz and Logistic models often provide more reliable fits in
#' these cases. By deriving VB-equivalent \eqn{k} from biological milestones,
#' users can select the growth model that best fits their data while still
#' using Chen-Watanabe mortality estimation and maintaining theoretical
#' coherence (since CW was derived under VB assumptions).
#'
#' @param Linf Numeric vector. Asymptotic length posterior draws.
#' @param L0 Numeric vector. Length at birth posterior draws.
#' @param Lmat Numeric vector. Length at maturity posterior draws.
#' @param tmat Numeric vector. Age at maturity posterior draws.
#' @param warn Logical. If \code{TRUE} (default), warns when draws produce
#'   invalid \eqn{k} values.
#'
#' @return Numeric vector of VB-equivalent \eqn{k} values. Invalid values
#'   (from non-positive arguments to log) are returned as \code{NA}.
#'
#' @examples
#' \dontrun{
#' # Extract from any growth model posterior
#' draws <- extract_growth_parameters(growth_fit, sex = 1)
#' k_vb <- compute_k_vb_equivalent(
#'   Linf = draws$Linf,
#'   L0   = draws$L0,
#'   Lmat = draws$Lmat,
#'   tmat = draws$tmat
#' )
#'
#' # Direct specification for sensitivity analysis
#' k_vb <- compute_k_vb_equivalent(
#'   Linf = rnorm(1000, 100, 5),
#'   L0   = rnorm(1000, 25, 2),
#'   Lmat = rnorm(1000, 70, 3),
#'   tmat = rnorm(1000, 10, 1)
#' )
#' }
#'
#' @seealso \code{\link{M_chen_watanabe_L0}} for the L0-parameterized CW model,
#'   \code{\link{extract_growth_parameters}} for posterior extraction.
#'
#' @export
compute_k_vb_equivalent <- function(Linf, L0, Lmat, tmat, warn = TRUE) {


  # Validate inputs

  n <- length(Linf)
  if (!all(c(length(L0), length(Lmat), length(tmat)) == n)) {
    stop("All input vectors must have the same length.", call. = FALSE)
  }

  # Compute VB-equivalent k: k = (1/tmat) * ln((Linf - L0) / (Linf - Lmat))
  numerator   <- Linf - L0
  denominator <- Linf - Lmat

  # Identify invalid values (would produce NaN or Inf)
  invalid <- numerator <= 0 | denominator <= 0 | tmat <= 0 |
    is.na(numerator) | is.na(denominator) | is.na(tmat)

  if (warn && any(invalid, na.rm = TRUE)) {
    n_invalid <- sum(invalid, na.rm = TRUE)
    warning(
      sprintf(
        "%d of %d draws (%.1f%%) produced invalid k values. Common causes:\n
        - Lmat >= Linf (maturity size exceeds asymptotic size)\n
        - L0 >= Linf (birth size exceeds asymptotic size)\n
        - tmat <= 0 (non-positive maturity age)\n
        These draws will be excluded from mortality calculations.",
        n_invalid, n, 100 * n_invalid / n
      ),
      call. = FALSE
    )
  }

  k <- (1 / tmat) * log(numerator / denominator)
  k[invalid] <- NA_real_

  k
}


# -----------------------------------------------------------------------------
# Core Helper: Extract Growth Parameters from Any Model Fit
# -----------------------------------------------------------------------------

#' Extract Life History Parameters from Growth Model Posterior
#'
#' @description
#' Extracts posterior draws of \eqn{(L_\infty, L_0, L_{mat}, t_{mat})} from a
#' vitalBayes growth model fit. Works identically for von Bertalanffy, Gompertz,
#' and Logistic models fitted via \code{\link{fit_bayesian_growth}}.
#'
#' @details
#' This function provides a unified interface for extracting the biological
#' parameters common to all growth models. These parameters — asymptotic length,
#' birth size, and maturity milestones — represent real biological quantities
#' that exist independently of the mathematical model used to describe growth.
#'
#' For maturity-based growth fits (\code{k_based = FALSE}), \eqn{L_{mat}} and
#' \eqn{t_{mat}} are directly estimated parameters. For k-based fits, these
#' must be supplied separately via \code{maturity_fit}.
#'
#' @param growth_fit A \code{CmdStanMCMC} object from
#'   \code{\link{fit_bayesian_growth}}.
#' @param maturity_fit Optional \code{CmdStanMCMC} object from
#'   \code{\link{fit_bayesian_maturity}} providing age-at-maturity. Required
#'   for k-based growth fits if \code{tmat} is needed.
#' @param sex Integer. Sex code (1 = female, 2 = male) for hierarchical models.
#'   If \code{NULL}, extracts from single-sex model or uses column 1.
#' @param n_draws Integer. Number of posterior draws to return. If \code{NULL},
#'   returns all available draws. If specified, draws are subsampled randomly.
#' @param seed Integer. Random seed for reproducible subsampling.
#'
#' @return A \code{data.table} with columns:
#' \describe{
#'   \item{draw}{Integer draw index}
#'   \item{Linf}{Asymptotic length}
#'   \item{L0}{Length at birth}
#'   \item{Lmat}{Length at maturity (if available)}
#'   \item{tmat}{Age at maturity (if available)}
#'   \item{k}{Growth coefficient (original model's k, not VB-equivalent)}
#'   \item{k_vb_equiv}{VB-equivalent k (if Lmat and tmat available)}
#' }
#'
#' @examples
#' \dontrun{
#' # From a Gompertz fit with maturity-based parameterization
#' params <- extract_growth_parameters(gomp_fit, sex = 1, n_draws = 2000)
#'
#' # Check VB-equivalent k distribution
#' hist(params$k_vb_equiv, main = "VB-Equivalent k from Gompertz Fit")
#'
#' # Compare to original Gompertz k
#' plot(params$k, params$k_vb_equiv,
#'      xlab = "Gompertz k", ylab = "VB-equivalent k")
#' }
#'
#' @import data.table
#' @export
extract_growth_parameters <- function(
    growth_fit,
    maturity_fit = NULL,
    sex = NULL,
    n_draws = NULL,
    seed = 1234
) {

  # Extract available parameters
  available_params <- growth_fit$metadata()$stan_variables

  # Determine if hierarchical (2-sex) model
  Linf_draws <- growth_fit$draws("Linf", format = "matrix")
  is_hierarchical <- ncol(Linf_draws) > 1

  # Set sex index
  if (is.null(sex)) {
    s <- 1L
    if (is_hierarchical) {
      message("Hierarchical model detected but sex not specified. Using sex = 1 (female).")
    }
  } else {
    s <- as.integer(sex)
    if (s < 1 || s > 2) stop("sex must be 1 (female) or 2 (male).", call. = FALSE)
  }

  # Extract core parameters
  if (is_hierarchical) {
    Linf <- Linf_draws[, s]
    L0   <- growth_fit$draws("L0", format = "matrix")[, s]
    k    <- growth_fit$draws("k", format = "matrix")[, s]
  } else {
    Linf <- as.vector(Linf_draws)
    L0   <- as.vector(growth_fit$draws("L0", format = "matrix"))
    k    <- as.vector(growth_fit$draws("k", format = "matrix"))
  }

  n_total <- length(Linf)

  # Extract maturity parameters if available (maturity-based models)
  has_Lmat <- "Lmat" %in% available_params
  has_tmat <- "tmat" %in% available_params

  if (has_Lmat) {
    Lmat_draws <- growth_fit$draws("Lmat", format = "matrix")
    Lmat <- if (is_hierarchical) Lmat_draws[, s] else as.vector(Lmat_draws)
  } else {
    Lmat <- rep(NA_real_, n_total)
  }

  if (has_tmat) {
    tmat_draws <- growth_fit$draws("tmat", format = "matrix")
    tmat <- if (is_hierarchical) tmat_draws[, s] else as.vector(tmat_draws)
  } else if (!is.null(maturity_fit)) {
    # Extract from separate maturity fit
    t50_draws <- maturity_fit$draws("t50", format = "matrix")
    mat_hierarchical <- ncol(t50_draws) > 1
    tmat_raw <- if (mat_hierarchical) t50_draws[, s] else as.vector(t50_draws)

    # Match lengths if different number of draws
    if (length(tmat_raw) != n_total) {
      set.seed(seed)
      tmat <- sample(tmat_raw, n_total, replace = TRUE)
    } else {
      tmat <- tmat_raw
    }
  } else {
    tmat <- rep(NA_real_, n_total)
  }

  # Compute VB-equivalent k if maturity parameters available
  if (has_Lmat && (has_tmat || !is.null(maturity_fit))) {
    k_vb_equiv <- compute_k_vb_equivalent(Linf, L0, Lmat, tmat, warn = FALSE)
  } else {
    k_vb_equiv <- rep(NA_real_, n_total)
  }

  # Build output data.table
  result <- data.table::data.table(
    draw      = seq_len(n_total),
    Linf      = Linf,
    L0        = L0,
    Lmat      = Lmat,
    tmat      = tmat,
    k         = k,
    k_vb_equiv = k_vb_equiv
  )

  # Subsample if requested
  if (!is.null(n_draws) && n_draws < n_total) {
    set.seed(seed)
    idx <- sample(n_total, n_draws, replace = FALSE)
    result <- result[idx]
    result[, draw := seq_len(.N)]
  }

  result
}


# -----------------------------------------------------------------------------
# L0-Parameterized Chen-Watanabe Model
# -----------------------------------------------------------------------------

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
#' where \eqn{t_0} is the theoretical age at length zero - a parameter with no
#' direct biological interpretation that can take implausible values,
#' particularly when growth data are sparse.
#'
#' We reparameterize using the relationship between \eqn{t_0} and \eqn{L_0}
#' (birth length) under von Bertalanffy dynamics:
#' \deqn{L_0 = L_\infty(1 - e^{kt_0})}
#'
#' After algebraic manipulation (see vignette), the \eqn{L_0}-parameterized form
#' becomes:
#' \deqn{M(t) = \frac{k \cdot L_\infty}{L(t)}}
#'
#' where \eqn{L(t) = L_\infty - (L_\infty - L_0)e^{-kt}} is the predicted length
#' at age \eqn{t}.
#'
#' This reformulation reveals that Chen-Watanabe mortality is inversely
#' proportional to body size - smaller (younger) individuals experience higher
#' mortality. The ratio \eqn{L_\infty / L(t)} represents how far an individual
#' is from asymptotic size, with mortality declining as this ratio approaches 1.
#'
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
#' @param k VB-equivalent growth coefficient. Can be computed from any growth
#'   model using \code{\link{compute_k_vb_equivalent}}.
#' @param tmax Maximum age. If \code{NULL}, estimated from growth parameters
#'   as age when \eqn{L(t) = } \code{Linf_factor} \eqn{\times L_\infty}.
#' @param Linf_factor Numeric in (0, 1). Fraction of \eqn{L_\infty} used to
#'   estimate \eqn{t_{max}}. Default 0.99 (age at 99% of asymptotic length).
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
#' # Single-phase CW
#' M_single <- M_chen_watanabe_L0(
#'   age = ages, Linf = 100, L0 = 25, k = 0.1,
#'   two_phase = FALSE
#' )
#'
#' # Two-phase with Gompertz senescence
#' M_two <- M_chen_watanabe_L0(
#'   age = ages, Linf = 100, L0 = 25, k = 0.1,
#'   two_phase = TRUE, tmat = 10, late_model = "gompertz"
#' )
#'
#' plot(ages, M_single, type = "l", ylim = c(0, 1))
#' lines(ages, M_two, col = "red")
#' }
#'
#' @seealso \code{\link{compute_k_vb_equivalent}} for deriving \eqn{k} from any
#'   growth model, \code{\link{get_stochastic_mortality}} for Monte Carlo
#'   mortality estimation with uncertainty.
#'
#' @export
M_chen_watanabe_L0 <- function(
    age,
    Linf,
    L0,
    k,
    tmax = NULL,
    Linf_factor = 0.99,
    two_phase = TRUE,
    tmat = NULL,
    late_model = c("gompertz", "logistic"),
    tm_factor = 2/3,
    M_mult = 2,
    smooth_factor = 1/3
) {

  late_model <- match.arg(late_model)

  # Validate inputs
  if (L0 >= Linf) {
    stop("L0 must be less than Linf.", call. = FALSE)
  }
  if (k <= 0) {
    stop("k must be positive.", call. = FALSE)
  }
  if (two_phase && is.null(tmat)) {
    stop("tmat required for two-phase model.", call. = FALSE)
  }

  # Estimate tmax if not provided
  if (is.null(tmax)) {
    # Age when L(t) = Linf_factor * Linf
    # L(t) = Linf - (Linf - L0) * exp(-k*t)
    # Linf_factor * Linf = Linf - (Linf - L0) * exp(-k*t)
    # exp(-k*t) = (Linf - Linf_factor * Linf) / (Linf - L0)
    # exp(-k*t) = Linf * (1 - Linf_factor) / (Linf - L0)
    tmax <- -log(Linf * (1 - Linf_factor) / (Linf - L0)) / k
  }

  # Predicted length at each age (VB equation)
  L_t <- Linf - (Linf - L0) * exp(-k * age)

  # Ensure L_t is positive (numerical safety for very young ages)
  L_t <- pmax(L_t, L0 * 0.01)

  # Core CW mortality: M(t) = k * Linf / L(t)
  M_cw <- k * Linf / L_t

  if (!two_phase) {
    return(M_cw)
  }

  # ----- Two-Phase Extension -----

  # Transition age (fraction of maturity age)
  tm <- tm_factor * tmat

  # Early-phase mortality (ages < tm): use CW
  M_early_at_tm <- k * Linf / (Linf - (Linf - L0) * exp(-k * tm))

  # Late-phase mortality setup
  if (late_model == "gompertz") {
    # Gompertz senescence: M(t) = M_s * exp(r * (t - tmax))
    # At t = tm, we want continuity: M(tm) = M_early_at_tm
    # At t = tmax, we want M(tmax) = M_mult * M_early_at_tm
    M_s <- M_mult * M_early_at_tm
    # Solve for r from M(tm) = M_s * exp(r * (tm - tmax)) = M_early_at_tm
    r <- log(M_early_at_tm / M_s) / (tm - tmax)

    M_late <- function(t) M_s * exp(r * (t - tmax))

  } else {
    # Logistic senescence: M(t) = K / (1 + exp(-r * (t - tmax)))
    # At t = tmax, M = K/2
    # We want M(tmax) ~ M_mult * M_early_at_tm, so K = 2 * M_mult * M_early_at_tm
    K <- 2 * M_mult * M_early_at_tm
    # Solve for r from continuity at tm
    # M(tm) = K / (1 + exp(-r * (tm - tmax))) = M_early_at_tm
    # 1 + exp(-r * (tm - tmax)) = K / M_early_at_tm
    # exp(-r * (tm - tmax)) = K / M_early_at_tm - 1
    r <- -log(K / M_early_at_tm - 1) / (tm - tmax)

    M_late <- function(t) K / (1 + exp(-r * (t - tmax)))
  }

  # Smooth transition between phases
  smooth_width <- smooth_factor * abs(tmax - tm)
  # Logistic weight: 0 at tm, 1 well past tm
  weight <- 1 / (1 + exp(-2 * (age - tm) / smooth_width))

  # Blend early and late phases
  M <- (1 - weight) * M_cw + weight * M_late(age)

  # For ages well below tm, use pure CW
  M[age < tm * 0.5] <- M_cw[age < tm * 0.5]

  # Cap any negative values (numerical artifacts)
  M <- pmax(M, 0)

  M
}


# -----------------------------------------------------------------------------
# Peterson-Wroblewski Model (unchanged, but documented for completeness)
# -----------------------------------------------------------------------------

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
#' relationship.
#'
#' @param age Numeric vector of ages at which to compute mortality.
#' @param Linf Asymptotic length.
#' @param L0 Length at birth.
#' @param k Growth coefficient (model-specific, not necessarily VB).
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

  # Predict length at age based on growth model
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

  # Convert to weight and compute mortality
  W_t <- lw_fun(L_t)
  1.92 * W_t^(-0.25)
}


# -----------------------------------------------------------------------------
# Lorenzen Model
# -----------------------------------------------------------------------------

#' Lorenzen Natural Mortality Model
#'
#' @description
#' Computes size-dependent natural mortality following Lorenzen (1996, 2022).
#' Supports both weight-based and growth-based formulations.
#'
#' @details
#' Two formulations are available:
#'
#' Weight-based (Lorenzen 1996):
#' \deqn{M(W) = \alpha \cdot W^{\beta}}
#' where \eqn{\alpha \sim N(3.69, 0.502)} and \eqn{\beta \sim N(-0.305, 0.029)}.
#'
#' Growth-based (Lorenzen 2022):
#' \deqn{\ln M = 0.28 - 1.30 \ln(L/L_\infty) + 1.08 \ln(k)}
#' This formulation was calibrated using von Bertalanffy parameters, so
#' \eqn{k} should be the VB-equivalent \eqn{k} when using fits from other
#' growth models.
#'
#' @param age Numeric vector of ages at which to compute mortality.
#' @param Linf Asymptotic length.
#' @param L0 Length at birth.
#' @param k Growth coefficient. For \code{weight_based = FALSE}, should be
#'   VB-equivalent \eqn{k} (use \code{\link{compute_k_vb_equivalent}}).
#' @param lw_fun Function mapping length to weight (required if
#'   \code{weight_based = TRUE}).
#' @param weight_based Logical. If \code{TRUE}, uses weight-based formulation.
#'   If \code{FALSE} (default), uses growth-based formulation.
#' @param growth_model Character. Growth model for length prediction (only
#'   used for weight-based formulation): \code{"vb"}, \code{"gompertz"}, or
#'   \code{"logistic"}.
#' @param sample_params Logical. If \code{TRUE}, samples allometric parameters
#'   from their distributions. If \code{FALSE}, uses mean values.
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
    k,
    lw_fun = NULL,
    weight_based = FALSE,
    growth_model = c("vb", "gompertz", "logistic"),
    sample_params = FALSE
) {

  growth_model <- match.arg(growth_model)

  if (weight_based) {
    # Weight-based formulation
    if (is.null(lw_fun) || !is.function(lw_fun)) {
      stop("Weight-based Lorenzen requires 'lw_fun'.", call. = FALSE)
    }

    # Predict length at age
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

    W_t <- lw_fun(L_t)

    # Allometric parameters
    if (sample_params) {
      alpha <- stats::rnorm(1, 3.69, 0.502)
      beta  <- stats::rnorm(1, -0.305, 0.029)
    } else {
      alpha <- 3.69
      beta  <- -0.305
    }

    M <- alpha * W_t^beta

  } else {
    # Growth-based formulation (Lorenzen 2022)
    # Note: k should be VB-equivalent k for theoretical consistency

    # Predict length at age (always VB for this formulation)
    L_t <- Linf - (Linf - L0) * exp(-k * age)
    L_ratio <- L_t / Linf

    # ln(M) = 0.28 - 1.30*ln(L/Linf) + 1.08*ln(k)
    if (sample_params) {
      intercept <- stats::rnorm(1, 0.28, 0.105)
      coef_L    <- stats::rnorm(1, -1.30, 0.059)
      coef_k    <- stats::rnorm(1, 1.08, 0.082)
    } else {
      intercept <- 0.28
      coef_L    <- -1.30
      coef_k    <- 1.08
    }

    log_M <- intercept + coef_L * log(L_ratio) + coef_k * log(k)
    M <- exp(log_M)
  }

  M
}


# -----------------------------------------------------------------------------
# Mortality Scaling Helper
# -----------------------------------------------------------------------------

#' Scale Mortality Schedule to Target Mean
#'
#' @description
#' Rescales an age-specific mortality schedule so its mean equals a target
#' value derived from empirical relationships (e.g., Hoenig, Then et al.) or
#' survival probability constraints.
#'
#' @details
#' The scaling applies:
#' \deqn{M_{scaled}(t) = M_{raw}(t) \times \frac{M_{target}}{\bar{M}_{raw}}}
#'
#' This preserves the \emph{shape} of the age-specific mortality curve while
#' adjusting its overall level. Scaling is useful because theoretical mortality
#' models often produce absolute levels that don't match empirical observations,
#' but the relative age pattern may still be informative.
#'
#' @param M Numeric vector of instantaneous mortality rates.
#' @param M_target Target mean mortality. Can be a numeric scalar (fixed target),
#'   a function of tmax (e.g., \code{function(tmax) 4.899 * tmax^(-0.916)}), or
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
#' M_raw <- M_chen_watanabe_L0(0:30, Linf = 100, L0 = 25, k = 0.1,
#'                              two_phase = FALSE)
#'
#' # Scale to fixed target
#' M_scaled <- scale_mortality(M_raw, M_target = 0.2)
#'
#' # Scale using Then et al. (2015) relationship
#' then_2015 <- function(tmax) 4.899 * tmax^(-0.916)
#' M_scaled <- scale_mortality(M_raw, M_target = then_2015, tmax = 30)
#'
#' # Scale to survival probability
#' M_scaled <- scale_mortality(M_raw, M_target = NULL, tmax = 30, p = 0.01)
#' }
#'
#' @export
scale_mortality <- function(M, M_target = NULL, tmax = NULL, p = 0.001) {

  # Derive target if not specified
  if (is.null(M_target)) {
    if (is.null(tmax)) {
      stop("tmax required when M_target is NULL.", call. = FALSE)
    }
    # M such that exp(-M * tmax) = p  =>  M = -log(p) / tmax
    M_target <- -log(p) / tmax

  } else if (is.function(M_target)) {
    if (is.null(tmax)) {
      stop("tmax required when M_target is a function.", call. = FALSE)
    }
    M_target <- M_target(tmax)
  }

  # Scale
  M_mean <- mean(M, na.rm = TRUE)
  if (M_mean <= 0) {
    warning("Mean mortality is non-positive; scaling not applied.", call. = FALSE)
    return(M)
  }

  M * (M_target / M_mean)
}


# -----------------------------------------------------------------------------
# Main Stochastic Mortality Function
# -----------------------------------------------------------------------------

#' Stochastic Estimation of Age-Specific Natural Mortality
#'
#' @description
#' Monte Carlo simulation of age-specific natural mortality schedules with
#' full uncertainty propagation from growth model posteriors. Supports
#' Chen-Watanabe, Peterson-Wroblewski, and Lorenzen models with automatic
#' derivation of VB-equivalent \eqn{k} from any growth model fit.
#'
#' @details
#' A key feature of this function is growth-model-agnostic mortality estimation.
#' When a growth fit from \code{\link{fit_bayesian_growth}} is provided, the
#' function extracts biological milestones (Linf, L0, Lmat, tmat)
#' and computes the VB-equivalent k needed for Chen-Watanabe and growth-based
#' Lorenzen models.
#'
#' This allows users to fit whichever growth model (von Bertalanffy, Gompertz,
#' or Logistic) best describes their data, then estimate mortality without
#' theoretical compromise.
#'
#' When \code{growth_fit} is provided, parameters are drawn from the joint
#' posterior distribution, preserving correlations. This yields mortality
#' estimates with appropriate (often narrower) uncertainty bounds compared to
#' independent sampling of each parameter.
#'
#' Three mortality models are available: CW (Chen-Watanabe 1989 with L0
#' parameterization and optional two-phase senescence), PW (Peterson-Wroblewski
#' 1984 weight-based allometric model), and L (Lorenzen 1996/2022 in weight-based
#' or growth-based form).
#'
#' @param method Character. Mortality model: \code{"CW"}, \code{"PW"}, or
#'   \code{"L"}.
#' @param growth_fit Optional \code{CmdStanMCMC} object from
#'   \code{\link{fit_bayesian_growth}}. If provided, parameters are extracted
#'   from the joint posterior.
#' @param maturity_fit Optional \code{CmdStanMCMC} object from
#'   \code{\link{fit_bayesian_maturity}} providing age-at-maturity for
#'   k-based growth fits or two-phase CW.
#' @param sex Integer. Sex code (1 = female, 2 = male) for hierarchical models.
#' @param Linf,L0,k,tmat Alternative to \code{growth_fit}: specify parameters
#'   directly as \code{c(mean, sd)} vectors for independent normal sampling.
#' @param Linf_factor Numeric in (0, 1). Fraction of \eqn{L_\infty} for
#'   \eqn{t_{max}} estimation. Default 0.99.
#' @param age_seq Function or numeric vector defining ages for mortality
#'   calculation. Default \code{function(tmax) seq(0, ceiling(tmax), length.out = 500)}.
#' @param iter Number of Monte Carlo iterations. Default 2000.
#' @param scaled Logical. If \code{TRUE} (default), scales mortality to
#'   \code{M_target} or survival probability \code{p}.
#' @param M_target Target mean mortality. Can be numeric scalar, function of
#'   tmax, or \code{NULL} for survival-probability-based scaling.
#' @param p Survival probability to \eqn{t_{max}} for scaling. Default 0.001.
#' @param two_phase Logical. For CW model, use two-phase senescence?
#'   Default \code{TRUE}.
#' @param late_model Character. Senescence model: \code{"gompertz"} or
#'   \code{"logistic"}. Default \code{"gompertz"}.
#' @param tm_factor,M_mult,smooth_factor Two-phase model parameters.
#' @param lw_fun Length-weight function for PW and weight-based Lorenzen.
#' @param weight_based Logical. For Lorenzen, use weight-based formulation?
#'   Default \code{FALSE}.
#' @param growth_model Character. Growth model type when using manual
#'   parameters: \code{"vb"}, \code{"gompertz"}, or \code{"logistic"}.
#' @param seed Random seed for reproducibility. Default 1234.
#' @param palette Color palette for plot: \code{"synthwave"}, \code{"viridis"},
#'   \code{"okabe"}, \code{"plasma"}, or \code{"inferno"}.
#' @param print_plot Logical. Print plot on completion? Default \code{TRUE}.
#' @param show_progress Logical. Show progress messages? Default \code{TRUE}.
#'
#' @return A list with components: Schedules (data.table of all mortality
#'   schedules with columns set_id, age, M, M_scaled), Parameters (data.table
#'   of sampled life history parameters), Summary (data.table with median and
#'   95 percent CI by age), and Plot (ggplot2 object).
#'
#' @examples
#' \dontrun{
#' # From a Gompertz growth fit (maturity-based parameterization)
#' mort <- get_stochastic_mortality(
#'   method     = "CW",
#'   growth_fit = gomp_fit,  # Any growth model works!
#'   sex        = 1,
#'   iter       = 2000,
#'   scaled     = TRUE,
#'   p          = 0.001
#' )
#'
#' # View plot
#' mort$Plot
#'
#' # Check VB-equivalent k distribution
#' hist(mort$Parameters$k_vb_equiv)
#'
#' # Manual specification for sensitivity analysis
#' mort <- get_stochastic_mortality(
#'   method = "CW",
#'   Linf = c(100, 5),
#'   L0   = c(25, 2),
#'   k    = c(0.1, 0.02),  # VB k or VB-equivalent k
#'   tmat = c(10, 1)
#' )
#' }
#'
#' @import data.table
#' @importFrom stats rnorm quantile median
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line labs theme_bw
#' @export
get_stochastic_mortality <- function(
    method = c("CW", "PW", "L"),
    growth_fit = NULL,
    maturity_fit = NULL,
    sex = NULL,
    Linf = NULL,
    L0 = NULL,
    k = NULL,
    tmat = NULL,
    Linf_factor = 0.99,
    age_seq = function(tmax) seq(0.1, ceiling(tmax), length.out = 500),
    iter = 2000,
    scaled = TRUE,
    M_target = NULL,
    p = 0.001,
    two_phase = TRUE,
    late_model = c("gompertz", "logistic"),
    tm_factor = 2/3,
    M_mult = 2,
    smooth_factor = 1/3,
    lw_fun = NULL,
    weight_based = FALSE,
    growth_model = c("vb", "gompertz", "logistic"),
    seed = 1234,
    palette = c("synthwave", "viridis", "okabe", "plasma", "inferno"),
    print_plot = TRUE,
    show_progress = TRUE
) {

  method <- match.arg(method)
  late_model <- match.arg(late_model)
  growth_model <- match.arg(growth_model)
  palette <- match.arg(palette)

  set.seed(seed)
  use_posterior <- !is.null(growth_fit)

  # -------------------------------------------------------------------------
  # Parameter Extraction / Generation
  # -------------------------------------------------------------------------

  if (use_posterior) {
    if (show_progress) message("Extracting parameters from growth model posterior...")

    # Extract from posterior (works for any growth model)
    params <- extract_growth_parameters(
      growth_fit   = growth_fit,
      maturity_fit = maturity_fit,
      sex          = sex,
      n_draws      = iter,
      seed         = seed
    )

    # For CW and growth-based Lorenzen, we need VB-equivalent k
    if (method == "CW" || (method == "L" && !weight_based)) {
      if (all(is.na(params$k_vb_equiv))) {
        stop(
          "Chen-Watanabe and growth-based Lorenzen require maturity parameters ",
          "(Lmat, tmat) for VB-equivalent k derivation.\n",
          "Use maturity-based growth fit (k_based = FALSE) or provide maturity_fit.",
          call. = FALSE
        )
      }
      k_for_mort <- params$k_vb_equiv
    } else {
      # PW and weight-based Lorenzen can use the native k
      k_for_mort <- params$k
    }

    par_draws <- data.table::data.table(
      set_id     = seq_len(nrow(params)),
      Linf       = params$Linf,
      L0         = params$L0,
      Lmat       = params$Lmat,
      tmat       = params$tmat,
      k_original = params$k,
      k_vb_equiv = params$k_vb_equiv,
      k_for_mort = k_for_mort
    )

  } else {
    # Manual parameter specification
    if (show_progress) message("Generating parameters from specified distributions...")

    if (is.null(Linf) || is.null(L0) || is.null(k)) {
      stop("Must provide growth_fit OR all of: Linf, L0, k", call. = FALSE)
    }
    if ((method == "CW" && two_phase) || (method == "L" && !weight_based)) {
      if (is.null(tmat)) {
        stop("tmat required for CW two-phase or growth-based Lorenzen.", call. = FALSE)
      }
    }

    # Sample from specified distributions
    Linf_draws <- stats::rnorm(iter, Linf[1], Linf[2])
    L0_draws   <- stats::rnorm(iter, L0[1], L0[2])
    k_draws    <- stats::rnorm(iter, k[1], k[2])

    if (!is.null(tmat)) {
      tmat_draws <- stats::rnorm(iter, tmat[1], tmat[2])
    } else {
      tmat_draws <- rep(NA_real_, iter)
    }

    # Ensure biological constraints
    Linf_draws <- pmax(Linf_draws, L0_draws + 1)
    L0_draws   <- pmax(L0_draws, 0.1)
    k_draws    <- pmax(k_draws, 0.001)
    tmat_draws <- pmax(tmat_draws, 0.1)

    par_draws <- data.table::data.table(
      set_id     = seq_len(iter),
      Linf       = Linf_draws,
      L0         = L0_draws,
      Lmat       = rep(NA_real_, iter),  # Not available for manual input
      tmat       = tmat_draws,
      k_original = k_draws,
      k_vb_equiv = k_draws,  # Assume user provides VB-equivalent k
      k_for_mort = k_draws
    )
  }

  # Estimate tmax for each parameter set
  par_draws[, tmax := -log(Linf * (1 - Linf_factor) / (Linf - L0)) / k_for_mort]

  # Remove invalid draws
  valid_mask <- !is.na(par_draws$k_for_mort) &
    par_draws$k_for_mort > 0 &
    par_draws$Linf > par_draws$L0 &
    par_draws$tmax > 0 &
    is.finite(par_draws$tmax)

  n_invalid <- sum(!valid_mask)
  if (n_invalid > 0) {
    if (show_progress) {
      message(sprintf("Removing %d invalid parameter sets (%.1f%%)",
                      n_invalid, 100 * n_invalid / nrow(par_draws)))
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
  ages <- ages[ages > 0]  # Ensure positive ages

  # Initialize storage
  schedules_list <- vector("list", nrow(par_draws))

  for (i in seq_len(nrow(par_draws))) {

    p_i <- par_draws[i]

    # Compute mortality based on method
    M_raw <- switch(
      method,

      "CW" = M_chen_watanabe_L0(
        age          = ages,
        Linf         = p_i$Linf,
        L0           = p_i$L0,
        k            = p_i$k_for_mort,
        tmax         = p_i$tmax,
        Linf_factor  = Linf_factor,
        two_phase    = two_phase,
        tmat         = p_i$tmat,
        late_model   = late_model,
        tm_factor    = tm_factor,
        M_mult       = M_mult,
        smooth_factor = smooth_factor
      ),

      "PW" = {
        if (is.null(lw_fun)) {
          stop("Peterson-Wroblewski requires 'lw_fun'.", call. = FALSE)
        }
        M_peterson_wroblewski(
          age          = ages,
          Linf         = p_i$Linf,
          L0           = p_i$L0,
          k            = p_i$k_original,  # PW uses native k
          lw_fun       = lw_fun,
          growth_model = if (use_posterior) "vb" else growth_model
        )
      },

      "L" = M_lorenzen(
        age          = ages,
        Linf         = p_i$Linf,
        L0           = p_i$L0,
        k            = if (weight_based) p_i$k_original else p_i$k_for_mort,
        lw_fun       = lw_fun,
        weight_based = weight_based,
        growth_model = if (use_posterior) "vb" else growth_model,
        sample_params = TRUE
      )
    )

    # Scale if requested
    if (scaled) {
      M_scaled <- scale_mortality(M_raw, M_target = M_target, tmax = p_i$tmax, p = p)
    } else {
      M_scaled <- M_raw
    }

    schedules_list[[i]] <- data.table::data.table(
      set_id   = p_i$set_id,
      age      = ages,
      M_raw    = M_raw,
      M_scaled = M_scaled
    )
  }

  # Combine all schedules
  schedules <- data.table::rbindlist(schedules_list)

  # -------------------------------------------------------------------------
  # Summary Statistics
  # -------------------------------------------------------------------------

  if (show_progress) message("Computing summary statistics...")

  # Round ages for grouping
  schedules[, age_round := round(age, 1)]

  summary_dt <- schedules[, .(
    M_median = stats::median(M_scaled, na.rm = TRUE),
    M_mean   = mean(M_scaled, na.rm = TRUE),
    M_lower  = stats::quantile(M_scaled, 0.025, na.rm = TRUE),
    M_upper  = stats::quantile(M_scaled, 0.975, na.rm = TRUE)
  ), by = age_round]

  # Parameter summary
  tmax_summary <- par_draws[, .(
    mean  = mean(tmax),
    lower = stats::quantile(tmax, 0.025),
    upper = stats::quantile(tmax, 0.975)
  )]

  # -------------------------------------------------------------------------
  # Plotting
  # -------------------------------------------------------------------------

  if (show_progress) message("Generating plot...")

  # Color palette
  pal <- switch(
    palette,
    "synthwave" = c("#FF6B9D", "#C490D1", "#9B6DFF", "#00D4AA"),
    "viridis"   = viridis::viridis(4),
    "okabe"     = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
    "plasma"    = viridis::plasma(4),
    "inferno"   = viridis::inferno(4)
  )

  fill_color <- pal[1]
  line_color <- pal[3]

  caption <- sprintf(
    "Estimated tmax: %.1f years (95%% CI: %.1f - %.1f) | Method: %s%s",
    tmax_summary$mean, tmax_summary$lower, tmax_summary$upper,
    method,
    if (method == "CW" && two_phase) paste0(" (two-phase, ", late_model, ")") else ""
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
      x        = "Age (years)",
      y        = "Instantaneous Mortality (M)",
      title    = "Age-Specific Natural Mortality",
      subtitle = "Median with 95% credible interval",
      caption  = caption
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

  # -------------------------------------------------------------------------
  # Return
  # -------------------------------------------------------------------------

  if (show_progress) message("Done.")

  list(
    Schedules  = schedules,
    Parameters = par_draws,
    Summary    = summary_dt,
    Plot       = mort_plot
  )
}

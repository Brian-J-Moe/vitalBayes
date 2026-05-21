#' Extract Life History Parameters from Growth Model Posterior
#'
#' @description
#' Extracts posterior draws of \eqn{(L_\infty, L_0, L_{mat}, t_{mat})} from a
#' vitalBayes growth model fit. Works identically for von Bertalanffy, Gompertz,
#' and Logistic models fitted via \code{\link{fit_bayesian_growth}}.
#'
#' @details
#' This function provides a unified interface for extracting the biological
#' parameters common to all growth models. These parameters --- asymptotic
#' length, birth size, and maturity milestones --- represent real biological
#' quantities that exist independently of the mathematical model used to
#' describe growth.
#'
#' For maturity-based growth fits (\code{k_based = FALSE}), \eqn{L_{mat}} and
#' \eqn{t_{mat}} are directly estimated parameters. For k-based fits, these
#' must be supplied separately via \code{maturity_fit}.
#'
#' The function also computes the VB-equivalent \eqn{k} needed for mortality
#' models that were derived under VB assumptions (Chen-Watanabe and growth-based
#' Lorenzen). Three derivation strategies are available, applied in order of
#' preference:
#' \enumerate{
#'   \item \strong{Milestone-based}: If \eqn{L_{mat}} and \eqn{t_{mat}} are
#'     available (from maturity-based fits or a separate \code{maturity_fit}),
#'     \eqn{k_{VB}^{equiv}} is derived by inverting the VB equation through the
#'     biological milestones \eqn{(0, L_0)} and \eqn{(t_{mat}, L_{mat})}.
#'   \item \strong{Derivative matching}: For k-based Gompertz or Logistic fits
#'     without maturity data, \eqn{k_{VB}^{equiv}} is derived by matching
#'     the instantaneous growth rate at birth (\eqn{t = 0}).
#'   \item \strong{Identity}: For VB fits, the native \eqn{k} is itself the
#'     VB growth coefficient.
#' }
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
#'   \item{k}{Growth coefficient (native to the fitted model)}
#'   \item{k_vb_equiv}{VB-equivalent k (always computed; see Details)}
#'   \item{growth_model}{Growth model type: \code{"vb"}, \code{"gompertz"},
#'     or \code{"logistic"}}
#' }
#'
#' The \code{growth_model} column is constant across all rows. The
#' \code{k_vb_equiv} derivation method is stored as an attribute
#' \code{"k_vb_method"}: one of \code{"milestone"}, \code{"derivative"},
#' or \code{"identity"}.
#'
#' @examples
#' \dontrun{
#' # From a Gompertz fit with maturity-based parameterization
#' params <- extract_growth_parameters(gomp_fit, sex = 1, n_draws = 2000)
#'
#' # Check VB-equivalent k distribution
#' hist(params$k_vb_equiv, main = "VB-Equivalent k from Gompertz Fit")
#'
#' # Verify the growth model was auto-detected
#' unique(params$growth_model)
#' #> [1] "gompertz"
#'
#' # Check which derivation method was used
#' attr(params, "k_vb_method")
#' #> [1] "milestone"
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

  # -----------------------------------------------------------------------
  # Auto-detect growth model type from fit attributes
  # -----------------------------------------------------------------------

  growth_model <- attr(growth_fit, "vb_growth_model")
  k_based      <- attr(growth_fit, "vb_k_based")

  # Fallback for fits created before metadata was attached
  if (is.null(growth_model)) {
    growth_model <- "vb"
    message(
      "Growth model type not found in fit attributes (pre-v0.X fit?).\n",
      "Assuming von Bertalanffy. Re-fit with current vitalBayes for ",
      "automatic detection."
    )
  }
  if (is.null(k_based)) {
    # Infer from whether Lmat is in the Stan variables
    available_params <- growth_fit$metadata()$stan_variables
    k_based <- !("Lmat" %in% available_params)
  }

  # -----------------------------------------------------------------------
  # Extract posterior draws
  # -----------------------------------------------------------------------

  available_params <- growth_fit$metadata()$stan_variables

  # Determine if hierarchical (2-sex) model
  Linf_draws <- growth_fit$draws("Linf", format = "matrix")
  is_hierarchical <- ncol(Linf_draws) > 1

  # Set sex index
  if (is.null(sex)) {
    s <- 1L
    if (is_hierarchical) {
      message(
        "Hierarchical model detected but sex not specified. ",
        "Using sex = 1 (female)."
      )
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

  # -----------------------------------------------------------------------
  # Extract maturity parameters if available
  # -----------------------------------------------------------------------

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

  # -----------------------------------------------------------------------
  # Compute VB-equivalent k (three-tier strategy)
  # -----------------------------------------------------------------------

  has_milestones <- has_Lmat && (has_tmat || !is.null(maturity_fit))

  if (has_milestones) {
    # Strategy 1: Milestone-based (preferred)
    k_vb_equiv <- compute_k_vb_equivalent(Linf, L0, Lmat, tmat, warn = FALSE)
    k_vb_method <- "milestone"

  } else if (growth_model == "vb") {
    # Strategy 3: Identity (VB k IS the VB-equivalent k)
    k_vb_equiv <- k
    k_vb_method <- "identity"

  } else {
    # Strategy 2: Derivative matching at birth (fallback for k-based non-VB)
    k_vb_equiv <- .derive_k_vb_derivative(k, Linf, L0, growth_model)
    k_vb_method <- "derivative"

    message(
      sprintf(
        "No maturity milestones available for %s fit. ",
        tools::toTitleCase(growth_model)
      ),
      "VB-equivalent k derived via derivative matching at birth.\n",
      "For more accurate derivation, provide maturity_fit or use ",
      "maturity-based growth parameterization (k_based = FALSE)."
    )
  }

  # -----------------------------------------------------------------------
  # Build output
  # -----------------------------------------------------------------------

  result <- data.table::data.table(
    draw         = seq_len(n_total),
    Linf         = Linf,
    L0           = L0,
    Lmat         = Lmat,
    tmat         = tmat,
    k            = k,
    k_vb_equiv   = k_vb_equiv,
    growth_model = growth_model
  )

  data.table::setnames(result, c("draws", "Linf", "L0",
                                 "Lmat", "tmat", "k",
                                 "k_vb_equiv", "growth_model"))
  # Subsample if requested
  if (!is.null(n_draws) && n_draws < n_total) {
    set.seed(seed)
    idx <- sample(n_total, n_draws, replace = FALSE)
    result <- result[idx]
    result[, draw := seq_len(.N)]
  }

  # Store derivation method as attribute
  attr(result, "k_vb_method") <- k_vb_method

  result
}

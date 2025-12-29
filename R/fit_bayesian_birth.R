# =============================================================================
# vitalBayes Birth Model Fitting Function
# =============================================================================
# Bayesian estimation of length-at-birth (b50) using binomial regression
# with probit link function. Uses instantiate for precompiled Stan models.
# =============================================================================

#' Fit Bayesian Length-at-Birth Model
#'
#' @description
#' Fits a binomial regression with probit link function for estimating the length
#' at which 50% of individuals are expected to transition from embryo to free-swimming
#' status (b50). The model is implemented in Stan via precompiled models using
#' the instantiate package.
#'
#' The probit link is chosen based on a threshold-crossing interpretation where latent
#' developmental readiness (influenced by maternal size, condition, temperature, etc.)
#' is assumed normally distributed. This contrasts with the logit link which assumes
#' log-odds, which are less interpretable in this biological context.
#'
#' @details
#' The linear predictor is parameterized directly in terms of b50:
#' \deqn{\eta_i = \text{slope} \times (\text{length}_i - b_{50})}
#'
#' Parameters are estimated using lognormal priors to ensure positivity.
#'
#' @section Prior Specification:
#' Priors are constructed using a coefficient of variation approach. The prior mean
#' for b50 is automatically computed as the midpoint between the maximum embryo
#' length and minimum free-swimming length (the transition zone). The prior SD
#' is then \code{mean_b50 * cv_b50}.
#'
#' @section Initialization:
#' Initial values are set to the data-derived midpoint estimate for b50 and a
#' reasonable default for slope, ensuring chains start near the high-probability
#' region of the posterior.
#'
#' @param embryo_lts Numeric vector of embryo lengths (internally coded as status = 0).
#' @param free_swimming_lts Numeric vector of free-swimming lengths (internally coded as status = 1).
#' @param mean_b50 Numeric. Prior mean for b50. If \code{NULL} (default), computed
#'   as the midpoint between min(free_swimming_lts) and max(embryo_lts).
#' @param cv_b50 Numeric. Coefficient of variation for b50 prior. Default 0.3.
#' @param mean_slope Numeric. Prior mean for slope on log scale. Default 0.
#' @param sd_slope Numeric. Prior SD for slope on log scale. Default 1.
#' @param parallel Logical. Run chains in parallel? Default \code{TRUE}.
#' @param chains Integer. Number of MCMC chains. Default 4.
#' @param iter_warmup Integer. Warmup iterations per chain. Default 1000.
#' @param iter_sampling Integer. Sampling iterations per chain. Default 1000.
#' @param refresh Integer. Progress update frequency. Set to 0 for no output. Default 500.
#' @param seed Integer. Random seed for reproducibility. Default 1234.
#' @param ... Additional arguments passed to \code{$sample()}.
#'
#' @return A \code{CmdStanMCMC} object containing:
#' \describe{
#'   \item{b50}{Length at 50% birth probability}
#'   \item{slope}{Transition steepness (probit scale)}
#'   \item{b05, b95}{Lengths at 5% and 95% birth probability}
#'   \item{transition_width}{Width of transition zone (b95 - b05)}
#'   \item{log_lik}{Log-likelihood values for LOO-CV}
#'   \item{p_pred}{Predicted probabilities}
#'   \item{status_rep}{Posterior predictive replications}
#'   \item{mean_p_embryo, mean_p_freeswim}{Mean predicted probabilities by group}
#'   \item{prop_correct_rep}{Classification accuracy in replications}
#' }
#'
#' @seealso
#' \code{vignette("fit_bayesian_birth")} for usage examples with gulper shark data.
#'
#' \code{vignette("complete_workflow")} for the full three-stage analysis pipeline.
#'
#' \href{../doc/vitalBayes_stats_explained.html#birth}{Statistical Methods: Birth Size Estimation}
#' for the mathematical derivation.
#'
#' \href{../doc/vitalBayes_stats_explained.html#probit}{Statistical Methods: Probit Link}
#' for justification of the probit over logit link.
#'
#' \code{\link{fit_bayesian_maturity}}, \code{\link{fit_bayesian_growth}} for
#' downstream models that use birth estimates as priors.
#'
#' @examples
#' \dontrun{
#' # Fit birth model with data-derived priors
#' birth_fit <- fit_bayesian_birth(
#'   embryo_lts        = sharks[embryo == TRUE, length],
#'   free_swimming_lts = sharks[embryo == FALSE, length]
#' )
#'
#' # View summary
#' birth_fit$summary(c("b50", "slope", "transition_width"))
#'
#' # Use b50 posterior as L0 prior for growth model
#' growth_fit <- fit_bayesian_growth(
#'   ...,
#'   birth_stanfit = birth_fit
#' )
#' }
#'
#' @import data.table
#' @export
fit_bayesian_birth <- function(
    embryo_lts,
    free_swimming_lts,
    mean_b50      = NULL,
    cv_b50        = 0.3,
    mean_slope    = 0,
    sd_slope      = 1,
    parallel      = TRUE,
    chains        = 4,
    iter_warmup   = 1000,
    iter_sampling = 1000,
    refresh       = 500,
    seed          = 1234,
    ...
) {
  
  # =========================================================================
  # Check Dependencies
  # =========================================================================
  
  if (!requireNamespace("instantiate", quietly = TRUE)) {
    stop(
      "Package 'instantiate' is required for precompiled Stan models.\n",
      "Install via: install.packages('instantiate')",
      call. = FALSE
    )
  }
  
  # =========================================================================
  # Validate Inputs
  # =========================================================================
  
  if (missing(embryo_lts) || missing(free_swimming_lts)) {
    stop("Both 'embryo_lts' and 'free_swimming_lts' are required.", call. = FALSE)
  }
  
  if (!is.numeric(embryo_lts) || !is.numeric(free_swimming_lts)) {
    stop("'embryo_lts' and 'free_swimming_lts' must be numeric vectors.", 
         call. = FALSE)
  }
  
  # Remove NAs
  embryo_lts <- embryo_lts[!is.na(embryo_lts)]
  free_swimming_lts <- free_swimming_lts[!is.na(free_swimming_lts)]
  
  if (length(embryo_lts) == 0 || length(free_swimming_lts) == 0) {
    stop("Both embryo and free-swimming samples must have non-NA values.", 
         call. = FALSE)
  }
  
  # =========================================================================
  # Build Combined Data
  # =========================================================================
  
  lengths <- c(embryo_lts, free_swimming_lts)
  status  <- c(rep(0L, length(embryo_lts)), rep(1L, length(free_swimming_lts)))
  
  newdat <- data.table::data.table(length = lengths, status = status)
  
  n_embryo <- sum(newdat$status == 0)
  n_free <- sum(newdat$status == 1)
  
  message("Data summary: ", n_embryo, " embryos, ", n_free, " free-swimming")
  
  # =========================================================================
  # Data-Driven Prior Center
  # =========================================================================
  
  embryo_max <- max(newdat[status == 0, length])
  free_min <- min(newdat[status == 1, length])
  
  # Check for overlap
  if (embryo_max < free_min) {
    message("Note: No overlap between embryo and free-swimming lengths.")
    message("  Max embryo: ", round(embryo_max, 1), " cm")
    message("  Min free-swimming: ", round(free_min, 1), " cm")
  }
  
  # Midpoint estimate for b50
  midpoint_b50 <- (embryo_max + free_min) / 2
  
  if (is.null(mean_b50)) {
    mean_b50 <- midpoint_b50
    message("Prior for b50 centered at data-derived midpoint: ", round(mean_b50, 2), " cm")
  } else {
    message("Using user-specified prior mean for b50: ", round(mean_b50, 2), " cm")
  }
  
  # Compute SD from CV
  sd_b50 <- mean_b50 * cv_b50
  message("Prior SD for b50 (CV=", cv_b50, "): ", round(sd_b50, 2), " cm")
  
  # Convert to log scale for lognormal prior
  b50_prior <- .natural_to_log_prior(mean_b50, sd_b50)
  
  # =========================================================================
  # Build Stan Data List
  # =========================================================================
  
  stan_data <- list(
    N               = nrow(newdat),
    length          = newdat$length,
    status          = newdat$status,
    prior_b50_mu    = b50_prior$log_mean,
    prior_b50_sigma = b50_prior$log_sd,
    prior_slope_mu  = mean_slope,
    prior_slope_sigma = sd_slope
  )
  
  # =========================================================================
  # Initial Values (data-driven)
  # =========================================================================
  
  init_fun <- function() {
    list(
      log_b50   = log(midpoint_b50),
      log_slope = 0
    )
  }
  
  # =========================================================================
  # Load Precompiled Model
  # =========================================================================
  
  message("\nLoading precompiled Stan model...")
  
  model <- instantiate::stan_package_model(
    name = "length_at_birth",
    package = "vitalBayes"
  )
  
  # =========================================================================
  # Fit Model
  # =========================================================================
  
  message("Fitting model with ", chains, " chains...")
  
  n_cores <- if (parallel) min(chains, parallel::detectCores() - 1) else 1
  
  fit <- model$sample(
    data            = stan_data,
    seed            = seed,
    chains          = chains,
    parallel_chains = n_cores,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    refresh         = refresh,
    init            = init_fun,
    ...
  )
  
  # =========================================================================
  # Print Summary
  # =========================================================================
  
  message("\n", paste(rep("=", 60), collapse = ""))
  message("Birth Model (b50) Posterior Summary")
  message(paste(rep("=", 60), collapse = ""))
  
  message("\nCore Parameters:")
  print(fit$summary(variables = c("b50", "slope")))
  
  message("\nDerived Quantities:")
  print(fit$summary(variables = c("b05", "b95", "transition_width")))
  
  message("\nPosterior Predictive Checks:")
  print(fit$summary(variables = c("mean_p_embryo", "mean_p_freeswim", "prop_correct_rep")))
  
  # =========================================================================
  # Check Convergence
  # =========================================================================
  
  diag <- fit$diagnostic_summary()
  if (any(diag$num_divergent > 0)) {
    warning("Divergent transitions detected. Consider more informative priors.", 
            call. = FALSE)
  }
  
  return(fit)
}

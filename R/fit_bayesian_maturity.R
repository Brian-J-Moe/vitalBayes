# =============================================================================
# vitalBayes Maturity Model Fitting Function
# =============================================================================
# Bayesian estimation of length-at-maturity (L50) and/or age-at-maturity (t50)
# using binomial regression with probit link function.
# Uses instantiate for precompiled Stan models.
# =============================================================================

#' Fit Bayesian Maturity Models
#'
#' @description
#' Fits binomial regression(s) with probit link function for estimating length-at-50%-maturity
#' (L50) and/or age-at-50%-maturity (t50). Supports both single-sex and hierarchical two-sex
#' models with optional partial pooling between sexes.
#'
#' The probit link is chosen for consistency with the birth model and for its threshold-crossing
#' interpretation: latent developmental readiness is assumed normally distributed, and an
#' individual matures when readiness exceeds a threshold. This interpretation is more
#' biologically intuitive than log-odds, which are not standard reporting metrics in
#' elasmobranch research.
#'
#' @details
#' The linear predictor is parameterized directly in terms of x50:
#' \deqn{p_i = \Phi[\beta \times (x_i - x_{50})]}
#'
#' where \eqn{\Phi} is the standard normal CDF, and x represents either length or age.
#'
#' @section Prior Specification:
#' Priors are constructed using a coefficient of variation approach. The prior mean
#' for L50/t50 is automatically computed as the midpoint between the smallest mature
#' individual and largest immature individual. The prior SD is \code{mean * cv}.
#'
#' @section Initialization:
#' Initial values are set to the data-derived midpoint estimates, ensuring chains
#' start near the high-probability region of the posterior.
#'
#' @param maturity Column name in \code{data} or numeric vector of binary maturity
#'   indicators (0 = immature, 1 = mature).
#' @param lt Optional column name or vector of observed lengths (> 0). Required for L50.
#' @param age Optional column name or vector of observed ages (>= 0). Required for t50.
#' @param sex Optional column name or vector for sex. If provided, fits two-sex model.
#'   Supports auto-detection of common coding conventions including: F/M, Female/Male,
#'   1/2 (numeric), and equivalents in Spanish, Portuguese, French, German, Italian,
#'   Japanese, Chinese, and Russian. See Details.
#' @param female Optional. Explicit specification of how females are coded in the data.
#'   If \code{NULL} (default), auto-detection is attempted. Use with \code{male} for
#'   non-standard coding schemes.
#' @param male Optional. Explicit specification of how males are coded in the data.
#'   Must be provided together with \code{female}.
#' @param data A data.frame or data.table containing referenced columns.
#' @param mean_L50,mean_t50 Prior means on natural scale. If \code{NULL} (default),
#'   computed from data as midpoint of transition zone. For two-sex models, can be
#'   length-2 vectors for sex-specific priors.
#' @param cv_L50,cv_t50 Coefficient of variation for priors. Default 0.3.
#' @param mean_slope,sd_slope Prior mean and SD for slope on log scale.
#'   Default \code{mean_slope = 0}, \code{sd_slope = 1}.
#' @param use_pooling Logical. For two-sex models, use partial pooling? Default \code{TRUE}.
#' @param prior_tau Scale for half-normal prior on between-sex SD. Default 0.5.
#' @param parallel Logical. Run chains in parallel? Default \code{TRUE}.
#' @param chains Integer. Number of MCMC chains. Default 4.
#' @param iter_warmup Integer. Warmup iterations per chain. Default 1000.
#' @param iter_sampling Integer. Sampling iterations per chain. Default 1000.
#' @param refresh Integer. Progress update frequency. Default 500.
#' @param seed Integer. Random seed. Default 1234.
#' @param ... Additional arguments passed to \code{$sample()}.
#'
#' @return For single model (length or age only): A \code{CmdStanMCMC} object.
#'   For both models: A named list with \code{$length} and \code{$age} fits.
#'
#' @seealso
#' \code{vignette("fit_bayesian_maturity")} for usage examples with gulper shark data.
#'
#' \code{vignette("partial_pooling")} for detailed explanation of hierarchical
#' modeling for imbalanced sex ratios.
#'
#' \href{../doc/vitalBayes_stats_explained.html#maturity}{Statistical Methods: Maturity Estimation}
#' for mathematical derivation including partial pooling.
#'
#' \href{../doc/vitalBayes_stats_explained.html#probit}{Statistical Methods: Probit Link}
#' for justification of probit over logit.
#'
#' \code{\link{fit_bayesian_birth}}, \code{\link{fit_bayesian_growth}},
#' \code{\link{plot_maturity_ogive}}
#'
#' @examples
#' \dontrun{
#' data(gulper_data)
#'
#' # Single-sex length-at-maturity (females only)
#' L50_fit <- fit_bayesian_maturity(
#'   maturity = "mat",
#'   lt       = "fl",
#'   data     = gulper_data[sex == 1 & embryo == FALSE]
#' )
#'
#' # Two-sex model with partial pooling (integer sex auto-detected)
#' L50_fit_2sex <- fit_bayesian_maturity(
#'   maturity    = "mat",
#'   lt          = "fl",
#'   sex         = "sex",
#'   data        = gulper_data[embryo == FALSE],
#'   use_pooling = TRUE
#' )
#'
#' # Explicit sex coding for non-standard data
#' L50_fit_explicit <- fit_bayesian_maturity(
#'   maturity = "mat",
#'   lt       = "fl",
#'   sex      = "sex",
#'   female   = 1,
#'   male     = 2,
#'   data     = gulper_data[embryo == FALSE]
#' )
#'
#' # Fit both length and age maturity models
#' mat_fits <- fit_bayesian_maturity(
#'   maturity = "mat",
#'   lt       = "fl",
#'   age      = "age1",
#'   sex      = "sex",
#'   data     = gulper_data[embryo == FALSE]
#' )
#' }
#'
#' @import data.table
#' @export
fit_bayesian_maturity <- function(
    maturity,
    lt            = NULL,
    age           = NULL,
    sex           = NULL,
    female        = NULL,
    male          = NULL,
    data          = NULL,
    mean_L50      = NULL,
    mean_t50      = NULL,
    cv_L50        = 0.3,
    cv_t50        = 0.3,
    mean_slope    = 0,
    sd_slope      = 1,
    use_pooling   = TRUE,
    prior_tau     = 0.5,
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
  # Parse Inputs
  # =========================================================================
  
  maturity_vec <- .resolve_or_vector(substitute(maturity), data, "numeric", "maturity")
  
  lt_vec <- if (!missing(lt) && !is.null(substitute(lt))) {
    tryCatch(
      .resolve_or_vector(substitute(lt), data, "numeric", "length"),
      error = function(e) NULL
    )
  } else NULL
  
  age_vec <- if (!missing(age) && !is.null(substitute(age))) {
    tryCatch(
      .resolve_or_vector(substitute(age), data, "numeric", "age"),
      error = function(e) NULL
    )
  } else NULL
  
  sex_vec <- if (!missing(sex) && !is.null(substitute(sex))) {
    # Try numeric first, then character
    tryCatch({
      .resolve_or_vector(substitute(sex), data, "numeric", "sex")
    }, error = function(e) {
      tryCatch({
        .resolve_or_vector(substitute(sex), data, "character", "sex")
      }, error = function(e2) NULL)
    })
  } else NULL
  
  # =========================================================================
  # Validate
  # =========================================================================
  
  if (is.null(lt_vec) && is.null(age_vec)) {
    stop("Must provide at least one of 'lt' or 'age'.", call. = FALSE)
  }
  
  if (!all(maturity_vec %in% c(0, 1, NA))) {
    stop("'maturity' must contain only 0, 1, or NA values.", call. = FALSE)
  }
  
  fit_length <- !is.null(lt_vec)
  fit_age <- !is.null(age_vec)
  is_twosex <- !is.null(sex_vec)
  
  # =========================================================================
  # Build Data Table
  # =========================================================================
  
  newdat <- data.table::data.table(maturity = as.integer(maturity_vec))
  
  if (fit_length) newdat[, length := lt_vec]
  if (fit_age) newdat[, age := age_vec]
  
  # Standardize sex coding using the helper function
  sex_info <- NULL
  if (is_twosex) {
    sex_info <- standardize_sex(sex_vec, female = female, male = male, silent = FALSE)
    newdat[, sex_code := sex_info$sex_int]
  }
  
  newdat <- stats::na.omit(newdat)
  
  if (nrow(newdat) == 0L) {
    stop("No valid observations after removing NAs.", call. = FALSE)
  }
  
  # Report sample sizes
  if (is_twosex) {
    n_by_sex <- newdat[, .N, by = sex_code]
    female_label <- sex_info$labels["1"]
    male_label <- sex_info$labels["2"]
    message("Sample sizes: ",
            "Females (", female_label, ") = ", sex_info$n_female,
            ", Males (", male_label, ") = ", sex_info$n_male)
    
    ratio <- max(n_by_sex$N) / min(n_by_sex$N)
    if (ratio > 2 && use_pooling) {
      message("Note: Imbalanced sex ratio (", round(ratio, 1), 
              ":1). Partial pooling recommended.")
    }
  } else {
    message("Sample size: n = ", nrow(newdat))
  }
  
  n_cores <- if (parallel) min(chains, parallel::detectCores() - 1) else 1
  
  results <- list()
  
  # =========================================================================
  # Length-at-Maturity Model
  # =========================================================================
  
  if (fit_length) {
    
    message("\n", paste(rep("-", 50), collapse = ""))
    message("Fitting Length-at-Maturity Model")
    message(paste(rep("-", 50), collapse = ""))
    
    if (is_twosex) {
      # --- Two-sex length model ---
      
      # Compute sex-specific midpoints
      midpoint_L50 <- newdat[, {
        min_mature <- min(length[maturity == 1], na.rm = TRUE)
        max_immature <- max(length[maturity == 0], na.rm = TRUE)
        (min_mature + max_immature) / 2
      }, by = sex_code][order(sex_code)]$V1
      
      if (is.null(mean_L50)) {
        mean_L50_vec <- midpoint_L50
        message("Prior L50 from midpoints: F = ", round(mean_L50_vec[1], 1), 
                " cm, M = ", round(mean_L50_vec[2], 1), " cm")
      } else {
        mean_L50_vec <- .standardize_scalar(mean_L50, 2)
      }
      
      sd_L50_vec <- mean_L50_vec * cv_L50
      
      # Convert to log scale
      L50_priors <- lapply(1:2, function(s) {
        .natural_to_log_prior(mean_L50_vec[s], sd_L50_vec[s])
      })
      
      stan_data <- list(
        N               = nrow(newdat),
        length          = newdat$length,
        mature          = newdat$maturity,
        sex             = newdat$sex_code,
        use_pooling     = as.integer(use_pooling),
        prior_L50_mu    = sapply(L50_priors, `[[`, "log_mean"),
        prior_L50_sigma = sapply(L50_priors, `[[`, "log_sd"),
        prior_slope_mu  = mean_slope,
        prior_slope_sigma = sd_slope,
        prior_tau_L50   = prior_tau
      )
      
      # Initialization - all params always present now
      init_fun <- function() {
        list(
          mu_L50    = mean(log(midpoint_L50)),
          tau_L50   = 0.1,
          raw_L50   = if (use_pooling) c(0, 0) else log(midpoint_L50),
          log_slope = c(0, 0)
        )
      }
      
      model <- instantiate::stan_package_model(
        name = "length_at_maturity_twosex",
        package = "vitalBayes"
      )
      
    } else {
      # --- Single-sex length model ---
      
      # Compute midpoint
      min_mature <- newdat[maturity == 1, min(length, na.rm = TRUE)]
      max_immature <- newdat[maturity == 0, max(length, na.rm = TRUE)]
      midpoint_L50 <- (min_mature + max_immature) / 2
      
      if (is.null(mean_L50)) {
        mean_L50_use <- midpoint_L50
        message("Prior L50 from midpoint: ", round(mean_L50_use, 1), " cm")
      } else {
        mean_L50_use <- mean_L50[1]
      }
      
      sd_L50_use <- mean_L50_use * cv_L50
      L50_prior <- .natural_to_log_prior(mean_L50_use, sd_L50_use)
      
      stan_data <- list(
        N               = nrow(newdat),
        length          = newdat$length,
        mature          = newdat$maturity,
        prior_L50_mu    = L50_prior$log_mean,
        prior_L50_sigma = L50_prior$log_sd,
        prior_slope_mu  = mean_slope,
        prior_slope_sigma = sd_slope
      )
      
      init_fun <- function() {
        list(
          log_L50   = log(midpoint_L50),
          log_slope = 0
        )
      }
      
      model <- instantiate::stan_package_model(
        name = "length_at_maturity_single",
        package = "vitalBayes"
      )
    }
    
    message("Fitting model...")
    fit_L <- model$sample(
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
    
    .print_maturity_summary(fit_L, "length", is_twosex, use_pooling)
    
    results$length <- fit_L
  }
  
  # =========================================================================
  # Age-at-Maturity Model
  # =========================================================================
  
  if (fit_age) {
    
    message("\n", paste(rep("-", 50), collapse = ""))
    message("Fitting Age-at-Maturity Model")
    message(paste(rep("-", 50), collapse = ""))
    
    if (is_twosex) {
      # --- Two-sex age model ---
      
      # Compute sex-specific midpoints
      midpoint_t50 <- newdat[, {
        min_mature <- min(age[maturity == 1], na.rm = TRUE)
        max_immature <- max(age[maturity == 0], na.rm = TRUE)
        (min_mature + max_immature) / 2
      }, by = sex_code][order(sex_code)]$V1
      
      if (is.null(mean_t50)) {
        mean_t50_vec <- midpoint_t50
        message("Prior t50 from midpoints: F = ", round(mean_t50_vec[1], 1), 
                " yrs, M = ", round(mean_t50_vec[2], 1), " yrs")
      } else {
        mean_t50_vec <- .standardize_scalar(mean_t50, 2)
      }
      
      sd_t50_vec <- mean_t50_vec * cv_t50
      
      t50_priors <- lapply(1:2, function(s) {
        .natural_to_log_prior(mean_t50_vec[s], sd_t50_vec[s])
      })
      
      stan_data <- list(
        N               = nrow(newdat),
        age             = newdat$age,
        mature          = newdat$maturity,
        sex             = newdat$sex_code,
        use_pooling     = as.integer(use_pooling),
        prior_t50_mu    = sapply(t50_priors, `[[`, "log_mean"),
        prior_t50_sigma = sapply(t50_priors, `[[`, "log_sd"),
        prior_slope_mu  = mean_slope,
        prior_slope_sigma = sd_slope,
        prior_tau_t50   = prior_tau
      )
      
      # Initialization - all params always present now
      init_fun <- function() {
        list(
          mu_t50    = mean(log(midpoint_t50)),
          tau_t50   = 0.1,
          raw_t50   = if (use_pooling) c(0, 0) else log(midpoint_t50),
          log_slope = c(0, 0)
        )
      }
      
      model <- instantiate::stan_package_model(
        name = "age_at_maturity_twosex",
        package = "vitalBayes"
      )
      
    } else {
      # --- Single-sex age model ---
      
      min_mature <- newdat[maturity == 1, min(age, na.rm = TRUE)]
      max_immature <- newdat[maturity == 0, max(age, na.rm = TRUE)]
      midpoint_t50 <- (min_mature + max_immature) / 2
      
      if (is.null(mean_t50)) {
        mean_t50_use <- midpoint_t50
        message("Prior t50 from midpoint: ", round(mean_t50_use, 1), " yrs")
      } else {
        mean_t50_use <- mean_t50[1]
      }
      
      sd_t50_use <- mean_t50_use * cv_t50
      t50_prior <- .natural_to_log_prior(mean_t50_use, sd_t50_use)
      
      stan_data <- list(
        N               = nrow(newdat),
        age             = newdat$age,
        mature          = newdat$maturity,
        prior_t50_mu    = t50_prior$log_mean,
        prior_t50_sigma = t50_prior$log_sd,
        prior_slope_mu  = mean_slope,
        prior_slope_sigma = sd_slope
      )
      
      init_fun <- function() {
        list(
          log_t50   = log(midpoint_t50),
          log_slope = 0
        )
      }
      
      model <- instantiate::stan_package_model(
        name = "age_at_maturity_single",
        package = "vitalBayes"
      )
    }
    
    message("Fitting model...")
    fit_t <- model$sample(
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
    
    .print_maturity_summary(fit_t, "age", is_twosex, use_pooling)
    
    results$age <- fit_t
  }
  
  # =========================================================================
  # Return
  # =========================================================================
  
  if (length(results) == 1L) {
    return(results[[1L]])
  } else {
    return(results)
  }
}

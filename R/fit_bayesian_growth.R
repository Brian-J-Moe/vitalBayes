# =============================================================================
# vitalBayes Growth Model Fitting Function
# =============================================================================
# Bayesian estimation of individual growth parameters using von Bertalanffy,
# Gompertz, or Logistic growth models. Supports single-sex and two-sex models
# with optional partial pooling. Uses instantiate for precompiled Stan models.
# =============================================================================

#' Fit Bayesian Individual Growth Models
#'
#' @description
#' Fits von Bertalanffy, Gompertz, or Logistic growth models using precompiled
#' Stan models via the instantiate package. Supports either traditional k-based
#' or maturity-based parameterization, and both single-sex and hierarchical
#' two-sex models with optional partial pooling between sexes.
#'
#' @details
#' ## Growth Model Equations
#'
#' **von Bertalanffy (model = "v"):**
#' \deqn{L(t) = L_\infty - (L_\infty - L_0) e^{-kt}}
#'
#' **Gompertz (model = "g"):**
#' \deqn{L(t) = L_\infty \exp\left(-\ln\left(\frac{L_\infty}{L_0}\right) e^{-kt}\right)}
#'
#' **Logistic (model = "l"):**
#' \deqn{L(t) = \frac{L_\infty}{1 + \left(\frac{L_\infty}{L_0} - 1\right) e^{-kt}}}
#'
#' ## Prior Specification
#'
#' Priors are constructed using coefficient of variation (CV) arguments:
#'
#' **Linf:** Prior mean defaults to \code{1.05 * Lmax}, where Lmax is the maximum
#' observed length (including incomplete cases without age data). The prior SD
#' is \code{mean * CV_Linf}. Users can override Lmax via the \code{Lmax} argument.
#'
#' **k:** Prior mean is estimated from the data by solving the growth equation
#' for k at each observation and averaging. Prior SD is \code{mean * CV_k}.
#'
#' **L0, Lmat, tmat:** Can be informed by prior birth/maturity model fits, or
#' specified manually.
#'
#' ## Initialization
#'
#' Initial values are set to data-derived estimates. For k, the Gulland-Holt
#' linearization method is used to provide a reasonable starting value.
#'
#' @param lt Column name in \code{data} or numeric vector of observed lengths (> 0).
#' @param age Column name in \code{data} or numeric vector of observed ages (>= 0).
#' @param sex Optional column name or vector for sex. If provided, fits two-sex model.
#'   Supports auto-detection of common coding conventions including: F/M, Female/Male,
#'   1/2 (numeric), and equivalents in multiple languages. See \code{\link{fit_bayesian_maturity}}
#'   for full list.
#' @param female Optional. Explicit specification of how females are coded in the data.
#'   If \code{NULL} (default), auto-detection is attempted.
#' @param male Optional. Explicit specification of how males are coded in the data.
#'   Must be provided together with \code{female}.
#' @param data A data.frame or data.table containing referenced columns. May include
#'   incomplete cases (length without age) which are used only to determine Lmax.
#'
#' @param model Character. Growth model type: \code{"v"} (von Bertalanffy),
#'   \code{"g"} (Gompertz), or \code{"l"} (Logistic). Default \code{"v"}.
#' @param k_based Logical. If \code{TRUE}, uses traditional k-based parameterization.
#'   If \code{FALSE} (default), derives k from maturity parameters (Lmat, tmat).
#'
#' @param birth_stanfit Optional CmdStanMCMC from \code{\link{fit_bayesian_birth}}.
#'   Used to set informative priors for L0.
#' @param length.mature_stanfit Optional CmdStanMCMC from \code{\link{fit_bayesian_maturity}}
#'   with L50 parameter. Used to set priors for Lmat.
#' @param age.mature_stanfit Optional CmdStanMCMC from \code{\link{fit_bayesian_maturity}}
#'   with t50 parameter. Used to set priors for tmat.
#'
#' @param Lmax Numeric. Maximum observed length by sex for setting Linf prior and
#'   lower bound. If \code{NULL} (default), computed from all length data in \code{data}
#'   (including rows with missing age). Can be a single value (applied to both sexes)
#'   or length-2 vector for sex-specific values c(female, male).
#' @param Linf_multiplier Numeric. Multiplier for Lmax to get Linf prior mean.
#'   Default 1.05 (i.e., 5% larger than observed max).
#' @param CV_Linf Coefficient of variation for Linf prior. Default 0.2.
#' @param CV_k Coefficient of variation for k prior. Default 0.5.
#' @param CV_L0 Coefficient of variation for L0 prior (if not from birth fit). Default 0.3.
#' @param CV_Lmat Coefficient of variation for Lmat prior (if not from maturity fit). Default 0.2.
#' @param CV_tmat Coefficient of variation for tmat prior (if not from maturity fit). Default 0.3.
#'
#' @param prior_L0 Manual prior for L0 on natural scale: \code{c(mean, sd)}.
#'   Ignored if \code{birth_stanfit} provided.
#' @param prior_Lmat Manual prior for Lmat on natural scale. Only for maturity-based.
#' @param prior_tmat Manual prior for tmat on natural scale. Only for maturity-based.
#' @param prior_k Manual prior for k on natural scale. Only for k-based models.
#'
#' @param use_pooling Logical. For two-sex models, use partial pooling on location
#'   parameters? Default \code{TRUE}. Slope parameters are never pooled.
#' @param prior_tau Scale for half-normal priors on between-sex SD. Default 0.2.
#'
#' @param robust Logical. Use Student-t observation model? Default \code{FALSE}.
#' @param loc_sig,scale_sig Location and scale for Cauchy prior on sigma.
#'   Defaults: \code{loc_sig = 0}, \code{scale_sig = 1}.
#'
#' @param parallel Logical. Run chains in parallel? Default \code{TRUE}.
#' @param chains Integer. Number of MCMC chains. Default 4.
#' @param iter_warmup Integer. Warmup iterations per chain. Default 1000.
#' @param iter_sampling Integer. Sampling iterations per chain. Default 1000.
#' @param refresh Integer. Progress update frequency. Default 500.
#' @param seed Integer. Random seed. Default 1234.
#' @param ... Additional arguments passed to \code{$sample()}.
#'
#' @return A \code{CmdStanMCMC} object.
#'
#' @seealso
#' \code{vignette("fit_bayesian_growth")} for usage examples with gulper shark data.
#'
#' \code{vignette("partial_pooling")} for hierarchical modeling of imbalanced sex ratios.
#'
#' \href{../doc/vitalBayes_stats_explained.html#growth}{Statistical Methods: Growth Models}
#' for equations and model comparison.
#'
#' \href{../doc/vitalBayes_stats_explained.html#maturity-param}{Statistical Methods: Maturity-Based Parameterization}
#' for the derivation of k from maturity parameters.
#'
#' \href{../doc/vitalBayes_stats_explained.html#linf}{Statistical Methods: L-infinity Constraint}
#' for why Linf must exceed maximum observed length.
#'
#' \code{\link{fit_bayesian_birth}}, \code{\link{fit_bayesian_maturity}},
#' \code{\link{plot_growth_curve}}, \code{\link{compare_growth_models}}
#'
#' @examples
#' \dontrun{
#' # Workflow: birth -> maturity -> growth
#'
#' # Step 1: Fit birth model
#' birth_fit <- fit_bayesian_birth(
#'   embryo_lts = sharks[embryo == TRUE, length],
#'   free_swimming_lts = sharks[embryo == FALSE, length]
#' )
#'
#' # Step 2: Fit maturity models
#' mat_fits <- fit_bayesian_maturity(
#'   maturity = mat, lt = length, age = age, sex = sex,
#'   data = sharks[embryo == FALSE]
#' )
#'
#' # Step 3: Fit growth model (maturity-based, pooled)
#' growth_fit <- fit_bayesian_growth(
#'   lt = length, age = age, sex = sex,
#'   data = sharks,  # Can include incomplete cases for Lmax
#'   model = "v",
#'   k_based = FALSE,
#'   birth_stanfit = birth_fit,
#'   length.mature_stanfit = mat_fits$length,
#'   age.mature_stanfit = mat_fits$age,
#'   use_pooling = TRUE,
#'   robust = TRUE
#' )
#'
#' # Alternative: Specify Lmax if you know it exceeds observed data
#' growth_fit <- fit_bayesian_growth(
#'   lt = length, age = age, sex = sex,
#'   data = sharks[!is.na(age)],  # Complete cases only
#'   Lmax = c(150, 140),  # Known max: females, males
#'   ...
#' )
#' }
#'
#' @import data.table
#' @export
fit_bayesian_growth <- function(
    lt,
    age,
    sex           = NULL,
    female        = NULL,
    male          = NULL,
    data          = NULL,
    model         = c("v", "g", "l"),
    k_based       = FALSE,

    # Prior sources from previous fits
    birth_stanfit         = NULL,
    length.mature_stanfit = NULL,
    age.mature_stanfit    = NULL,

    # Lmax specification
    Lmax              = NULL,
    Linf_multiplier   = 1.05,

    # CV arguments for priors
    CV_Linf       = 0.2,
    CV_k          = 0.5,
    CV_L0         = 0.3,
    CV_Lmat       = 0.2,
    CV_tmat       = 0.3,

    # Manual priors (override defaults)
    prior_L0      = NULL,
    prior_Lmat    = NULL,
    prior_tmat    = NULL,
    prior_k       = NULL,

    # Pooling options
    use_pooling   = TRUE,
    prior_tau     = 0.2,

    # Observation error
    robust        = FALSE,
    loc_sig       = 0,
    scale_sig     = 1,

    # Sampling settings
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
  # Input Validation
  # =========================================================================

  if (missing(lt) || missing(age)) {
    stop("Both 'lt' (length) and 'age' are required for growth models.",
         call. = FALSE)
  }

  model <- match.arg(model)
  which_model <- switch(model, v = 1L, g = 2L, l = 3L)
  model_name <- switch(model,
                       v = "von Bertalanffy",
                       g = "Gompertz",
                       l = "Logistic")

  # Check maturity-based requirements
  if (!k_based) {
    if (is.null(length.mature_stanfit) && is.null(prior_Lmat)) {
      stop(
        "Maturity-based models require either 'length.mature_stanfit' or 'prior_Lmat'.\n",
        "Alternatively, set k_based = TRUE for traditional parameterization.",
        call. = FALSE
      )
    }
    if (is.null(age.mature_stanfit) && is.null(prior_tmat)) {
      stop(
        "Maturity-based models require either 'age.mature_stanfit' or 'prior_tmat'.\n",
        "Alternatively, set k_based = TRUE for traditional parameterization.",
        call. = FALSE
      )
    }
  }

  # =========================================================================
  # Parse Input Data
  # =========================================================================

  lt_vec <- .resolve_or_vector(substitute(lt), data, "numeric", "length")
  age_vec <- .resolve_or_vector(substitute(age), data, "numeric", "age")

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

  is_twosex <- !is.null(sex_vec)
  n_sex <- if (is_twosex) 2L else 1L

  # Standardize sex coding if two-sex model
  sex_info <- NULL
  if (is_twosex) {
    sex_info <- standardize_sex(sex_vec, female = female, male = male, silent = FALSE)
  }

  # =========================================================================
  # Determine Lmax (including incomplete cases)
  # =========================================================================

  # Build data including all observations for Lmax determination
  all_dat <- data.table::data.table(length = lt_vec)
  if (is_twosex) {
    all_dat[, sex_code := sex_info$sex_int]
  }

  if (is.null(Lmax)) {
    if (is_twosex) {
      # Compute Lmax by sex from ALL data (including incomplete cases)
      Lmax_vec <- all_dat[, .(max_lt = max(length, na.rm = TRUE)),
                          by = sex_code][order(sex_code)]$max_lt
      message("Lmax from data (including incomplete cases): ",
              "F = ", round(Lmax_vec[1], 1), " cm, ",
              "M = ", round(Lmax_vec[2], 1), " cm")
    } else {
      Lmax_vec <- max(lt_vec, na.rm = TRUE)
      message("Lmax from data (including incomplete cases): ", round(Lmax_vec, 1), " cm")
    }
  } else {
    Lmax_vec <- .standardize_scalar(Lmax, n_sex)
    message("Using user-specified Lmax: ", paste(round(Lmax_vec, 1), collapse = ", "), " cm")
  }

  # =========================================================================
  # Filter to Complete Cases for Fitting
  # =========================================================================

  newdat <- data.table::data.table(
    length = lt_vec,
    age = age_vec,
    row_id = seq_along(lt_vec)
  )

  if (is_twosex) {
    newdat[, sex_code := sex_info$sex_int]
  }

  # Filter to complete cases (both length and age present)
  n_before <- nrow(newdat)
  newdat <- newdat[!is.na(length) & !is.na(age)]
  n_after <- nrow(newdat)

  if (n_before > n_after) {
    message("Filtered to complete cases: ", n_before, " -> ", n_after, " observations")
  }

  if (nrow(newdat) == 0L) {
    stop("No complete observations (both length and age) after filtering.", call. = FALSE)
  }

  # Validate data quality
  .validate_growth_data(newdat$length, newdat$age,
                        if (is_twosex) newdat$sex else NULL)

  # Sample size checks
  .check_sample_size(newdat, if (is_twosex) newdat$sex else NULL, k_based)

  # =========================================================================
  # Report Model Configuration
  # =========================================================================

  message("\n", paste(rep("=", 60), collapse = ""))
  message("Fitting ", model_name, " Growth Model")
  message(paste(rep("=", 60), collapse = ""))
  message("Parameterization: ", if (k_based) "k-based" else "maturity-based")
  message("Sex structure: ", if (is_twosex) "two-sex" else "single-sex")
  if (is_twosex) {
    message("Pooling: ", if (use_pooling) "partial pooling" else "independent")
  }
  message("Observation error: ", if (robust) "Student-t (robust)" else "lognormal")
  message("Sample size: ", nrow(newdat))

  # =========================================================================
  # Build Priors
  # =========================================================================

  # --- Linf priors ---
  mean_Linf_nat <- Lmax_vec * Linf_multiplier
  sd_Linf_nat <- mean_Linf_nat * CV_Linf
  message("Linf prior: mean = ", paste(round(mean_Linf_nat, 1), collapse = ", "),
          " cm (", Linf_multiplier, " x Lmax, CV = ", CV_Linf, ")")

  Linf_priors <- lapply(1:n_sex, function(s) {
    .natural_to_log_prior(mean_Linf_nat[s], sd_Linf_nat[s])
  })

  # Linf lower bounds = observed Lmax
  Linf_lower <- Lmax_vec

  # --- L0 priors ---
  if (!is.null(birth_stanfit)) {
    # Extract log-scale parameters directly from posterior samples
    L0_prior_info <- .extract_birth_prior(birth_stanfit, multiplier = 1.0)
    mean_L0_nat <- rep(L0_prior_info$natural["mean"], n_sex)
    # Use log-scale values directly (no double conversion!)
    L0_priors <- lapply(1:n_sex, function(s) {
      list(log_mean = L0_prior_info$log["mean"], log_sd = L0_prior_info$log["sd"])
    })
    message("L0 prior from birth fit: mean = ", round(mean_L0_nat[1], 1), " cm")
  } else if (!is.null(prior_L0)) {
    prior_L0 <- .standardize_prior(prior_L0, n_sex)
    mean_L0_nat <- prior_L0[, 1]
    sd_L0_nat <- prior_L0[, 2]
    L0_priors <- lapply(1:n_sex, function(s) {
      .natural_to_log_prior(mean_L0_nat[s], sd_L0_nat[s])
    })
  } else {
    mean_L0_nat <- rep(min(newdat$length) * 0.8, n_sex)
    sd_L0_nat <- mean_L0_nat * CV_L0
    message("L0 prior (default): mean = ", round(mean_L0_nat[1], 1), " cm (CV = ", CV_L0, ")")
    L0_priors <- lapply(1:n_sex, function(s) {
      .natural_to_log_prior(mean_L0_nat[s], sd_L0_nat[s])
    })
  }

  # --- k priors (for k-based models) ---
  if (k_based) {
    if (!is.null(prior_k)) {
      prior_k <- .standardize_prior(prior_k, n_sex)
      mean_k_nat <- prior_k[, 1]
      sd_k_nat <- prior_k[, 2]
    } else {
      # Estimate k from data
      mean_k_nat <- sapply(1:n_sex, function(s) {
        if (is_twosex) {
          dat_s <- newdat[sex_code == s]
        } else {
          dat_s <- newdat
        }
        .estimate_k_from_data(dat_s$length, dat_s$age,
                              mean_Linf_nat[s], mean_L0_nat[s], model)
      })
      sd_k_nat <- mean_k_nat * CV_k
      message("k prior (data-derived): mean = ",
              paste(round(mean_k_nat, 3), collapse = ", "), " (CV = ", CV_k, ")")
    }

    k_priors <- lapply(1:n_sex, function(s) {
      .natural_to_log_prior(mean_k_nat[s], sd_k_nat[s])
    })
  }

  # --- Lmat priors (for maturity-based models) ---
  if (!k_based) {
    if (!is.null(length.mature_stanfit)) {
      # Extract log-scale parameters directly from posterior samples
      Lmat_info <- .extract_maturity_prior(length.mature_stanfit, "L50", multiplier = 1.0)
      mean_Lmat_nat <- sapply(Lmat_info, function(x) x$natural["mean"])
      # Use log-scale values directly (no double conversion!)
      Lmat_priors <- lapply(Lmat_info, function(x) {
        list(log_mean = x$log["mean"], log_sd = x$log["sd"])
      })
      if (length(mean_Lmat_nat) < n_sex) {
        mean_Lmat_nat <- rep(mean_Lmat_nat, n_sex)
        Lmat_priors <- rep(Lmat_priors, n_sex)
      }
      message("Lmat prior from maturity fit: mean = ",
              paste(round(mean_Lmat_nat, 1), collapse = ", "), " cm")
    } else {
      prior_Lmat <- .standardize_prior(prior_Lmat, n_sex)
      mean_Lmat_nat <- prior_Lmat[, 1]
      sd_Lmat_nat <- prior_Lmat[, 2]
      Lmat_priors <- lapply(1:n_sex, function(s) {
        .natural_to_log_prior(mean_Lmat_nat[s], sd_Lmat_nat[s])
      })
    }

    # --- tmat priors ---
    if (!is.null(age.mature_stanfit)) {
      # Extract log-scale parameters directly from posterior samples
      tmat_info <- .extract_maturity_prior(age.mature_stanfit, "t50", multiplier = 1.0)
      mean_tmat_nat <- sapply(tmat_info, function(x) x$natural["mean"])
      # Use log-scale values directly (no double conversion!)
      tmat_priors <- lapply(tmat_info, function(x) {
        list(log_mean = x$log["mean"], log_sd = x$log["sd"])
      })
      if (length(mean_tmat_nat) < n_sex) {
        mean_tmat_nat <- rep(mean_tmat_nat, n_sex)
        tmat_priors <- rep(tmat_priors, n_sex)
      }
      message("tmat prior from maturity fit: mean = ",
              paste(round(mean_tmat_nat, 1), collapse = ", "), " yrs")
    } else {
      prior_tmat <- .standardize_prior(prior_tmat, n_sex)
      mean_tmat_nat <- prior_tmat[, 1]
      sd_tmat_nat <- prior_tmat[, 2]
      tmat_priors <- lapply(1:n_sex, function(s) {
        .natural_to_log_prior(mean_tmat_nat[s], sd_tmat_nat[s])
      })
    }
  }

  # =========================================================================
  # Build Stan Data List
  # =========================================================================

  if (k_based) {
    # ------ k-based models ------

    if (is_twosex) {
      stan_data <- list(
        N           = nrow(newdat),
        length      = newdat$length,
        age         = newdat$age,
        sex         = newdat$sex_code,
        row_id      = newdat$row_id,
        which_model = which_model,
        robust      = as.integer(robust),
        use_pooling = as.integer(use_pooling),

        prior_Linf_mu    = sapply(Linf_priors, `[[`, "log_mean"),
        prior_Linf_sigma = sapply(Linf_priors, `[[`, "log_sd"),
        Linf_lower       = Linf_lower,

        prior_L0_mu    = sapply(L0_priors, `[[`, "log_mean"),
        prior_L0_sigma = sapply(L0_priors, `[[`, "log_sd"),

        prior_k_mu    = sapply(k_priors, `[[`, "log_mean"),
        prior_k_sigma = sapply(k_priors, `[[`, "log_sd"),

        prior_tau_Linf = prior_tau,
        prior_tau_L0   = prior_tau,
        prior_tau_k    = prior_tau,

        loc_sig   = loc_sig,
        scale_sig = scale_sig
      )

      stan_model_name <- "growth_twosex_k"

    } else {
      stan_data <- list(
        N           = nrow(newdat),
        length      = newdat$length,
        age         = newdat$age,
        row_id      = newdat$row_id,
        which_model = which_model,
        robust      = as.integer(robust),

        prior_Linf_mu    = Linf_priors[[1]]$log_mean,
        prior_Linf_sigma = Linf_priors[[1]]$log_sd,
        Linf_lower       = Linf_lower[1],

        prior_L0_mu    = L0_priors[[1]]$log_mean,
        prior_L0_sigma = L0_priors[[1]]$log_sd,

        prior_k_mu    = k_priors[[1]]$log_mean,
        prior_k_sigma = k_priors[[1]]$log_sd,

        loc_sig   = loc_sig,
        scale_sig = scale_sig
      )

      stan_model_name <- "growth_single_k"
    }

  } else {
    # ------ Maturity-based models ------

    if (is_twosex) {
      stan_data <- list(
        N           = nrow(newdat),
        length      = newdat$length,
        age         = newdat$age,
        sex         = newdat$sex_code,
        row_id      = newdat$row_id,
        which_model = which_model,
        robust      = as.integer(robust),
        use_pooling = as.integer(use_pooling),

        prior_Linf_mu    = sapply(Linf_priors, `[[`, "log_mean"),
        prior_Linf_sigma = sapply(Linf_priors, `[[`, "log_sd"),
        Linf_lower       = Linf_lower,

        prior_L0_mu    = sapply(L0_priors, `[[`, "log_mean"),
        prior_L0_sigma = sapply(L0_priors, `[[`, "log_sd"),

        prior_Lmat_mu    = sapply(Lmat_priors, `[[`, "log_mean"),
        prior_Lmat_sigma = sapply(Lmat_priors, `[[`, "log_sd"),

        prior_tmat_mu    = sapply(tmat_priors, `[[`, "log_mean"),
        prior_tmat_sigma = sapply(tmat_priors, `[[`, "log_sd"),

        prior_tau_Linf = prior_tau,
        prior_tau_L0   = prior_tau,
        prior_tau_Lmat = prior_tau,
        prior_tau_tmat = prior_tau,

        loc_sig   = loc_sig,
        scale_sig = scale_sig
      )

      stan_model_name <- "growth_twosex_maturity"

    } else {
      stan_data <- list(
        N           = nrow(newdat),
        length      = newdat$length,
        age         = newdat$age,
        row_id      = newdat$row_id,
        which_model = which_model,
        robust      = as.integer(robust),

        prior_Linf_mu    = Linf_priors[[1]]$log_mean,
        prior_Linf_sigma = Linf_priors[[1]]$log_sd,
        Linf_lower       = Linf_lower[1],

        prior_L0_mu    = L0_priors[[1]]$log_mean,
        prior_L0_sigma = L0_priors[[1]]$log_sd,

        prior_Lmat_mu    = Lmat_priors[[1]]$log_mean,
        prior_Lmat_sigma = Lmat_priors[[1]]$log_sd,

        prior_tmat_mu    = tmat_priors[[1]]$log_mean,
        prior_tmat_sigma = tmat_priors[[1]]$log_sd,

        loc_sig   = loc_sig,
        scale_sig = scale_sig
      )

      stan_model_name <- "growth_single_maturity"
    }
  }

  # =========================================================================
  # Build Initialization Function
  # =========================================================================

  # Use Gulland-Holt for k initialization
  if (is_twosex) {
    k_init <- sapply(1:2, function(s) {
      dat_s <- newdat[sex_code == s]
      .gulland_holt_k(dat_s$length, dat_s$age)
    })
  } else {
    k_init <- .gulland_holt_k(newdat$length, newdat$age)
  }

  init_fun <- function() {
    base_init <- if (k_based) {
      if (is_twosex) {
        list(
          mu_Linf  = mean(log(mean_Linf_nat)),
          mu_L0    = mean(log(mean_L0_nat)),
          mu_k     = mean(log(k_init)),
          tau_Linf = 0.1,
          tau_L0   = 0.1,
          tau_k    = 0.1,
          raw_Linf = if (use_pooling) c(0, 0) else log(mean_Linf_nat),
          raw_L0   = if (use_pooling) c(0, 0) else log(mean_L0_nat),
          raw_k = if (use_pooling) c(0, 0) else log(as.numeric(k_init)),
          sigma    = c(0.1, 0.1)
        )
      } else {
        list(
          log_Linf = log(mean_Linf_nat[1]),
          log_L0   = log(mean_L0_nat[1]),
          log_k    = log(k_init),
          sigma    = 0.1
        )
      }
    } else {
      # Maturity-based
      if (is_twosex) {

        .init_noncentered_from_priors <- function(prior_mu, tau_init_floor = 0.25) {
          mu0  <- mean(prior_mu)
          sdmu <- stats::sd(prior_mu)
          tau0 <- max(tau_init_floor, sdmu)
          raw0 <- (prior_mu - mu0) / tau0
          list(mu0 = mu0, tau0 = tau0, raw0 = raw0)
        }

        if(use_pooling) {
          Lmat_init <- .init_noncentered_from_priors(stan_data$prior_Lmat_mu)
          tmat_init <- .init_noncentered_from_priors(stan_data$prior_tmat_mu)
        }

        list(
          mu_Linf  = mean(log(mean_Linf_nat)),
          mu_L0    = mean(log(mean_L0_nat)),
          mu_Lmat = if (use_pooling) Lmat_init$mu0 else mean(log(mean_Lmat_nat)),
          mu_tmat = if (use_pooling) tmat_init$mu0 else mean(log(mean_tmat_nat)),
          tau_Linf = 0.1,
          tau_L0   = 0.1,
          tau_Lmat = if (use_pooling) Lmat_init$tau0 else 0.1,
          tau_tmat = if (use_pooling) tmat_init$tau0 else 0.1,
          raw_Linf = if (use_pooling) c(0, 0) else log(mean_Linf_nat),
          raw_L0   = if (use_pooling) c(0, 0) else log(mean_L0_nat),
          raw_Lmat = if (use_pooling) as.numeric(Lmat_init$raw0) else log(mean_Lmat_nat),
          raw_tmat = if (use_pooling) as.numeric(tmat_init$raw0) else log(mean_tmat_nat),
          sigma    = c(0.1, 0.1)
        )
      } else {
        list(
          log_Linf = log(mean_Linf_nat[1]),
          log_L0   = log(mean_L0_nat[1]),
          log_Lmat = log(mean_Lmat_nat[1]),
          log_tmat = log(mean_tmat_nat[1]),
          sigma    = 0.1
        )
      }
    }

    # Add nu_raw initialization for robust models
    if (robust) {
      base_init$nu_raw <- 10  # Moderate tails as starting point
    }

    base_init
  }


  # =========================================================================
  # Load Precompiled Model and Fit
  # =========================================================================

  message("\nLoading precompiled Stan model: ", stan_model_name)

  model_obj <- instantiate::stan_package_model(
    name = stan_model_name,
    package = "vitalBayes"
  )

  n_cores <- if (parallel) min(chains, parallel::detectCores() - 1) else 1

  message("Fitting model with ", chains, " chains...")

  fit <- model_obj$sample(
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

  .print_growth_summary(fit, k_based, is_twosex, use_pooling, robust)

  # =========================================================================
  # Check Convergence
  # =========================================================================

  diag <- fit$diagnostic_summary()
  if (any(diag$num_divergent > 0)) {
    warning(
      "Divergent transitions detected (", sum(diag$num_divergent), " total).\n",
      "Consider: more informative priors, or increasing adapt_delta.",
      call. = FALSE
    )
  }

  return(fit)
}





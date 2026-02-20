# =============================================================================
# vitalBayes Mortality Estimation Functions (Updated)
# =============================================================================
#
# This module implements natural mortality estimation using biological
# milestones (Linf, L0, Lmat, tmat) with model-dependent G(t).
#
# UPDATED: Now accepts vitalBayes stanfit objects OR manual c(mean, sd)
# specifications. When using stanfit objects, posterior correlations are
# preserved. When maturity parameters come from separate fits, user can
# specify assumed correlation to maintain biological plausibility.
#
# Core mortality formula:
#   M(t) = M_inf / G(t)
#
# where:
#   M_inf = (1/tmat) * ln[(Linf - L0)/(Linf - Lmat)]  (VB-derived, unified anchor)
#   G(t) = L(t)/Linf  (computed using native growth model trajectory)
#
# =============================================================================


# =============================================================================
# SECTION 1: CORE HELPER FUNCTIONS
# =============================================================================

#' Compute Native Growth Coefficient for a Specific Model
#'
#' @description
#' Derives the growth coefficient \eqn{k} for a specified growth model from
#' biological milestones \eqn{(L_\infty, L_0, L_{mat}, t_{mat})}.
#'
#' @param Linf Numeric vector. Asymptotic length.
#' @param L0 Numeric vector. Length at birth.
#' @param Lmat Numeric vector. Length at maturity.
#' @param tmat Numeric vector. Age at maturity.
#' @param growth_model Character. Growth model: \code{"vb"}, \code{"gompertz"},
#'   or \code{"logistic"}.
#' @param warn Logical. Warn on invalid values?
#'
#' @return Numeric vector of native k values. Invalid values returned as NA.
#'
#' @export
compute_k_native <- function(Linf, L0, Lmat, tmat,
                             growth_model = c("vb", "gompertz", "logistic"),
                             warn = TRUE) {

  growth_model <- match.arg(growth_model)

  n <- max(length(Linf), length(L0), length(Lmat), length(tmat))
  Linf <- rep_len(Linf, n)
  L0   <- rep_len(L0, n)
  Lmat <- rep_len(Lmat, n)
  tmat <- rep_len(tmat, n)

  k <- switch(
    growth_model,
    "vb" = {
      numerator   <- Linf - L0
      denominator <- Linf - Lmat
      invalid <- numerator <= 0 | denominator <= 0 | tmat <= 0
      k_out <- (1 / tmat) * log(numerator / denominator)
      k_out[invalid] <- NA_real_
      k_out
    },
    "gompertz" = {
      ratio_mat <- Linf / Lmat
      ratio_0   <- Linf / L0
      invalid <- ratio_mat <= 1 | ratio_0 <= 1 | tmat <= 0
      k_out <- -(1 / tmat) * log(log(ratio_mat) / log(ratio_0))
      k_out[invalid] <- NA_real_
      k_out
    },
    "logistic" = {
      term_mat <- Linf / Lmat - 1
      term_0   <- Linf / L0 - 1
      invalid <- term_mat <= 0 | term_0 <= 0 | tmat <= 0
      k_out <- -(1 / tmat) * log(term_mat / term_0)
      k_out[invalid] <- NA_real_
      k_out
    }
  )

  if (warn && any(is.na(k))) {
    n_invalid <- sum(is.na(k))
    warning(sprintf("%d of %d values produced invalid k.", n_invalid, n), call. = FALSE)
  }

  k
}


#' Compute Asymptotic Mortality Rate (M_inf)
#'
#' @description
#' Computes the asymptotic mortality rate \eqn{M_\infty} from biological
#' milestones using the VB formula. This serves as the unified mortality
#' scaling factor in the Chen-Watanabe framework.
#'
#' @param Linf Numeric vector. Asymptotic length.
#' @param L0 Numeric vector. Length at birth.
#' @param Lmat Numeric vector. Length at maturity.
#' @param tmat Numeric vector. Age at maturity.
#' @param warn Logical. Warn on invalid values?
#'
#' @return Numeric vector of asymptotic mortality rates.
#'
#' @export
compute_Minf <- function(Linf, L0, Lmat, tmat, warn = TRUE) {

  n <- max(length(Linf), length(L0), length(Lmat), length(tmat))
  Linf <- rep_len(Linf, n)
  L0   <- rep_len(L0, n)
  Lmat <- rep_len(Lmat, n)
  tmat <- rep_len(tmat, n)

  numerator   <- Linf - L0
  denominator <- Linf - Lmat

  invalid <- numerator <= 0 | denominator <= 0 | tmat <= 0 |
    is.na(numerator) | is.na(denominator) | is.na(tmat)

  if (warn && any(invalid, na.rm = TRUE)) {
    n_invalid <- sum(invalid, na.rm = TRUE)
    warning(sprintf("%d of %d draws produced invalid M_inf.", n_invalid, n), call. = FALSE)
  }

  Minf <- (1 / tmat) * log(numerator / denominator)
  Minf[invalid] <- NA_real_

  Minf
}


#' Compute Relative Size G(t) = L(t)/Linf
#'
#' @description
#' Computes the relative size at specified ages using a given growth model.
#'
#' @param age Numeric vector. Ages at which to compute G(t).
#' @param Linf Numeric. Asymptotic length.
#' @param L0 Numeric. Length at birth.
#' @param k Numeric. Native growth coefficient.
#' @param growth_model Character. Growth model.
#'
#' @return Numeric vector of G(t) values bounded in (0, 1].
#'
#' @export
compute_G <- function(age, Linf, L0, k,
                      growth_model = c("vb", "gompertz", "logistic")) {

  growth_model <- match.arg(growth_model)

  G <- switch(
    growth_model,
    "vb" = 1 - ((Linf - L0) / Linf) * exp(-k * age),
    "gompertz" = {
      r0 <- log(Linf / L0)
      exp(-r0 * exp(-k * age))
    },
    "logistic" = {
      c <- Linf / L0 - 1
      1 / (1 + c * exp(-k * age))
    }
  )

  pmax(G, 1e-6)
}


#' Compute Length at Age L(t)
#'
#' @param age Numeric vector. Ages.
#' @param Linf Numeric. Asymptotic length.
#' @param L0 Numeric. Length at birth.
#' @param k Numeric. Native growth coefficient.
#' @param growth_model Character. Growth model.
#'
#' @return Numeric vector of predicted lengths.
#'
#' @export
compute_L <- function(age, Linf, L0, k,
                      growth_model = c("vb", "gompertz", "logistic")) {

  growth_model <- match.arg(growth_model)

  switch(
    growth_model,
    "vb" = Linf - (Linf - L0) * exp(-k * age),
    "gompertz" = {
      r0 <- log(Linf / L0)
      Linf * exp(-r0 * exp(-k * age))
    },
    "logistic" = {
      c <- Linf / L0 - 1
      Linf / (1 + c * exp(-k * age))
    }
  )
}


#' Compute Maximum Age from Growth Parameters
#'
#' @param Linf Numeric. Asymptotic length.
#' @param L0 Numeric. Length at birth.
#' @param k Numeric. Native growth coefficient.
#' @param growth_model Character. Growth model.
#' @param Linf_factor Numeric in (0,1). Fraction of Linf (default 0.99).
#'
#' @return Numeric tmax value.
#'
#' @export
compute_tmax <- function(Linf, L0, k,
                         growth_model = c("vb", "gompertz", "logistic"),
                         Linf_factor = 0.99) {

  growth_model <- match.arg(growth_model)

  if (Linf_factor <= 0 || Linf_factor >= 1) {
    stop("Linf_factor must be in (0, 1).", call. = FALSE)
  }

  switch(
    growth_model,
    "vb" = -log(Linf * (1 - Linf_factor) / (Linf - L0)) / k,
    "gompertz" = {
      r0 <- log(Linf / L0)
      -log(-log(Linf_factor) / r0) / k
    },
    "logistic" = {
      c <- Linf / L0 - 1
      -log((1 / Linf_factor - 1) / c) / k
    }
  )
}


# =============================================================================
# SECTION 2: PARAMETER EXTRACTION
# =============================================================================


#' Check if Object is a CmdStanMCMC Fit
#'
#' @param x Object to check.
#' @return Logical.
#' @keywords internal
is_stanfit <- function(x) {
  if (is.null(x)) return(FALSE)
  inherits(x, "CmdStanMCMC")
}


#' Extract Posterior Draws from a Stanfit Object
#'
#' @description
#' Extracts posterior draws for a parameter, handling both single-sex and
#' two-sex models.
#'
#' @param fit CmdStanMCMC object.
#' @param param Character. Parameter name.
#' @param sex Integer or NULL. Sex index for two-sex models.
#'
#' @return Numeric vector of draws.
#' @keywords internal
extract_param_draws <- function(fit, param, sex = NULL) {

  if (!param %in% fit$metadata()$stan_variables) {
    return(NULL)
  }

  draws <- fit$draws(param, format = "matrix")

  if (is.null(sex) || ncol(draws) == 1) {
    return(as.vector(draws[, 1]))
  }

  if (sex > ncol(draws)) {
    warning(sprintf("Sex index %d exceeds available columns for %s. Using column 1.",
                    sex, param), call. = FALSE)
    return(as.vector(draws[, 1]))
  }

  as.vector(draws[, sex])
}


#' Sample from Bivariate Normal with Specified Correlation
#'
#' @description
#' Generates correlated samples from two normal distributions.
#'
#' @param n Number of samples.
#' @param mu1,mu2 Means.
#' @param sd1,sd2 Standard deviations.
#' @param rho Correlation coefficient.
#' @param seed Random seed.
#'
#' @return List with components x1 and x2.
#' @keywords internal
sample_bivariate_normal <- function(n, mu1, sd1, mu2, sd2, rho, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Generate independent standard normals
  z1 <- stats::rnorm(n)
  z2 <- stats::rnorm(n)

  # Apply Cholesky factorization for correlation
  x1 <- mu1 + sd1 * z1
  x2 <- mu2 + sd2 * (rho * z1 + sqrt(1 - rho^2) * z2)

  list(x1 = x1, x2 = x2)
}


#' Extract Life History Parameters from vitalBayes Fits or Manual Specification
#'
#' @description
#' Extracts life history parameters from a combination of vitalBayes stanfit
#' objects and/or manual \code{c(mean, sd)} specifications. When stanfit objects
#' are provided, posterior correlations are preserved. When maturity parameters
#' come from separate fits, optional correlation can be specified.
#'
#' @param growth_fit CmdStanMCMC object from \code{fit_bayesian_growth}, or NULL.
#' @param birth_fit CmdStanMCMC object from \code{fit_bayesian_birth}, or NULL.
#' @param length_maturity_fit CmdStanMCMC object from \code{fit_bayesian_maturity}
#'   for length-at-maturity, or NULL.
#' @param age_maturity_fit CmdStanMCMC object from \code{fit_bayesian_maturity}
#'   for age-at-maturity, or NULL.
#' @param Linf,L0,Lmat,tmat Manual specification as \code{c(mean, sd)}, or NULL.
#' @param maturity_cor Correlation between Lmat and tmat when both are from
#'   manual specification or separate fits. Default 0.5 (positive correlation
#'   reflecting that larger individuals tend to mature later). Set to 0 for
#'   independent sampling, or NA to attempt estimation from data.
#' @param sex Integer. Sex code (1 = female, 2 = male) for extracting from
#'   two-sex models. Default 1.
#' @param iter Number of parameter sets to generate. Default 2000.
#' @param seed Random seed. Default 1234.
#' @param show_progress Logical. Show messages? Default TRUE.
#'
#' @return A data.table with columns: set_id, Linf, L0, Lmat, tmat, plus
#'   attribute "sources" indicating parameter origins.
#'
#' @keywords internal
extract_lh_params <- function(
    growth_fit = NULL,
    birth_fit = NULL,
    length_maturity_fit = NULL,
    age_maturity_fit = NULL,
    Linf = NULL,
    L0 = NULL,
    Lmat = NULL,
    tmat = NULL,
    maturity_cor = 0.5,
    sex = 1L,
    iter = 2000,
    seed = 1234,
    show_progress = TRUE
) {

  set.seed(seed)
  sex <- as.integer(sex)

  # ---------------------------------------------------------------------------
  # Determine sources and check availability
  # ---------------------------------------------------------------------------

  sources <- list()
  has_growth <- is_stanfit(growth_fit)
  has_birth <- is_stanfit(birth_fit)
  has_Lmat_fit <- is_stanfit(length_maturity_fit)
  has_tmat_fit <- is_stanfit(age_maturity_fit)

  # Check what's available from growth_fit
  gf_has_Linf <- has_growth && "Linf" %in% growth_fit$metadata()$stan_variables
  gf_has_L0 <- has_growth && "L0" %in% growth_fit$metadata()$stan_variables
  gf_has_Lmat <- has_growth && "Lmat" %in% growth_fit$metadata()$stan_variables
  gf_has_tmat <- has_growth && "tmat" %in% growth_fit$metadata()$stan_variables

  # Determine source for each parameter
  # Priority: growth_fit > specialized fit > manual

  # Linf
  if (gf_has_Linf) {
    sources$Linf <- "growth_fit"
  } else if (!is.null(Linf) && length(Linf) == 2) {
    sources$Linf <- "manual"
  } else {
    stop("Linf must be provided via growth_fit or as c(mean, sd).", call. = FALSE)
  }

  # L0
  if (gf_has_L0) {
    sources$L0 <- "growth_fit"
  } else if (has_birth && "b50" %in% birth_fit$metadata()$stan_variables) {
    sources$L0 <- "birth_fit"
  } else if (!is.null(L0) && length(L0) == 2) {
    sources$L0 <- "manual"
  } else {
    stop("L0 must be provided via growth_fit, birth_fit, or as c(mean, sd).", call. = FALSE)
  }

  # Lmat
  if (gf_has_Lmat) {
    sources$Lmat <- "growth_fit"
  } else if (has_Lmat_fit && "L50" %in% length_maturity_fit$metadata()$stan_variables) {
    sources$Lmat <- "length_maturity_fit"
  } else if (!is.null(Lmat) && length(Lmat) == 2) {
    sources$Lmat <- "manual"
  } else {
    stop("Lmat must be provided via growth_fit, length_maturity_fit, or as c(mean, sd).", call. = FALSE)
  }

  # tmat
  if (gf_has_tmat) {
    sources$tmat <- "growth_fit"
  } else if (has_tmat_fit && "t50" %in% age_maturity_fit$metadata()$stan_variables) {
    sources$tmat <- "age_maturity_fit"
  } else if (!is.null(tmat) && length(tmat) == 2) {
    sources$tmat <- "manual"
  } else {
    stop("tmat must be provided via growth_fit, age_maturity_fit, or as c(mean, sd).", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Report sources
  # ---------------------------------------------------------------------------

  if (show_progress) {
    posterior_params <- names(sources)[sources != "manual"]
    manual_params <- names(sources)[sources == "manual"]

    if (length(posterior_params) > 0) {
      message(sprintf("Parameters from posterior draws: %s",
                      paste(posterior_params, collapse = ", ")))
    }
    if (length(manual_params) > 0) {
      message(sprintf("Parameters from manual specification: %s",
                      paste(manual_params, collapse = ", ")))
    }
  }

  # ---------------------------------------------------------------------------
  # Extract draws with correlation preservation
  # ---------------------------------------------------------------------------

  # Case 1: All four parameters from growth_fit (full correlation preservation)
  if (all(sapply(sources, function(x) x == "growth_fit"))) {

    if (show_progress) message("All parameters from growth_fit - correlations fully preserved.")

    all_draws <- growth_fit$draws(format = "draws_df")
    n_draws <- nrow(all_draws)

    # Sample indices
    if (n_draws >= iter) {
      idx <- sample.int(n_draws, iter, replace = FALSE)
    } else {
      idx <- sample.int(n_draws, iter, replace = TRUE)
    }

    # Extract with sex handling
    get_col <- function(param) {
      col_1 <- paste0(param, "[1]")
      col_2 <- paste0(param, "[2]")
      if (sex == 2 && col_2 %in% names(all_draws)) {
        return(all_draws[[col_2]][idx])
      } else if (col_1 %in% names(all_draws)) {
        return(all_draws[[col_1]][idx])
      } else if (param %in% names(all_draws)) {
        return(all_draws[[param]][idx])
      } else {
        stop(sprintf("Could not find %s in growth_fit.", param), call. = FALSE)
      }
    }

    par_draws <- data.table::data.table(
      set_id = seq_len(iter),
      Linf   = get_col("Linf"),
      L0     = get_col("L0"),
      Lmat   = get_col("Lmat"),
      tmat   = get_col("tmat")
    )

    attr(par_draws, "sources") <- sources
    attr(par_draws, "correlation_status") <- "full"

    return(.enforce_constraints(par_draws))
  }

  # Case 2: Mixed sources - extract/generate each parameter
  draws_list <- list()

  # --- Linf ---
  if (sources$Linf == "growth_fit") {
    draws_list$Linf <- extract_param_draws(growth_fit, "Linf", sex)
  } else {
    draws_list$Linf <- stats::rnorm(iter, Linf[1], Linf[2])
  }

  # --- L0 ---
  if (sources$L0 == "growth_fit") {
    draws_list$L0 <- extract_param_draws(growth_fit, "L0", sex)
  } else if (sources$L0 == "birth_fit") {
    draws_list$L0 <- extract_param_draws(birth_fit, "b50", NULL)
  } else {
    draws_list$L0 <- stats::rnorm(iter, L0[1], L0[2])
  }

  # --- Lmat and tmat (handle correlation) ---

  # Check if Lmat and tmat come from the same fit (correlation preserved)
  lmat_tmat_same_source <- sources$Lmat == sources$tmat &&
    sources$Lmat %in% c("growth_fit", "length_maturity_fit", "age_maturity_fit")

  # Special case: both from separate maturity fits - check if they can be aligned
  lmat_tmat_separate_fits <- sources$Lmat == "length_maturity_fit" &&
    sources$tmat == "age_maturity_fit"

  if (lmat_tmat_same_source && sources$Lmat == "growth_fit") {
    # Already handled above in full correlation case, but for partial cases:
    draws_list$Lmat <- extract_param_draws(growth_fit, "Lmat", sex)
    draws_list$tmat <- extract_param_draws(growth_fit, "tmat", sex)

    if (show_progress) message("Lmat and tmat from same growth_fit - correlation preserved.")

  } else if (lmat_tmat_separate_fits) {
    # Both from maturity fits - sample with assumed correlation

    # Get marginal posteriors
    Lmat_draws_raw <- extract_param_draws(length_maturity_fit, "L50", sex)
    tmat_draws_raw <- extract_param_draws(age_maturity_fit, "t50", sex)

    if (is.na(maturity_cor)) {
      # Attempt to estimate correlation from raw draws (if same length)
      if (length(Lmat_draws_raw) == length(tmat_draws_raw)) {
        maturity_cor <- stats::cor(Lmat_draws_raw, tmat_draws_raw)
        if (show_progress) {
          message(sprintf("Estimated Lmat-tmat correlation from aligned draws: %.3f", maturity_cor))
        }
      } else {
        maturity_cor <- 0.5
        if (show_progress) {
          message("Could not estimate correlation; using default rho = 0.5")
        }
      }
    }

    # Generate correlated samples matching marginal statistics
    mu_Lmat <- mean(Lmat_draws_raw)
    sd_Lmat <- stats::sd(Lmat_draws_raw)
    mu_tmat <- mean(tmat_draws_raw)
    sd_tmat <- stats::sd(tmat_draws_raw)

    corr_samples <- sample_bivariate_normal(
      n = iter,
      mu1 = mu_Lmat, sd1 = sd_Lmat,
      mu2 = mu_tmat, sd2 = sd_tmat,
      rho = maturity_cor,
      seed = seed + 1
    )

    draws_list$Lmat <- corr_samples$x1
    draws_list$tmat <- corr_samples$x2

    if (show_progress) {
      message(sprintf("Lmat and tmat from separate fits - using specified correlation (rho = %.2f).",
                      maturity_cor))
    }

  } else {
    # Handle Lmat and tmat independently or with specified correlation

    # Extract/generate Lmat
    if (sources$Lmat == "growth_fit") {
      draws_list$Lmat <- extract_param_draws(growth_fit, "Lmat", sex)
    } else if (sources$Lmat == "length_maturity_fit") {
      draws_list$Lmat <- extract_param_draws(length_maturity_fit, "L50", sex)
    } else {
      # Manual - will handle correlation below
      draws_list$Lmat <- NULL
    }

    # Extract/generate tmat
    if (sources$tmat == "growth_fit") {
      draws_list$tmat <- extract_param_draws(growth_fit, "tmat", sex)
    } else if (sources$tmat == "age_maturity_fit") {
      draws_list$tmat <- extract_param_draws(age_maturity_fit, "t50", sex)
    } else {
      # Manual - will handle correlation below
      draws_list$tmat <- NULL
    }

    # If both are manual, generate with correlation
    if (sources$Lmat == "manual" && sources$tmat == "manual") {

      if (abs(maturity_cor) > 0.001) {
        corr_samples <- sample_bivariate_normal(
          n = iter,
          mu1 = Lmat[1], sd1 = Lmat[2],
          mu2 = tmat[1], sd2 = tmat[2],
          rho = maturity_cor,
          seed = seed + 1
        )
        draws_list$Lmat <- corr_samples$x1
        draws_list$tmat <- corr_samples$x2

        if (show_progress) {
          message(sprintf("Lmat and tmat manual - sampled with correlation (rho = %.2f).", maturity_cor))
        }
      } else {
        draws_list$Lmat <- stats::rnorm(iter, Lmat[1], Lmat[2])
        draws_list$tmat <- stats::rnorm(iter, tmat[1], tmat[2])

        if (show_progress) {
          message("Lmat and tmat manual - sampled independently (rho = 0).")
        }
      }

    } else {
      # One from posterior, one manual - sample independently
      if (is.null(draws_list$Lmat)) {
        draws_list$Lmat <- stats::rnorm(iter, Lmat[1], Lmat[2])
      }
      if (is.null(draws_list$tmat)) {
        draws_list$tmat <- stats::rnorm(iter, tmat[1], tmat[2])
      }

      if (show_progress) {
        message("Mixed Lmat/tmat sources - correlation not preserved.")
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Align lengths and build output
  # ---------------------------------------------------------------------------

  # Ensure all draws have same length
  target_n <- iter
  for (nm in names(draws_list)) {
    n_i <- length(draws_list[[nm]])
    if (n_i > target_n) {
      idx <- sample.int(n_i, target_n, replace = FALSE)
      draws_list[[nm]] <- draws_list[[nm]][idx]
    } else if (n_i < target_n) {
      idx <- sample.int(n_i, target_n, replace = TRUE)
      draws_list[[nm]] <- draws_list[[nm]][idx]
    }
  }

  par_draws <- data.table::data.table(
    set_id = seq_len(iter),
    Linf   = draws_list$Linf,
    L0     = draws_list$L0,
    Lmat   = draws_list$Lmat,
    tmat   = draws_list$tmat
  )

  attr(par_draws, "sources") <- sources
  attr(par_draws, "correlation_status") <- "partial"
  attr(par_draws, "maturity_cor") <- maturity_cor

  .enforce_constraints(par_draws)
}


#' Enforce Biological Constraints on Parameter Draws
#'
#' @param par_draws data.table with Linf, L0, Lmat, tmat columns.
#' @return data.table with constraints enforced.
#' @keywords internal
.enforce_constraints <- function(par_draws) {
  par_draws[, Linf := pmax(Linf, Lmat + 1)]
  par_draws[, Lmat := pmax(Lmat, L0 + 1)]
  par_draws[, L0 := pmax(L0, 0.1)]
  par_draws[, tmat := pmax(tmat, 0.1)]
  par_draws
}


# =============================================================================
# SECTION 3: MORTALITY MODELS
# =============================================================================

#' Chen-Watanabe Natural Mortality (Model-Dependent)
#'
#' @description
#' Computes age-specific natural mortality using a generalized Chen-Watanabe
#' framework where G(t) is derived from the native growth model trajectory.
#'
#' @param age Numeric vector of ages.
#' @param Linf Asymptotic length.
#' @param L0 Length at birth.
#' @param Lmat Length at maturity.
#' @param tmat Age at maturity.
#' @param growth_model Character. Growth model for G(t).
#' @param tmax Maximum age (computed if NULL).
#' @param Linf_factor Fraction of Linf for tmax estimation.
#' @param two_phase Use two-phase senescence model?
#' @param late_model Senescence model: \code{"gompertz"} or \code{"logistic"}.
#' @param tm_factor Transition age as fraction of tmat.
#' @param M_mult Mortality multiplier for senescence.
#' @param smooth_factor Transition smoothness.
#'
#' @return Numeric vector of instantaneous mortality rates.
#'
#' @export
M_chen_watanabe <- function(
    age,
    Linf,
    L0,
    Lmat,
    tmat,
    growth_model = c("vb", "gompertz", "logistic"),
    tmax = NULL,
    Linf_factor = 0.99,
    two_phase = FALSE,
    late_model = c("gompertz", "logistic"),
    tm_factor = 2/3,
    M_mult = 2,
    smooth_factor = 1/3
) {

  growth_model <- match.arg(growth_model)
  late_model   <- match.arg(late_model)

  if (L0 >= Linf) stop("L0 must be less than Linf.", call. = FALSE)
  if (Lmat >= Linf) stop("Lmat must be less than Linf.", call. = FALSE)
  if (L0 >= Lmat) stop("L0 must be less than Lmat.", call. = FALSE)
  if (tmat <= 0) stop("tmat must be positive.", call. = FALSE)

  Minf <- compute_Minf(Linf, L0, Lmat, tmat, warn = FALSE)
  k_native <- compute_k_native(Linf, L0, Lmat, tmat, growth_model, warn = FALSE)

  if (is.na(Minf) || Minf <= 0) stop("Invalid M_inf.", call. = FALSE)
  if (is.na(k_native) || k_native <= 0) stop("Invalid native k.", call. = FALSE)

  if (is.null(tmax)) {
    tmax <- compute_tmax(Linf, L0, k_native, growth_model, Linf_factor)
  }

  G_t <- compute_G(age, Linf, L0, k_native, growth_model)
  M_cw <- Minf / G_t

  if (!two_phase) return(M_cw)

  # Two-phase extension
  tm <- tm_factor * tmat
  G_at_tm <- compute_G(tm, Linf, L0, k_native, growth_model)
  M_early_at_tm <- Minf / G_at_tm

  if (late_model == "gompertz") {
    M_s <- M_mult * M_early_at_tm
    r <- log(M_early_at_tm / M_s) / (tm - tmax)
    M_late <- function(t) M_s * exp(r * (t - tmax))
  } else {
    K <- 2 * M_mult * M_early_at_tm
    r <- -log(K / M_early_at_tm - 1) / (tm - tmax)
    M_late <- function(t) K / (1 + exp(-r * (t - tmax)))
  }

  smooth_width <- smooth_factor * abs(tmax - tm)
  weight <- 1 / (1 + exp(-2 * (age - tm) / smooth_width))

  M <- (1 - weight) * M_cw + weight * M_late(age)
  M[age < tm * 0.5] <- M_cw[age < tm * 0.5]
  pmax(M, 0)
}


#' Peterson-Wroblewski Natural Mortality Model
#'
#' @param age Numeric vector of ages.
#' @param Linf Asymptotic length.
#' @param L0 Length at birth.
#' @param Lmat Length at maturity.
#' @param tmat Age at maturity.
#' @param lw_fun Length-weight function.
#' @param growth_model Character. Growth model for L(t).
#'
#' @return Numeric vector of instantaneous mortality rates.
#'
#' @export
M_peterson_wroblewski <- function(
    age,
    Linf,
    L0,
    Lmat,
    tmat,
    lw_fun,
    growth_model = c("vb", "gompertz", "logistic")
) {

  growth_model <- match.arg(growth_model)

  if (is.null(lw_fun) || !is.function(lw_fun)) {
    stop("PW model requires a length-weight function 'lw_fun(L)'.", call. = FALSE)
  }

  k_native <- compute_k_native(Linf, L0, Lmat, tmat, growth_model, warn = FALSE)
  L_t <- compute_L(age, Linf, L0, k_native, growth_model)
  W_t <- lw_fun(L_t)
  W_t <- pmax(W_t, 1)

  1.92 * W_t^(-0.25)
}


#' Lorenzen Natural Mortality Model
#'
#' @param age Numeric vector of ages.
#' @param Linf Asymptotic length.
#' @param L0 Length at birth.
#' @param Lmat Length at maturity.
#' @param tmat Age at maturity.
#' @param lw_fun Length-weight function (required if weight_based = TRUE).
#' @param weight_based Use weight-based formulation?
#' @param growth_model Growth model for L(t).
#' @param sample_params Sample parameters from their distributions?
#'
#' @return Numeric vector of instantaneous mortality rates.
#'
#' @export
M_lorenzen <- function(
    age,
    Linf,
    L0,
    Lmat,
    tmat,
    lw_fun = NULL,
    weight_based = FALSE,
    growth_model = c("vb", "gompertz", "logistic"),
    sample_params = TRUE
) {

  growth_model <- match.arg(growth_model)

  k_native <- compute_k_native(Linf, L0, Lmat, tmat, growth_model, warn = FALSE)
  Minf <- compute_Minf(Linf, L0, Lmat, tmat, warn = FALSE)

  if (weight_based) {
    if (is.null(lw_fun) || !is.function(lw_fun)) {
      stop("Weight-based Lorenzen requires 'lw_fun'.", call. = FALSE)
    }

    L_t <- compute_L(age, Linf, L0, k_native, growth_model)
    W_t <- lw_fun(L_t)
    W_t <- pmax(W_t, 1)

    if (sample_params) {
      alpha <- stats::rnorm(1, 3.69, 0.502)
      beta  <- stats::rnorm(1, -0.305, 0.029)
    } else {
      alpha <- 3.69
      beta  <- -0.305
    }
    M <- alpha * W_t^beta

  } else {
    L_t <- compute_L(age, Linf, L0, k_native, growth_model)
    L_ratio <- L_t / Linf

    if (sample_params) {
      intercept <- stats::rnorm(1, 0.28, 0.105)
      coef_L    <- stats::rnorm(1, -1.30, 0.059)
      coef_k    <- stats::rnorm(1, 1.08, 0.082)
    } else {
      intercept <- 0.28
      coef_L    <- -1.30
      coef_k    <- 1.08
    }

    log_M <- intercept + coef_L * log(L_ratio) + coef_k * log(Minf)
    M <- exp(log_M)
  }

  M
}


#' Scale Mortality Schedule to Target Mean
#'
#' @param M Numeric vector of mortality rates.
#' @param M_target Target mean mortality. Can be scalar, function of tmax, or NULL.
#' @param tmax Maximum age (required if M_target is function or NULL).
#' @param p Survival probability (used if M_target NULL).
#'
#' @return Scaled mortality vector.
#'
#' @export
scale_mortality <- function(M, M_target = NULL, tmax = NULL, p = 0.001) {

  if (is.null(M_target)) {
    if (is.null(tmax)) stop("tmax required when M_target is NULL.", call. = FALSE)
    M_target <- -log(p) / tmax
  } else if (is.function(M_target)) {
    if (is.null(tmax)) stop("tmax required when M_target is a function.", call. = FALSE)
    M_target <- M_target(tmax)
  }

  M_mean <- mean(M, na.rm = TRUE)
  if (M_mean <= 0) {
    warning("Mean mortality <= 0; scaling not applied.", call. = FALSE)
    return(M)
  }

  M * (M_target / M_mean)
}


# =============================================================================
# SECTION 4: STOCHASTIC MORTALITY ESTIMATION
# =============================================================================


#' Stochastic Estimation of Age-Specific Natural Mortality
#'
#' @description
#' Monte Carlo simulation of age-specific natural mortality schedules with
#' uncertainty propagation from growth and maturity parameters. Accepts either
#' vitalBayes stanfit objects (preserving posterior correlations) or manual
#' \code{c(mean, sd)} specifications.
#'
#' @details
#' This function samples life history parameters and computes mortality
#' schedules using the chosen method. When vitalBayes fit objects are provided,
#' posterior correlations between parameters are preserved, yielding more
#' realistic uncertainty bounds than independent sampling.
#'
#' When maturity parameters (Lmat, tmat) come from separate model fits, the
#' \code{maturity_cor} argument allows specification of assumed correlation
#' between these parameters. Larger individuals typically mature at older ages,
#' so a positive correlation (default 0.5) is biologically reasonable.
#'
#' @section Parameter Sources:
#' Parameters can be supplied from multiple sources with the following priority:
#' \describe{
#'   \item{Linf}{growth_fit > manual}
#'   \item{L0}{growth_fit > birth_fit > manual}
#'   \item{Lmat}{growth_fit > length_maturity_fit > manual}
#'   \item{tmat}{growth_fit > age_maturity_fit > manual}
#' }
#'
#' @param method Character. Mortality model: \code{"CW"}, \code{"PW"}, or \code{"L"}.
#'
#' @param growth_fit CmdStanMCMC object from \code{fit_bayesian_growth}. If provided,
#'   extracts Linf, L0, and (for maturity-based fits) Lmat, tmat with preserved
#'   correlations. Default NULL.
#' @param birth_fit CmdStanMCMC object from \code{fit_bayesian_birth}. Provides
#'   L0 (birth size) if not available from growth_fit. Default NULL.
#' @param length_maturity_fit CmdStanMCMC object from \code{fit_bayesian_maturity}
#'   for length-at-maturity. Provides Lmat if not in growth_fit. Default NULL.
#' @param age_maturity_fit CmdStanMCMC object from \code{fit_bayesian_maturity}
#'   for age-at-maturity. Provides tmat if not in growth_fit. Default NULL.
#' @param sex Integer. Sex code (1 = female, 2 = male) for extracting from
#'   two-sex model fits. Default 1.
#'
#' @param Linf,L0,Lmat,tmat Numeric vectors of length 2: \code{c(mean, sd)}.
#'   Used when corresponding stanfit is not provided.
#' @param maturity_cor Numeric. Correlation between Lmat and tmat when both
#'   come from manual specification or separate fits. Default 0.5 (positive
#'   correlation is biologically reasonable). Set to 0 for independent sampling.
#'
#' @param growth_model Character. Growth model for G(t) and L(t): \code{"vb"},
#'   \code{"gompertz"}, or \code{"logistic"}. Default \code{"vb"}.
#' @param Linf_factor Fraction of Linf used to estimate tmax. Default 0.99.
#' @param age_seq Function or numeric vector for age grid.
#' @param iter Number of Monte Carlo iterations. Default 2000.
#'
#' @param scaled Logical. Scale mortality to M_target? Default TRUE.
#' @param M_target Target mean mortality. Can be scalar, function of tmax, or NULL.
#' @param p Survival probability for scaling when M_target is NULL. Default 0.001.
#'
#' @param two_phase Logical. Use CW two-phase senescence? Default FALSE.
#' @param late_model Character. Senescence model. Default \code{"gompertz"}.
#' @param tm_factor Numeric. Transition age as fraction of tmat. Default 2/3.
#' @param M_mult Numeric. Mortality multiplier for senescence. Default 2.
#' @param smooth_factor Numeric. Transition smoothness. Default 1/3.
#'
#' @param weight_based Logical. For Lorenzen: use weight-based? Default FALSE.
#' @param lw_fun Function. Length-weight relationship for PW and weight-based Lorenzen.
#'
#' @param seed Integer. Random seed. Default 1234.
#' @param palette Character. Color palette for plot.
#' @param print_plot Logical. Print the plot? Default TRUE.
#' @param show_progress Logical. Show progress messages? Default TRUE.
#'
#' @return A list with components:
#' \describe{
#'   \item{Schedules}{data.table with columns: set_id, age, M}
#'   \item{Parameters}{data.table with life history parameters and derived quantities}
#'   \item{Summary}{data.table with age-wise summary statistics}
#'   \item{Plot}{ggplot2 object}
#' }
#'
#' @examples
#' \dontrun{
#' # Using vitalBayes fits (full correlation preservation)
#' mort <- get_stochastic_mortality(
#'   method = "CW",
#'   growth_fit = growth_fit,  # From maturity-based fit_bayesian_growth
#'   sex = 1,
#'   growth_model = "vb"
#' )
#'
#' # Using separate maturity fits with specified correlation
#' mort <- get_stochastic_mortality(
#'   method = "CW",
#'   growth_fit = growth_fit,           # For Linf, L0
#'   length_maturity_fit = L50_fit,     # For Lmat
#'   age_maturity_fit = t50_fit,        # For tmat
#'   maturity_cor = 0.6,                # Assumed Lmat-tmat correlation
#'   sex = 1
#' )
#'
#' # Using manual parameters with correlation
#' mort <- get_stochastic_mortality(
#'   method = "CW",
#'   Linf = c(126, 10),
#'   L0 = c(35, 3),
#'   Lmat = c(83, 5),
#'   tmat = c(47, 4),
#'   maturity_cor = 0.5,
#'   growth_model = "gompertz"
#' )
#' }
#'
#' @import data.table
#' @importFrom stats rnorm quantile median cor sd
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line labs theme_bw theme
#' @importFrom ggplot2 element_text scale_x_continuous scale_y_continuous expansion
#' @export
get_stochastic_mortality <- function(
    method = c("CW", "PW", "L"),
    # vitalBayes fit objects
    growth_fit = NULL,
    birth_fit = NULL,
    length_maturity_fit = NULL,
    age_maturity_fit = NULL,
    sex = 1L,
    # Manual parameter specification
    Linf = NULL,
    L0 = NULL,
    Lmat = NULL,
    tmat = NULL,
    maturity_cor = 0.5,
    # Growth model
    growth_model = c("vb", "gompertz", "logistic"),
    Linf_factor = 0.99,
    age_seq = function(tmax) seq(0.1, ceiling(tmax), length.out = 500),
    iter = 2000,
    # Scaling
    scaled = TRUE,
    M_target = NULL,
    p = 0.001,
    # CW two-phase parameters
    two_phase = FALSE,
    late_model = c("gompertz", "logistic"),
    tm_factor = 2/3,
    M_mult = 2,
    smooth_factor = 1/3,
    # Lorenzen parameters
    weight_based = FALSE,
    lw_fun = NULL,
    # Other
    seed = 1234,
    palette = c("synthwave", "viridis", "okabe", "plasma", "inferno"),
    print_plot = TRUE,
    show_progress = TRUE
) {

  method <- match.arg(method)
  growth_model <- match.arg(growth_model)
  late_model <- match.arg(late_model)
  palette <- match.arg(palette)
  sex <- as.integer(sex)

  # -------------------------------------------------------------------------
  # Parameter Extraction
  # -------------------------------------------------------------------------

  if (show_progress) message("Extracting life history parameters...")

  par_draws <- extract_lh_params(
    growth_fit = growth_fit,
    birth_fit = birth_fit,
    length_maturity_fit = length_maturity_fit,
    age_maturity_fit = age_maturity_fit,
    Linf = Linf,
    L0 = L0,
    Lmat = Lmat,
    tmat = tmat,
    maturity_cor = maturity_cor,
    sex = sex,
    iter = iter,
    seed = seed,
    show_progress = show_progress
  )

  # -------------------------------------------------------------------------
  # Compute Derived Quantities
  # -------------------------------------------------------------------------

  if (show_progress) message("Computing derived quantities...")

  par_draws[, Minf := compute_Minf(Linf, L0, Lmat, tmat, warn = FALSE)]
  par_draws[, k_native := compute_k_native(Linf, L0, Lmat, tmat,
                                           growth_model = growth_model,
                                           warn = FALSE)]

  par_draws[, tmax := mapply(
    compute_tmax,
    Linf = Linf,
    L0 = L0,
    k = k_native,
    MoreArgs = list(growth_model = growth_model, Linf_factor = Linf_factor)
  )]

  # -------------------------------------------------------------------------
  # Remove Invalid Draws
  # -------------------------------------------------------------------------

  valid_mask <- !is.na(par_draws$Minf) &
    par_draws$Minf > 0 &
    !is.na(par_draws$k_native) &
    par_draws$k_native > 0 &
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

  median_tmax <- stats::median(par_draws$tmax)
  if (is.function(age_seq)) {
    ages <- age_seq(median_tmax)
  } else {
    ages <- age_seq
  }
  ages <- ages[ages > 0]

  n_sets <- nrow(par_draws)
  step <- max(1L, as.integer(n_sets * 0.1))
  notify <- unique(c(seq(step, n_sets, by = step), n_sets))

  schedules_list <- vector("list", n_sets)

  for (i in seq_len(n_sets)) {

    p_i <- par_draws[i]

    M_raw <- switch(
      method,

      "CW" = {
        G_t <- compute_G(ages, p_i$Linf, p_i$L0, p_i$k_native, growth_model)
        M_cw <- p_i$Minf / G_t

        if (two_phase && !is.na(p_i$tmat)) {
          tm <- tm_factor * p_i$tmat
          G_at_tm <- compute_G(tm, p_i$Linf, p_i$L0, p_i$k_native, growth_model)
          M_early_at_tm <- p_i$Minf / G_at_tm

          if (late_model == "gompertz") {
            M_s <- M_mult * M_early_at_tm
            r <- log(M_early_at_tm / M_s) / (tm - p_i$tmax)
            M_late <- M_s * exp(r * (ages - p_i$tmax))
          } else {
            K <- 2 * M_mult * M_early_at_tm
            r <- -log(K / M_early_at_tm - 1) / (tm - p_i$tmax)
            M_late <- K / (1 + exp(-r * (ages - p_i$tmax)))
          }

          smooth_width <- smooth_factor * abs(p_i$tmax - tm)
          weight <- 1 / (1 + exp(-2 * (ages - tm) / smooth_width))

          M_cw <- (1 - weight) * M_cw + weight * M_late
          early_idx <- ages < tm * 0.5
          G_early <- compute_G(ages[early_idx], p_i$Linf, p_i$L0, p_i$k_native, growth_model)
          M_cw[early_idx] <- p_i$Minf / G_early
          M_cw <- pmax(M_cw, 0)
        }
        M_cw
      },

      "PW" = {
        if (is.null(lw_fun)) {
          stop("Peterson-Wroblewski requires 'lw_fun'.", call. = FALSE)
        }
        L_t <- compute_L(ages, p_i$Linf, p_i$L0, p_i$k_native, growth_model)
        W_t <- lw_fun(L_t)
        W_t <- pmax(W_t, 1)
        1.92 * W_t^(-0.25)
      },

      "L" = {
        if (weight_based) {
          if (is.null(lw_fun)) {
            stop("Weight-based Lorenzen requires 'lw_fun'.", call. = FALSE)
          }
          L_t <- compute_L(ages, p_i$Linf, p_i$L0, p_i$k_native, growth_model)
          W_t <- lw_fun(L_t)
          W_t <- pmax(W_t, 1)
          alpha <- stats::rnorm(1, 3.69, 0.502)
          beta  <- stats::rnorm(1, -0.305, 0.029)
          alpha * W_t^beta
        } else {
          L_t <- compute_L(ages, p_i$Linf, p_i$L0, p_i$k_native, growth_model)
          L_ratio <- L_t / p_i$Linf
          intercept <- stats::rnorm(1, 0.28, 0.105)
          coef_L    <- stats::rnorm(1, -1.30, 0.059)
          coef_k    <- stats::rnorm(1, 1.08, 0.082)
          log_M <- intercept + coef_L * log(L_ratio) + coef_k * log(p_i$Minf)
          exp(log_M)
        }
      }
    )

    if (scaled) {
      if (is.null(M_target)) {
        M_target_i <- -log(p) / p_i$tmax
      } else if (is.function(M_target)) {
        M_target_i <- M_target(p_i$tmax)
      } else {
        M_target_i <- M_target
      }
      M_final <- M_raw / mean(M_raw, na.rm = TRUE) * M_target_i
    } else {
      M_final <- M_raw
    }

    schedules_list[[i]] <- data.table::data.table(
      set_id = par_draws$set_id[i],
      age    = ages,
      M      = M_final
    )

    if (show_progress && i %in% notify) {
      message(sprintf("  Progress: %d/%d (%.0f%%)", i, n_sets, 100 * i / n_sets))
    }
  }

  schedules <- data.table::rbindlist(schedules_list)

  # -------------------------------------------------------------------------
  # Summary Statistics
  # -------------------------------------------------------------------------

  if (show_progress) message("Computing summary statistics...")

  schedules[, age_round := round(age, 2)]

  summary_dt <- schedules[, .(
    M_median = stats::median(M, na.rm = TRUE),
    M_mean   = mean(M, na.rm = TRUE),
    M_lower  = stats::quantile(M, 0.025, na.rm = TRUE),
    M_upper  = stats::quantile(M, 0.975, na.rm = TRUE)
  ), by = age_round]

  tmax_summary <- par_draws[, .(
    mean  = mean(tmax),
    lower = stats::quantile(tmax, 0.025),
    upper = stats::quantile(tmax, 0.975)
  )]

  # -------------------------------------------------------------------------
  # Plotting
  # -------------------------------------------------------------------------

  if (show_progress) message("Generating plot...")

  pal <- switch(
    palette,
    "synthwave" = c("#FF6B9D", "#C490D1", "#9B6DFF", "#00D4AA"),
    "viridis"   = c("#440154", "#31688E", "#35B779", "#FDE725"),
    "okabe"     = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
    "plasma"    = c("#0D0887", "#7E03A8", "#CC4678", "#F89441"),
    "inferno"   = c("#000004", "#57106E", "#BC3754", "#F98E09")
  )

  fill_color <- pal[1]
  line_color <- pal[3]

  # Determine correlation status for caption
  cor_status <- attr(par_draws, "correlation_status")
  cor_str <- if (!is.null(cor_status)) {
    switch(cor_status,
           "full" = "posterior (full correlation)",
           "partial" = sprintf("mixed (maturity rho=%.2f)",
                               attr(par_draws, "maturity_cor")),
           "")
  } else {
    ""
  }

  caption <- sprintf(
    "tmax: %.1f yrs (95%% CI: %.1f-%.1f) | Method: %s | Growth: %s | %s",
    tmax_summary$mean, tmax_summary$lower, tmax_summary$upper,
    method, growth_model, cor_str
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
      x       = "Age (years)",
      y       = "Instantaneous Mortality (M)",
      title   = "Age-Specific Natural Mortality",
      subtitle = "Median with 95% credible interval",
      caption = caption
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

  if (show_progress) message("Done.")

  # -------------------------------------------------------------------------
  # Return
  # -------------------------------------------------------------------------

  schedules[, age_round := NULL]

  list(
    Schedules  = schedules,
    Parameters = par_draws,
    Summary    = summary_dt,
    Plot       = mort_plot
  )
}

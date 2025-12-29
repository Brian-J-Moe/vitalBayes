# =============================================================================
# vitalBayes Posterior Predictive Check Functions
# =============================================================================
# Functions for model diagnostics, validation, and comparison using
# posterior predictive checks and LOO cross-validation.
# =============================================================================

#' Posterior Predictive Summary
#'
#' @description
#' Generates a comprehensive summary of posterior predictive checks for any
#' vitalBayes model fit.
#'
#' @param fit A CmdStanMCMC object from any vitalBayes fitting function.
#' @param model_type Character. Type of model: \code{"growth"}, \code{"maturity"},
#'   or \code{"birth"}. If \code{NULL}, auto-detected.
#' @param data Optional data.table with observed data for additional checks.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list with class \code{"vitalBayes_ppc"} containing diagnostic summaries.
#'
#' @seealso
#' \code{vignette("model_diagnostics")} for comprehensive diagnostic workflow.
#'
#' \href{../doc/vitalBayes_stats_explained.html#assessment}{Statistical Methods: Model Assessment}
#' for posterior predictive check theory.
#'
#' \code{\link{compute_loo}}, \code{\link{compare_loo}}, \code{\link{plot_residuals}}
#'
#' @examples
#' \dontrun{
#' ppc <- ppc_summary(growth_fit, model_type = "growth")
#' print(ppc)
#' }
#'
#' @export
ppc_summary <- function(fit, model_type = NULL, data = NULL, ...) {
  
  # Auto-detect model type
  if (is.null(model_type)) {
    vars <- fit$metadata()$stan_variables
    
    if ("b50" %in% vars) {
      model_type <- "birth"
    } else if ("L50" %in% vars || "t50" %in% vars) {
      model_type <- "maturity"
    } else if ("Linf" %in% vars) {
      model_type <- "growth"
    } else {
      stop("Could not auto-detect model type. Please specify 'model_type'.",
           call. = FALSE)
    }
  }
  
  result <- switch(model_type,
    growth = .ppc_growth(fit),
    maturity = .ppc_maturity(fit),
    birth = .ppc_birth(fit),
    stop("Unknown model_type: ", model_type, call. = FALSE)
  )
  
  result$model_type <- model_type
  class(result) <- c("vitalBayes_ppc", "list")
  
  return(result)
}


#' @noRd
.ppc_growth <- function(fit) {
  
  vars <- fit$metadata()$stan_variables
  is_twosex <- "Linf_diff" %in% vars
  
  result <- list()
  
  # Overall fit statistics
  if ("rmse" %in% vars) {
    result$rmse <- fit$summary("rmse")
    result$mae <- fit$summary("mae")
  } else if ("rmse_f" %in% vars) {
    result$rmse_f <- fit$summary("rmse_f")
    result$rmse_m <- fit$summary("rmse_m")
  }
  
  # Coverage
  if ("n_in_CI" %in% vars) {
    result$coverage <- fit$summary("n_in_CI")
  } else if ("n_in_CI_f" %in% vars) {
    result$coverage_f <- fit$summary("n_in_CI_f")
    result$coverage_m <- fit$summary("n_in_CI_m")
  }
  
  # Residuals
  if ("mean_residual" %in% vars) {
    result$mean_residual <- fit$summary("mean_residual")
    result$sd_residual <- fit$summary("sd_residual")
  }
  
  result$is_twosex <- is_twosex
  
  return(result)
}


#' @noRd
.ppc_maturity <- function(fit) {
  
  vars <- fit$metadata()$stan_variables
  is_length <- "L50" %in% vars
  is_twosex <- if (is_length) "L50_diff" %in% vars else "t50_diff" %in% vars
  
  result <- list(
    type = if (is_length) "length" else "age",
    is_twosex = is_twosex
  )
  
  # Classification accuracy
  if ("prop_correct_rep" %in% vars) {
    result$accuracy <- fit$summary("prop_correct_rep")
  }
  
  # Calibration by group
  if (is_twosex) {
    result$calib_mature_f <- fit$summary("mean_p_mature_f")
    result$calib_mature_m <- fit$summary("mean_p_mature_m")
    result$calib_immature_f <- fit$summary("mean_p_immature_f")
    result$calib_immature_m <- fit$summary("mean_p_immature_m")
  } else {
    result$calib_mature <- fit$summary("mean_p_mature")
    result$calib_immature <- fit$summary("mean_p_immature")
  }
  
  return(result)
}


#' @noRd
.ppc_birth <- function(fit) {
  
  result <- list()
  
  result$accuracy <- fit$summary("prop_correct_rep")
  result$calib_embryo <- fit$summary("mean_p_embryo")
  result$calib_freeswim <- fit$summary("mean_p_freeswim")
  
  return(result)
}


#' Print PPC Summary
#'
#' @param x A vitalBayes_ppc object.
#' @param ... Additional arguments (unused).
#'
#' @export
print.vitalBayes_ppc <- function(x, ...) {
  
  cat("\n", paste(rep("=", 50), collapse = ""), "\n")
  cat("Posterior Predictive Check Summary\n")
  cat("Model type:", x$model_type, "\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
  
  if (x$model_type == "growth") {
    cat("Fit Statistics:\n")
    if (!is.null(x$rmse)) {
      cat("  RMSE: ", round(x$rmse$median, 2), " [",
          round(x$rmse$q5, 2), ", ", round(x$rmse$q95, 2), "]\n", sep = "")
      cat("  MAE:  ", round(x$mae$median, 2), " [",
          round(x$mae$q5, 2), ", ", round(x$mae$q95, 2), "]\n", sep = "")
    } else if (!is.null(x$rmse_f)) {
      cat("  RMSE (Female): ", round(x$rmse_f$median, 2), "\n", sep = "")
      cat("  RMSE (Male):   ", round(x$rmse_m$median, 2), "\n", sep = "")
    }
    
    cat("\n95% CI Coverage:\n")
    if (!is.null(x$coverage)) {
      cat("  Proportion: ", round(x$coverage$median, 0), " observations\n", sep = "")
    } else if (!is.null(x$coverage_f)) {
      cat("  Female: ", round(x$coverage_f$median, 0), " observations\n", sep = "")
      cat("  Male:   ", round(x$coverage_m$median, 0), " observations\n", sep = "")
    }
    
  } else if (x$model_type == "maturity") {
    cat("Maturity type:", x$type, "\n\n")
    
    cat("Classification Accuracy:\n")
    cat("  ", round(x$accuracy$median * 100, 1), "% [",
        round(x$accuracy$q5 * 100, 1), "%, ",
        round(x$accuracy$q95 * 100, 1), "%]\n\n", sep = "")
    
    cat("Calibration (mean predicted probability):\n")
    if (x$is_twosex) {
      cat("  Mature females:   ", round(x$calib_mature_f$median, 3), "\n", sep = "")
      cat("  Mature males:     ", round(x$calib_mature_m$median, 3), "\n", sep = "")
      cat("  Immature females: ", round(x$calib_immature_f$median, 3), "\n", sep = "")
      cat("  Immature males:   ", round(x$calib_immature_m$median, 3), "\n", sep = "")
    } else {
      cat("  Mature:   ", round(x$calib_mature$median, 3), 
          " (ideal: close to 1)\n", sep = "")
      cat("  Immature: ", round(x$calib_immature$median, 3), 
          " (ideal: close to 0)\n", sep = "")
    }
    
  } else if (x$model_type == "birth") {
    cat("Classification Accuracy:\n")
    cat("  ", round(x$accuracy$median * 100, 1), "% [",
        round(x$accuracy$q5 * 100, 1), "%, ",
        round(x$accuracy$q95 * 100, 1), "%]\n\n", sep = "")
    
    cat("Calibration (mean predicted probability):\n")
    cat("  Free-swimming: ", round(x$calib_freeswim$median, 3), 
        " (ideal: close to 1)\n", sep = "")
    cat("  Embryos:       ", round(x$calib_embryo$median, 3), 
        " (ideal: close to 0)\n", sep = "")
  }
  
  invisible(x)
}


# -----------------------------------------------------------------------------
# LOO Cross-Validation
# -----------------------------------------------------------------------------

#' LOO Cross-Validation for vitalBayes Models
#'
#' @description
#' Computes approximate leave-one-out cross-validation using PSIS-LOO.
#' Requires the \code{loo} package.
#'
#' @param fit A CmdStanMCMC object with log_lik generated quantities.
#' @param ... Additional arguments passed to \code{\link[loo]{loo}}.
#'
#' @return A \code{loo} object.
#'
#' @seealso
#' \code{vignette("model_diagnostics")} for comprehensive model comparison workflow.
#'
#' \href{../doc/vitalBayes_stats_explained.html#assessment}{Statistical Methods: Model Assessment}
#' for LOO-CV theory and interpretation.
#'
#' \code{\link{compare_loo}}, \code{\link{create_loo_table}}, \code{\link{ppc_summary}}
#'
#' @examples
#' \dontrun{
#' loo_growth <- compute_loo(growth_fit)
#' print(loo_growth)
#' }
#'
#' @export
compute_loo <- function(fit, ...) {
  
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("Package 'loo' is required. Install with: install.packages('loo')",
         call. = FALSE)
  }
  
  log_lik <- fit$draws("log_lik", format = "matrix")
  
  loo::loo(log_lik, ...)
}


#' Compare Models Using LOO
#'
#' @description
#' Compares multiple models using LOO-CV. Returns a comparison table
#' with expected log predictive density (elpd) differences.
#'
#' @param ... Named CmdStanMCMC objects or loo objects to compare.
#' @param criterion Character. Comparison criterion: \code{"loo"} (default)
#'   or \code{"waic"}.
#'
#' @return A data.table with model comparison statistics.
#'
#' @seealso
#' \code{vignette("model_diagnostics")} for comprehensive model comparison workflow.
#'
#' \href{../doc/vitalBayes_stats_explained.html#assessment}{Statistical Methods: Model Assessment}
#' for interpreting elpd differences.
#'
#' \code{\link{compute_loo}}, \code{\link{create_loo_table}}, \code{\link{compare_growth_models}}
#'
#' @examples
#' \dontrun{
#' comp <- compare_loo(
#'   "von Bertalanffy" = vb_fit,
#'   "Gompertz" = gomp_fit,
#'   "Logistic" = log_fit
#' )
#' print(comp)
#' }
#'
#' @import data.table
#' @export
compare_loo <- function(..., criterion = c("loo", "waic")) {
  
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("Package 'loo' is required. Install with: install.packages('loo')",
         call. = FALSE)
  }
  
  criterion <- match.arg(criterion)
  
  fits <- list(...)
  n_models <- length(fits)
  model_names <- names(fits)
  
  if (is.null(model_names) || any(model_names == "")) {
    model_names <- paste0("Model_", seq_len(n_models))
  }
  
  # Compute LOO for each model if needed
  loo_list <- lapply(fits, function(x) {
    if (inherits(x, "loo")) {
      return(x)
    } else if (inherits(x, "CmdStanMCMC")) {
      log_lik <- x$draws("log_lik", format = "matrix")
      return(loo::loo(log_lik))
    } else {
      stop("Each input must be a CmdStanMCMC or loo object.", call. = FALSE)
    }
  })
  
  names(loo_list) <- model_names
  
  # Create comparison
  comp <- loo::loo_compare(loo_list)
  
  # Convert to data.table
  comp_dt <- data.table::as.data.table(comp, keep.rownames = "model")
  
  # Add model weights
  weights <- loo::loo_model_weights(loo_list, method = "stacking")
  comp_dt[, weight := weights[match(model, names(weights))]]
  
  return(comp_dt)
}


#' Create LOO Comparison Table
#'
#' @description
#' Creates a publication-ready table from LOO comparison results.
#'
#' @param comp_dt A data.table from \code{\link{compare_loo}}.
#' @param format Character. Output format: \code{"data.table"} (default),
#'   \code{"kable"}, or \code{"gt"}.
#' @param digits Integer. Number of decimal places. Default 2.
#'
#' @return Formatted table in specified format.
#'
#' @export
create_loo_table <- function(comp_dt, format = "data.table", digits = 2) {
  
  # Create clean table
  out <- data.table::copy(comp_dt)
  
  # Round numeric columns
  num_cols <- names(out)[sapply(out, is.numeric)]
  for (col in num_cols) {
    data.table::set(out, j = col, value = round(out[[col]], digits))
  }
  
  # Rename columns for clarity
  data.table::setnames(out, 
    old = c("elpd_diff", "se_diff", "elpd_loo", "se_elpd_loo", "p_loo", "looic"),
    new = c("ELPD Diff", "SE Diff", "ELPD", "SE", "p_LOO", "LOOIC"),
    skip_absent = TRUE
  )
  
  if (format == "kable") {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      warning("Package 'knitr' not available. Returning data.table.",
              call. = FALSE)
      return(out)
    }
    return(knitr::kable(out, format = "pipe", digits = digits))
    
  } else if (format == "gt") {
    if (!requireNamespace("gt", quietly = TRUE)) {
      warning("Package 'gt' not available. Returning data.table.",
              call. = FALSE)
      return(out)
    }
    return(gt::gt(out))
  }
  
  return(out)
}


# -----------------------------------------------------------------------------
# Residual Diagnostics
# -----------------------------------------------------------------------------

#' Plot Residual Diagnostics
#'
#' @description
#' Creates diagnostic plots for growth model residuals.
#'
#' @param fit A CmdStanMCMC object from \code{\link{fit_bayesian_growth}}.
#' @param data Optional data.table with observed data.
#' @param age_col Character. Column name for age.
#' @param length_col Character. Column name for length.
#' @param sex_col Character. Column name for sex.
#' @param type Character. Type of plot: \code{"all"} (default), \code{"fitted"},
#'   \code{"qq"}, \code{"histogram"}, or \code{"age"}.
#' @param colors Color palette.
#' @param colorblind Logical. Use colorblind-safe palette? Default \code{FALSE}.
#' @param base_size Base font size.
#'
#' @return A ggplot object or list of ggplot objects.
#'
#' @examples
#' \dontrun{
#' plot_residuals(growth_fit, data = shark_data)
#' plot_residuals(growth_fit, type = "qq")
#' }
#'
#' @import ggplot2
#' @import data.table
#' @export
plot_residuals <- function(
    fit,
    data        = NULL,
    age_col     = "age",
    length_col  = "length",
    sex_col     = "sex",
    type        = c("all", "fitted", "qq", "histogram", "age"),
    colors      = NULL,
    colorblind  = FALSE,
    base_size   = 12
) {
  
  type <- match.arg(type)
  
  # Extract residuals
  resid_draws <- fit$draws("std_residual", format = "matrix")
  y_pred_draws <- fit$draws("y_pred", format = "matrix")
  
  # Use posterior means
  std_resid <- colMeans(resid_draws)
  y_pred <- colMeans(y_pred_draws)
  
  # Build data for plotting
  resid_dt <- data.table::data.table(
    fitted = y_pred,
    std_residual = std_resid,
    obs_id = seq_along(std_resid)
  )
  
  # Add observed data if provided
  if (!is.null(data)) {
    data_dt <- data.table::as.data.table(data)
    
    if (age_col %in% names(data_dt)) {
      resid_dt[, age := data_dt[[age_col]][obs_id]]
    }
    
    if (length_col %in% names(data_dt)) {
      resid_dt[, observed := data_dt[[length_col]][obs_id]]
    }
    
    if (sex_col %in% names(data_dt)) {
      resid_dt[, sex := ifelse(
        grepl("f", tolower(data_dt[[sex_col]][obs_id])), "Female", "Male"
      )]
    }
  }
  
  # Set colors from vitalBayes palette - handle colorblind option
  if (colorblind) {
    colors <- vital_colors(2, "sex_cb")
  } else if (is.null(colors)) {
    colors <- vital_colors(2, "sex")
  }
  
  has_sex <- "sex" %in% names(resid_dt)
  
  # Helper function for individual plots
  make_fitted_plot <- function() {
    p <- ggplot2::ggplot(resid_dt, ggplot2::aes(x = fitted, y = std_residual))
    
    if (has_sex) {
      p <- p + ggplot2::geom_point(ggplot2::aes(color = sex), alpha = 0.5)
      p <- p + ggplot2::scale_color_manual(values = colors, name = "Sex")
    } else {
      p <- p + ggplot2::geom_point(alpha = 0.5, color = colors[1])
    }
    
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#7B2D8E")
    p <- p + ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dotted", 
                                  color = "#7B2D8E", alpha = 0.6)
    p <- p + ggplot2::labs(x = "Fitted Values", y = "Standardized Residuals")
    p <- p + theme_vital(base_size = base_size)
    p
  }
  
  make_qq_plot <- function() {
    p <- ggplot2::ggplot(resid_dt, ggplot2::aes(sample = std_residual))
    
    if (has_sex) {
      p <- p + ggplot2::stat_qq(ggplot2::aes(color = sex), alpha = 0.5)
      p <- p + ggplot2::stat_qq_line(ggplot2::aes(color = sex))
      p <- p + ggplot2::scale_color_manual(values = colors, name = "Sex")
    } else {
      p <- p + ggplot2::stat_qq(alpha = 0.5, color = colors[1])
      p <- p + ggplot2::stat_qq_line(color = colors[1])
    }
    
    p <- p + ggplot2::labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
    p <- p + theme_vital(base_size = base_size)
    p
  }
  
  make_histogram_plot <- function() {
    p <- ggplot2::ggplot(resid_dt, ggplot2::aes(x = std_residual))
    
    if (has_sex) {
      p <- p + ggplot2::geom_histogram(
        ggplot2::aes(fill = sex), 
        bins = 30, alpha = 0.6, position = "identity"
      )
      p <- p + ggplot2::scale_fill_manual(values = colors, name = "Sex")
    } else {
      p <- p + ggplot2::geom_histogram(bins = 30, fill = colors[1], alpha = 0.7)
    }
    
    p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "#7B2D8E")
    p <- p + ggplot2::labs(x = "Standardized Residuals", y = "Count")
    p <- p + theme_vital(base_size = base_size)
    p
  }
  
  make_age_plot <- function() {
    if (!"age" %in% names(resid_dt)) {
      message("Age data not available. Returning fitted plot instead.")
      return(make_fitted_plot())
    }
    
    p <- ggplot2::ggplot(resid_dt, ggplot2::aes(x = age, y = std_residual))
    
    if (has_sex) {
      p <- p + ggplot2::geom_point(ggplot2::aes(color = sex), alpha = 0.5)
      p <- p + ggplot2::geom_smooth(ggplot2::aes(color = sex), method = "loess", 
                                     se = FALSE, linewidth = 0.8)
      p <- p + ggplot2::scale_color_manual(values = colors, name = "Sex")
    } else {
      p <- p + ggplot2::geom_point(alpha = 0.5, color = colors[1])
      p <- p + ggplot2::geom_smooth(method = "loess", se = FALSE, 
                                     color = colors[1], linewidth = 0.8)
    }
    
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#7B2D8E")
    p <- p + ggplot2::labs(x = "Age (years)", y = "Standardized Residuals")
    p <- p + theme_vital(base_size = base_size)
    p
  }
  
  if (type == "all") {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      message("Package 'patchwork' not available. Returning list of plots.")
      return(list(
        fitted = make_fitted_plot(),
        qq = make_qq_plot(),
        histogram = make_histogram_plot(),
        age = make_age_plot()
      ))
    }
    
    p1 <- make_fitted_plot() + ggplot2::ggtitle("Residuals vs Fitted")
    p2 <- make_qq_plot() + ggplot2::ggtitle("Q-Q Plot")
    p3 <- make_histogram_plot() + ggplot2::ggtitle("Residual Distribution")
    p4 <- make_age_plot() + ggplot2::ggtitle("Residuals vs Age")
    
    return((p1 + p2) / (p3 + p4))
    
  } else {
    switch(type,
      fitted = make_fitted_plot(),
      qq = make_qq_plot(),
      histogram = make_histogram_plot(),
      age = make_age_plot()
    )
  }
}


# -----------------------------------------------------------------------------
# Pooling Comparison (Two-Sex Models)
# -----------------------------------------------------------------------------

#' Compare Pooling Effects
#'
#' @description
#' Compares parameter estimates between pooled and unpooled two-sex models
#' to assess the effect of partial pooling.
#'
#' @param fit_pooled CmdStanMCMC object from a pooled two-sex model.
#' @param fit_unpooled CmdStanMCMC object from an unpooled two-sex model.
#' @param params Character vector of parameter names to compare.
#'   If \code{NULL}, uses c("Linf", "L0", "k").
#'
#' @return A data.table with parameter comparisons.
#'
#' @examples
#' \dontrun{
#' comp <- compare_pooling(growth_pooled, growth_unpooled)
#' print(comp)
#' }
#'
#' @import data.table
#' @export
compare_pooling <- function(fit_pooled, fit_unpooled, params = NULL) {
  
  if (is.null(params)) {
    params <- c("Linf", "L0", "k")
  }
  
  results <- list()
  
  for (param in params) {
    pooled_draws <- fit_pooled$draws(param, format = "matrix")
    unpooled_draws <- fit_unpooled$draws(param, format = "matrix")
    
    for (s in 1:ncol(pooled_draws)) {
      sex_label <- if (s == 1) "Female" else "Male"
      
      pooled_mean <- mean(pooled_draws[, s])
      pooled_sd <- stats::sd(pooled_draws[, s])
      
      unpooled_mean <- mean(unpooled_draws[, s])
      unpooled_sd <- stats::sd(unpooled_draws[, s])
      
      shrinkage <- 1 - (pooled_sd / unpooled_sd)
      
      results[[length(results) + 1]] <- data.table::data.table(
        parameter = param,
        sex = sex_label,
        pooled_mean = pooled_mean,
        pooled_sd = pooled_sd,
        unpooled_mean = unpooled_mean,
        unpooled_sd = unpooled_sd,
        shrinkage_pct = shrinkage * 100
      )
    }
  }
  
  result_dt <- data.table::rbindlist(results)
  
  message("\nPooling Effect Summary:")
  message("  Positive shrinkage = uncertainty reduction from pooling")
  message("  Negative shrinkage = uncertainty increase (rare)\n")
  
  return(result_dt)
}


# -----------------------------------------------------------------------------
# Parameter Table Export
# -----------------------------------------------------------------------------

#' Create Parameter Summary Table
#'
#' @description
#' Creates a publication-ready parameter summary table from a vitalBayes fit.
#'
#' @param fit A CmdStanMCMC object from any vitalBayes fitting function.
#' @param params Character vector of parameter names. If \code{NULL}, auto-selects
#'   key parameters based on model type.
#' @param probs Numeric vector of quantiles. Default c(0.025, 0.5, 0.975).
#' @param digits Integer. Decimal places. Default 3.
#' @param format Character. Output format: \code{"data.table"}, \code{"kable"}, 
#'   or \code{"gt"}.
#'
#' @return Formatted parameter table.
#'
#' @examples
#' \dontrun{
#' create_parameter_table(growth_fit)
#' create_parameter_table(growth_fit, format = "kable")
#' }
#'
#' @import data.table
#' @export
create_parameter_table <- function(
    fit,
    params = NULL,
    probs  = c(0.025, 0.5, 0.975),
    digits = 3,
    format = "data.table"
) {
  
  # Auto-select parameters if not specified
  if (is.null(params)) {
    vars <- fit$metadata()$stan_variables
    
    if ("b50" %in% vars) {
      params <- c("b50", "slope", "transition_width")
    } else if ("L50" %in% vars) {
      params <- c("L50", "slope", "L05", "L95", "transition_width")
    } else if ("t50" %in% vars) {
      params <- c("t50", "slope", "t05", "t95", "transition_width")
    } else if ("Linf" %in% vars) {
      params <- c("Linf", "L0", "k", "sigma")
      if ("Lmat" %in% vars) {
        params <- c("Linf", "L0", "Lmat", "tmat", "k", "sigma")
      }
    }
  }
  
  # Get summary
  summ <- fit$summary(variables = params, ~stats::quantile(.x, probs = probs))
  
  # Convert to data.table
  summ_dt <- data.table::as.data.table(summ)
  
  # Round
  num_cols <- names(summ_dt)[sapply(summ_dt, is.numeric)]
  for (col in num_cols) {
    data.table::set(summ_dt, j = col, value = round(summ_dt[[col]], digits))
  }
  
  # Clean column names
  data.table::setnames(summ_dt, 
    old = c("variable", "2.5%", "50%", "97.5%"),
    new = c("Parameter", "Lower 95%", "Median", "Upper 95%"),
    skip_absent = TRUE
  )
  
  if (format == "kable") {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      warning("Package 'knitr' not available. Returning data.table.")
      return(summ_dt)
    }
    return(knitr::kable(summ_dt, format = "pipe", digits = digits))
    
  } else if (format == "gt") {
    if (!requireNamespace("gt", quietly = TRUE)) {
      warning("Package 'gt' not available. Returning data.table.")
      return(summ_dt)
    }
    return(gt::gt(summ_dt))
  }
  
  return(summ_dt)
}

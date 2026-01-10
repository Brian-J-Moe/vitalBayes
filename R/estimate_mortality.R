# -----------------------------------------------------------------------------
# Individual Mortality Model Functions
# -----------------------------------------------------------------------------

#' Chen-Watanabe Natural Mortality Model
#'
#' @description
#' Computes age-specific natural mortality using the Chen & Watanabe (1989) model,
#' with optional late-life mortality phase following Gompertz or logistic senescence.
#'
#' @details
#' The Chen-Watanabe model expresses mortality as a function of the von Bertalanffy
#' growth parameters:
#' \deqn{M(t) = k / (1 - e^{-k(t - t_0)})}
#'
#' For the two-phase model, mortality transitions from the stable CW formulation to
#' a senescence phase (Gompertz or logistic) after a transition age \eqn{t_m}.
#'
#' @param age Numeric vector of ages at which to compute mortality.
#' @param Linf Asymptotic length.
#' @param k von Bertalanffy growth coefficient.
#' @param t0 Theoretical age at length zero.
#' @param tmax Maximum age (used for late-phase parameterization).
#' @param tmat Optional age at maturity (used to estimate transition age).
#' @param two_phase Logical. If TRUE, includes late-life senescence phase.
#' @param late_model Character. Senescence model: \code{"gompertz"} or \code{"logistic"}.
#' @param tm_factor Proportion of adult lifespan where late-phase transition occurs.
#'   Default 2/3.
#' @param M_mult Multiplier for mortality at tmax relative to M at transition age.
#' @param M_cap_factor Cap for M_tmax as proportion of maximum early-phase M.
#' @param smooth_factor Controls width of transition zone between phases.
#' @param mode Logistic mode: \code{"K_mult"} or \code{"r_target"}.
#' @param alpha Logistic K multiplier (if mode = "K_mult").
#' @param r_given Fixed logistic rate (if mode = "r_target").
#'
#' @return Numeric vector of instantaneous mortality rates.
#'
#' @references
#' Chen, S., & Watanabe, S. (1989). Age dependence of natural mortality coefficient
#' in fish population dynamics. Nippon Suisan Gakkaishi, 55(2), 205-208.
#'
#' @examples
#' \dontrun{
#' ages <- seq(0.1, 30, length.out = 100)
#' M <- M_chen_watanabe(ages, Linf = 200, k = 0.15, t0 = -1, tmax = 30)
#' plot(ages, M, type = "l")
#' }
#'
#' @export
M_chen_watanabe <- function(age,
                            Linf,
                            k,
                            t0,
                            tmax,
                            tmat = NULL,
                            two_phase = TRUE,
                            late_model = c("gompertz", "logistic"),
                            tm_factor = 2/3,
                            M_mult = 2,
                            M_cap_factor = 1/4,
                            smooth_factor = 1/3,
                            mode = c("K_mult", "r_target"),
                            alpha = 3,
                            r_given = NULL) {

  late_model <- match.arg(late_model)
  mode <- match.arg(mode)

  t <- age

  if (two_phase) {
    # Transition age estimation
    tm_pred <- -1/k * log(abs(1 - exp(k * t0))) + t0
    if (!is.null(tmat)) {
      tm_est <- tmat + (diff(c(tmat, tmax)) * tm_factor)
      tm <- max(tm_est, tm_pred)
    } else {
      tm <- tm_pred
    }

    # Early phase
    M_early <- k / (1 - exp(-k * (t[t <= tm] - t0)))
    M_s <- utils::tail(M_early, 1)

    if (late_model == "gompertz") {
      M_cap <- max(M_early) * M_cap_factor
      M_tmax <- min(M_s * M_mult, M_cap)
      B <- log(M_tmax / M_s) / (tmax - tm)
      A <- M_s * exp(-B * tm)

      smooth_w <- abs(tmax - tm) * smooth_factor
      weight <- 1 / (1 + exp(-2 * (t - tm) / smooth_w))

      a0 <- 1 - exp(-k * (tm - t0))
      a1 <- k * exp(-k * (tm - t0))
      a2 <- -0.5 * k^2 * exp(-k * (tm - t0))
      M_cw_late <- k / (a0 + a1 * (t - tm) + a2 * (t - tm)^2)

      M_raw <- (1 - weight) * M_cw_late + weight * (A * exp(B * t))
      M_raw[t < tm] <- M_early

    } else {
      # Logistic late phase
      compute_logistic <- function(tm, t0, k, tmax, mode, alpha, r_given) {
        x <- pmax(tm - t0, 1e-12)
        M_s <- k / (1 - exp(-k * x))
        if (mode == "K_mult") {
          K <- alpha * M_s
          r <- -(1 / (tm - tmax)) * log((K / M_s) - 1)
        } else {
          r <- if (is.null(r_given)) log(3) / (tmax - tm) else r_given
          K <- M_s * (1 + exp(-r * (tm - tmax)))
        }
        function(tt) K / (1 + exp(-r * (tt - tmax)))
      }
      M_log <- compute_logistic(tm, t0, k, tmax, mode, alpha, r_given)

      smooth_w <- abs(tmax - tm) * smooth_factor
      weight <- 1 / (1 + exp(-2 * (t - tm) / smooth_w))

      a0 <- 1 - exp(-k * (tm - t0))
      a1 <- k * exp(-k * (tm - t0))
      a2 <- -0.5 * k^2 * exp(-k * (tm - t0))
      M_cw_late <- k / (a0 + a1 * (t - tm) + a2 * (t - tm)^2)

      M_raw <- (1 - weight) * M_cw_late + weight * M_log(t)
      M_raw[t < tm] <- M_early
    }

    # Handle negative values
    neg_idx <- which(M_raw < 0)
    if (length(neg_idx) > 0) {
      first_neg <- min(neg_idx)
      if (first_neg > 1) {
        last_pos <- M_raw[first_neg - 1]
        M_raw[first_neg:length(M_raw)] <- last_pos
      }
    }

  } else {
    # Single-phase CW
    M_raw <- k / (1 - exp(-k * (t - t0)))
  }

  M_raw
}


#' Peterson-Wroblewski Natural Mortality Model
#'
#' @description
#' Computes weight-based natural mortality following Peterson & Wroblewski (1984).
#'
#' @details
#' The model expresses mortality as an allometric function of body weight:
#' \deqn{M(W) = 1.92 W^{-0.25}}
#' where \eqn{W} is body weight in grams.
#'
#' @param age Numeric vector of ages at which to compute mortality.
#' @param Linf Asymptotic length.
#' @param k von Bertalanffy growth coefficient.
#' @param t0 Theoretical age at length zero.
#' @param lw_fun Function mapping length to weight in grams: \code{lw_fun(length)}.
#'
#' @return Numeric vector of instantaneous mortality rates.
#'
#' @references
#' Peterson, I., & Wroblewski, J. S. (1984). Mortality rate of fishes in the
#' pelagic ecosystem. Canadian Journal of Fisheries and Aquatic Sciences, 41(7),
#' 1117-1120.
#'
#' @examples
#' \dontrun{
#' lw <- function(lt) 1e-5 * lt^3  # Weight in grams
#' ages <- seq(0.1, 30, length.out = 100)
#' M <- M_peterson_wroblewski(ages, Linf = 200, k = 0.15, t0 = -1, lw_fun = lw)
#' plot(ages, M, type = "l")
#' }
#'
#' @export
M_peterson_wroblewski <- function(age, Linf, k, t0, lw_fun) {

  if (is.null(lw_fun) || !is.function(lw_fun)) {
    stop("PW model requires a length-weight function 'lw_fun(length)'.", call. = FALSE)
  }

  lt <- Linf * (1 - exp(-k * (age - t0)))
  Wt <- lw_fun(lt)
  1.92 * (Wt ^ -0.25)
}


#' Lorenzen Natural Mortality Model
#'
#' @description
#' Computes size-dependent natural mortality following Lorenzen (1996, 2022).
#' Supports both weight-based and von Bertalanffy growth-based formulations.
#'
#' @details
#' \strong{Weight-based formulation:}
#' \deqn{M(W) = \alpha W^{\beta}}
#' where \eqn{\alpha \sim N(3.69, 0.502)} and \eqn{\beta \sim N(-0.305, 0.029)}.
#'
#' \strong{Growth-based formulation (Lorenzen 2022):}
#' \deqn{\ln M = 0.28 - 1.30 \ln(L/L_\infty) + 1.08 \ln(k)}
#' with uncertainty incorporated via parameter distributions.
#'
#' @param age Numeric vector of ages at which to compute mortality.
#' @param Linf Asymptotic length.
#' @param k von Bertalanffy growth coefficient.
#' @param t0 Theoretical age at length zero.
#' @param weight_based Logical. If TRUE, uses weight-based formulation requiring
#'   \code{lw_fun}. If FALSE, uses growth-based formulation.
#' @param lw_fun Function mapping length to weight in grams (required if weight_based = TRUE).
#' @param sample_params Logical. If TRUE, samples allometric parameters from their
#'   distributions (for Monte Carlo simulation). If FALSE, uses mean values.
#'
#' @return Numeric vector of instantaneous mortality rates.
#'
#' @references
#' Lorenzen, K. (1996). The relationship between body weight and natural mortality
#' in juvenile and adult fish: a comparison of natural ecosystems and aquaculture.
#' Journal of Fish Biology, 49(4), 627-642.
#'
#' Lorenzen, K. (2022). Size- and age-dependent natural mortality in fish populations:
#' Biology, models, implications, and a generalized length-inverse model.
#' Fisheries Research, 255, 106454.
#'
#' @examples
#' \dontrun{
#' # Growth-based formulation
#' ages <- seq(0.1, 30, length.out = 100)
#' M <- M_lorenzen(ages, Linf = 200, k = 0.15, t0 = -1, weight_based = FALSE)
#' plot(ages, M, type = "l")
#' }
#'
#' @importFrom stats rnorm
#' @export
M_lorenzen <- function(age,
                       Linf,
                       k,
                       t0,
                       weight_based = FALSE,
                       lw_fun = NULL,
                       sample_params = TRUE) {

  lt <- Linf * (1 - exp(-k * (age - t0)))

  if (weight_based) {
    if (is.null(lw_fun) || !is.function(lw_fun)) {
      stop("Lorenzen weight-based model requires 'lw_fun'.", call. = FALSE)
    }
    Wt <- lw_fun(lt)

    if (sample_params) {
      alpha <- stats::rnorm(1, 3.69, 0.502)
      beta <- stats::rnorm(1, -0.305, 0.029)
    } else {
      alpha <- 3.69
      beta <- -0.305
    }

    M_raw <- alpha * Wt^beta

  } else {
    # Growth-based formulation
    if (sample_params) {
      lnM <- stats::rnorm(1, 0.28, 0.105) +
        stats::rnorm(1, -1.30, 0.059) * log(lt / Linf) +
        stats::rnorm(1, 1.08, 0.082) * log(k)
    } else {
      lnM <- 0.28 - 1.30 * log(lt / Linf) + 1.08 * log(k)
    }
    M_raw <- exp(lnM)
  }

  M_raw
}


#' Scale Mortality Schedule to Target
#'
#' @description
#' Rescales an age-specific mortality schedule so the mean mortality equals a
#' target value.
#'
#' @param M Numeric vector of instantaneous mortality rates.
#' @param M_target Target mean mortality. Can be a numeric scalar, a function
#'   of tmax (e.g., \code{function(tmax) ...}), or NULL to derive from survival
#'   probability \code{p}.
#' @param tmax Maximum age (required if M_target is a function or NULL).
#' @param p Probability of surviving to tmax. Used only if M_target is NULL.
#'   Default 0.001 (0.1 percent survival to tmax).
#'
#' @return Numeric vector of scaled mortality rates.
#'
#' @details
#' The scaling applies: \eqn{M_{scaled} = M_{raw} / \bar{M}_{raw} \times M_{target}}
#'
#' If \code{M_target = NULL}, it is derived as \eqn{-\ln(p) / t_{max}}.
#'
#' @examples
#' \dontrun{
#' M_raw <- M_chen_watanabe(0:30, Linf = 200, k = 0.15, t0 = -1, tmax = 30)
#' M_scaled <- scale_mortality(M_raw, M_target = 0.2)
#'
#' # Using Hoenig-style target
#' hoenig_target <- function(tmax) 4.899 * tmax^(-0.916)
#' M_scaled <- scale_mortality(M_raw, M_target = hoenig_target, tmax = 30)
#' }
#'
#' @export
scale_mortality <- function(M, M_target = NULL, tmax = NULL, p = 0.001) {

  if (is.null(M_target)) {
    if (is.null(tmax)) stop("tmax required when M_target is NULL.")
    M_target <- -log(p) / tmax
  } else if (is.function(M_target)) {
    if (is.null(tmax)) stop("tmax required when M_target is a function.")
    M_target <- M_target(tmax)
  }

  M / mean(M, na.rm = TRUE) * M_target
}


# -----------------------------------------------------------------------------
# Main Stochastic Mortality Function
# -----------------------------------------------------------------------------

#' Stochastic Simulation of Age-Specific Natural Mortality
#'
#' @description
#' Monte Carlo simulation of age-specific natural mortality schedules under
#' parameter uncertainty. Life-history parameters can be drawn from vitalBayes
#' model posteriors (preserving correlations) or specified manually.
#'
#' @details
#' Three mortality models are available:
#' \itemize{
#'   \item \strong{CW} - Chen & Watanabe (1989), with optional late-life senescence
#'   \item \strong{PW} - Peterson & Wroblewski (1984), weight-based
#'   \item \strong{L} - Lorenzen (1996/2022), weight- or growth-based
#' }
#'
#' Schedules can be scaled to match a target mean mortality derived from:
#' empirical relationships (e.g., Hoenig 1983, Then et al. 2015), a fixed value,
#' or probability of survival to maximum age.
#'
#' @param method Character. One of \code{"CW"}, \code{"PW"}, \code{"L"}.
#' @param growth_fit Optional CmdStanMCMC object from \code{fit_bayesian_growth()}.
#'   If provided, life-history parameters are drawn from the joint posterior.
#' @param maturity_fit Optional CmdStanMCMC object from \code{fit_bayesian_maturity()}.
#'   Provides age-at-maturity (tmat) for CW model.
#' @param sex Integer. Sex code (1 = female, 2 = male) for hierarchical models.
#'   Required if fits are hierarchical.
#' @param Linf,k,t0,tmat Manual parameter specification as \code{c(mean, sd)}.
#'   Ignored if corresponding fit objects are provided.
#' @param Linf_factor Proportion of asymptotic length for tmax estimation.
#'   \code{tmax = -(1/k) * log(1 - Linf_factor) + t0}. Default 0.999.
#' @param age_seq Function returning evaluation ages given tmax.
#'   Default: \code{function(tmax) seq(0, ceiling(tmax), length.out = 2000)}.
#' @param iter Integer. Number of Monte Carlo iterations. Default 2000.
#' @param scaled Logical. Scale schedules to target mortality? Default TRUE.
#' @param M_target Target mean mortality. Can be NULL (derive from p),
#'   numeric, or function of tmax (e.g., Hoenig). Default NULL.
#' @param p Probability of survival to tmax (if M_target = NULL). Default 0.001.
#' @param two_phase Logical. Use two-phase CW model? Default TRUE.
#' @param late_model CW late-phase: \code{"gompertz"} or \code{"logistic"}.
#' @param tm_factor,M_mult,M_cap_factor,smooth_factor CW late-phase parameters.
#' @param mode,alpha,r_given CW logistic late-phase parameters.
#' @param weight_based For Lorenzen: TRUE = weight-based, FALSE = growth-based.
#' @param lw_fun Length-weight function \code{lw_fun(length)} returning grams.
#'   Required for PW and Lorenzen weight-based.
#' @param seed Random seed. Default 1234.
#' @param palette Color palette: \code{"synthwave"}, \code{"viridis"}, \code{"okabe"}.
#' @param fill_color,line_color Override specific colors. If NULL, from palette.
#' @param ribbon_alpha Ribbon transparency. Default 0.6.
#' @param linewidth Line width. Default 1.
#' @param plot_round Decimal places for age rounding in plot summary. Default 2.
#' @param xlab,ylab,title,subtitle Plot labels.
#' @param base_size Base font size. Default 11.
#' @param font_family Font family. Default "serif".
#' @param print_plot Print plot on completion? Default TRUE.
#' @param show_progress Show progress messages? Default TRUE.
#'
#' @return A list with:
#' \describe{
#'   \item{Schedules}{data.table of all mortality schedules (set_id, age, M)}
#'   \item{Parameters}{data.table of sampled life-history parameters}
#'   \item{Summary}{data.table with mean, 2.5 percent, 97.5 percent by age}
#'   \item{Plot}{ggplot2 object}
#' }
#'
#' @examples
#' \dontrun{
#' # Using vitalBayes fits
#' mort <- get_stochastic_mortality(
#'   method = "CW",
#'   growth_fit = growth_fit,
#'   maturity_fit = maturity_fit,
#'   sex = 2,
#'   scaled = TRUE,
#'   M_target = function(tmax) 4.899 * tmax^(-0.916)
#' )
#'
#' # Manual specification
#' mort <- get_stochastic_mortality(
#'   method = "L",
#'   Linf = c(484, 10),
#'   k = c(0.17, 0.02),
#'   t0 = c(-0.97, 0.05),
#'   weight_based = FALSE,
#'   palette = "okabe"
#' )
#' }
#'
#' @import data.table
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line labs
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous theme_bw theme
#' @importFrom ggplot2 element_text expansion
#' @importFrom stats rnorm quantile sd
#' @export
get_stochastic_mortality <- function(
    method = c("CW", "PW", "L"),
    # vitalBayes fit inputs (preferred)
    growth_fit = NULL,
    maturity_fit = NULL,
    sex = NULL,
    # Manual parameter inputs (alternative)
    Linf = NULL,
    k = NULL,
    t0 = NULL,
    tmat = NULL,
    # Model configuration
    Linf_factor = 0.999,
    age_seq = function(tmax) seq(0, ceiling(tmax), length.out = 2000),
    iter = 2000,
    scaled = TRUE,
    M_target = NULL,
    p = 0.001,
    # CW-specific
    two_phase = TRUE,
    late_model = "gompertz",
    tm_factor = 2/3,
    M_mult = 2,
    M_cap_factor = 1/4,
    smooth_factor = 1/3,
    mode = "K_mult",
    alpha = 3,
    r_given = NULL,
    # Lorenzen-specific
    weight_based = FALSE,
    # Shared
    lw_fun = NULL,
    seed = 1234,
    # Plotting
    palette = c("synthwave", "viridis", "okabe", "plasma", "inferno"),
    fill_color = NULL,
    line_color = NULL,
    ribbon_alpha = 0.6,
    linewidth = 1,
    plot_round = 2,
    xlab = "Age (yrs)",
    ylab = "Instantaneous Mortality (M)",
    title = "Natural Mortality Rates",
    subtitle = "Median with 95% credible interval",
    base_size = 11,
    font_family = "serif",
    print_plot = TRUE,
    show_progress = TRUE
) {

  method <- match.arg(method)
  palette <- match.arg(palette)

  # Determine parameter source
  use_posterior <- !is.null(growth_fit)

  if (use_posterior) {
    # Extract from fits
    if (show_progress) message("Extracting parameters from posterior...")

    par_draws <- extract_lh_params(
      growth_fit = growth_fit,
      maturity_fit = maturity_fit,
      sex = sex,
      format = "draws",
      n_draws = iter
    )

    # Ensure we have required columns
    req_cols <- c("Linf", "k", "t0")
    missing <- setdiff(req_cols, names(par_draws))
    if (length(missing) > 0) {
      stop("Missing required parameters from fit: ", paste(missing, collapse = ", "))
    }

    # Add tmat from manual if not in fit
    if (!"tmat" %in% names(par_draws) && !is.null(tmat)) {
      if (length(tmat) == 2) {
        set.seed(seed)
        par_draws[, tmat := stats::rnorm(.N, tmat[1], tmat[2])]
      }
    }

  } else {
    # Manual parameters
    if (is.null(Linf) || is.null(k) || is.null(t0)) {
      stop("Must provide either growth_fit or manual Linf, k, t0 parameters.")
    }
    if (any(sapply(list(Linf, k, t0), function(x) length(x) != 2))) {
      stop("Manual parameters must each be c(mean, sd).")
    }

    if (show_progress) message("Sampling parameters from specified distributions...")

    set.seed(seed)
    par_draws <- data.table::data.table(
      Linf = stats::rnorm(iter, Linf[1], Linf[2]),
      k = pmax(stats::rnorm(iter, k[1], k[2]), 1e-6),
      t0 = stats::rnorm(iter, t0[1], t0[2])
    )

    if (!is.null(tmat) && length(tmat) == 2) {
      par_draws[, tmat := stats::rnorm(iter, tmat[1], tmat[2])]
    }
  }

  # Compute tmax for each parameter set
  par_draws[, tmax := -(1 / k) * log(1 - Linf_factor) + t0]
  par_draws[, set_id := .I]

  # Build schedules
  if (show_progress) message("Building mortality schedules...")

  n_iter <- nrow(par_draws)
  step <- max(1L, as.integer(n_iter * 0.1))
  notify <- unique(c(seq(step, n_iter, by = step), n_iter))

  schedule_list <- vector("list", n_iter)

  for (i in seq_len(n_iter)) {
    prm <- par_draws[i]
    ages <- age_seq(prm$tmax)

    M_raw <- switch(
      method,
      CW = M_chen_watanabe(
        age = ages, Linf = prm$Linf, k = prm$k, t0 = prm$t0,
        tmax = prm$tmax, tmat = prm$tmat,
        two_phase = two_phase, late_model = late_model,
        tm_factor = tm_factor, M_mult = M_mult,
        M_cap_factor = M_cap_factor, smooth_factor = smooth_factor,
        mode = mode, alpha = alpha, r_given = r_given
      ),
      PW = M_peterson_wroblewski(
        age = ages, Linf = prm$Linf, k = prm$k, t0 = prm$t0,
        lw_fun = lw_fun
      ),
      L = M_lorenzen(
        age = ages, Linf = prm$Linf, k = prm$k, t0 = prm$t0,
        weight_based = weight_based, lw_fun = lw_fun
      )
    )

    if (scaled) {
      M_raw <- scale_mortality(M_raw, M_target = M_target, tmax = prm$tmax, p = p)
    }

    schedule_list[[i]] <- data.table::data.table(
      set_id = i,
      age = ages,
      M = M_raw
    )

    if (show_progress && i %in% notify) {
      message(sprintf("  Progress: %d/%d (%.0f%%)", i, n_iter, 100 * i / n_iter))
    }
  }

  schedules <- data.table::rbindlist(schedule_list)

  # Summarize
  if (show_progress) message("Summarizing schedules...")

  schedules[, age_round := round(age, plot_round)]
  summary_dt <- schedules[, .(
    M_median = stats::median(M, na.rm = TRUE),
    M_mean = mean(M, na.rm = TRUE),
    M_lower = stats::quantile(M, 0.025, na.rm = TRUE),
    M_upper = stats::quantile(M, 0.975, na.rm = TRUE)
  ), by = age_round]

  tmax_summary <- par_draws[, .(
    mean = mean(tmax),
    lower = stats::quantile(tmax, 0.025),
    upper = stats::quantile(tmax, 0.975)
  )]

  # Plot
  if (show_progress) message("Generating plot...")

  pal <- vital_palette(n = 4, type = palette)
  if (is.null(fill_color)) fill_color <- pal[1]
  if (is.null(line_color)) line_color <- pal[3]

  caption <- sprintf(
    "Estimated maximum age: %.1f yrs (95%% CI: %.1f - %.1f)",
    tmax_summary$mean, tmax_summary$lower, tmax_summary$upper
  )

  mort_plot <- ggplot2::ggplot(summary_dt, ggplot2::aes(x = age_round)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = M_lower, ymax = M_upper),
      fill = fill_color, alpha = ribbon_alpha
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = M_median),
      color = line_color, linewidth = linewidth
    ) +
    ggplot2::labs(
      x = xlab, y = ylab,
      title = title, subtitle = subtitle,
      caption = caption
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = font_family),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(face = "italic"),
      plot.caption = ggplot2::element_text(hjust = 0),
      axis.title = ggplot2::element_text(face = "bold")
    )

  if (print_plot) print(mort_plot)

  list(
    Schedules = schedules,
    Parameters = par_draws,
    Summary = summary_dt,
    Plot = mort_plot
  )
}



#' Get Mortality Plot Colors
#'
#' Internal helper to generate consistent colors for mortality/survival plots.
#'
#' @param type Character. Palette type passed to \code{vital_palette()}.
#' @param fill_idx Index for fill color in palette.
#' @param line_idx Index for line color in palette.
#'
#' @return Named list with fill and line colors.
#'
#' @noRd
.get_mort_colors <- function(type = "synthwave", fill_idx = 1, line_idx = 3) {
  pal <- vital_palette(n = max(fill_idx, line_idx), type = type)
  list(fill = pal[fill_idx], line = pal[line_idx])
}


# -----------------------------------------------------------------------------
# Utility Functions
# -----------------------------------------------------------------------------

#' Approximate a Standard Deviation from a Confidence Interval
#'
#' @description
#' Uses a normal approximation to convert confidence interval bounds to a
#' standard deviation: SD = (upper - lower) / (2 * z), where z = qnorm(1 - (1-CI)/2).
#' Useful for extracting uncertainty from published confidence intervals.
#'
#' @param lower Numeric vector of lower CI bounds.
#' @param upper Numeric vector of upper CI bounds (same length as lower).
#' @param CI Scalar or numeric vector in (0, 1); default 0.95.
#'
#' @return A numeric vector of approximate standard deviations.
#'
#' @examples
#' \dontrun{
#' approx_sd(8, 12)
#' approx_sd(c(8, 10), c(12, 15), 0.90)
#' }
#'
#' @importFrom stats qnorm
#' @export
approx_sd <- function(lower, upper, CI = 0.95) {

  if (!identical(length(lower), length(upper))) {
    stop("`lower` and `upper` must have the same length.", call. = FALSE)
  }
  if (any(CI <= 0 | CI >= 1, na.rm = TRUE)) {
    stop("`CI` must be strictly between 0 and 1.", call. = FALSE)
  }

  # Recycle CI if needed
  if (length(CI) == 1L) CI <- rep(CI, length(lower))

  tails <- (1 - CI) / 2
  z <- stats::qnorm(1 - tails)
  (upper - lower) / (2 * z)
}


#' Extract Life History Parameters from vitalBayes Fits
#'
#' @description
#' Extracts von Bertalanffy growth parameters and maturity information from
#' vitalBayes model fits. Can return full posterior draws (preserving correlations)
#' or summary statistics.
#'
#' @param growth_fit Optional CmdStanMCMC object from \code{fit_bayesian_growth()}.
#' @param maturity_fit Optional CmdStanMCMC object from \code{fit_bayesian_maturity()}
#'   for age-at-maturity (t50).
#' @param sex Integer. Sex code (1 = female, 2 = male) for hierarchical models.
#'   If NULL and model is hierarchical, extracts both sexes.
#' @param format Character. Either \code{"draws"} for full posterior draws as
#'   data.table (default), or \code{"summary"} for mean and SD for each parameter.
#' @param n_draws Integer. Number of posterior draws to sample (if format = "draws").
#'   If NULL, uses all draws. Default NULL.
#'
#' @return If format = "draws": data.table with columns for each parameter and
#'   optional \code{.draw} index and \code{sex} indicator.
#'   If format = "summary": data.table with columns \code{parameter}, \code{mean},
#'   \code{sd}, and optional \code{sex}.
#'
#' @details
#' This function preserves the joint posterior distribution including correlations
#' between parameters. This is statistically superior to independent resampling from
#' marginal \code{c(mean, sd)} specifications, particularly when correlations are
#' substantial (as is often the case for Linf, k, and t0).
#'
#' @examples
#' \dontrun{
#' # Extract full draws for males
#' pars <- extract_lh_params(growth_fit, maturity_fit, sex = 2, format = "draws")
#'
#' # Extract summaries for both sexes
#' pars <- extract_lh_params(growth_fit, maturity_fit, sex = NULL, format = "summary")
#' }
#'
#' @import data.table
#' @export
extract_lh_params <- function(growth_fit = NULL,
                              maturity_fit = NULL,
                              sex = NULL,
                              format = c("draws", "summary"),
                              n_draws = NULL) {

  format <- match.arg(format)

  if (is.null(growth_fit) && is.null(maturity_fit)) {
    stop("At least one of growth_fit or maturity_fit must be provided.")
  }

  # Helper to extract parameter draws
  .extract_draws <- function(fit, param, sex_code = NULL) {
    if (is.null(fit)) return(NULL)

    draws <- tryCatch(
      fit$draws(param, format = "matrix"),
      error = function(e) NULL
    )
    if (is.null(draws)) return(NULL)

    is_hier <- ncol(draws) > 1

    if (is_hier && !is.null(sex_code)) {
      return(draws[, sex_code, drop = TRUE])
    } else if (is_hier && is.null(sex_code)) {
      # Return both as list
      return(list(F = draws[, 1], M = draws[, 2]))
    } else {
      return(draws[, 1, drop = TRUE])
    }
  }

  # Check if hierarchical
  is_hierarchical <- FALSE
  if (!is.null(growth_fit)) {
    test_draws <- growth_fit$draws("Linf", format = "matrix")
    is_hierarchical <- ncol(test_draws) > 1
  }

  # Build parameter list
  params <- c("Linf", "k", "L0")
  t0_param <- if (!is.null(growth_fit)) {
    if ("t0" %in% growth_fit$metadata()$stan_variables) "t0" else NULL
  } else NULL

  # Extract draws
  extract_for_sex <- function(s) {
    dt_list <- list()

    if (!is.null(growth_fit)) {
      dt_list$Linf <- .extract_draws(growth_fit, "Linf", s)
      dt_list$k <- .extract_draws(growth_fit, "k", s)
      dt_list$L0 <- .extract_draws(growth_fit, "L0", s)

      # t0 derivation: for VBGM, t0 = -1/k * log(1 - L0/Linf)
      if (!is.null(dt_list$Linf) && !is.null(dt_list$L0) && !is.null(dt_list$k)) {
        dt_list$t0 <- -1/dt_list$k * log(pmax(1e-10, 1 - dt_list$L0/dt_list$Linf))
      }
    }

    if (!is.null(maturity_fit)) {
      dt_list$tmat <- .extract_draws(maturity_fit, "t50", s)
    }

    # Remove NULLs
    dt_list <- dt_list[!sapply(dt_list, is.null)]

    if (length(dt_list) == 0) return(NULL)

    # Convert to data.table
    dt <- data.table::as.data.table(dt_list)
    dt[, .draw := .I]

    # Subsample if requested
    if (!is.null(n_draws) && n_draws < nrow(dt)) {
      idx <- sample.int(nrow(dt), n_draws, replace = FALSE)
      dt <- dt[idx]
      dt[, .draw := .I]
    }

    dt
  }

  # Extract based on sex specification
  if (is_hierarchical && is.null(sex)) {
    # Both sexes
    dt_f <- extract_for_sex(1)
    dt_m <- extract_for_sex(2)
    if (!is.null(dt_f)) dt_f[, sex := "F"]
    if (!is.null(dt_m)) dt_m[, sex := "M"]
    result <- data.table::rbindlist(list(dt_f, dt_m), fill = TRUE)
  } else {
    result <- extract_for_sex(sex)
    if (!is.null(result) && is_hierarchical) {
      result[, sex := ifelse(sex == 1, "F", "M")]
    }
  }

  # Convert to summary if requested
  if (format == "summary") {
    if ("sex" %in% names(result)) {
      param_cols <- setdiff(names(result), c(".draw", "sex"))
      result <- result[, lapply(.SD, function(x) list(mean = mean(x), sd = stats::sd(x))),
                       .SDcols = param_cols, by = sex]
      result <- data.table::melt(result, id.vars = "sex",
                                 variable.name = "parameter",
                                 value.name = "stats")
      result[, `:=`(mean = sapply(stats, `[[`, "mean"),
                    sd = sapply(stats, `[[`, "sd"))]
      result[, stats := NULL]
    } else {
      param_cols <- setdiff(names(result), ".draw")
      result <- data.table::data.table(
        parameter = param_cols,
        mean = sapply(param_cols, function(p) mean(result[[p]])),
        sd = sapply(param_cols, function(p) stats::sd(result[[p]]))
      )
    }
  }

  result
}

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
#' function automatically detects the growth model type and extracts posterior
#' draws of the native growth coefficient and biological milestones. All
#' mortality models use the native growth trajectory for body-size predictions,
#' while CW and growth-based Lorenzen additionally receive the VB-equivalent
#' \eqn{k} for their model-specific parameters that were derived or calibrated
#' under VB assumptions.
#'
#' @section Growth Coefficient Roles:
#' The three mortality models use \eqn{k} in fundamentally different ways:
#' \describe{
#'   \item{CW}{\eqn{M(t) = k_{VB} \cdot L_\infty / L(t)}: the \eqn{k_{VB}}
#'     appears as the asymptotic mortality rate constant (derived under VB
#'     assumptions), while \eqn{L(t)} is predicted by the native growth model.}
#'   \item{Lorenzen growth-based}{\eqn{\ln M = a_0 + a_1 \ln(L/L_\infty) +
#'     a_2 \ln(k_{VB})}: the \eqn{k_{VB}} enters as a calibration coefficient
#'     (fitted against VB parameters), while \eqn{L(t)/L_\infty} uses the
#'     native trajectory.}
#'   \item{PW and Lorenzen weight-based}{Operate entirely on predicted body
#'     weight: \eqn{M(W(t))}. No \eqn{k} appears in the mortality equation
#'     itself; only the native growth trajectory matters.}
#' }
#'
#' @section Manual Parameter Specification:
#' When a growth fit is not available, parameters can be specified manually.
#' Two sampling modes are supported:
#' \describe{
#'   \item{Bivariate (preferred)}{When both \code{Lmat} and \code{tmat} are
#'     provided, they are sampled jointly from a bivariate normal with
#'     correlation \code{rho_Lmat_tmat}. The VB-equivalent \eqn{k} is then
#'     derived from the sampled milestones.}
#'   \item{Independent (fallback)}{When only \code{k} is provided, it is
#'     sampled independently from a normal distribution and assumed to be
#'     VB-equivalent for CW and growth-based Lorenzen. For non-VB growth
#'     models, the derivative-matching approach provides a VB-equivalent.}
#' }
#'
#' @param method Character. Mortality model: \code{"CW"}, \code{"PW"}, or
#'   \code{"L"}.
#' @param growth_fit Optional \code{CmdStanMCMC} object from
#'   \code{\link{fit_bayesian_growth}}. If provided, growth model type is
#'   auto-detected and parameters are extracted from the joint posterior.
#' @param maturity_fit Optional \code{CmdStanMCMC} object from
#'   \code{\link{fit_bayesian_maturity}} providing age-at-maturity.
#' @param sex Integer. Sex code (1 = female, 2 = male) for hierarchical models.
#' @param Linf,L0 Required for manual mode: \code{c(mean, sd)} vectors.
#' @param k Optional \code{c(mean, sd)} for native growth coefficient (manual
#'   mode). Required for PW and weight-based Lorenzen. Optional for CW and
#'   growth-based Lorenzen when \code{Lmat} and \code{tmat} are provided (the
#'   VB-equivalent will be derived from milestones).
#' @param Lmat,tmat Optional \code{c(mean, sd)} vectors for maturity milestones.
#'   When both are provided, enables bivariate sampling and milestone-based
#'   \eqn{k_{VB}^{equiv}} derivation.
#' @param rho_Lmat_tmat Numeric in (-1, 1). Correlation between \eqn{L_{mat}}
#'   and \eqn{t_{mat}} for bivariate sampling. Default 0.5.
#' @param Linf_factor Numeric in (0, 1). Fraction of \eqn{L_\infty} for
#'   \eqn{t_{max}} estimation. Default 0.99.
#' @param age_seq Function or numeric vector defining age grid. Default
#'   \code{function(tmax) seq(0.1, ceiling(tmax), length.out = 500)}.
#' @param iter Number of Monte Carlo iterations. Default 2000.
#' @param scaled Logical. Scale mortality to target? Default \code{TRUE}.
#' @param M_target Target mean mortality (numeric, function of tmax, or
#'   \code{NULL} for survival-based scaling).
#' @param p Survival probability to \eqn{t_{max}}. Default 0.001.
#' @param two_phase Logical. For CW, use two-phase senescence? Default
#'   \code{TRUE}.
#' @param late_model Character. Senescence model: \code{"gompertz"} or
#'   \code{"logistic"}.
#' @param tm_factor,M_mult,smooth_factor Two-phase model parameters.
#' @param lw_fun Length-weight function for PW and weight-based Lorenzen.
#' @param weight_based Logical. For Lorenzen, use weight-based? Default
#'   \code{FALSE}.
#' @param growth_model Character. Growth model type for manual specification:
#'   \code{"vb"}, \code{"gompertz"}, or \code{"logistic"}. Ignored when
#'   \code{growth_fit} is provided (auto-detected).
#' @param seed Random seed. Default 1234.
#' @param palette Color palette for plot.
#' @param print_plot Logical. Print plot on completion? Default \code{TRUE}.
#' @param show_progress Logical. Show messages? Default \code{TRUE}.
#'
#' @return A list with components:
#' \describe{
#'   \item{Schedules}{data.table of all mortality schedules (set_id, age,
#'     M_raw, M_scaled)}
#'   \item{Parameters}{data.table of sampled life history parameters including
#'     k_original (native), k_vb_equiv (VB-equivalent), and growth_model}
#'   \item{Summary}{data.table with median and 95% CI by age}
#'   \item{Plot}{ggplot2 object}
#' }
#'
#' @examples
#' \dontrun{
#' # From a Gompertz growth fit (auto-detects model, uses native trajectory)
#' mort <- get_stochastic_mortality(
#'   method     = "CW",
#'   growth_fit = gomp_fit,
#'   sex        = 1,
#'   iter       = 2000
#' )
#'
#' # Manual: maturity-based (bivariate Lmat-tmat sampling)
#' mort_mat <- get_stochastic_mortality(
#'   method        = "CW",
#'   Linf          = c(100, 8),
#'   L0            = c(25, 2),
#'   Lmat          = c(72, 4),
#'   tmat          = c(12, 1.5),
#'   rho_Lmat_tmat = 0.5,
#'   growth_model  = "gompertz",
#'   k             = c(0.15, 0.02),  # native Gompertz k for L(t)
#'   iter          = 2000
#' )
#'
#' # Manual: k-based (backward-compatible, VB assumed)
#' mort_k <- get_stochastic_mortality(
#'   method = "CW",
#'   Linf   = c(100, 8),
#'   L0     = c(25, 2),
#'   k      = c(0.08, 0.01),
#'   tmat   = c(12, 1.5),
#'   iter   = 2000
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
    Lmat = NULL,
    tmat = NULL,
    rho_Lmat_tmat = 0.5,
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

  # =========================================================================
  # Parameter Extraction / Generation
  # =========================================================================

  if (use_posterior) {

    # ----- Posterior path: extract from growth fit -----

    if (show_progress) message("Extracting parameters from growth model posterior...")

    params <- extract_growth_parameters(
      growth_fit   = growth_fit,
      maturity_fit = maturity_fit,
      sex          = sex,
      n_draws      = iter,
      seed         = seed
    )

    # Auto-detect growth model from fit metadata or extraction result
    gm <- unique(params$growth_model)
    if (length(gm) != 1) {
      stop("Unexpected: multiple growth models in parameter draws.", call. = FALSE)
    }

    # Report auto-detection
    if (show_progress && gm != "vb") {
      k_vb_method <- attr(params, "k_vb_method") %||% "unknown"
      message(
        sprintf("  Growth model: %s (k_vb derived via %s)",
                tools::toTitleCase(gm), k_vb_method)
      )
    }

    # Validate that k_vb_equiv is available for methods that need it
    needs_k_vb <- method == "CW" || (method == "L" && !weight_based)
    if (needs_k_vb && all(is.na(params$k_vb_equiv))) {
      stop(
        "CW and growth-based Lorenzen require VB-equivalent k.\n",
        "Provide maturity_fit for milestone-based derivation, or ",
        "use maturity-based growth parameterization (k_based = FALSE).",
        call. = FALSE
      )
    }

    par_draws <- data.table::data.table(
      set_id       = seq_len(nrow(params)),
      Linf         = params$Linf,
      L0           = params$L0,
      Lmat         = params$Lmat,
      tmat         = params$tmat,
      k_native     = params$k,
      k_vb_equiv   = params$k_vb_equiv,
      growth_model = gm
    )

  } else {

    # ----- Manual parameter specification -----

    if (show_progress) message("Generating parameters from specified distributions...")

    has_Lmat <- !is.null(Lmat)
    has_tmat <- !is.null(tmat)
    has_k    <- !is.null(k)
    use_bivariate <- has_Lmat && has_tmat

    # --- Validate required inputs ---

    if (is.null(Linf) || is.null(L0)) {
      stop("Must provide growth_fit OR at minimum: Linf, L0.", call. = FALSE)
    }

    # k is always needed for non-VB growth models (for L(t) prediction)
    if (growth_model != "vb" && !has_k) {
      stop(
        "Non-VB growth models require 'k' (the native growth coefficient) ",
        "for length-at-age prediction.",
        call. = FALSE
      )
    }

    # k is required unless Lmat+tmat can derive k_vb (for VB, k = k_vb)
    if (!has_k && !use_bivariate) {
      stop(
        "Must provide 'k' (as c(mean, sd)) OR both 'Lmat' and 'tmat' ",
        "(to derive VB-equivalent k from maturity milestones).",
        call. = FALSE
      )
    }

    # CW two-phase always needs tmat
    if (method == "CW" && two_phase && !has_tmat) {
      stop("tmat required for CW two-phase model.", call. = FALSE)
    }

    # Validate rho
    if (use_bivariate && abs(rho_Lmat_tmat) >= 1) {
      stop("rho_Lmat_tmat must be strictly between -1 and 1.", call. = FALSE)
    }

    # --- Sample parameters ---

    Linf_draws <- stats::rnorm(iter, Linf[1], Linf[2])
    L0_draws   <- stats::rnorm(iter, L0[1], L0[2])

    # Sample k (native growth coefficient) if provided
    if (has_k) {
      k_draws <- stats::rnorm(iter, k[1], k[2])
    }

    if (use_bivariate) {
      # ----- Bivariate Lmat-tmat sampling -----

      if (show_progress) {
        message(
          sprintf("  Sampling Lmat and tmat jointly (rho = %.2f)...",
                  rho_Lmat_tmat)
        )
      }

      mat_draws <- .sample_bivariate_normal(
        n   = iter,
        mu1 = Lmat[1], sd1 = Lmat[2],
        mu2 = tmat[1], sd2 = tmat[2],
        rho = rho_Lmat_tmat
      )

      Lmat_draws <- mat_draws[, 1]
      tmat_draws <- mat_draws[, 2]

      # Enforce biological constraints
      Lmat_draws <- pmax(Lmat_draws, L0_draws + 1)
      tmat_draws <- pmax(tmat_draws, 0.1)

      # Derive VB-equivalent k from milestones (milestone-based strategy)
      k_vb_draws <- compute_k_vb_equivalent(
        Linf = Linf_draws,
        L0   = L0_draws,
        Lmat = Lmat_draws,
        tmat = tmat_draws,
        warn = TRUE
      )

      # If user didn't provide k, derive native k too
      # For VB: native k = k_vb. For non-VB: user must provide k.
      if (!has_k) {
        # Only reachable when growth_model == "vb" (validated above)
        k_draws <- k_vb_draws
      }

    } else {
      # ----- Independent sampling (backward-compatible path) -----

      Lmat_draws <- rep(NA_real_, iter)

      if (has_tmat) {
        tmat_draws <- stats::rnorm(iter, tmat[1], tmat[2])
      } else {
        tmat_draws <- rep(NA_real_, iter)
      }

      # Derive VB-equivalent k based on growth model
      if (growth_model == "vb") {
        # VB: native k IS the VB-equivalent k
        k_vb_draws <- k_draws
      } else {
        # Non-VB: derivative matching at birth
        k_vb_draws <- .derive_k_vb_derivative(
          k_draws, Linf_draws, L0_draws, growth_model
        )
        if (show_progress) {
          message(
            sprintf(
              "  %s growth model: VB-equivalent k derived via derivative ",
              tools::toTitleCase(growth_model)
            ),
            "matching at birth."
          )
        }
      }
    }

    # --- Enforce biological constraints ---
    Linf_draws <- pmax(Linf_draws, L0_draws + 1)
    L0_draws   <- pmax(L0_draws, 0.1)
    k_draws    <- pmax(k_draws, 0.001)
    k_vb_draws <- pmax(k_vb_draws, 0.001)
    if (has_tmat || use_bivariate) {
      tmat_draws <- pmax(tmat_draws, 0.1)
    }

    par_draws <- data.table::data.table(
      set_id       = seq_len(iter),
      Linf         = Linf_draws,
      L0           = L0_draws,
      Lmat         = Lmat_draws,
      tmat         = tmat_draws,
      k_native     = k_draws,
      k_vb_equiv   = k_vb_draws,
      growth_model = growth_model
    )
  }

  # =========================================================================
  # Estimate tmax
  # =========================================================================
  # tmax is computed using the VB equation with k_vb_equiv, because it
  # represents the lifespan implied by the CW/scaling framework. This keeps
  # tmax consistent across growth model choices for the same biological
  # milestones.

  par_draws[, tmax := -log(Linf * (1 - Linf_factor) / (Linf - L0)) / k_vb_equiv]

  # Remove invalid draws
  valid_mask <- !is.na(par_draws$k_vb_equiv) &
    !is.na(par_draws$k_native) &
    par_draws$k_vb_equiv > 0 &
    par_draws$k_native > 0 &
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

  # =========================================================================
  # Mortality Calculation
  # =========================================================================

  if (show_progress) message(sprintf("Computing %s mortality schedules...", method))

  # Common age grid (based on median tmax)
  median_tmax <- stats::median(par_draws$tmax)
  if (is.function(age_seq)) {
    ages <- age_seq(median_tmax)
  } else {
    ages <- age_seq
  }
  ages <- ages[ages > 0]

  # Growth model (constant across draws)
  gm <- par_draws$growth_model[1]

  schedules_list <- vector("list", nrow(par_draws))

  for (i in seq_len(nrow(par_draws))) {

    p_i <- par_draws[i]

    M_raw <- switch(
      method,

      # ---------------------------------------------------------------
      # Chen-Watanabe: k_vb for M_inf, native k for L(t)
      # ---------------------------------------------------------------
      "CW" = M_chen_watanabe_L0(
        age          = ages,
        Linf         = p_i$Linf,
        L0           = p_i$L0,
        k            = p_i$k_vb_equiv,   # M_inf = k_vb * Linf / L(t)
        k_native     = p_i$k_native,     # For native L(t) prediction
        growth_model = gm,
        tmax         = p_i$tmax,
        Linf_factor  = Linf_factor,
        two_phase    = two_phase,
        tmat         = p_i$tmat,
        late_model   = late_model,
        tm_factor    = tm_factor,
        M_mult       = M_mult,
        smooth_factor = smooth_factor
      ),

      # ---------------------------------------------------------------
      # Peterson-Wroblewski: native k + native growth model for L(t)
      # No k enters the mortality equation itself.
      # ---------------------------------------------------------------
      "PW" = {
        if (is.null(lw_fun)) {
          stop("Peterson-Wroblewski requires 'lw_fun'.", call. = FALSE)
        }
        M_peterson_wroblewski(
          age          = ages,
          Linf         = p_i$Linf,
          L0           = p_i$L0,
          k            = p_i$k_native,    # Native k for L(t)
          lw_fun       = lw_fun,
          growth_model = gm                # Native growth model
        )
      },

      # ---------------------------------------------------------------
      # Lorenzen: native k for L(t), k_vb for growth-based coefficient
      # ---------------------------------------------------------------
      "L" = M_lorenzen(
        age          = ages,
        Linf         = p_i$Linf,
        L0           = p_i$L0,
        k            = p_i$k_native,      # Native k for L(t)
        k_vb         = p_i$k_vb_equiv,    # VB-equiv for ln(k) coefficient
        lw_fun       = lw_fun,
        weight_based = weight_based,
        growth_model = gm,                 # Native growth model
        sample_params = TRUE
      )
    )

    # Scale if requested
    if (scaled) {
      M_scaled <- scale_mortality(M_raw, age = ages, M_target = M_target, tmax = p_i$tmax, p = p)
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

  # =========================================================================
  # Summary Statistics
  # =========================================================================

  if (show_progress) message("Computing summary statistics...")

  schedules[, age_round := round(age, 2)]

  summary_dt <- schedules[, .(
    M_median = stats::median(M_scaled, na.rm = TRUE),
    M_mean   = mean(M_scaled, na.rm = TRUE),
    M_lower  = stats::quantile(M_scaled, 0.025, na.rm = TRUE),
    M_upper  = stats::quantile(M_scaled, 0.975, na.rm = TRUE)
  ), by = age_round]

  tmax_summary <- par_draws[, .(
    mean  = mean(tmax),
    lower = stats::quantile(tmax, 0.025),
    upper = stats::quantile(tmax, 0.975)
  )]

  # =========================================================================
  # Plotting
  # =========================================================================

  if (show_progress) message("Generating plot...")

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
    "Estimated tmax: %.1f years (95%% CI: %.1f - %.1f) | Method: %s%s | Growth: %s",
    tmax_summary$mean, tmax_summary$lower, tmax_summary$upper,
    method,
    if (method == "CW" && two_phase) paste0(" (two-phase, ", late_model, ")") else "",
    tools::toTitleCase(gm)
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

  # =========================================================================
  # Return
  # =========================================================================

  if (show_progress) message("Done.")

  list(
    Schedules  = schedules,
    Parameters = par_draws,
    Summary    = summary_dt,
    Plot       = mort_plot
  )
}

#' Monte Carlo Simulation of Survivorship
#'
#' @description
#' Simulates cohort survivorship trajectories using age-specific mortality schedules.
#' Computes cumulative survival, mean age at death, and probabilities of reaching
#' maturity and maximum age.
#'
#' @details
#' For each iteration, the function:
#' \enumerate{
#'   \item Integrates the mortality schedule to get cumulative survival: \eqn{S(t) = \exp(-\int_0^t M(a) da)}
#'   \item Simulates discrete cohort transitions using binomial survival
#'   \item Tracks deaths to compute age-at-death distribution
#' }
#'
#' @param mc_object Output from \code{get_stochastic_mortality()}, containing
#'   \code{$Schedules} and \code{$Parameters}.
#' @param n Starting cohort size. Default 50000.
#' @param n_iter Number of simulation iterations. Default 5000.
#' @param mode Character. Either:
#'   \itemize{
#'     \item \code{"random"} - Each iteration uses a randomly sampled parameter set
#'     \item \code{"per_set"} - Run n_iter simulations for each parameter set
#'   }
#' @param seed Random seed. Default 1234.
#' @param palette Color palette: \code{"synthwave"}, \code{"viridis"}, \code{"okabe"}.
#' @param bar_fill,bar_color Colors for survivorship bars.
#' @param bar_alpha Transparency for bars. Default 0.7.
#' @param error_color Color for error bars. Default "black".
#' @param error_width Width of error bars. Default 0.4.
#' @param vline_color,vline_width,vline_type Mean age-at-death line aesthetics.
#' @param rect_fill,rect_alpha Uncertainty band for mean age-at-death.
#' @param x_breaks Function returning x-axis breaks given max tmax.
#' @param title,subtitle Plot title and subtitle.
#' @param base_size Base font size. Default 11.
#' @param font_family Font family. Default "serif".
#' @param print_plot Print plot? Default TRUE.
#' @param show_progress Show progress messages? Default TRUE.
#'
#' @return A list with:
#' \describe{
#'   \item{Aggregate}{List with overall survival curve, age-at-death, and
#'     survival probabilities to tmat/tmax}
#'   \item{Per_Set}{List with per-parameter-set summaries}
#'   \item{Raw}{data.table of all trajectory data}
#'   \item{Plot}{ggplot2 object}
#' }
#'
#' @examples
#' \dontrun{
#' mort <- get_stochastic_mortality(method = "CW", ...)
#' surv <- simulate_survivorship(mort, n = 1e5, n_iter = 2000)
#' surv$Plot
#' }
#'
#' @import data.table
#' @importFrom pracma trapz
#' @importFrom stats rbinom quantile sd approx
#' @importFrom ggplot2 ggplot aes geom_col geom_errorbar annotate geom_vline
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous labs theme_bw theme
#' @export
simulate_survivorship <- function(
    mc_object,
    n = 50000,
    n_iter = 5000,
    mode = c("random", "per_set"),
    seed = 1234,
    # Plot aesthetics
    palette = c("synthwave", "viridis", "okabe", "plasma", "inferno"),
    bar_fill = NULL,
    bar_color = NULL,
    bar_alpha = 0.7,
    error_color = "black",
    error_width = 0.4,
    vline_color = NULL,
    vline_width = 1,
    vline_type = "solid",
    rect_fill = NULL,
    rect_alpha = 0.4,
    x_breaks = function(tmax) seq(0, ceiling(tmax), by = 5),
    title = NULL,
    subtitle = "Aggregated across simulations",
    base_size = 11,
    font_family = "serif",
    print_plot = TRUE,
    show_progress = TRUE
) {

  mode <- match.arg(mode)
  palette <- match.arg(palette)

  # Validate input
  if (is.null(mc_object$Schedules) || is.null(mc_object$Parameters)) {
    stop("`mc_object` must contain $Schedules and $Parameters from get_stochastic_mortality().")
  }

  Schedules <- data.table::as.data.table(mc_object$Schedules)
  Params <- data.table::as.data.table(mc_object$Parameters)

  if (!"set_id" %in% names(Params)) Params[, set_id := .I]
  if (!"set_id" %in% names(Schedules)) stop("`mc_object$Schedules` must include `set_id`.")

  set_ids <- sort(unique(Params$set_id))
  tmax_max <- max(Params$tmax, na.rm = TRUE)

  # Set title
  if (is.null(title)) {
    title <- sprintf("Cumulative Survivorship (Monte Carlo, %d iterations)", n_iter)
  }

  # Survival simulator (lean, vectorized where possible)
  .simulate_one <- function(mc_df, n) {
    max_age <- max(mc_df$age)
    ages <- 0:floor(max_age)
    n_ages <- length(ages)

    # Cumulative survival via trapezoidal integration
    cum_survival <- numeric(n_ages)
    cum_survival[1] <- 1

    for (k in 2:n_ages) {
      idx <- mc_df$age <= ages[k]
      integral_M <- pracma::trapz(mc_df$age[idx], mc_df$M[idx])
      cum_survival[k] <- exp(-integral_M)
    }

    # Discrete survival simulation
    n_alive <- integer(n_ages)
    n_alive[1] <- n

    for (k in 2:n_ages) {
      prev <- n_alive[k - 1]
      if (prev == 0) break
      p_surv <- if (cum_survival[k - 1] > 0) {
        cum_survival[k] / cum_survival[k - 1]
      } else 0
      n_alive[k] <- stats::rbinom(1, size = prev, prob = max(0, min(1, p_surv)))
    }

    data.table::data.table(
      age = ages,
      survivors = n_alive,
      prop_surv = n_alive / n
    )
  }

  # Build task list
  if (mode == "per_set") {
    task_set <- rep(set_ids, each = n_iter)
    task_iter <- rep(seq_len(n_iter), times = length(set_ids))
  } else {
    set.seed(seed)
    task_set <- sample(set_ids, size = n_iter, replace = TRUE)
    task_iter <- seq_len(n_iter)
  }
  n_tasks <- length(task_set)

  # Split schedules by set_id for fast lookup
  S_split <- split(Schedules, by = "set_id")

  if (show_progress) message("Running survivorship simulations...")

  step <- max(1L, as.integer(n_tasks * 0.1))
  notify <- unique(c(seq(step, n_tasks, by = step), n_tasks))

  results <- vector("list", n_tasks)

  for (i in seq_len(n_tasks)) {
    sid <- task_set[i]
    mc_df <- S_split[[as.character(sid)]]

    out <- .simulate_one(mc_df, n)
    out[, `:=`(iter = task_iter[i], set_id = sid)]
    results[[i]] <- out

    if (show_progress && i %in% notify) {
      message(sprintf("  Progress: %d/%d (%.0f%%)", i, n_tasks, 100 * i / n_tasks))
    }
  }

  sim_df <- data.table::rbindlist(results)

  # Compute summaries
  if (show_progress) message("Computing summaries...")

  # Per-set survival curves
  per_set_summary <- sim_df[, .(
    mean_surv = mean(prop_surv),
    se_surv = stats::sd(prop_surv) / sqrt(.N)
  ), by = .(set_id, age)]
  per_set_summary <- merge(per_set_summary, Params, by = "set_id", all.x = TRUE)

  # Age at death
  deaths_per_iter <- sim_df[, {
    d <- c(survivors[-length(survivors)] - survivors[-1], survivors[length(survivors)])
    d <- pmax(d, 0)
    .(mean_age_death = if (sum(d) > 0) sum(age * d) / sum(d) else NA_real_)
  }, by = .(set_id, iter)]

  age_death_per_set <- deaths_per_iter[, .(
    mean_age = mean(mean_age_death, na.rm = TRUE),
    sd_age = stats::sd(mean_age_death, na.rm = TRUE)
  ), by = set_id]
  age_death_per_set <- merge(age_death_per_set, Params, by = "set_id", all.x = TRUE)

  age_death_agg <- deaths_per_iter[, .(
    mean = mean(mean_age_death, na.rm = TRUE),
    sd = stats::sd(mean_age_death, na.rm = TRUE),
    lower = stats::quantile(mean_age_death, 0.025, na.rm = TRUE),
    upper = stats::quantile(mean_age_death, 0.975, na.rm = TRUE)
  )]

  # Aggregate survival curve
  agg_survival <- sim_df[, .(
    mean_surv = mean(prop_surv),
    se_surv = stats::sd(prop_surv) / sqrt(.N),
    lower = stats::quantile(prop_surv, 0.025),
    upper = stats::quantile(prop_surv, 0.975)
  ), by = age]

  # Survival to tmat and tmax
  if ("tmat" %in% names(Params)) {
    surv_tmat <- per_set_summary[, {
      tmat_val <- unique(tmat)
      if (length(tmat_val) == 1 && !is.na(tmat_val)) {
        .(prob_surv_tmat = stats::approx(x = age, y = mean_surv, xout = tmat_val)$y)
      } else {
        .(prob_surv_tmat = NA_real_)
      }
    }, by = set_id]
    surv_tmat <- merge(surv_tmat, Params, by = "set_id", all.x = TRUE)

    surv_tmat_agg <- surv_tmat[, .(
      mean = mean(prob_surv_tmat, na.rm = TRUE),
      sd = stats::sd(prob_surv_tmat, na.rm = TRUE),
      lower = stats::quantile(prob_surv_tmat, 0.025, na.rm = TRUE),
      upper = stats::quantile(prob_surv_tmat, 0.975, na.rm = TRUE)
    )]
  } else {
    surv_tmat <- NULL
    surv_tmat_agg <- NULL
  }

  surv_tmax <- per_set_summary[, {
    tmax_val <- unique(tmax)
    if (length(tmax_val) == 1 && !is.na(tmax_val)) {
      .(prob_surv_tmax = stats::approx(x = age, y = mean_surv, xout = tmax_val)$y)
    } else {
      .(prob_surv_tmax = NA_real_)
    }
  }, by = set_id]
  surv_tmax <- merge(surv_tmax, Params, by = "set_id", all.x = TRUE)

  surv_tmax_agg <- surv_tmax[, .(
    mean = mean(prob_surv_tmax, na.rm = TRUE),
    sd = stats::sd(prob_surv_tmax, na.rm = TRUE),
    lower = stats::quantile(prob_surv_tmax, 0.025, na.rm = TRUE),
    upper = stats::quantile(prob_surv_tmax, 0.975, na.rm = TRUE)
  )]

  # Plot
  if (show_progress) message("Generating plot...")

  # Get colors
  pal <- vital_palette(n = 6, type = palette)
  if (is.null(bar_fill)) bar_fill <- pal[5]
  if (is.null(bar_color)) bar_color <- pal[4]
  if (is.null(vline_color)) vline_color <- pal[1]
  if (is.null(rect_fill)) rect_fill <- pal[2]

  # Format function for caption
  fmt <- function(x) ifelse(x < 0.01, sprintf("%.2e", x), round(x, 4))

  caption <- sprintf(
    "Mean Age-at-Death: %.1f yrs (95%% CI: %.1f - %.1f)\nProbability of Reaching Max Age: %s (%s - %s)",
    age_death_agg$mean, age_death_agg$lower, age_death_agg$upper,
    fmt(surv_tmax_agg$mean), fmt(surv_tmax_agg$lower), fmt(surv_tmax_agg$upper)
  )

  x_brk <- x_breaks(tmax_max)

  plot_obj <- ggplot2::ggplot(agg_survival, ggplot2::aes(x = age, y = mean_surv)) +
    ggplot2::geom_col(fill = bar_fill, color = bar_color, alpha = bar_alpha) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      width = error_width, color = error_color
    ) +
    ggplot2::annotate(
      "rect",
      xmin = age_death_agg$lower, xmax = age_death_agg$upper,
      ymin = -Inf, ymax = Inf,
      fill = rect_fill, alpha = rect_alpha, color = NA
    ) +
    ggplot2::geom_vline(
      xintercept = age_death_agg$mean,
      color = vline_color, linewidth = vline_width, linetype = vline_type
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_brk,
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::labs(
      title = title, subtitle = subtitle,
      caption = caption,
      x = "Age", y = "Proportion Surviving"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = font_family),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(face = "italic"),
      plot.caption = ggplot2::element_text(hjust = 0),
      axis.title = ggplot2::element_text(face = "bold")
    )

  if (print_plot) print(plot_obj)

  list(
    Aggregate = list(
      Survival = agg_survival,
      Age_of_Death = age_death_agg,
      Survival_to_tmat = surv_tmat_agg,
      Survival_to_tmax = surv_tmax_agg
    ),
    Per_Set = list(
      Survival = per_set_summary,
      Age_of_Death = age_death_per_set,
      Survival_to_tmat = surv_tmat,
      Survival_to_tmax = surv_tmax
    ),
    Raw = sim_df,
    Plot = plot_obj
  )
}

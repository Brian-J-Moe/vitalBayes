#' Precompile Stan models shipped with vitalBayes
#'
#' Compiles one or more Stan models located in `inst/stan/` and caches the
#' executables in the user's cache directory. This avoids compilation during
#' package installation and speeds up first use.
#'
#' @param models Character vector of model base names (no `.stan` extension).
#'   Use `NULL` to compile all models in `inst/stan/`.
#' @param force_recompile Logical. Force recompilation even if cached exe exists.
#' @param cpp_options,stanc_options Lists forwarded to cmdstanr::cmdstan_model().
#' @param quiet Logical. If `FALSE`, prints CmdStan compilation output.
#' @param stop_on_error Logical. If `TRUE`, stops at the first error.
#'
#' @return A `data.table` with columns: `model`, `ok`, `exe_file`, `elapsed_sec`, `error`.
#' @export
precompile_models <- function(models = NULL,
                                 force_recompile = FALSE,
                                 cpp_options = list(),
                                 stanc_options = list(),
                                 quiet = TRUE,
                                 stop_on_error = FALSE) {

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for precompile_models().", call. = FALSE)
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("Package 'cli' is required for precompile_models().", call. = FALSE)
  }

  # Your known model set (from inst/stan)
  all_models <- c(
    "age_at_maturity_single",
    "age_at_maturity_twosex",
    "growth_single_k",
    "growth_single_maturity",
    "growth_twosex_k",
    "growth_twosex_maturity",
    "length_at_birth",
    "length_at_maturity_single",
    "length_at_maturity_twosex"
  )

  if (is.null(models)) {
    models <- all_models
  } else {
    stopifnot(is.character(models))
    bad <- setdiff(models, all_models)
    if (length(bad) > 0L) {
      stop(sprintf(
        "Unknown model(s): %s\nAvailable: %s",
        paste(bad, collapse = ", "),
        paste(all_models, collapse = ", ")
      ), call. = FALSE)
    }
  }

  pb <- cli::cli_progress_bar(
    total = length(models),
    format = "Compiling Stan models {cli::pb_bar} {cli::pb_percent} | {cur}/{total} | {model}",
    clear = FALSE
  )
  on.exit(cli::cli_progress_done(pb), add = TRUE)

  res <- vector("list", length(models))

  for (i in seq_along(models)) {
    m <- models[[i]]
    cli::cli_progress_update(pb, cur = i, total = length(models), model = m)

    t0 <- proc.time()[["elapsed"]]
    out <- tryCatch(
      {
        mod <- .vb_cmdstan_model_cached(
          model = m,
          force_recompile = force_recompile,
          cpp_options = cpp_options,
          stanc_options = stanc_options,
          quiet = quiet
        )
        exe <- tryCatch(mod$exe_file(), error = function(e) NA_character_)

        data.table::data.table(
          model = m,
          ok = TRUE,
          exe_file = exe,
          elapsed_sec = proc.time()[["elapsed"]] - t0,
          error = NA_character_
        )
      },
      error = function(e) {
        dt <- data.table::data.table(
          model = m,
          ok = FALSE,
          exe_file = NA_character_,
          elapsed_sec = proc.time()[["elapsed"]] - t0,
          error = conditionMessage(e)
        )
        if (isTRUE(stop_on_error)) stop(e)
        dt
      }
    )

    res[[i]] <- out
  }

  data.table::rbindlist(res, use.names = TRUE, fill = TRUE)
}

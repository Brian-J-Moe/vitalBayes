

#' Precompile Stan models shipped with vitalBayes
#'
#' Compiles one or more Stan models located in `inst/stan/` and caches the
#' executables in the user's cache directory. This avoids compilation during
#' package installation and speeds up first use.
#'
#' @param models Character vector of model base names (no `.stan` extension).
#'   Use `NULL` to compile all models shipped in `inst/stan/`.
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
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Package 'cmdstanr' is required to compile Stan models.", call. = FALSE)
  }

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

  n <- length(models)
  if (n == 0L) {
    return(data.table::data.table(
      model = character(), ok = logical(), exe_file = character(),
      elapsed_sec = numeric(), error = character()
    ))
  }

  # ---- progress

  n <- length(models)
  pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
  on.exit(close(pb), add = TRUE)


  res <- vector("list", n)

  for (i in seq_len(n)) {
    m <- models[[i]]

    # update progress
    utils::setTxtProgressBar(pb, i)
    message(sprintf("→ [%d/%d] %s", i, n, m))


    t0 <- proc.time()[["elapsed"]]

    out <- tryCatch(
      {
        # Stan source lives in inst/stan/, accessible via system.file()
        stan_file <- system.file("stan", paste0(m, ".stan"), package = "vitalBayes")
        if (!nzchar(stan_file)) {
          stop(sprintf("Stan file not found for model '%s'.", m), call. = FALSE)
        }

        # Choose a stable cache location (per user) so we don't recompile every session.
        cache_dir <- file.path(tools::R_user_dir("vitalBayes", "cache"), "cmdstan_models")
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

        # Cache key: md5 of stan file + options + cmdstanr + CmdStan version
        stan_hash <- unname(tools::md5sum(stan_file))
        cmdstanr_ver <- as.character(utils::packageVersion("cmdstanr"))
        cmdstan_ver <- tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) NA_character_)

        opt_hash <- .vb_hash_object(list(cpp_options = cpp_options, stanc_options = stanc_options))
        key <- paste(stan_hash, cmdstanr_ver, cmdstan_ver, opt_hash, sep = "-")
        exe_ext <- if (.Platform$OS.type == "windows") ".exe" else ""
        exe_dir <- file.path(cache_dir, m, key)
        dir.create(exe_dir, recursive = TRUE, showWarnings = FALSE)
        exe_file <- file.path(exe_dir, paste0(m, exe_ext))

        mod <- cmdstanr::cmdstan_model(
          stan_file = stan_file,
          exe_file = exe_file,
          force_recompile = force_recompile,
          cpp_options = cpp_options,
          stanc_options = stanc_options,
          quiet = quiet
        )

        data.table::data.table(
          model = m,
          ok = TRUE,
          exe_file = tryCatch(mod$exe_file(), error = function(e) exe_file),
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


#' Hash an R object deterministically
#' @noRd
.vb_hash_object <- function(x) {
  raw <- base::serialize(x, connection = NULL, xdr = FALSE)
  tf <- tempfile(fileext = ".bin")
  base::writeBin(raw, tf)
  on.exit(base::unlink(tf), add = TRUE)
  unname(tools::md5sum(tf))
}

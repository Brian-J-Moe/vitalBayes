#' Stan model cache directory (per-user)
#'
#' @return Character scalar path
#' @noRd
.vb_stan_cache_dir <- function() {
  dir <- file.path(tools::R_user_dir("vitalBayes", "cache"), "cmdstan_models")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

#' Locate a Stan program shipped with the package
#'
#' @param model Base name (no extension)
#' @return Path to .stan file
#' @noRd
.vb_stan_file <- function(model) {
  stopifnot(is.character(model), length(model) == 1L, nzchar(model))
  path <- system.file("stan", paste0(model, ".stan"), package = "vitalBayes")
  if (!nzchar(path)) {
    stop(sprintf(
      "Stan file not found for model '%s'. Expected inst/stan/%s.stan",
      model, model
    ), call. = FALSE)
  }
  path
}

#' Hash an R object using base tools (no extra deps)
#'
#' @param x Any R object
#' @return Character scalar hash
#' @noRd
.vb_hash_object <- function(x) {
  tf <- tempfile(fileext = ".bin")
  con <- base::file(tf, open = "wb")
  on.exit({
    base::close(con)
    base::unlink(tf)
  }, add = TRUE)

  base::serialize(x, con, xdr = FALSE)
  unname(tools::md5sum(tf))
}

#' Build a stable cache key for a Stan model build
#'
#' @param stan_file Path to .stan
#' @param cpp_options,stanc_options Lists passed to cmdstanr::cmdstan_model()
#' @return Character scalar cache key
#' @noRd
.vb_stan_cache_key <- function(stan_file, cpp_options, stanc_options) {
  stan_hash <- unname(tools::md5sum(stan_file))

  cmdstanr_ver <- if (requireNamespace("cmdstanr", quietly = TRUE)) {
    as.character(utils::packageVersion("cmdstanr"))
  } else {
    NA_character_
  }

  cmdstan_ver <- tryCatch(
    as.character(cmdstanr::cmdstan_version()),
    error = function(e) NA_character_
  )

  opt_hash <- .vb_hash_object(list(cpp_options = cpp_options, stanc_options = stanc_options))

  paste(stan_hash, cmdstanr_ver, cmdstan_ver, opt_hash, sep = "-")
}

#' Compute the executable path for a cached Stan model
#'
#' @param model Base name (no extension)
#' @param cpp_options,stanc_options Lists passed to cmdstanr::cmdstan_model()
#' @return Character scalar path to exe_file
#' @noRd
.vb_stan_exe_file <- function(model, cpp_options, stanc_options) {
  stan_file <- .vb_stan_file(model)
  key <- .vb_stan_cache_key(stan_file, cpp_options, stanc_options)

  exe_ext <- if (.Platform$OS.type == "windows") ".exe" else ""
  exe_dir <- file.path(.vb_stan_cache_dir(), model, key)
  dir.create(exe_dir, recursive = TRUE, showWarnings = FALSE)

  file.path(exe_dir, paste0(model, exe_ext))
}

#' Get a CmdStanModel with per-user caching (uses cmdstanr::cmdstan_model)
#'
#' @param model Base name (no extension)
#' @param force_recompile Force recompilation even if cached exe exists
#' @param cpp_options,stanc_options Lists passed to cmdstanr::cmdstan_model()
#' @param quiet Passed to cmdstanr::cmdstan_model()
#' @return cmdstanr::CmdStanModel
#' @noRd
.vb_cmdstan_model_cached <- function(model,
                                     force_recompile = FALSE,
                                     cpp_options = list(),
                                     stanc_options = list(),
                                     quiet = TRUE) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("cmdstanr is required to compile Stan models. Install it first.", call. = FALSE)
  }

  stan_file <- .vb_stan_file(model)
  exe_file  <- .vb_stan_exe_file(model, cpp_options, stanc_options)

  cmdstanr::cmdstan_model(
    stan_file = stan_file,
    exe_file = exe_file,
    force_recompile = force_recompile,
    cpp_options = cpp_options,
    stanc_options = stanc_options,
    quiet = quiet
  )
}

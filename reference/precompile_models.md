# Precompile Stan models shipped with vitalBayes

Compiles one or more Stan models located in `inst/stan/` and caches the
executables in the user's cache directory. This avoids compilation
during package installation and speeds up first use.

## Usage

``` r
precompile_models(
  models = NULL,
  force_recompile = FALSE,
  cpp_options = list(),
  stanc_options = list(),
  quiet = TRUE,
  stop_on_error = FALSE
)
```

## Arguments

- models:

  Character vector of model base names (no `.stan` extension). Use
  `NULL` to compile all models in `inst/stan/`.

- force_recompile:

  Logical. Force recompilation even if cached exe exists.

- cpp_options, stanc_options:

  Lists forwarded to cmdstanr::cmdstan_model().

- quiet:

  Logical. If `FALSE`, prints CmdStan compilation output.

- stop_on_error:

  Logical. If `TRUE`, stops at the first error.

## Value

A `data.table` with columns: `model`, `ok`, `exe_file`, `elapsed_sec`,
`error`.

#' Simulated von Bertalanffy growth datasets
#'
#' Three example datasets generated with [simulate_vb_growth_data()] to illustrate
#' common sample-size scenarios (imbalanced, larger balanced, and limited sample sizes)
#' for fitting and testing von Bertalanffy growth and maturity relationships.
#'
#' @format A [data.table::data.table] with 5 columns:
#' \describe{
#'   \item{sex}{Character. Sex for non-embryo observations (\code{"female"} or \code{"male"}); \code{NA} for embryo rows.}
#'   \item{mat}{Integer. Maturity indicator for non-embryo observations (\code{0} = immature, \code{1} = mature); \code{NA} for embryo rows.}
#'   \item{fl}{Numeric. Fork length (cm). For embryos, this is the embryo length measurement.}
#'   \item{age}{Numeric. Age (years) for non-embryo observations; \code{NA} for embryo rows.}
#'   \item{embryo}{Logical. \code{TRUE} for embryo observations; \code{FALSE} otherwise.}
#' }
#'
#' @details
#' All three datasets share the same structure and are intended as lightweight,
#' reproducible examples for package vignettes, unit tests, and demonstrations.
#' Non-embryo observations are simulated from a von Bertalanffy growth model with
#' sex-specific parameter distributions, observation error on length, a truncated-exponential
#' age distribution (older ages increasingly rare), and a sex-specific logistic maturity ogive
#' on observed length. Embryo rows include lengths only and are flagged with \code{embryo = TRUE}.
#'
#' Dataset-specific sizes:
#' \describe{
#'   \item{\code{imbalanced_data}}{Default call to [simulate_vb_growth_data()] (defaults: 150 females, 34 males, 13 embryos).}
#'   \item{\code{growth_data}}{Larger dataset (189 females, 176 males, 26 embryos; \code{seed = 1234}).}
#'   \item{\code{limited_data}}{Small dataset (24 females, 18 males, 5 embryos; \code{seed = 12345}).}
#' }
#'
#' @name simulated_datasets
#' @aliases imbalanced_data growth_data limited_data
#' @docType data
#' @keywords datasets
NULL

#' @rdname simulated_vb_growth_datasets
"imbalanced_data"

#' @rdname simulated_vb_growth_datasets
"growth_data"

#' @rdname simulated_vb_growth_datasets
"limited_data"

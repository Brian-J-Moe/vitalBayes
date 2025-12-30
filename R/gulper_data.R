#' Gulper Shark Life History Data
#'
#' Length, age, and maturity data for the little gulper shark, *Centrophorus uyato*.
#' Contains `NA`s to illustrate that [fit_bayesian_growth()] and
#' [fit_bayesian_maturity()] can handle missing data.
#'
#' @docType data
#'
#' @usage data(gulper_data)
#'
#' @format A `data.table` with 668 rows and 6 variables:
#' \describe{
#'   \item{sex}{Integer. Sex indicator: `1` = female, `2` = male.}
#'   \item{mat}{Integer. Maturity status: `0` = immature, `1` = mature.}
#'   \item{fl}{Numeric. Fork length in cm.}
#'   \item{age1}{Numeric. Age in years estimated from the dorsal fin spine "inner"
#'     growth band region.}
#'   \item{age2}{Numeric. Age in years estimated from the sum of dorsal fin spine
#'     "inner" and "outer" growth band regions.}
#'   \item{embryo}{Logical. If `TRUE`, individual was an embryo rather than
#'     free-swimming.}
#' }
#'
#' @details
#' The sex and maturity variables use integer coding for compatibility with
#' Stan models:
#' \itemize{
#'   \item \code{sex}: 1 = female, 2 = male (auto-detected by vitalBayes functions)
#'   \item \code{mat}: 0 = immature, 1 = mature (standard binary coding)
#' }
#'
#' Age estimates from two band-counting methodologies are provided to allow
#' comparison of growth models under different aging assumptions.
#'
#' @section Data Subsets:
#' Common data subsets for analysis:
#' \itemize{
#'   \item Embryos: \code{gulper_data[embryo == TRUE]}
#'   \item Free-swimming: \code{gulper_data[embryo == FALSE]}
#'   \item Females: \code{gulper_data[sex == 1]}
#'   \item Males: \code{gulper_data[sex == 2]}
#'   \item Mature females: \code{gulper_data[sex == 1 & mat == 1]}
#' }
#'
#' @keywords dataset
#'
#' @references
#' Moe, B.J. (unpublished). Life history of the little gulper shark
#' (*Centrophorus uyato*) in the Gulf of Mexico.
#'
#' @examples
#' data("gulper_data")
#'
#' # Data overview
#' str(gulper_data)
#'
#' # Sample sizes by sex
#' gulper_data[embryo == FALSE, .N, by = sex]
#'
#' # Maturity by sex (free-swimming only)
#' gulper_data[embryo == FALSE, .(
#'   n = .N,
#'   n_mature = sum(mat, na.rm = TRUE),
#'   prop_mature = mean(mat, na.rm = TRUE)
#' ), by = sex]
#'
#' # Extract embryo and free-swimming lengths for birth model
#' embryo_fl <- gulper_data[embryo == TRUE, fl]
#' freeswim_fl <- gulper_data[embryo == FALSE, fl]
#'
#' # Subset females for single-sex analysis
#' female_data <- gulper_data[sex == 1]
#'
"gulper_data"

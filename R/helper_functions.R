# =============================================================================
# vitalBayes Helper Functions
# =============================================================================
# Internal utility functions for data processing, prior conversion, initialization,
# and parameter extraction. These functions are not exported but are used
# throughout the package.
# =============================================================================

# -----------------------------------------------------------------------------
# Sex Coding Standardization
# -----------------------------------------------------------------------------

#' Standardize Sex Coding to Integer Format
#'
#' @description
#' Detects and standardizes sex coding from various formats to integer (1 = female,
#' 2 = male). Supports auto-detection of common conventions in multiple languages
#' used in elasmobranch research (English, Spanish, Portuguese, French, German,
#' Italian). Users can also specify explicit mappings for non-standard codings.
#'
#' @param sex_values Vector of sex values from the data (character, factor, or numeric).
#' @param female Character or numeric. How females are coded in the data. If \code{NULL},
#'   attempts auto-detection.
#' @param male Character or numeric. How males are coded in the data. If \code{NULL},
#'   attempts auto-detection.
#' @param silent Logical. Suppress informational messages? Default \code{FALSE}.
#'
#' @return A list containing:
#'   \item{sex_int}{Integer vector with 1 = female, 2 = male, NA for missing/unknown}
#'   \item{labels}{Named character vector mapping internal codes to original labels}
#'   \item{method}{Character describing how mapping was determined}
#'   \item{n_female}{Count of female observations}
#'   \item{n_male}{Count of male observations}
#'   \item{n_missing}{Count of NA values}
#'
#' @section Supported Auto-Detection Patterns:
#' The function recognizes the following sex coding conventions:
#'
#' \describe{
#'   \item{English}{"F"/"M", "f"/"m", "Female"/"Male", "female"/"male", "FEMALE"/"MALE"}
#'   \item{Spanish}{"H"/"M", "h"/"m", "Hembra"/"Macho", "hembra"/"macho"}
#'   \item{Portuguese}{"F"/"M", "Fêmea"/"Macho", "fêmea"/"macho", "Femea"/"Macho"}
#'   \item{French}{"F"/"M", "Femelle"/"Mâle", "femelle"/"mâle", "Femelle"/"Male"}
#'   \item{German}{"W"/"M", "w"/"m", "Weibchen"/"Männchen", "weibchen"/"männchen"}
#'   \item{Italian}{"F"/"M", "Femmina"/"Maschio", "femmina"/"maschio"}
#'   \item{Japanese}{"メス"/"オス" (mesu/osu), "雌"/"雄" (kanji)}
#'   \item{Chinese/Taiwanese}{"母"/"公", "雌"/"雄", "雌性"/"雄性"}
#'   \item{Russian}{"Ж"/"М", "Самка"/"Самец", "самка"/"самец"}
#'   \item{Symbols}{"♀"/"♂" (Unicode symbols)}
#'   \item{Numeric}{1/2 (female/male)}
#' }
#'
#' @section Note on Ambiguous Patterns:
#' Some patterns like "F"/"M" are shared across languages but always follow the

#' same female/male assignment. The pattern "H"/"M" (Spanish) is distinguished from
#' numeric patterns by character type detection.
#'
#' @examples
#' \dontrun{
#' # Auto-detection examples
#' standardize_sex(c("F", "M", "F", "F", "M"))
#' standardize_sex(c("Hembra", "Macho", "Hembra"))
#' standardize_sex(c("Femelle", "Mâle", "Femelle"))
#' standardize_sex(c(0, 1, 0, 1, 1))
#'
#' # Explicit mapping for non-standard coding
#' standardize_sex(c("A", "B", "A"), female = "A", male = "B")
#' }
#'
#' @noRd
standardize_sex <- function(sex_values, female = NULL, male = NULL, silent = FALSE) {


  # Convert factors to character
  if (is.factor(sex_values)) {
    sex_values <- as.character(sex_values)
  }

  # Get unique non-NA values
  unique_vals <- unique(sex_values[!is.na(sex_values)])

  # Handle case where all values are NA

  if (length(unique_vals) == 0) {
    stop("All sex values are NA. Cannot determine sex coding.", call. = FALSE)
  }

  # Check for more than 2 unique values (excluding NA)
  if (length(unique_vals) > 2) {
    stop(
      "Found more than 2 unique sex values: ",
      paste(unique_vals, collapse = ", "), "\n",
      "Expected exactly 2 values representing female and male.\n",
      "Please clean your data or specify `female` and `male` arguments.",
      call. = FALSE
    )
  }

  # If only one sex present, warn but continue
  if (length(unique_vals) == 1) {
    warning(
      "Only one sex value found: '", unique_vals, "'. ",
      "Two-sex models require both sexes.",
      call. = FALSE
    )
  }

  # ---- User-Specified Mapping ----
  if (!is.null(female) && !is.null(male)) {

    # Validate that specified values exist in data
    if (!(female %in% unique_vals)) {
      stop("Specified female code '", female, "' not found in data. ",
           "Found values: ", paste(unique_vals, collapse = ", "), call. = FALSE)
    }
    if (!(male %in% unique_vals)) {
      stop("Specified male code '", male, "' not found in data. ",
           "Found values: ", paste(unique_vals, collapse = ", "), call. = FALSE)
    }

    sex_int <- ifelse(sex_values == female, 1L,
                      ifelse(sex_values == male, 2L, NA_integer_))

    result <- list(
      sex_int = sex_int,
      labels = c("1" = as.character(female), "2" = as.character(male)),
      method = "user-specified",
      n_female = sum(sex_int == 1L, na.rm = TRUE),
      n_male = sum(sex_int == 2L, na.rm = TRUE),
      n_missing = sum(is.na(sex_int))
    )

    if (!silent) {
      message("Sex coding (user-specified): Female = '", female, "', Male = '", male, "'")
    }

    return(result)
  }

  # If only one argument provided, require both
  if (xor(is.null(female), is.null(male))) {
    stop("Both `female` and `male` must be specified together, or neither (for auto-detection).",
         call. = FALSE)
  }

  # ---- Auto-Detection ----

  # Define patterns: each is c(female_code, male_code)
  # Order matters for matching - more specific patterns first
  patterns <- list(
    # === ENGLISH ===
    c("Female", "Male"),
    c("female", "male"),
    c("FEMALE", "MALE"),
    c("F", "M"),
    c("f", "m"),

    # === SPANISH ===
    # "Hembra" = female, "Macho" = male
    c("Hembra", "Macho"),
    c("hembra", "macho"),
    c("HEMBRA", "MACHO"),
    c("H", "M"),   # Spanish abbreviation (note: same M for male as English)
    c("h", "m"),

    # === PORTUGUESE ===
    # "Fêmea" = female, "Macho" = male
    c("Fêmea", "Macho"),
    c("fêmea", "macho"),
    c("FÊMEA", "MACHO"),
    c("Femea", "Macho"),   # Without accent (common in datasets)
    c("femea", "macho"),

    # === FRENCH ===
    # "Femelle" = female, "Mâle" = male
    c("Femelle", "Mâle"),
    c("femelle", "mâle"),
    c("FEMELLE", "MÂLE"),
    c("Femelle", "Male"),  # Without accent on male (common)
    c("femelle", "male"),

    # === GERMAN ===
    # "Weibchen" = female, "Männchen" = male
    c("Weibchen", "Männchen"),
    c("weibchen", "männchen"),
    c("WEIBCHEN", "MÄNNCHEN"),
    c("Weibchen", "Mannchen"),  # Without umlaut
    c("weibchen", "mannchen"),
    c("W", "M"),   # German abbreviation
    c("w", "m"),

    # === ITALIAN ===
    # "Femmina" = female, "Maschio" = male
    c("Femmina", "Maschio"),
    c("femmina", "maschio"),
    c("FEMMINA", "MASCHIO"),

    # === JAPANESE ===
    # "メス" (mesu) = female, "オス" (osu) = male
    c("メス", "オス"),
    c("雌", "雄"),        # Kanji: female/male
    c("♀", "♂"),          # Unicode symbols (used across languages)

    # === TAIWANESE/MANDARIN CHINESE ===
    # "雌" (cí) = female, "雄" (xióng) = male (same kanji as Japanese)
    # "母" (mǔ) = female/mother, "公" (gōng) = male (common for animals)
    c("母", "公"),
    c("雌性", "雄性"),    # With -性 suffix meaning "sex/gender"
    c("女", "男"),        # Less common for animals but sometimes used

    # === RUSSIAN ===
    # "Самка" (samka) = female, "Самец" (samets) = male
    c("Самка", "Самец"),
    c("самка", "самец"),
    c("САМКА", "САМЕЦ"),
    c("С", "С"),          # Abbreviation problematic (both start with С)
    c("Ж", "М"),          # Женский/Мужской abbreviation
    c("ж", "м"),

    # === NUMERIC PATTERNS ===
    # Convention: 1 = female, 2 = male (common in fisheries)
    c(1, 2)
  )

  # Add language labels for informative messages
  pattern_languages <- c(
    rep("English", 5),
    rep("Spanish", 5),
    rep("Portuguese", 5),
    rep("French", 5),
    rep("German", 7),
    rep("Italian", 3),
    rep("Japanese", 3),
    rep("Chinese/Taiwanese", 3),
    rep("Russian", 5),
    "Numeric"
  )

  # Try to match patterns
  for (i in seq_along(patterns)) {
    pat <- patterns[[i]]

    # Check if unique values match this pattern (order-independent)
    if (setequal(unique_vals, pat)) {

      female_code <- pat[1]
      male_code <- pat[2]

      sex_int <- ifelse(sex_values == female_code, 1L,
                        ifelse(sex_values == male_code, 2L, NA_integer_))

      result <- list(
        sex_int = sex_int,
        labels = c("1" = as.character(female_code), "2" = as.character(male_code)),
        method = paste0("auto-detected (", pattern_languages[i], ")"),
        n_female = sum(sex_int == 1L, na.rm = TRUE),
        n_male = sum(sex_int == 2L, na.rm = TRUE),
        n_missing = sum(is.na(sex_int))
      )

      if (!silent) {
        message(
          "Sex coding auto-detected (", pattern_languages[i], "): ",
          "Female = '", female_code, "', Male = '", male_code, "'"
        )
      }

      return(result)
    }
  }

  # ---- No Match Found ----
  stop(
    "Could not auto-detect sex coding.\n",
    "Found values: ", paste(unique_vals, collapse = ", "), "\n\n",
    "Please specify the `female` and `male` arguments explicitly.\n",
    "Example: fit_bayesian_maturity(..., female = '", unique_vals[1],
    "', male = '", unique_vals[2], "')\n\n",
    "Supported auto-detection patterns include:\n",
    "  English:    F/M, Female/Male\n",
    "  Spanish:    H/M, Hembra/Macho\n",
    "  Portuguese: F/M, Fêmea/Macho, Femea/Macho\n",
    "  French:     F/M, Femelle/Mâle, Femelle/Male\n",
    "  German:     W/M, Weibchen/Männchen\n",
    "  Italian:    F/M, Femmina/Maschio\n",
    "  Japanese:   メス/オス, 雌/雄\n",
    "  Chinese:    母/公, 雌/雄, 雌性/雄性\n",
    "  Russian:    Ж/М, Самка/Самец\n",
    "  Symbols:    ♀/♂\n",
    "  Numeric:    1/2",
    call. = FALSE
  )
}


#' Get Sex Label from Standardized Coding
#'
#' @description
#' Returns human-readable sex labels from the standardized sex result.
#'
#' @param sex_result Output from \code{standardize_sex()}.
#' @param code Integer code (1 or 2) or "both".
#' @param style Character. Output style: "original" uses the original data coding,
#'   "english" uses "Female"/"Male".
#'
#' @return Character string or vector of labels.
#'
#' @noRd
.get_sex_label <- function(sex_result, code = "both", style = c("english", "original")) {

  style <- match.arg(style)

  if (style == "original") {
    labels <- sex_result$labels
  } else {
    labels <- c("1" = "Female", "2" = "Male")
  }

  if (code == "both") {
    return(labels)
  } else if (code == 1) {
    return(unname(labels["1"]))
  } else if (code == 2) {
    return(unname(labels["2"]))
  } else {
    stop("Invalid code. Use 1, 2, or 'both'.", call. = FALSE)
  }
}


# -----------------------------------------------------------------------------
# Truncated Normal Sampling
# -----------------------------------------------------------------------------

#' Sample from Truncated Normal Distribution
#'
#' @description
#' Generates samples from a truncated normal distribution using inverse CDF
#' method. Supports vectorized parameters.
#'
#' @param n Integer. Number of samples per distribution.
#' @param mu Numeric vector. Mean(s).
#' @param sigma Numeric vector. SD(s). Must be positive.
#' @param lower Numeric vector. Lower bound(s). Default 0.
#' @param upper Numeric vector. Upper bound(s). Default Inf.
#'
#' @return Matrix of samples.
#'
#' @noRd
.rtruncnorm <- function(n, mu, sigma, lower = 0, upper = Inf) {

  stopifnot(n > 0)

  k <- max(length(mu), length(sigma), length(lower), length(upper))
  mu    <- rep_len(mu, k)
  sigma <- rep_len(sigma, k)
  lower <- rep_len(lower, k)
  upper <- rep_len(upper, k)

  if (any(sigma <= 0)) stop("All sigma values must be > 0")
  if (any(lower >= upper)) stop("Each 'lower' must be < corresponding 'upper'")

  p_lower <- stats::pnorm(lower, mean = mu, sd = sigma)
  p_upper <- stats::pnorm(upper, mean = mu, sd = sigma)

  U <- matrix(stats::runif(n * k), nrow = n, ncol = k)
  P <- sweep(U, 2, (p_upper - p_lower), `*`)
  P <- sweep(P, 2, p_lower, `+`)

  stats::qnorm(P,
               mean = matrix(mu, n, k, byrow = TRUE),
               sd = matrix(sigma, n, k, byrow = TRUE))
}


# -----------------------------------------------------------------------------
# Argument Resolution
# -----------------------------------------------------------------------------

#' Resolve Column Name from Data
#'
#' @param arg The substituted argument.
#' @param data A data.frame/data.table.
#' @param type Expected type: "numeric" or "character".
#'
#' @return The resolved vector.
#'
#' @noRd
.resolve_arg <- function(arg, data, type = c("numeric", "character")) {

  type <- match.arg(type)
  arg_name <- as.character(arg)

  if (is.null(data)) {
    stop("Must provide 'data' when using column names.", call. = FALSE)
  }

  if (!(arg_name %in% names(data))) {
    stop("Column '", arg_name, "' not found in data.", call. = FALSE)
  }

  vec <- data[[arg_name]]

  if (type == "character") {
    if (!is.character(vec) && !is.factor(vec)) {
      stop("Column '", arg_name, "' must be character or factor type.", call. = FALSE)
    }
    vec <- as.character(vec)
  } else {
    if (!is.numeric(vec)) {
      stop("Column '", arg_name, "' must be numeric type.", call. = FALSE)
    }
  }

  return(vec)
}


#' Resolve Argument as Vector or Column Reference
#'
#' @param arg_expr The substituted expression.
#' @param data Optional data containing columns.
#' @param type Expected type.
#' @param name Argument name for errors.
#'
#' @return The resolved vector.
#'
#' @noRd
.resolve_or_vector <- function(arg_expr, data, type, name) {

  arg_value <- tryCatch(
    eval(arg_expr, envir = parent.frame(2)),
    error = function(e) NULL
  )

  if (!is.null(arg_value) && length(arg_value) > 1) {
    if (type == "numeric" && is.numeric(arg_value)) {
      return(arg_value)
    } else if (type == "character" && (is.character(arg_value) || is.factor(arg_value))) {
      return(as.character(arg_value))
    }
  }

  .resolve_arg(arg_expr, data, type)
}


# -----------------------------------------------------------------------------
# Prior Conversion Utilities
# -----------------------------------------------------------------------------

#' Convert Natural Scale Prior to Log Scale
#'
#' @param mean_nat Prior mean on natural scale.
#' @param sd_nat Prior SD on natural scale.
#' @param n_sim Number of simulations. Default 1000.
#'
#' @return Named list with log_mean and log_sd.
#'
#' @noRd
.natural_to_log_prior <- function(mean_nat, sd_nat, n_sim = 1000) {

  samples <- .rtruncnorm(n_sim, mean_nat, sd_nat, lower = 0)
  log_samples <- log(samples)

  list(
    log_mean = mean(log_samples),
    log_sd = stats::sd(log_samples)
  )
}


#' Compute Prior from Mean and CV
#'
#' @param mean_nat Prior mean on natural scale.
#' @param cv Coefficient of variation.
#'
#' @return Named list with log_mean and log_sd.
#'
#' @noRd
.cv_to_log_prior <- function(mean_nat, cv) {
  sd_nat <- mean_nat * cv
  .natural_to_log_prior(mean_nat, sd_nat)
}


#' Standardize Prior Input Format
#'
#' @param prior Prior specification.
#' @param n_groups Number of groups (1 or 2).
#'
#' @return Matrix with n_groups rows and 2 columns (mean, sd).
#'
#' @noRd
.standardize_prior <- function(prior, n_groups) {

  if (is.null(prior)) return(NULL)

  if (is.matrix(prior) || is.data.frame(prior)) {
    prior <- as.matrix(prior)
    if (ncol(prior) != 2) {
      stop("Prior matrix must have 2 columns (mean, sd).", call. = FALSE)
    }
    if (nrow(prior) == 1 && n_groups > 1) {
      prior <- matrix(rep(prior, n_groups), nrow = n_groups, byrow = TRUE)
    }
    return(prior)
  }

  if (is.list(prior) && all(c("mean", "sd") %in% names(prior))) {
    prior <- c(prior$mean, prior$sd)
  }

  if (is.numeric(prior)) {
    if (length(prior) == 2) {
      return(matrix(rep(prior, n_groups), nrow = n_groups, byrow = TRUE))
    } else if (length(prior) == n_groups * 2) {
      return(matrix(prior, nrow = n_groups, ncol = 2, byrow = TRUE))
    }
  }

  stop("Invalid prior format. Provide c(mean, sd) or matrix with columns (mean, sd).",
       call. = FALSE)
}


#' Standardize Scalar to Vector
#'
#' @param x Numeric scalar or vector.
#' @param n Required length.
#'
#' @return Numeric vector of length n.
#'
#' @noRd
.standardize_scalar <- function(x, n) {
  if (length(x) == 1) {
    return(rep(x, n))
  } else if (length(x) == n) {
    return(x)
  } else {
    stop("Value must be length 1 or ", n, ".", call. = FALSE)
  }
}


# -----------------------------------------------------------------------------
# Data-Driven Prior and Initialization Estimation
# -----------------------------------------------------------------------------

#' Estimate Midpoint for Maturity/Birth Models
#'
#' @description
#' Computes the midpoint between the smallest "positive" individual and the
#' largest "negative" individual for binary response models.
#'
#' @param x_vec Numeric vector. Predictor values.
#' @param y_vec Integer vector. Binary response (0/1).
#'
#' @return Numeric. Midpoint estimate.
#'
#' @noRd
.estimate_midpoint <- function(x_vec, y_vec) {

  complete <- !is.na(x_vec) & !is.na(y_vec)
  x_vec <- x_vec[complete]
  y_vec <- y_vec[complete]

  if (length(x_vec) == 0) {
    stop("No complete cases for midpoint estimation.", call. = FALSE)
  }

  min_positive <- min(x_vec[y_vec == 1], na.rm = TRUE)
  max_negative <- max(x_vec[y_vec == 0], na.rm = TRUE)

  if (is.infinite(min_positive) || is.infinite(max_negative)) {
    warning("One class missing. Using median as midpoint estimate.", call. = FALSE)
    return(stats::median(x_vec))
  }

  (min_positive + max_negative) / 2
}


#' Estimate k from Observed Length-at-Age Data
#'
#' @description
#' Estimates the growth rate k by solving the growth equation for each
#' observation and averaging. Excludes endpoints.
#'
#' @param lt_vec Numeric vector. Observed lengths.
#' @param age_vec Numeric vector. Observed ages.
#' @param Linf Numeric. Estimated asymptotic length.
#' @param L0 Numeric. Estimated length at birth.
#' @param model Character. Growth model: "v", "g", or "l".
#'
#' @return Numeric. Estimated k value.
#'
#' @noRd
.estimate_k_from_data <- function(lt_vec, age_vec, Linf, L0, model = "v") {

  valid <- age_vec > 0 & lt_vec > L0 & lt_vec < Linf
  lt_use <- lt_vec[valid]
  age_use <- age_vec[valid]

  if (length(lt_use) < 3) {
    warning("Insufficient valid observations for k estimation. Using default 0.1.",
            call. = FALSE)
    return(0.1)
  }

  if (model == "v") {
    k_vec <- -log((Linf - lt_use) / (Linf - L0)) / age_use
  } else if (model == "g") {
    r0 <- log(L0 / Linf)
    k_vec <- -log(log(lt_use / Linf) / r0) / age_use
  } else if (model == "l") {
    k_vec <- -log((Linf / lt_use - 1) / (Linf / L0 - 1)) / age_use
  } else {
    stop("Unknown model: ", model, call. = FALSE)
  }

  k_vec <- k_vec[is.finite(k_vec) & k_vec > 0]

  if (length(k_vec) == 0) {
    warning("Could not estimate k from data. Using default 0.1.", call. = FALSE)
    return(0.1)
  }

  stats::median(k_vec)
}


#' Gulland-Holt k Estimation for Initialization
#'
#' @description
#' Estimates k using the Gulland-Holt linearization method.
#'
#' @param lt_vec Numeric vector. Observed lengths.
#' @param age_vec Numeric vector. Observed ages.
#'
#' @return Numeric. Estimated k value, or 0.1 if estimation fails.
#'
#' @noRd
.gulland_holt_k <- function(lt_vec, age_vec) {

  ord <- order(age_vec)
  lt_sorted <- lt_vec[ord]
  age_sorted <- age_vec[ord]

  n <- length(lt_sorted)
  if (n < 5) return(0.1)

  dL <- diff(lt_sorted)
  dt <- diff(age_sorted)
  L_mean <- (lt_sorted[-n] + lt_sorted[-1]) / 2

  valid <- dt > 0 & dL > 0
  if (sum(valid) < 3) return(0.1)

  growth_rate <- dL[valid] / dt[valid]
  L_mid <- L_mean[valid]

  fit <- tryCatch(
    stats::lm(growth_rate ~ L_mid),
    error = function(e) NULL
  )

  if (is.null(fit)) return(0.1)

  k_est <- -stats::coef(fit)[2]

  if (is.na(k_est) || k_est <= 0) return(0.1)

  as.numeric(k_est)
}


# -----------------------------------------------------------------------------
# Parameter Extraction from Stan Fits
# -----------------------------------------------------------------------------

#' Extract Prior Parameters from Maturity Fit
#'
#' @param fit A CmdStanMCMC object from fit_bayesian_maturity().
#' @param param Parameter name: "L50" or "t50".
#' @param multiplier Multiplier for SD. Default 1.5.
#'
#' @return List with natural and log scale parameters for each sex.
#'
#' @noRd
.extract_maturity_prior <- function(fit, param, multiplier = 1.5) {

  draws <- fit$draws(param, format = "matrix")
  n_sex <- ncol(draws)

  result <- list()

  for (s in 1:n_sex) {
    samples <- draws[, s]

    nat_mean <- stats::median(samples)
    nat_sd <- stats::sd(samples) * multiplier

    log_mean <- mean(log(samples))
    log_sd <- stats::sd(log(samples)) * multiplier

    result[[s]] <- list(
      natural = c(mean = nat_mean, sd = nat_sd),
      log = c(mean = log_mean, sd = log_sd)
    )
  }

  class(result) <- c("vitalBayes_prior", "list")
  attr(result, "param") <- param
  attr(result, "n_sex") <- n_sex

  return(result)
}


#' Extract Prior Parameters from Birth Fit
#'
#' @param fit A CmdStanMCMC object from fit_bayesian_birth().
#' @param multiplier Multiplier for SD. Default 1.5.
#'
#' @return List with natural and log scale parameters.
#'
#' @noRd
.extract_birth_prior <- function(fit, multiplier = 1.5) {

  samples <- as.vector(fit$draws("b50", format = "matrix"))

  nat_mean <- stats::median(samples)
  nat_sd <- stats::sd(samples) * multiplier

  log_mean <- mean(log(samples))
  log_sd <- stats::sd(log(samples)) * multiplier

  list(
    natural = c(mean = nat_mean, sd = nat_sd),
    log = c(mean = log_mean, sd = log_sd)
  )
}


# -----------------------------------------------------------------------------
# Data Validation
# -----------------------------------------------------------------------------

#' Validate Growth Model Data
#'
#' @param lt_vec Numeric vector of lengths.
#' @param age_vec Numeric vector of ages.
#' @param sex_vec Optional character vector of sex.
#'
#' @return NULL (invisibly).
#'
#' @noRd
.validate_growth_data <- function(lt_vec, age_vec, sex_vec = NULL) {

  if (any(lt_vec <= 0, na.rm = TRUE)) {
    stop("All length values must be positive.", call. = FALSE)
  }

  if (any(age_vec < 0, na.rm = TRUE)) {
    stop("All age values must be non-negative.", call. = FALSE)
  }

  if (length(lt_vec) != length(age_vec)) {
    stop("Length and age vectors must have the same length.", call. = FALSE)
  }

  if (!is.null(sex_vec)) {
    if (length(sex_vec) != length(lt_vec)) {
      stop("Sex vector must have the same length as length/age vectors.",
           call. = FALSE)
    }

    unique_sex <- unique(tolower(as.character(sex_vec[!is.na(sex_vec)])))

    if (all(unique_sex %in% c("1", "2"))) {
      has_female <- any(unique_sex == "1")
      has_male   <- any(unique_sex == "2")
    } else {
      has_female <- any(grepl("f", unique_sex))
      has_male   <- any(grepl("m", unique_sex))
    }

    if (!has_female || !has_male) {
      warning("Expected both sexes. Found: ", paste(unique_sex, collapse = ", "),
              call. = FALSE)
    }
  }

  invisible(NULL)
}


#' Check Sample Size and Warn
#'
#' @param newdat A data.table with processed data.
#' @param sex_vec Original sex vector (or NULL).
#' @param k_based Whether using k-based parameterization.
#'
#' @return NULL (invisibly).
#'
#' @noRd
.check_sample_size <- function(newdat, sex_vec, k_based) {

  n_total <- nrow(newdat)
  n_params <- if (k_based) 4 else 5

  if (n_total < n_params * 10) {
    warning(
      "Small sample size (n=", n_total, ") relative to number of parameters (",
      n_params, ").\n",
      "Consider using informative priors or partial pooling.",
      call. = FALSE
    )
  }

  if (!is.null(sex_vec)) {
    n_by_sex <- newdat[, .N, by = sex_code]

    for (i in 1:nrow(n_by_sex)) {
      if (n_by_sex$N[i] < n_params * 5) {
        sex_label <- if (n_by_sex$sex_code[i] == 1) "females" else "males"
        warning(
          "Small sample size for ", sex_label, " (n=", n_by_sex$N[i], ").\n",
          "Consider using partial pooling (use_pooling = TRUE).",
          call. = FALSE
        )
      }
    }

    ratio <- max(n_by_sex$N) / min(n_by_sex$N)
    if (ratio > 3) {
      warning(
        "Imbalanced sex ratio (", round(ratio, 1), ":1). Partial pooling recommended.",
        call. = FALSE
      )
    }
  }

  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Growth Model Computations
# -----------------------------------------------------------------------------

#' Compute Growth Curve Predictions
#'
#' @param Linf Asymptotic length.
#' @param L0 Length at birth.
#' @param k Growth rate.
#' @param age_vec Ages at which to predict.
#' @param which_model 1=VBGM, 2=Gompertz, 3=Logistic.
#'
#' @return Numeric vector of predicted lengths.
#'
#' @noRd
.growth_curve <- function(Linf, L0, k, age_vec, which_model = 1L) {

  if (which_model == 1L) {
    mu <- Linf - (Linf - L0) * exp(-k * age_vec)
  } else if (which_model == 2L) {
    r0 <- log(Linf / L0)
    mu <- Linf * exp(-r0 * exp(-k * age_vec))
  } else {
    A <- exp(-k * age_vec)
    mu <- Linf / (1 + A * (Linf / L0 - 1))
  }

  return(mu)
}


#' Compute k from Maturity Parameters
#'
#' @param Linf Asymptotic length.
#' @param L0 Length at birth.
#' @param Lmat Length at maturity.
#' @param tmat Age at maturity.
#' @param which_model 1=VBGM, 2=Gompertz, 3=Logistic.
#'
#' @return Derived growth rate k.
#'
#' @noRd
.derive_k_from_maturity <- function(Linf, L0, Lmat, tmat, which_model = 1L) {

  if (which_model == 1L) {
    k <- (1 / tmat) * log((Linf - L0) / (Linf - Lmat))
  } else if (which_model == 2L) {
    r0 <- log(Linf / L0)
    rmt <- log(Lmat / Linf)
    k <- (1 / tmat) * log(-r0 / rmt)
  } else {
    num <- Lmat * (Linf - L0)
    den <- L0 * (Linf - Lmat)
    k <- (1 / tmat) * log(num / den)
  }

  return(k)
}


# -----------------------------------------------------------------------------
# Output Formatting
# -----------------------------------------------------------------------------

#' Print Growth Model Summary
#'
#' @param fit A CmdStanMCMC object.
#' @param k_based Whether model is k-based.
#' @param is_twosex Whether model is two-sex.
#' @param use_pooling Whether pooling was used.
#' @param robust Whether robust errors were used.
#'
#' @return NULL (invisibly).
#'
#' @noRd
.print_growth_summary <- function(fit, k_based, is_twosex, use_pooling, robust) {

  message("\n", paste(rep("=", 60), collapse = ""))
  message("Growth Model Posterior Summary")
  message(paste(rep("=", 60), collapse = ""))

  core_params <- c("Linf", "L0", "k", "sigma")
  if (!k_based) {
    core_params <- c("Linf", "L0", "Lmat", "tmat", "k", "sigma")
  }
  if (robust) {
    core_params <- c(core_params, "nu")
  }

  message("\nCore Parameters:")
  print(fit$summary(variables = core_params))

  if (is_twosex) {
    message("\nSex Differences (Female - Male):")
    diff_params <- if (k_based) {
      c("Linf_diff", "L0_diff", "k_diff")
    } else {
      c("Linf_diff", "L0_diff", "Lmat_diff", "tmat_diff", "k_diff")
    }
    print(fit$summary(variables = diff_params))
  }

  message("\nModel Fit Statistics:")
  ppc_vars <- c("rmse_f", "rmse_m", "n_in_CI_f", "n_in_CI_m")
  if (!is_twosex) {
    ppc_vars <- c("rmse", "mae", "n_in_CI")
  }

  available <- fit$metadata()$stan_variables
  ppc_vars <- ppc_vars[ppc_vars %in% available]

  if (length(ppc_vars) > 0) {
    print(fit$summary(variables = ppc_vars))
  }

  invisible(NULL)
}


#' Print Maturity Model Summary
#'
#' @param fit A CmdStanMCMC object.
#' @param type "length" or "age".
#' @param is_twosex Whether two-sex.
#' @param use_pooling Whether pooling was used.
#'
#' @return NULL (invisibly).
#'
#' @noRd
.print_maturity_summary <- function(fit, type, is_twosex, use_pooling) {

  param_name <- if (type == "length") "L50" else "t50"

  message("\n", paste(rep("=", 60), collapse = ""))
  message(tools::toTitleCase(type), "-at-Maturity Posterior Summary")
  message(paste(rep("=", 60), collapse = ""))

  message("\nCore Parameters:")
  print(fit$summary(variables = c(param_name, "slope")))

  if (is_twosex) {
    diff_param <- paste0(param_name, "_diff")
    message("\nSex Difference (Female - Male):")
    print(fit$summary(variables = diff_param))

    message("\nTransition Zones by Sex:")
    print(fit$summary(variables = "transition_width"))
  } else {
    message("\nDerived Quantities:")
    x05 <- if (type == "length") "L05" else "t05"
    x95 <- if (type == "length") "L95" else "t95"
    print(fit$summary(variables = c(x05, x95, "transition_width")))
  }

  message("\nPosterior Predictive Checks:")
  ppc_vars <- c("prop_correct_rep")
  if (is_twosex) {
    ppc_vars <- c(ppc_vars, "mean_p_mature_f", "mean_p_mature_m",
                  "mean_p_immature_f", "mean_p_immature_m")
  } else {
    ppc_vars <- c(ppc_vars, "mean_p_mature", "mean_p_immature")
  }
  print(fit$summary(variables = ppc_vars))

  invisible(NULL)
}


#' Compute non-centered raw values that reproduce target log-scale values
#'
#' @param target_log Numeric vector of target values on log scale (e.g., log(mean_Linf_nat)).
#' @param mu Numeric scalar population mean on log scale.
#' @param tau Numeric scalar population SD (must be > 0).
#' @param eps Small positive number to avoid division by 0.
#'
#' @return Numeric vector of same length as target_log.
#' @noRd
.safe_raw_from_target <- function(target_log, mu, tau, eps = 1e-6) {
  tau <- max(tau, eps)
  (target_log - mu) / tau
}

#' Enforce L0 < Lmat < Linf and Linf >= Linf_lower for stable maturity-based inits
#'
#' @param L0,Lmat,Linf Numeric vectors (length 1 or 2).
#' @param Linf_lower Numeric vector (length 1 or 2) lower bound(s) for Linf.
#' @param eps Margin to keep parameters away from boundaries.
#'
#' @return List with L0, Lmat, Linf adjusted.
#' @noRd
.safe_maturity_order <- function(L0, Lmat, Linf, Linf_lower, eps = 0.5) {
  L0 <- pmax(L0, eps)
  Linf <- pmax(Linf, Linf_lower + eps)

  # keep Lmat between L0 and Linf
  Lmat <- pmax(Lmat, L0 + eps)
  Lmat <- pmin(Lmat, Linf - eps)

  list(L0 = L0, Lmat = Lmat, Linf = Linf)
}


#' Clamp values to be positive and finite
#'
#' @param x Numeric vector.
#' @param eps Small positive lower bound.
#' @param name Character label used in warnings.
#'
#' @return Numeric vector with non-finite values replaced by eps, and values < eps set to eps.
#' @noRd
.safe_pos <- function(x, eps = 1e-6, name = "value") {
  x <- as.numeric(x)
  bad <- !is.finite(x) | is.na(x) | (x < eps)
  if (any(bad)) {
    warning("Non-finite or non-positive ", name, " in init; clamping/replacing with ", eps, call. = FALSE)
    x[bad] <- eps
  }
  x
}

#' Make Linf and L0 consistent with constraints
#'
#' Ensures: Linf >= Linf_lower + margin, L0 >= eps, and Linf >= L0 + margin.
#'
#' @param L0,Linf Numeric vectors (length 1 or 2).
#' @param Linf_lower Numeric vector of Linf lower bounds (length 1 or 2).
#' @param margin Positive margin to keep away from boundaries.
#' @param eps Small lower bound for positivity.
#'
#' @return List with adjusted L0 and Linf.
#' @noRd
.safe_L0_Linf <- function(L0, Linf, Linf_lower, margin = 0.5, eps = 1e-6) {
  L0   <- .safe_pos(L0,   eps = eps, name = "L0")
  Linf <- .safe_pos(Linf, eps = eps, name = "Linf")

  # enforce Linf above lower bound + margin
  Linf <- pmax(Linf, Linf_lower + margin)

  # enforce Linf above L0 + margin
  Linf <- pmax(Linf, L0 + margin)

  list(L0 = L0, Linf = Linf)
}

#' Compute non-centered raw values to reproduce target log-scale values
#'
#' @param target_log Numeric vector of targets on log scale.
#' @param mu Population mean on log scale.
#' @param tau Population SD (> 0).
#' @param eps Small positive number to avoid division by 0.
#'
#' @return Numeric vector of raws.
#' @noRd
.safe_raw_from_target <- function(target_log, mu, tau, eps = 1e-6) {
  tau <- max(as.numeric(tau), eps)
  (as.numeric(target_log) - as.numeric(mu)) / tau
}


#' Compute Von Bertalanffy-Equivalent Growth Coefficient
#'
#' @description
#' Derives the von Bertalanffy growth coefficient \eqn{k} from biological
#' milestones \eqn{(L_\infty, L_0, L_{mat}, t_{mat})}. This enables Chen-Watanabe
#' mortality estimation using posteriors from \emph{any} growth model (von
#' Bertalanffy, Gompertz, or Logistic).
#'
#' @details
#' The key insight is that while the three growth models use different functional
#' forms and produce different numerical \eqn{k} values, they all estimate the
#' same underlying biological quantities: asymptotic length, birth size, and
#' maturity milestones. The VB-equivalent \eqn{k} is computed as:
#'
#' \deqn{k_{VB}^{equiv} = \frac{1}{t_{mat}} \ln\left(\frac{L_\infty - L_0}{L_\infty - L_{mat}}\right)}
#'
#' This is the growth coefficient that would produce a von Bertalanffy curve
#' passing through the biological milestones \eqn{(0, L_0)} and
#' \eqn{(t_{mat}, L_{mat})} with asymptote \eqn{L_\infty}.
#'
#' When the input growth fit is von Bertalanffy with maturity-based
#' parameterization, this computation exactly reproduces the fitted \eqn{k}.
#' When the input is Gompertz or Logistic, it produces the VB-equivalent
#' \eqn{k} that encodes the same biological growth information.
#'
#' The von Bertalanffy model tends to produce unstable \eqn{L_\infty} estimates
#' when data are sparse at older ages - a common situation in elasmobranch
#' research. Gompertz and Logistic models often provide more reliable fits in
#' these cases. By deriving VB-equivalent \eqn{k} from biological milestones,
#' users can select the growth model that best fits their data while still
#' using Chen-Watanabe mortality estimation and maintaining theoretical
#' coherence (since CW was derived under VB assumptions).
#'
#' @param Linf Numeric vector. Asymptotic length posterior draws.
#' @param L0 Numeric vector. Length at birth posterior draws.
#' @param Lmat Numeric vector. Length at maturity posterior draws.
#' @param tmat Numeric vector. Age at maturity posterior draws.
#' @param warn Logical. If \code{TRUE} (default), warns when draws produce
#'   invalid \eqn{k} values.
#'
#' @return Numeric vector of VB-equivalent \eqn{k} values. Invalid values
#'   (from non-positive arguments to log) are returned as \code{NA}.
#'
#' @examples
#' \dontrun{
#' # Extract from any growth model posterior
#' draws <- extract_growth_parameters(growth_fit, sex = 1)
#' k_vb <- compute_k_vb_equivalent(
#'   Linf = draws$Linf,
#'   L0   = draws$L0,
#'   Lmat = draws$Lmat,
#'   tmat = draws$tmat
#' )
#'
#' # Direct specification for sensitivity analysis
#' k_vb <- compute_k_vb_equivalent(
#'   Linf = rnorm(1000, 100, 5),
#'   L0   = rnorm(1000, 25, 2),
#'   Lmat = rnorm(1000, 70, 3),
#'   tmat = rnorm(1000, 10, 1)
#' )
#' }
#'
#' @seealso \code{\link{M_chen_watanabe_L0}} for the L0-parameterized CW model,
#'   \code{\link{extract_growth_parameters}} for posterior extraction.
#'
#' @export
compute_k_vb_equivalent <- function(Linf, L0, Lmat, tmat, warn = TRUE) {


  # Validate inputs

  n <- length(Linf)
  if (!all(c(length(L0), length(Lmat), length(tmat)) == n)) {
    stop("All input vectors must have the same length.", call. = FALSE)
  }

  # Compute VB-equivalent k: k = (1/tmat) * ln((Linf - L0) / (Linf - Lmat))
  numerator   <- Linf - L0
  denominator <- Linf - Lmat

  # Identify invalid values (would produce NaN or Inf)
  invalid <- numerator <= 0 | denominator <= 0 | tmat <= 0 |
    is.na(numerator) | is.na(denominator) | is.na(tmat)

  if (warn && any(invalid, na.rm = TRUE)) {
    n_invalid <- sum(invalid, na.rm = TRUE)
    warning(
      sprintf(
        "%d of %d draws (%.1f%%) produced invalid k values. Common causes:\n
        - Lmat >= Linf (maturity size exceeds asymptotic size)\n
        - L0 >= Linf (birth size exceeds asymptotic size)\n
        - tmat <= 0 (non-positive maturity age)\n
        These draws will be excluded from mortality calculations.",
        n_invalid, n, 100 * n_invalid / n
      ),
      call. = FALSE
    )
  }

  k <- (1 / tmat) * log(numerator / denominator)
  k[invalid] <- NA_real_

  k
}
